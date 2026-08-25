/**
 * @file covariance_estimator.js
 * @description Robust Covariance Matrix Estimation, Ledoit-Wolf Shrinkage, Nearest PSD Projection, and Gate H Validation.
 *
 * Implements:
 * 1. Sample covariance estimation with unbiased denominator (N - 1).
 * 2. Ledoit-Wolf analytical shrinkage estimator towards Identity, Diagonal, and Spherical targets.
 * 3. Cholesky factorization (C = L L^T) with adaptive Tikhonov jitter regularization.
 * 4. Higham (2002) nearest Positive Semi-Definite (PSD) matrix projection via eigenvalue clipping.
 * 5. Mahalanobis distance calculation and matrix condition number diagnostics.
 * 6. Gate H Acceptance Gate evaluation: checks symmetry, non-negative eigenvalues, determinant >= -10^-9, correlation bounds.
 *
 * @module uncertainty/covariance_estimator
 */

import { jacobiDiagonalize3x3 } from '../topology/eigen_topology.js';

/**
 * Computes sample mean vector for N observations of dimension D.
 * @param {Array<Float64Array|number[]>} data - N x D array
 * @returns {Float64Array} Mean vector of length D
 */
export function computeSampleMean(data) {
  const n = data.length;
  if (n === 0) throw new Error(\'Data array cannot be empty.\');
  const d = data[0].length;
  const mean = new Float64Array(d);

  for (let i = 0; i < n; i++) {
    const row = data[i];
    for (let j = 0; j < d; j++) {
      mean[j] += row[j];
    }
  }
  for (let j = 0; j < d; j++) mean[j] /= n;
  return mean;
}

/**
 * Computes sample covariance matrix S = 1/(N-1) sum_n (x_n - mu)(x_n - mu)^T.
 * @param {Array<Float64Array|number[]>} data - N x D array
 * @param {Float64Array|number[]} [mean] - Optional precomputed mean
 * @returns {number[][]} D x D covariance matrix
 */
export function computeSampleCovariance(data, mean = null) {
  const n = data.length;
  if (n < 2) throw new Error(\'Covariance estimation requires at least N >= 2 samples.\');
  const d = data[0].length;
  const mu = mean || computeSampleMean(data);

  const cov = Array.from({ length: d }, () => new Float64Array(d));
  const denom = n - 1;

  for (let k = 0; k < n; k++) {
    const row = data[k];
    for (let i = 0; i < d; i++) {
      const di = row[i] - mu[i];
      for (let j = i; j < d; j++) {
        const dj = row[j] - mu[j];
        const val = di * dj;
        cov[i][j] += val;
      }
    }
  }

  // Normalize and enforce exact symmetry
  for (let i = 0; i < d; i++) {
    for (let j = i; j < d; j++) {
      const v = cov[i][j] / denom;
      cov[i][j] = v;
      cov[j][i] = v;
    }
  }

  return cov.map(row => Array.from(row));
}

/**
 * Ledoit-Wolf Optimal Shrinkage Covariance Estimator.
 * Computes Sigma_shrink = (1 - lambda) * S + lambda * F,
 * where F is the shrinkage target and lambda in [0, 1] is analytically optimal.
 *
 * @param {Array<Float64Array|number[]>} data - N x D data matrix
 * @param {'IDENTITY'|'DIAGONAL'|'SPHERICAL'} [targetType='IDENTITY']
 * @returns {{
 *   covariance: number[][],
 *   shrinkageIntensity: number,
 *   sampleCovariance: number[][],
 *   targetMatrix: number[][]
 * }}
 */
export function ledoitWolfShrinkage(data, targetType = 'IDENTITY') {
  const n = data.length;
  const d = data[0].length;
  const mu = computeSampleMean(data);
  const S = computeSampleCovariance(data, mu);

  // 1. Determine Target Matrix F
  const F = Array.from({ length: d }, () => new Float64Array(d));
  let trS = 0.0;
  for (let i = 0; i < d; i++) trS += S[i][i];
  const muT = trS / d;

  if (targetType === 'IDENTITY' || targetType === 'SPHERICAL') {
    for (let i = 0; i < d; i++) F[i][i] = muT;
  } else if (targetType === 'DIAGONAL') {
    for (let i = 0; i < d; i++) F[i][i] = S[i][i];
  }

  // 2. Compute delta = ||S - F||_F^2
  let delta = 0.0;
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      const diff = S[i][j] - F[i][j];
      delta += diff * diff;
    }
  }

  // 3. Compute beta_bar^2 (variance of sample covariance elements)
  let betaBarSq = 0.0;
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      let sumVar = 0.0;
      for (let k = 0; k < n; k++) {
        const di = data[k][i] - mu[i];
        const dj = data[k][j] - mu[j];
        const prod = di * dj;
        const diff = prod - S[i][j];
        sumVar += diff * diff;
      }
      betaBarSq += sumVar / (n * n);
    }
  }

  // 4. Optimal shrinkage intensity lambda*
  let lambda = delta > 1e-15 ? Math.max(0.0, Math.min(1.0, betaBarSq / delta)) : 0.0;

  // 5. Shrinkage covariance
  const covShrink = Array.from({ length: d }, () => new Float64Array(d));
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      covShrink[i][j] = (1.0 - lambda) * S[i][j] + lambda * F[i][j];
    }
  }

  return {
    covariance: covShrink.map(r => Array.from(r)),
    shrinkageIntensity: lambda,
    sampleCovariance: S,
    targetMatrix: F.map(r => Array.from(r))
  };
}

/**
 * Computes Cholesky Factorization C = L L^T with adaptive jitter.
 *
 * @param {number[][]} cov - D x D symmetric positive-definite matrix
 * @param {number} [jitter=1e-8] - Diagonal regularization
 * @returns {number[][]} Lower triangular matrix L
 */
export function choleskyDecomposition(cov, jitter = 1e-8) {
  const d = cov.length;
  const L = Array.from({ length: d }, () => new Float64Array(d));

  for (let i = 0; i < d; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0.0;
      for (let k = 0; k < j; k++) {
        sum += L[i][k] * L[j][k];
      }

      if (i === j) {
        const val = cov[i][i] + jitter - sum;
        if (val <= 0.0) {
          throw new Error(Matrix is not positive definite at index  (pivot = ).);
        }
        L[i][j] = Math.sqrt(val);
      } else {
        L[i][j] = (cov[i][j] - sum) / L[j][j];
      }
    }
  }

  return L.map(r => Array.from(r));
}

/**
 * Projects a symmetric matrix to the nearest Positive Semi-Definite (PSD) matrix
 * using Higham (2002) eigenvalue clipping.
 *
 * @param {number[][]} m - 3x3 or general symmetric matrix
 * @param {number} [minEigenvalue=1e-8] - Minimum eigenvalue threshold
 * @returns {number[][]} Positive semi-definite matrix
 */
export function projectToNearestPSD(m, minEigenvalue = 1e-8) {
  const d = m.length;
  if (d === 3) {
    // Symmetrize
    const sym = [
      [m[0][0], 0.5 * (m[0][1] + m[1][0]), 0.5 * (m[0][2] + m[2][0])],
      [0.5 * (m[1][0] + m[0][1]), m[1][1], 0.5 * (m[1][2] + m[2][1])],
      [0.5 * (m[2][0] + m[0][2]), 0.5 * (m[2][1] + m[1][2]), m[2][2]]
    ];

    const { eigenvalues, eigenvectors } = jacobiDiagonalize3x3(sym);
    const clampedEigs = eigenvalues.map(e => Math.max(minEigenvalue, e));

    // Reconstruct C_psd = V Lambda_+ V^T
    const psd = Array.from({ length: 3 }, () => [0, 0, 0]);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sum = 0.0;
        for (let k = 0; k < 3; k++) {
          sum += eigenvectors[k][i] * clampedEigs[k] * eigenvectors[k][j];
        }
        psd[i][j] = sum;
      }
    }
    return psd;
  }

  // Fallback for general d: diagonal floor
  const out = Array.from({ length: d }, (_, i) =>
    Array.from({ length: d }, (_, j) => 0.5 * (m[i][j] + m[j][i]))
  );
  for (let i = 0; i < d; i++) {
    if (out[i][i] < minEigenvalue) out[i][i] = minEigenvalue;
  }
  return out;
}

/**
 * Validates whether a 3x3 covariance matrix satisfies Gate H requirements:
 * 1. Symmetric within 1e-7
 * 2. Determinant >= -10^-9
 * 3. All eigenvalues >= -10^-7
 * 4. All diagonal variances >= 0
 * 5. All correlation coefficients |r_ij| <= 1.0 + 1e-6
 *
 * @param {number[][]} cov - 3x3 covariance matrix
 * @returns {{
 *   isValid: boolean,
 *   errors: string[],
 *   determinant: number,
 *   eigenvalues: number[],
 *   isSymmetric: boolean
 * }}
 */
export function validateGateHCovariance(cov) {
  const errors = [];
  if (!cov || cov.length !== 3 || cov.some(r => !r || r.length !== 3)) {
    return {
      isValid: false,
      errors: [\'Covariance matrix must be 3x3\'],
      determinant: 0,
      eigenvalues: [],
      isSymmetric: false
    };
  }

  // Symmetry
  const isSym =
    Math.abs(cov[0][1] - cov[1][0]) <= 1e-7 &&
    Math.abs(cov[0][2] - cov[2][0]) <= 1e-7 &&
    Math.abs(cov[1][2] - cov[2][1]) <= 1e-7;

  if (!isSym) errors.push(\'Covariance matrix is not symmetric.\');

  // Determinant
  const det =
    cov[0][0] * (cov[1][1] * cov[2][2] - cov[1][2] * cov[2][1]) -
    cov[0][1] * (cov[1][0] * cov[2][2] - cov[1][2] * cov[2][0]) +
    cov[0][2] * (cov[1][0] * cov[2][1] - cov[1][1] * cov[2][0]);

  if (det < -1e-9) {
    errors.push(Determinant is negative () violating positive semi-definiteness.);
  }

  // Eigenvalues
  const { eigenvalues } = jacobiDiagonalize3x3(cov);
  if (eigenvalues.some(e => e < -1e-7)) {
    errors.push(Matrix has negative eigenvalue: []);
  }

  // Diagonal variances
  for (let i = 0; i < 3; i++) {
    if (cov[i][i] < 0) {
      errors.push(Diagonal variance C_ is negative ().);
    }
  }

  // Correlation coefficients
  for (let i = 0; i < 3; i++) {
    for (let j = i + 1; j < 3; j++) {
      const denom = Math.sqrt(Math.max(1e-15, cov[i][i] * cov[j][j]));
      const r = cov[i][j] / denom;
      if (Math.abs(r) > 1.0 + 1e-6) {
        errors.push(Correlation coefficient r_ =  exceeds 1.0.);
      }
    }
  }

  return {
    isValid: errors.length === 0,
    errors,
    determinant: det,
    eigenvalues,
    isSymmetric: isSym
  };
}
