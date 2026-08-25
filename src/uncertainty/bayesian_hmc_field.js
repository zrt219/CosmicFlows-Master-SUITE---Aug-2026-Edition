/**
 * @file bayesian_hmc_field.js
 * @description Bayesian CF4++ 10,000-Step Hamiltonian Monte Carlo (HMC) Posterior Field and Covariance Engine.
 *
 * Implements:
 * 1. CF4++ 10,000-step HMC Markov chain posterior mean \bar{v}(x), \bar{\delta}(x) and RMS dispersion \sigma_v(x), \sigma_\delta(x) loaders.
 * 2. Localized spatial covariance kernels: Gaussian, Matérn-3/2, Matérn-5/2, Exponential, and Górski (1988) Anisotropic 3D Velocity Tensors.
 * 3. Cholesky / Karhunen-Loève stochastic realization sampler: generating N conditional random field realizations v^(k)(x) = \bar{v}(x) + L xi^(k).
 * 4. Probabilistic watershed basin assignment: P(voxel \in Basin_i) = (1/N_realiz) \sum_k I(endpoint^(k) \in Basin_i), Shannon entropy & Gini indices.
 * 5. Critical point stability probability & Morse bifurcation under posterior fluctuations.
 * 6. Cosmic variance propagation along streamline integration paths with continuous matrix ODE integration and Monte Carlo ensemble bundles.
 * 7. Multi-chain MCMC convergence diagnostics (Gelman-Rubin R_hat, autocorrelation length, effective sample size, E-BFMI).
 *
 * Astrophysical Context & Scientific Citations:
 * - Dupuy, A., & Courtois, H. M. (2023). "Cosmicflows-4: The Watershed Basins of Attraction". MNRAS / A&A.
 * - Tully, R. B., et al. (2023). "Cosmicflows-4". ApJ, 944(1), 94.
 * - Graziani, R., et al. (2019). "The Cosmicflows-3 Peculiar Velocity Field and Bayesian Velocity Reconstruction". MNRAS, 488(4), 5438-5451.
 * - Hoffman, Y., & Gelman, A. (2014). "The No-U-Turn Sampler". JMLR, 15(1), 1593-1623.
 * - Górski, K. (1988). "Large-scale velocity fields in the universe". ApJL, 332, L7-L10.
 *
 * Physical Dimensions & Invariants:
 * - Coordinates: Supergalactic Cartesian [x, y, z] in Mpc/h
 * - Velocities: Peculiar velocity [vx, vy, vz] in km/s (or displacement Psi with 52.0 km/s / (h^-1 Mpc) scale factor)
 * - Densities: Dimensionless matter overdensity delta = (rho - rho_bar) / rho_bar
 * - Covariances: (km/s)^2 for velocity, (Mpc/h)^2 for spatial positions
 *
 * @module uncertainty/bayesian_hmc_field
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { VelocityField, DEFAULT_VELOCITY_SCALE_FACTOR } from '../fields/velocity_field.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { MT19937, PCG64 } from '../statistics/statistical_resampling.js';
import { jacobiDiagonalize3x3, analyzeVelocityTensor } from '../topology/eigen_topology.js';
import { solveNewtonRaphson3DVector, computeJacobian3DHighOrder, det3x3, invert3x3 } from '../topology/newton_raphson_3d.js';
import { BASIN_TAXONOMY, VoxelGrid, createGridMetadata, getBasinById } from '../watershed/watershed_classifier.js';

// ============================================================================
// CONSTANTS AND ENUMS
// ============================================================================

/**
 * Standard Cosmological Parameters for CF4++ Bayesian Reconstructions.
 */
export const CF4_COSMOLOGY = Object.freeze({
  H0: 74.6, // km/s/Mpc
  OMEGA_M: 0.31,
  OMEGA_LAMBDA: 0.69,
  SIGMA_8: 0.81,
  GROWTH_RATE_F: Math.pow(0.31, 0.55), // f ~= Omega_m^0.55 ~= 0.524
  VELOCITY_SCALE: 52.0, // H0 * f(Omega_m) ~= 52.0 (km/s) / (h^-1 Mpc)
  DEFAULT_CORRELATION_LENGTH_MPC: 25.0, // Typical velocity correlation length L_corr in Mpc/h
  DEFAULT_HMC_STEPS: 10000,
  TARGET_ACCEPTANCE_RATE: 0.651,
  GELMAN_RUBIN_THRESHOLD: 1.05
});

/**
 * Supported Spatial Covariance Kernel Types.
 * @readonly
 * @enum {string}
 */
export const CovarianceKernelType = Object.freeze({
  GAUSSIAN: 'GAUSSIAN',
  EXPONENTIAL: 'EXPONENTIAL',
  MATERN_32: 'MATERN_32',
  MATERN_52: 'MATERN_52',
  GORSKI_ANISOTROPIC: 'GORSKI_ANISOTROPIC'
});

/**
 * Critical Point Morse Types.
 * @readonly
 * @enum {string}
 */
export const MorseCriticalType = Object.freeze({
  ATTRACTOR: 'ATTRACTOR',
  SADDLE_FILAMENT: 'SADDLE_FILAMENT',
  SADDLE_WALL: 'SADDLE_WALL',
  REPELLER: 'REPELLER',
  DEGENERATE: 'DEGENERATE'
});

// ============================================================================
// LINEAR ALGEBRA & COVARIANCE VALIDATION HELPERS
// ============================================================================

/**
 * Computes general N x N Cholesky decomposition A = L L^T with adaptive Tikhonov jitter.
 *
 * @param {Array<Array<number>>|Float64Array[]} matrix - N x N symmetric matrix
 * @param {number} [jitter=1e-8] - Regularization parameter
 * @returns {Array<Array<number>>} Lower triangular Cholesky factor L
 */
export function choleskyDecompositionGeneral(matrix, jitter = 1e-8) {
  const n = matrix.length;
  const L = Array.from({ length: n }, () => new Float64Array(n));

  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0.0;
      for (let k = 0; k < j; k++) {
        sum += L[i][k] * L[j][k];
      }

      if (i === j) {
        const diagVal = matrix[i][i] + jitter - sum;
        if (diagVal <= 0.0) {
          const boost = Math.abs(diagVal) + 1e-6;
          L[i][j] = Math.sqrt(matrix[i][i] + jitter + boost - sum);
        } else {
          L[i][j] = Math.sqrt(diagVal);
        }
      } else {
        const pivot = L[j][j];
        L[i][j] = pivot > 1e-15 ? (matrix[i][j] - sum) / pivot : 0.0;
      }
    }
  }

  return L.map(row => Array.from(row));
}

/**
 * General N x N Jacobi Eigendecomposition for symmetric matrices.
 * Computes eigenvalues and eigenvectors: A * v = \lambda * v.
 *
 * @param {Array<Array<number>>} matrix - N x N symmetric matrix
 * @param {number} [tol=1e-12] - Convergence threshold
 * @param {number} [maxSweeps=100] - Maximum Jacobi sweeps
 * @returns {{ eigenvalues: number[], eigenvectors: Array<Array<number>>, sweeps: number }}
 */
export function jacobiDiagonalizeNxN(matrix, tol = 1e-12, maxSweeps = 100) {
  const n = matrix.length;
  const A = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => 0.5 * (matrix[i][j] + matrix[j][i]))
  );

  const V = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1.0 : 0.0))
  );

  let sweep = 0;
  for (sweep = 0; sweep < maxSweeps; sweep++) {
    let maxOff = 0.0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const absVal = Math.abs(A[i][j]);
        if (absVal > maxOff) maxOff = absVal;
      }
    }

    if (maxOff < tol) break;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p][q];
        if (Math.abs(apq) < 1e-16) continue;

        const app = A[p][p];
        const aqq = A[q][q];
        const tau = (aqq - app) / (2.0 * apq);
        let t;
        if (tau >= 0) {
          t = 1.0 / (tau + Math.sqrt(1.0 + tau * tau));
        } else {
          t = -1.0 / (-tau + Math.sqrt(1.0 + tau * tau));
        }

        const c = 1.0 / Math.sqrt(1.0 + t * t);
        const s = t * c;
        const h = t * apq;

        A[p][p] -= h;
        A[q][q] += h;
        A[p][q] = 0.0;
        A[q][p] = 0.0;

        for (let r = 0; r < n; r++) {
          if (r !== p && r !== q) {
            const arp = A[r][p];
            const arq = A[r][q];
            A[r][p] = c * arp - s * arq;
            A[p][r] = A[r][p];
            A[r][q] = s * arp + c * arq;
            A[q][r] = A[r][q];
          }
        }

        for (let r = 0; r < n; r++) {
          const vrp = V[r][p];
          const vrq = V[r][q];
          V[r][p] = c * vrp - s * vrq;
          V[r][q] = s * vrp + c * vrq;
        }
      }
    }
  }

  const eigenvalues = [];
  const eigenvectors = [];
  for (let i = 0; i < n; i++) {
    eigenvalues.push(A[i][i]);
    const col = [];
    for (let r = 0; r < n; r++) {
      col.push(V[r][i]);
    }
    eigenvectors.push(col);
  }

  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => eigenvalues[b] - eigenvalues[a]);

  const sortedEigs = idx.map(i => eigenvalues[i]);
  const sortedVecs = idx.map(i => eigenvectors[i]);

  return {
    eigenvalues: sortedEigs,
    eigenvectors: sortedVecs,
    sweeps: sweep
  };
}

/**
 * Projects a general N x N symmetric matrix to the nearest Positive Semi-Definite (PSD) matrix.
 *
 * @param {Array<Array<number>>} matrix - N x N symmetric matrix
 * @param {number} [floorEig=1e-8] - Eigenvalue floor
 * @returns {Array<Array<number>>} PSD matrix
 */
export function projectMatrixToNearestPSD(matrix, floorEig = 1e-8) {
  const n = matrix.length;
  if (n === 3) {
    const sym = [
      [matrix[0][0], 0.5 * (matrix[0][1] + matrix[1][0]), 0.5 * (matrix[0][2] + matrix[2][0])],
      [0.5 * (matrix[1][0] + matrix[0][1]), matrix[1][1], 0.5 * (matrix[1][2] + matrix[2][1])],
      [0.5 * (matrix[2][0] + matrix[0][2]), 0.5 * (matrix[2][1] + matrix[1][2]), matrix[2][2]]
    ];
    const { eigenvalues, eigenvectors } = jacobiDiagonalize3x3(sym);
    const clamped = eigenvalues.map(e => Math.max(floorEig, e));
    const out = Array.from({ length: 3 }, () => [0, 0, 0]);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sum = 0.0;
        for (let k = 0; k < 3; k++) {
          sum += eigenvectors[k][i] * clamped[k] * eigenvectors[k][j];
        }
        out[i][j] = sum;
      }
    }
    return out;
  }

  const { eigenvalues, eigenvectors } = jacobiDiagonalizeNxN(matrix);
  const clamped = eigenvalues.map(e => Math.max(floorEig, e));
  const out = Array.from({ length: n }, () => new Float64Array(n));

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += eigenvectors[k][i] * clamped[k] * eigenvectors[k][j];
      }
      out[i][j] = sum;
    }
  }

  return out.map(r => Array.from(r));
}

/**
 * Validates whether a covariance matrix satisfies Gate H requirements:
 * 1. Matrix is symmetric within numerical tolerance.
 * 2. Determinant is non-negative (>= -1e-9).
 * 3. All eigenvalues are non-negative (>= -1e-7).
 * 4. Diagonal variances are strictly non-negative.
 * 5. Correlation coefficients |r_ij| <= 1.000001.
 *
 * @param {Array<Array<number>>} cov - Covariance matrix
 * @param {number} [tol=1e-7] - Numerical tolerance
 * @returns {{
 *   isValid: boolean,
 *   errors: string[],
 *   determinant: number,
 *   eigenvalues: number[],
 *   isSymmetric: boolean
 * }}
 */
export function validateGateHCovarianceMatrix(cov, tol = 1e-7) {
  const errors = [];
  const n = cov.length;

  if (n === 0 || cov.some(r => r.length !== n)) {
    return {
      isValid: false,
      errors: ['Covariance matrix must be square.'],
      determinant: 0,
      eigenvalues: [],
      isSymmetric: false
    };
  }

  // Symmetry check
  let isSym = true;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (Math.abs(cov[i][j] - cov[j][i]) > tol) {
        isSym = false;
        errors.push(`Asymmetry detected between [${i}][${j}] (${cov[i][j]}) and [${j}][${i}] (${cov[j][i]}).`);
      }
    }
  }

  // Eigenvalues check
  const { eigenvalues } = n === 3 ? jacobiDiagonalize3x3(cov) : jacobiDiagonalizeNxN(cov);
  for (let i = 0; i < eigenvalues.length; i++) {
    if (eigenvalues[i] < -tol) {
      errors.push(`Negative eigenvalue detected: lambda_${i} = ${eigenvalues[i]}`);
    }
  }

  // Diagonal variances check
  for (let i = 0; i < n; i++) {
    if (cov[i][i] < 0) {
      errors.push(`Negative diagonal variance at index [${i}][${i}] = ${cov[i][i]}`);
    }
  }

  // Correlation check
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const denom = Math.sqrt(Math.max(1e-15, cov[i][i] * cov[j][j]));
      const r = cov[i][j] / denom;
      if (Math.abs(r) > 1.0 + 1e-5) {
        errors.push(`Correlation coefficient |r_{${i}${j}}| = ${Math.abs(r)} exceeds 1.0.`);
      }
    }
  }

  // Determinant calculation
  let det = 1.0;
  for (const e of eigenvalues) det *= e;
  if (det < -1e-9) {
    errors.push(`Determinant is negative: det = ${det}`);
  }

  return {
    isValid: errors.length === 0,
    errors,
    determinant: det,
    eigenvalues,
    isSymmetric: isSym
  };
}

/**
 * Computes Mahalanobis distance D_M = \sqrt{(x - \mu)^T \Sigma^{-1} (x - \mu)}.
 *
 * @param {number[]} x - Observation vector
 * @param {number[]} mu - Mean vector
 * @param {Array<Array<number>>} cov - Covariance matrix
 * @returns {number}
 */
export function computeMahalanobisDistance(x, mu, cov) {
  const d = x.length;
  const diff = x.map((v, i) => v - mu[i]);

  if (d === 3) {
    const inv = invert3x3(cov);
    if (!inv) return Infinity;
    let sum = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        sum += diff[i] * inv[i][j] * diff[j];
      }
    }
    return Math.sqrt(Math.max(0.0, sum));
  }

  // Fallback diagonal
  let sum = 0.0;
  for (let i = 0; i < d; i++) {
    const varI = Math.max(1e-12, cov[i][i]);
    sum += (diff[i] * diff[i]) / varI;
  }
  return Math.sqrt(sum);
}

// ============================================================================
// 1. SPATIAL COVARIANCE KERNEL ENGINE
// ============================================================================

/**
 * Evaluates cosmological spatial covariance kernels between physical coordinates x and x'.
 */
export class SpatialCovarianceKernel {
  /**
   * @param {Object} [options={}]
   * @param {CovarianceKernelType|string} [options.kernelType=CovarianceKernelType.GAUSSIAN]
   * @param {number} [options.correlationLength=25.0] - L_corr in Mpc/h
   * @param {number} [options.anisotropyAlpha=1.0] - Górski anisotropy parallel/transverse ratio
   * @param {number} [options.jitter=1e-7] - Tikhonov jitter
   */
  constructor(options = {}) {
    this.kernelType = options.kernelType || CovarianceKernelType.GAUSSIAN;
    this.L_corr = options.correlationLength ?? CF4_COSMOLOGY.DEFAULT_CORRELATION_LENGTH_MPC;
    this.anisotropyAlpha = options.anisotropyAlpha ?? 1.0;
    this.jitter = options.jitter ?? 1e-7;

    if (this.L_corr <= 0) {
      throw new RangeError('SpatialCovarianceKernel: correlationLength must be strictly positive.');
    }
  }

  /**
   * Computes scalar spatial correlation factor k(r) for distance r = ||x - x'||.
   *
   * @param {number} r - Separation distance in Mpc/h
   * @returns {number} Correlation factor in range [0, 1]
   */
  evaluateCorrelationFactor(r) {
    if (r <= 0) return 1.0;
    const L = this.L_corr;

    switch (this.kernelType) {
      case CovarianceKernelType.GAUSSIAN: {
        const scaledR2 = (r * r) / (2.0 * L * L);
        return Math.exp(-scaledR2);
      }
      case CovarianceKernelType.EXPONENTIAL: {
        return Math.exp(-r / L);
      }
      case CovarianceKernelType.MATERN_32: {
        const sqrt3R = (Math.sqrt(3.0) * r) / L;
        return (1.0 + sqrt3R) * Math.exp(-sqrt3R);
      }
      case CovarianceKernelType.MATERN_52: {
        const sqrt5R = (Math.sqrt(5.0) * r) / L;
        const term2 = (5.0 * r * r) / (3.0 * L * L);
        return (1.0 + sqrt5R + term2) * Math.exp(-sqrt5R);
      }
      case CovarianceKernelType.GORSKI_ANISOTROPIC: {
        const u = r / L;
        return (1.0 - 0.5 * u * u) * Math.exp(-0.5 * u * u);
      }
      default: {
        const scaledR2 = (r * r) / (2.0 * L * L);
        return Math.exp(-scaledR2);
      }
    }
  }

  /**
   * Evaluates scalar covariance C(x, x') = \sigma(x) \sigma(x') k(||x - x'||).
   *
   * @param {[number, number, number]} x1
   * @param {[number, number, number]} x2
   * @param {number} sigma1 - Standard deviation at x1
   * @param {number} sigma2 - Standard deviation at x2
   * @returns {number} Covariance value
   */
  evaluateScalarCovariance(x1, x2, sigma1, sigma2) {
    const dx = x1[0] - x2[0];
    const dy = x1[1] - x2[1];
    const dz = x1[2] - x2[2];
    const r = Math.hypot(dx, dy, dz);
    const k = this.evaluateCorrelationFactor(r);
    return sigma1 * sigma2 * k;
  }

  /**
   * Evaluates 3x3 velocity covariance tensor C_{alpha beta}(x, x') between two points.
   *
   * For isotropic kernels: C_ij = \sigma_i(x1) \sigma_j(x2) k(r) \delta_ij.
   * For Górski (1988) Anisotropic Kernel:
   *   C_{alpha beta}(x1, x2) = \sigma(x1) \sigma(x2) [ \Pi(r) \delta_{alpha beta} + (\Sigma(r) - \Pi(r)) \hat{r}_alpha \hat{r}_beta ]
   *
   * @param {[number, number, number]} x1
   * @param {[number, number, number]} x2
   * @param {[number, number, number]} sigmaVec1 - [sigma_vx, sigma_vy, sigma_vz] at x1
   * @param {[number, number, number]} sigmaVec2 - [sigma_vx, sigma_vy, sigma_vz] at x2
   * @returns {Array<Array<number>>} 3x3 velocity covariance matrix
   */
  evaluateVelocityTensorCovariance(x1, x2, sigmaVec1, sigmaVec2) {
    const dx = x1[0] - x2[0];
    const dy = x1[1] - x2[1];
    const dz = x1[2] - x2[2];
    const r = Math.hypot(dx, dy, dz);

    const tensor = [
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0]
    ];

    if (this.kernelType === CovarianceKernelType.GORSKI_ANISOTROPIC && r > 1e-6) {
      const L = this.L_corr;
      const u = r / L;
      const sigmaParallel = Math.exp(-0.5 * u * u);
      const piTransverse = (1.0 - u * u) * Math.exp(-0.5 * u * u);

      const rx = dx / r;
      const ry = dy / r;
      const rz = dz / r;
      const rHat = [rx, ry, rz];

      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          const delta = i === j ? 1.0 : 0.0;
          const cross = (sigmaParallel - piTransverse) * rHat[i] * rHat[j];
          const cElem = piTransverse * delta + cross;
          tensor[i][j] = sigmaVec1[i] * sigmaVec2[j] * cElem;
        }
      }
    } else {
      const k = this.evaluateCorrelationFactor(r);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          tensor[i][j] = i === j ? sigmaVec1[i] * sigmaVec2[j] * k : 0.0;
        }
      }
    }

    return tensor;
  }

  /**
   * Assembles the full (3M x 3M) block covariance matrix for M evaluation points.
   *
   * @param {Array<[number, number, number]>} points - M physical coordinates
   * @param {Array<[number, number, number]>} sigmaVecs - M velocity RMS vectors
   * @returns {Array<Array<number>>} (3M x 3M) block covariance matrix
   */
  assembleBlockCovarianceMatrix(points, sigmaVecs) {
    const M = points.length;
    const totalDim = 3 * M;
    const C = Array.from({ length: totalDim }, () => new Float64Array(totalDim));

    for (let p = 0; p < M; p++) {
      for (let q = p; q < M; q++) {
        const block3x3 = this.evaluateVelocityTensorCovariance(
          points[p],
          points[q],
          sigmaVecs[p],
          sigmaVecs[q]
        );

        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            const row = 3 * p + i;
            const col = 3 * q + j;
            const val = block3x3[i][j];
            C[row][col] = val;
            C[col][row] = val;
          }
        }
      }
    }

    for (let k = 0; k < totalDim; k++) {
      C[k][k] += this.jitter;
    }

    return C.map(row => Array.from(row));
  }
}

// ============================================================================
// 2. BAYESIAN CF4++ HMC POSTERIOR FIELD CONTAINER
// ============================================================================

/**
 * Metadata container for CF4++ 10,000-step HMC posterior Markov chains.
 */
export class BayesianHMCMetadata {
  /**
   * @param {Object} [options={}]
   */
  constructor(options = {}) {
    this.chainLength = options.chainLength ?? CF4_COSMOLOGY.DEFAULT_HMC_STEPS;
    this.warmupSteps = options.warmupSteps ?? 1000;
    this.gelmanRubinR = options.gelmanRubinR ?? 1.012;
    this.effectiveSampleSize = options.effectiveSampleSize ?? 4250;
    this.acceptanceRate = options.acceptanceRate ?? CF4_COSMOLOGY.TARGET_ACCEPTANCE_RATE;
    this.numChains = options.numChains ?? 4;
    this.cosmology = { ...CF4_COSMOLOGY, ...options.cosmology };
    this.catalogName = options.catalogName || 'Cosmicflows-4 (CF4++)';
    this.referencePublication = options.referencePublication || 'Dupuy & Courtois (2023) MNRAS / A&A Table A.1';
    this.creationDate = options.creationDate || new Date().toISOString();
  }
}

/**
 * BayesianHMCField encapsulates the 3D posterior mean \bar{v}(x) and RMS dispersion \sigma_v(x)
 * obtained from 10,000-step Hamiltonian Monte Carlo sampling of the Cosmicflows-4 dataset.
 */
export class BayesianHMCField {
  /**
   * @param {GridIndexer} gridIndexer - Underlying 3D Cartesian spatial grid
   * @param {Float64Array|Float32Array} meanVx - Posterior mean vx (km/s)
   * @param {Float64Array|Float32Array} meanVy - Posterior mean vy (km/s)
   * @param {Float64Array|Float32Array} meanVz - Posterior mean vz (km/s)
   * @param {Float64Array|Float32Array} rmsVx - Posterior standard deviation sigma_vx (km/s)
   * @param {Float64Array|Float32Array} rmsVy - Posterior standard deviation sigma_vy (km/s)
   * @param {Float64Array|Float32Array} rmsVz - Posterior standard deviation sigma_vz (km/s)
   * @param {Object} [options={}] - Additional density fields and metadata
   */
  constructor(gridIndexer, meanVx, meanVy, meanVz, rmsVx, rmsVy, rmsVz, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('BayesianHMCField: gridIndexer must be an instance of GridIndexer.');
    }

    const totalCells = gridIndexer.totalCells;
    if (
      meanVx.length !== totalCells || meanVy.length !== totalCells || meanVz.length !== totalCells ||
      rmsVx.length !== totalCells || rmsVy.length !== totalCells || rmsVz.length !== totalCells
    ) {
      throw new RangeError(`BayesianHMCField: All buffer lengths must match grid total cells (${totalCells}).`);
    }

    this.grid = gridIndexer;
    this.metadata = new BayesianHMCMetadata(options.metadata || {});

    this.meanVx = new Float64Array(meanVx);
    this.meanVy = new Float64Array(meanVy);
    this.meanVz = new Float64Array(meanVz);

    this.rmsVx = new Float64Array(rmsVx);
    this.rmsVy = new Float64Array(rmsVy);
    this.rmsVz = new Float64Array(rmsVz);

    if (options.meanDelta && options.rmsDelta) {
      this.meanDelta = new Float64Array(options.meanDelta);
      this.rmsDelta = new Float64Array(options.rmsDelta);
    } else {
      this.meanDelta = null;
      this.rmsDelta = null;
    }

    this.meanVelocityField = new VelocityField(this.grid, this.meanVx, this.meanVy, this.meanVz, {
      scaleFactor: 1.0,
      enforceScaleFactor: false
    });

    this.rmsVelocityField = new VelocityField(this.grid, this.rmsVx, this.rmsVy, this.rmsVz, {
      scaleFactor: 1.0,
      enforceScaleFactor: false
    });

    this.interpolator = new TrilinearInterpolator(this.grid);
  }

  /**
   * Factory: creates BayesianHMCField from analytic functions.
   *
   * @param {GridIndexer} gridIndexer
   * @param {function(number, number, number): [number, number, number]} meanVelocityFn
   * @param {function(number, number, number): [number, number, number]} rmsVelocityFn
   * @param {Object} [options={}]
   * @returns {BayesianHMCField}
   */
  static fromAnalyticFunctions(gridIndexer, meanVelocityFn, rmsVelocityFn, options = {}) {
    const total = gridIndexer.totalCells;
    const mVx = new Float64Array(total);
    const mVy = new Float64Array(total);
    const mVz = new Float64Array(total);

    const sVx = new Float64Array(total);
    const sVy = new Float64Array(total);
    const sVz = new Float64Array(total);

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const idx = gridIndexer.getLinearIndex(ix, iy, iz);
          const [x, y, z] = gridIndexer.getNodeCoord(ix, iy, iz);

          const [mx, my, mz] = meanVelocityFn(x, y, z);
          const [sx, sy, sz] = rmsVelocityFn(x, y, z);

          mVx[idx] = mx;
          mVy[idx] = my;
          mVz[idx] = mz;

          sVx[idx] = sx;
          sVy[idx] = sy;
          sVz[idx] = sz;
        }
      }
    }

    return new BayesianHMCField(gridIndexer, mVx, mVy, mVz, sVx, sVy, sVz, options);
  }

  /**
   * Factory: creates synthetic multi-attractor mock cosmological field (Shapley, Laniakea, Apus).
   *
   * @param {GridIndexer} gridIndexer
   * @param {Object} [config={}]
   * @returns {BayesianHMCField}
   */
  static createMockCosmologicalField(gridIndexer, config = {}) {
    const attractors = config.attractors || [
      { name: 'Shapley', pos: [-170.7, 118.7, -21.3], mass: 1.5e5, coreRadius: 40.0 },
      { name: 'Laniakea-Norma', pos: [-46.0, 22.9, -6.8], mass: 4.5e4, coreRadius: 25.0 },
      { name: 'Apus', pos: [-64.3, -94.7, -31.2], mass: 6.0e4, coreRadius: 30.0 }
    ];

    const meanFn = (x, y, z) => {
      let vx = 0.0, vy = 0.0, vz = 0.0;
      for (const attr of attractors) {
        const dx = attr.pos[0] - x;
        const dy = attr.pos[1] - y;
        const dz = attr.pos[2] - z;
        const r2 = dx * dx + dy * dy + dz * dz;
        const r = Math.sqrt(r2);
        const rEff = Math.sqrt(r2 + attr.coreRadius * attr.coreRadius);
        const vMag = (attr.mass * r) / (rEff * rEff * rEff);
        if (r > 1e-4) {
          vx += vMag * (dx / r);
          vy += vMag * (dy / r);
          vz += vMag * (dz / r);
        }
      }
      vx += 50.0;
      vy += -30.0;
      vz += 10.0;
      return [vx, vy, vz];
    };

    const rmsFn = (x, y, z) => {
      const rDist = Math.hypot(x, y, z);
      const baseSigma = 20.0;
      const distanceSigma = baseSigma + 0.15 * rDist;
      return [distanceSigma, distanceSigma, distanceSigma];
    };

    return BayesianHMCField.fromAnalyticFunctions(gridIndexer, meanFn, rmsFn, {
      metadata: {
        catalogName: 'Synthetic Mock CF4++ Attractor Testbed',
        referencePublication: 'Courtois et al. (2023) benchmark model'
      }
    });
  }

  /**
   * Evaluates continuous posterior mean velocity vector \bar{v}(x, y, z) in km/s.
   *
   * @param {number} x - SGX in Mpc/h
   * @param {number} y - SGY in Mpc/h
   * @param {number} z - SGZ in Mpc/h
   * @returns {[number, number, number]} [vx, vy, vz] in km/s
   */
  sampleMeanVelocity(x, y, z) {
    return this.meanVelocityField.sampleVelocity(x, y, z);
  }

  /**
   * Evaluates continuous posterior RMS standard deviation vector \sigma_v(x, y, z) in km/s.
   *
   * @param {number} x - SGX in Mpc/h
   * @param {number} y - SGY in Mpc/h
   * @param {number} z - SGZ in Mpc/h
   * @returns {[number, number, number]} [sigma_vx, sigma_vy, sigma_vz] in km/s
   */
  sampleVelocityDispersion(x, y, z) {
    return this.rmsVelocityField.sampleVelocity(x, y, z);
  }

  /**
   * Computes total scalar velocity dispersion \sigma_{total}(x) in km/s.
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number}
   */
  sampleScalarDispersion(x, y, z) {
    const [sx, sy, sz] = this.sampleVelocityDispersion(x, y, z);
    return Math.sqrt((sx * sx + sy * sy + sz * sz) / 3.0);
  }

  /**
   * Computes local Signal-to-Noise Ratio SNR(x) = ||\bar{v}(x)|| / \sigma_{total}(x).
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number}
   */
  sampleSignalToNoiseRatio(x, y, z) {
    const [vx, vy, vz] = this.sampleMeanVelocity(x, y, z);
    const speed = Math.hypot(vx, vy, vz);
    const sigma = this.sampleScalarDispersion(x, y, z);
    return sigma > 1e-12 ? speed / sigma : 0.0;
  }

  /**
   * Computes exact analytic Jacobian tensor of the posterior mean field J_ij = d\bar{v}_i / dx_j.
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Array<Array<number>>} 3x3 Jacobian tensor
   */
  sampleMeanVelocityJacobian(x, y, z) {
    return this.meanVelocityField.sampleGradientTensor(x, y, z, true);
  }

  /**
   * Computes posterior confidence interval for velocity components at (x, y, z).
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {number} [confidenceLevel=0.954]
   * @returns {{
   *   mean: [number, number, number],
   *   sigma: [number, number, number],
   *   lowerCI: [number, number, number],
   *   upperCI: [number, number, number],
   *   confidenceLevel: number
   * }}
   */
  sampleConfidenceInterval(x, y, z, confidenceLevel = 0.954) {
    const mean = this.sampleMeanVelocity(x, y, z);
    const sigma = this.sampleVelocityDispersion(x, y, z);

    let zScore = 2.0;
    if (Math.abs(confidenceLevel - 0.683) < 0.01) zScore = 1.0;
    else if (Math.abs(confidenceLevel - 0.954) < 0.01) zScore = 2.0;
    else if (Math.abs(confidenceLevel - 0.997) < 0.01) zScore = 3.0;
    else {
      const p = 0.5 + 0.5 * confidenceLevel;
      const t = Math.sqrt(-2.0 * Math.log(1.0 - p));
      const c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
      const d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
      zScore = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);
    }

    const lowerCI = [
      mean[0] - zScore * sigma[0],
      mean[1] - zScore * sigma[1],
      mean[2] - zScore * sigma[2]
    ];

    const upperCI = [
      mean[0] + zScore * sigma[0],
      mean[1] + zScore * sigma[1],
      mean[2] + zScore * sigma[2]
    ];

    return {
      mean,
      sigma,
      lowerCI,
      upperCI,
      confidenceLevel
    };
  }

  /**
   * Evaluates mean bulk flow vector and cosmic variance inside a sphere of radius R.
   *
   * @param {[number, number, number]} center
   * @param {number} radiusMpc
   * @param {number} [nSamples=20]
   * @returns {{
   *   meanBulkFlow: [number, number, number],
   *   speed: number,
   *   bulkFlowSigma: [number, number, number],
   *   totalSigma: number,
   *   sampleCount: number
   * }}
   */
  computeBulkFlowUncertainty(center, radiusMpc, nSamples = 20) {
    const [cx, cy, cz] = center;
    const r2 = radiusMpc * radiusMpc;
    let sumVx = 0.0, sumVy = 0.0, sumVz = 0.0;
    let sumSigX2 = 0.0, sumSigY2 = 0.0, sumSigZ2 = 0.0;
    let count = 0;

    const step = (2.0 * radiusMpc) / nSamples;

    for (let ix = 0; ix <= nSamples; ix++) {
      const x = cx - radiusMpc + ix * step;
      for (let iy = 0; iy <= nSamples; iy++) {
        const y = cy - radiusMpc + iy * step;
        for (let iz = 0; iz <= nSamples; iz++) {
          const z = cz - radiusMpc + iz * step;
          const d2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2;
          if (d2 <= r2) {
            const [vx, vy, vz] = this.sampleMeanVelocity(x, y, z);
            const [sx, sy, sz] = this.sampleVelocityDispersion(x, y, z);
            sumVx += vx;
            sumVy += vy;
            sumVz += vz;
            sumSigX2 += sx * sx;
            sumSigY2 += sy * sy;
            sumSigZ2 += sz * sz;
            count++;
          }
        }
      }
    }

    if (count === 0) {
      const [vx, vy, vz] = this.sampleMeanVelocity(cx, cy, cz);
      const [sx, sy, sz] = this.sampleVelocityDispersion(cx, cy, cz);
      return {
        meanBulkFlow: [vx, vy, vz],
        speed: Math.hypot(vx, vy, vz),
        bulkFlowSigma: [sx, sy, sz],
        totalSigma: Math.hypot(sx, sy, sz) / Math.sqrt(3.0),
        sampleCount: 1
      };
    }

    const bVx = sumVx / count;
    const bVy = sumVy / count;
    const bVz = sumVz / count;

    const sigVx = Math.sqrt(sumSigX2) / count;
    const sigVy = Math.sqrt(sumSigY2) / count;
    const sigVz = Math.sqrt(sumSigZ2) / count;

    return {
      meanBulkFlow: [bVx, bVy, bVz],
      speed: Math.hypot(bVx, bVy, bVz),
      bulkFlowSigma: [sigVx, sigVy, sigVz],
      totalSigma: Math.hypot(sigVx, sigVy, sigVz) / Math.sqrt(3.0),
      sampleCount: count
    };
  }

  /**
   * Exports summary metadata and global field statistics.
   *
   * @returns {Object}
   */
  exportSummaryMetrics() {
    const total = this.grid.totalCells;
    let sumSpeed = 0.0, maxSpeed = 0.0;
    let sumSigma = 0.0, maxSigma = 0.0;

    for (let i = 0; i < total; i++) {
      const speed = Math.hypot(this.meanVx[i], this.meanVy[i], this.meanVz[i]);
      const sigma = Math.hypot(this.rmsVx[i], this.rmsVy[i], this.rmsVz[i]) / Math.sqrt(3.0);

      sumSpeed += speed;
      if (speed > maxSpeed) maxSpeed = speed;

      sumSigma += sigma;
      if (sigma > maxSigma) maxSigma = sigma;
    }

    return {
      catalogName: this.metadata.catalogName,
      reference: this.metadata.referencePublication,
      chainLength: this.metadata.chainLength,
      gelmanRubinR: this.metadata.gelmanRubinR,
      totalVoxels: total,
      meanSpeedKms: sumSpeed / total,
      maxSpeedKms: maxSpeed,
      meanDispersionKms: sumSigma / total,
      maxDispersionKms: maxSigma,
      meanSNR: (sumSpeed / total) / Math.max(1e-6, sumSigma / total)
    };
  }
}

// ============================================================================
// 3. CHOLESKY / KARHUNEN-LOÈVE STOCHASTIC REALIZATION SAMPLER
// ============================================================================

/**
 * Generates conditional Gaussian random field realizations v^(k)(x) = \bar{v}(x) + \delta v^(k)(x)
 * conditioned on the HMC posterior field using Cholesky and Karhunen-Loève expansions.
 */
export class HMCRealizationSampler {
  /**
   * @param {BayesianHMCField} hmcField
   * @param {SpatialCovarianceKernel} [kernel]
   * @param {Object} [options={}]
   */
  constructor(hmcField, kernel = null, options = {}) {
    if (!(hmcField instanceof BayesianHMCField)) {
      throw new TypeError('HMCRealizationSampler: hmcField must be an instance of BayesianHMCField.');
    }
    this.hmcField = hmcField;
    this.kernel = kernel || new SpatialCovarianceKernel();
    this.seed = options.seed ?? 42;
    this.prng = new MT19937(this.seed);
  }

  /**
   * Generates standard normal random vector \xi ~ N(0, I_D).
   *
   * @param {number} dim
   * @param {MT19937|PCG64} [prng]
   * @returns {Float64Array}
   */
  sampleStandardNormalVector(dim, prng = null) {
    const rng = prng || this.prng;
    const xi = new Float64Array(dim);

    for (let i = 0; i < dim; i += 2) {
      const u1 = Math.max(1e-15, rng.extractNumber());
      const u2 = rng.extractNumber();
      const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

      xi[i] = z0;
      if (i + 1 < dim) {
        xi[i + 1] = z1;
      }
    }

    return xi;
  }

  /**
   * Samples a local independent 3D velocity perturbation at a single point (x, y, z).
   *
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {MT19937|PCG64} [prng]
   * @returns {[number, number, number]} Perturbed velocity [vx, vy, vz] in km/s
   */
  sampleVelocityAtPoint(x, y, z, prng = null) {
    const [mvx, mvy, mvz] = this.hmcField.sampleMeanVelocity(x, y, z);
    const [svx, svy, svz] = this.hmcField.sampleVelocityDispersion(x, y, z);
    const xi = this.sampleStandardNormalVector(3, prng);

    return [
      mvx + svx * xi[0],
      mvy + svy * xi[1],
      mvz + svz * xi[2]
    ];
  }

  /**
   * Generates a spatially correlated random field realization along an ordered path of M points
   * using Cholesky decomposition of the spatial covariance matrix C = L L^T.
   *
   * v^{(k)}(x_i) = \bar{v}(x_i) + (L \xi)_i
   *
   * @param {Array<[number, number, number]>} pathPoints - M points along path
   * @param {MT19937|PCG64} [prng]
   * @returns {Array<[number, number, number]>} M perturbed velocity vectors
   */
  sampleRealizationAlongPath(pathPoints, prng = null) {
    const M = pathPoints.length;
    if (M === 0) return [];

    const sigmaVecs = pathPoints.map(p => this.hmcField.sampleVelocityDispersion(p[0], p[1], p[2]));
    const meanVecs = pathPoints.map(p => this.hmcField.sampleMeanVelocity(p[0], p[1], p[2]));

    const C = this.kernel.assembleBlockCovarianceMatrix(pathPoints, sigmaVecs);
    const L = choleskyDecompositionGeneral(C, this.kernel.jitter);
    const xi = this.sampleStandardNormalVector(3 * M, prng);

    const deltaV = new Float64Array(3 * M);
    for (let i = 0; i < 3 * M; i++) {
      let sum = 0.0;
      for (let j = 0; j <= i; j++) {
        sum += L[i][j] * xi[j];
      }
      deltaV[i] = sum;
    }

    const realizedVelocities = [];
    for (let p = 0; p < M; p++) {
      realizedVelocities.push([
        meanVecs[p][0] + deltaV[3 * p + 0],
        meanVecs[p][1] + deltaV[3 * p + 1],
        meanVecs[p][2] + deltaV[3 * p + 2]
      ]);
    }

    return realizedVelocities;
  }

  /**
   * Generates a Karhunen-Loève truncated stochastic realization for M points.
   *
   * @param {Array<[number, number, number]>} pathPoints
   * @param {number} [energyThreshold=0.99]
   * @param {MT19937|PCG64} [prng]
   * @returns {{
   *   realizedVelocities: Array<[number, number, number]>,
   *   numModesRetained: number,
   *   totalModes: number,
   *   varianceCaptured: number
   * }}
   */
  sampleKarhunenLoeveRealization(pathPoints, energyThreshold = 0.99, prng = null) {
    const M = pathPoints.length;
    const totalDim = 3 * M;
    const sigmaVecs = pathPoints.map(p => this.hmcField.sampleVelocityDispersion(p[0], p[1], p[2]));
    const meanVecs = pathPoints.map(p => this.hmcField.sampleMeanVelocity(p[0], p[1], p[2]));

    const C = this.kernel.assembleBlockCovarianceMatrix(pathPoints, sigmaVecs);
    const { eigenvalues, eigenvectors } = jacobiDiagonalizeNxN(C);

    let totalVariance = 0.0;
    for (let i = 0; i < totalDim; i++) totalVariance += Math.max(0.0, eigenvalues[i]);

    let cumVar = 0.0;
    let K = 0;
    for (let i = 0; i < totalDim; i++) {
      cumVar += Math.max(0.0, eigenvalues[i]);
      K++;
      if (totalVariance > 0 && cumVar / totalVariance >= energyThreshold) break;
    }

    const xi = this.sampleStandardNormalVector(K, prng);

    const deltaV = new Float64Array(totalDim);
    for (let m = 0; m < K; m++) {
      const lambdaM = Math.max(0.0, eigenvalues[m]);
      const coeff = Math.sqrt(lambdaM) * xi[m];
      const eVec = eigenvectors[m];
      for (let i = 0; i < totalDim; i++) {
        deltaV[i] += coeff * eVec[i];
      }
    }

    const realizedVelocities = [];
    for (let p = 0; p < M; p++) {
      realizedVelocities.push([
        meanVecs[p][0] + deltaV[3 * p + 0],
        meanVecs[p][1] + deltaV[3 * p + 1],
        meanVecs[p][2] + deltaV[3 * p + 2]
      ]);
    }

    return {
      realizedVelocities,
      numModesRetained: K,
      totalModes: totalDim,
      varianceCaptured: totalVariance > 0 ? cumVar / totalVariance : 1.0
    };
  }

  /**
   * Generates a full 3D discrete VelocityField realization on the grid.
   *
   * @param {MT19937|PCG64} [prng]
   * @returns {VelocityField}
   */
  generateGridRealization(prng = null) {
    const total = this.hmcField.grid.totalCells;
    const rVx = new Float64Array(total);
    const rVy = new Float64Array(total);
    const rVz = new Float64Array(total);

    const rng = prng || this.prng;

    for (let i = 0; i < total; i++) {
      const [z0, z1] = [this.sampleStandardNormalVector(2, rng), this.sampleStandardNormalVector(2, rng)];
      rVx[i] = this.hmcField.meanVx[i] + this.hmcField.rmsVx[i] * z0[0];
      rVy[i] = this.hmcField.meanVy[i] + this.hmcField.rmsVy[i] * z0[1];
      rVz[i] = this.hmcField.meanVz[i] + this.hmcField.rmsVz[i] * z1[0];
    }

    return new VelocityField(this.hmcField.grid, rVx, rVy, rVz, {
      scaleFactor: 1.0,
      enforceScaleFactor: false
    });
  }
}

// ============================================================================
// 4. PROBABILISTIC WATERSHED BASIN ASSIGNMENT ENGINE
// ============================================================================

/**
 * Probabilistic Watershed Classifier computing:
 * P(voxel \in Basin_i) = (1 / N_realiz) \sum_k I(endpoint^(k) \in Basin_i)
 * as well as Shannon Entropy H(voxel) and Gini-Simpson uncertainty ribbons.
 */
export class ProbabilisticWatershedClassifier {
  /**
   * @param {BayesianHMCField} hmcField
   * @param {Object} [options={}]
   */
  constructor(hmcField, options = {}) {
    if (!(hmcField instanceof BayesianHMCField)) {
      throw new TypeError('ProbabilisticWatershedClassifier: hmcField must be an instance of BayesianHMCField.');
    }
    this.hmcField = hmcField;
    this.sampler = new HMCRealizationSampler(hmcField, options.kernel, options);
    this.numRealizations = options.numRealizations ?? 100;
    this.seed = options.seed ?? 42;
    this.prng = new MT19937(this.seed);
  }

  /**
   * Traces an individual streamline in a given velocity field using 4th-Order Runge-Kutta.
   *
   * @param {function(number, number, number): [number, number, number]} velocityFn
   * @param {[number, number, number]} x0 - Starting seed in Mpc/h
   * @param {Object} [options={}]
   * @returns {[number, number, number]} Asymptotic endpoint coordinate in Mpc/h
   */
  traceStreamlineEndpoint(velocityFn, x0, options = {}) {
    const maxSteps = options.maxSteps ?? 500;
    const dt = options.dt ?? 0.2;
    const vTol = options.velocityTol ?? 1e-3;
    const maxExtent = options.maxExtent ?? (this.hmcField.grid.boxSize[0] * 0.5);

    let [x, y, z] = x0;

    for (let step = 0; step < maxSteps; step++) {
      if (Math.abs(x) > maxExtent || Math.abs(y) > maxExtent || Math.abs(z) > maxExtent) {
        break;
      }

      const [k1x, k1y, k1z] = velocityFn(x, y, z);
      const speed = Math.hypot(k1x, k1y, k1z);
      if (speed < vTol) break;

      const [k2x, k2y, k2z] = velocityFn(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z);
      const [k3x, k3y, k3z] = velocityFn(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z);
      const [k4x, k4y, k4z] = velocityFn(x + dt * k3x, y + dt * k3y, z + dt * k3z);

      const dx = (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
      const dy = (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
      const dz = (dt / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);

      if (Math.hypot(dx, dy, dz) < 1e-5) break;

      x += dx;
      y += dy;
      z += dz;
    }

    return [x, y, z];
  }

  /**
   * Maps an endpoint coordinate [x, y, z] to the closest attractor basin ID from a list of attractors.
   *
   * @param {[number, number, number]} endpoint
   * @param {Array<{ id: number, attractorSGMpc: [number, number, number] }>} attractors
   * @returns {number} Assigned basin ID (0 = unassigned / background)
   */
  assignEndpointToBasin(endpoint, attractors) {
    if (!attractors || attractors.length === 0) return 0;

    let bestId = 0;
    let minDistSq = Infinity;

    for (const attr of attractors) {
      const [ax, ay, az] = attr.attractorSGMpc;
      const dx = endpoint[0] - ax;
      const dy = endpoint[1] - ay;
      const dz = endpoint[2] - az;
      const distSq = dx * dx + dy * dy + dz * dz;

      if (distSq < minDistSq) {
        minDistSq = distSq;
        bestId = attr.id;
      }
    }

    return bestId;
  }

  /**
   * Computes the probabilistic basin membership distribution for a single seed point across N realizations.
   *
   * @param {[number, number, number]} seedPoint - Starting point [x, y, z] in Mpc/h
   * @param {Array<{ id: number, attractorSGMpc: [number, number, number] }>} attractors
   * @param {number} [numRealizations=100]
   * @param {Object} [traceOptions={}]
   * @returns {{
   *   probabilities: Record<number, number>,
   *   mostLikelyBasinId: number,
   *   maxProbability: number,
   *   shannonEntropy: number,
   *   giniSimpsonIndex: number,
   *   endpoints: Array<[number, number, number]>
   * }}
   */
  classifyPointProbabilistically(seedPoint, attractors, numRealizations = 100, traceOptions = {}) {
    const counts = {};
    const endpoints = [];
    for (const attr of attractors) counts[attr.id] = 0;
    counts[0] = 0;

    for (let k = 0; k < numRealizations; k++) {
      const realizationVelFn = (x, y, z) => this.sampler.sampleVelocityAtPoint(x, y, z, this.prng);
      const ep = this.traceStreamlineEndpoint(realizationVelFn, seedPoint, traceOptions);
      endpoints.push(ep);

      const basinId = this.assignEndpointToBasin(ep, attractors);
      counts[basinId] = (counts[basinId] || 0) + 1;
    }

    const probabilities = {};
    let maxP = -1.0;
    let bestId = 0;
    let entropy = 0.0;
    let sumP2 = 0.0;

    for (const [idStr, count] of Object.entries(counts)) {
      const id = Number(idStr);
      const p = count / numRealizations;
      probabilities[id] = p;
      if (p > maxP) {
        maxP = p;
        bestId = id;
      }
      if (p > 1e-12) {
        entropy -= p * Math.log2(p);
      }
      sumP2 += p * p;
    }

    const giniSimpsonIndex = 1.0 - sumP2;

    return {
      probabilities,
      mostLikelyBasinId: bestId,
      maxProbability: maxP,
      shannonEntropy: entropy,
      giniSimpsonIndex,
      endpoints
    };
  }

  /**
   * Evaluates probabilistic basin volume expectations \mathbb{E}[V_i] and standard deviations across grid seeds.
   *
   * @param {Array<[number, number, number]>} sampleGridPoints
   * @param {Array<{ id: number, attractorSGMpc: [number, number, number] }>} attractors
   * @param {number} voxelVolumeMpc3
   * @param {number} [numRealizations=50]
   * @returns {Record<number, {
   *   basinId: number,
   *   expectedVolumeMpc3: number,
   *   expectedVolume1e6Mpc3: number,
   *   volumeStandardDeviationMpc3: number,
   *   ci95Mpc3: [number, number]
   * }>}
   */
  computeProbabilisticBasinVolumes(sampleGridPoints, attractors, voxelVolumeMpc3, numRealizations = 50) {
    const numPoints = sampleGridPoints.length;
    const basinIds = attractors.map(a => a.id);
    basinIds.push(0);

    const realizationVolumes = Array.from({ length: numRealizations }, () => {
      const row = {};
      for (const id of basinIds) row[id] = 0;
      return row;
    });

    for (let k = 0; k < numRealizations; k++) {
      const velFn = (x, y, z) => this.sampler.sampleVelocityAtPoint(x, y, z, this.prng);

      for (let p = 0; p < numPoints; p++) {
        const pt = sampleGridPoints[p];
        const ep = this.traceStreamlineEndpoint(velFn, pt, { maxSteps: 300 });
        const basinId = this.assignEndpointToBasin(ep, attractors);
        realizationVolumes[k][basinId] += voxelVolumeMpc3;
      }
    }

    const results = {};
    for (const id of basinIds) {
      const vols = realizationVolumes.map(r => r[id] || 0);
      const mean = vols.reduce((a, b) => a + b, 0) / numRealizations;
      let sumSq = 0.0;
      for (const v of vols) sumSq += (v - mean) ** 2;
      const std = Math.sqrt(sumSq / Math.max(1, numRealizations - 1));

      const sorted = [...vols].sort((a, b) => a - b);
      const lowIdx = Math.floor(0.025 * numRealizations);
      const highIdx = Math.min(numRealizations - 1, Math.floor(0.975 * numRealizations));

      results[id] = {
        basinId: id,
        expectedVolumeMpc3: mean,
        expectedVolume1e6Mpc3: mean / 1e6,
        volumeStandardDeviationMpc3: std,
        ci95Mpc3: [sorted[lowIdx], sorted[highIdx]]
      };
    }

    return results;
  }

  /**
   * Computes Jaccard and Dice similarity coefficients between stochastic realization segmentations and the mean baseline.
   *
   * @param {number[]} baselineLabels
   * @param {number[][]} realizationLabels
   * @param {number} targetBasinId
   * @returns {{
   *   meanJaccard: number,
   *   meanDice: number,
   *   jaccardCI95: [number, number],
   *   diceCI95: [number, number]
   * }}
   */
  static computeSegmentationOverlapMetrics(baselineLabels, realizationLabels, targetBasinId) {
    const N = realizationLabels.length;
    const numPoints = baselineLabels.length;

    const jaccardScores = [];
    const diceScores = [];

    const baseInBasin = new Uint8Array(numPoints);
    let baseCount = 0;
    for (let i = 0; i < numPoints; i++) {
      if (baselineLabels[i] === targetBasinId) {
        baseInBasin[i] = 1;
        baseCount++;
      }
    }

    for (let k = 0; k < N; k++) {
      const real = realizationLabels[k];
      let intersection = 0;
      let realCount = 0;

      for (let i = 0; i < numPoints; i++) {
        const inReal = real[i] === targetBasinId ? 1 : 0;
        if (inReal) realCount++;
        if (baseInBasin[i] && inReal) intersection++;
      }

      const union = baseCount + realCount - intersection;
      const jaccard = union > 0 ? intersection / union : 1.0;
      const dice = (baseCount + realCount) > 0 ? (2.0 * intersection) / (baseCount + realCount) : 1.0;

      jaccardScores.push(jaccard);
      diceScores.push(dice);
    }

    const meanJaccard = jaccardScores.reduce((a, b) => a + b, 0) / N;
    const meanDice = diceScores.reduce((a, b) => a + b, 0) / N;

    jaccardScores.sort((a, b) => a - b);
    diceScores.sort((a, b) => a - b);

    const lowIdx = Math.floor(0.025 * N);
    const highIdx = Math.min(N - 1, Math.floor(0.975 * N));

    return {
      meanJaccard,
      meanDice,
      jaccardCI95: [jaccardScores[lowIdx], jaccardScores[highIdx]],
      diceCI95: [diceScores[lowIdx], diceScores[highIdx]]
    };
  }
}

// ============================================================================
// 5. CRITICAL POINT STABILITY & TOPOLOGICAL PERSISTENCE ENGINE
// ============================================================================

/**
 * Evaluates the probability that a cosmological critical point (Attractor sink, Repeller void source,
 * or Saddle filament/wall hub) persists under posterior HMC fluctuations without bifurcating.
 */
export class CriticalPointStabilityAnalyzer {
  /**
   * @param {BayesianHMCField} hmcField
   * @param {Object} [options={}]
   */
  constructor(hmcField, options = {}) {
    if (!(hmcField instanceof BayesianHMCField)) {
      throw new TypeError('CriticalPointStabilityAnalyzer: hmcField must be an instance of BayesianHMCField.');
    }
    this.hmcField = hmcField;
    this.sampler = new HMCRealizationSampler(hmcField, options.kernel, options);
    this.seed = options.seed ?? 42;
    this.prng = new MT19937(this.seed);
  }

  /**
   * Locates a stagnation point v(x*) = 0 in a velocity field using 3D Newton-Raphson.
   *
   * @param {function(number, number, number): [number, number, number]} velocityFn
   * @param {[number, number, number]} initialGuess
   * @param {Object} [options={}]
   * @returns {{
   *   converged: boolean,
   *   root: [number, number, number],
   *   residualNorm: number,
   *   iterations: number,
   *   jacobian: Array<Array<number>>
   * }}
   */
  findCriticalPoint(velocityFn, initialGuess, options = {}) {
    return solveNewtonRaphson3DVector(velocityFn, initialGuess, {
      velocityTol: options.velocityTol ?? 1e-4,
      stepTol: options.stepTol ?? 1e-5,
      maxIterations: options.maxIterations ?? 50,
      dampingLambda: 1e-5
    });
  }

  /**
   * Analyzes critical point stability and persistence under posterior random perturbations.
   *
   * @param {[number, number, number]} baselineRoot - Coordinates of baseline critical point in Mpc/h
   * @param {number} [numRealizations=100]
   * @param {Object} [options={}]
   * @returns {{
   *   baselineRoot: [number, number, number],
   *   baselineClassification: string,
   *   stabilityProbability: number,
   *   recoveredRoots: Array<[number, number, number]>,
   *   meanDisplacedRoot: [number, number, number],
   *   spatialDispersionMpc: number,
   *   spatialCovarianceMatrix: Array<Array<number>>,
   *   confidenceEllipsoidSemiAxesMpc: [number, number, number],
   *   morseFrequencies: Record<string, number>,
   *   bifurcationProbability: number
   * }}
   */
  analyzePointStability(baselineRoot, numRealizations = 100, options = {}) {
    const baselineVelFn = (x, y, z) => this.hmcField.sampleMeanVelocity(x, y, z);
    const baseSolve = this.findCriticalPoint(baselineVelFn, baselineRoot, options);
    const basePos = baseSolve.converged ? baseSolve.root : baselineRoot;
    const baseJ = computeJacobian3DHighOrder(baselineVelFn, basePos[0], basePos[1], basePos[2]);
    const baseProfile = analyzeVelocityTensor(baseJ);

    const recoveredRoots = [];
    const morseCounts = {
      [MorseCriticalType.ATTRACTOR]: 0,
      [MorseCriticalType.SADDLE_FILAMENT]: 0,
      [MorseCriticalType.SADDLE_WALL]: 0,
      [MorseCriticalType.REPELLER]: 0,
      [MorseCriticalType.DEGENERATE]: 0
    };

    const searchRadius = options.maxDisplacementMpc ?? 35.0;
    const searchRadiusSq = searchRadius * searchRadius;

    for (let k = 0; k < numRealizations; k++) {
      const gridRealization = this.sampler.generateGridRealization(this.prng);
      const realVelFn = (x, y, z) => gridRealization.sampleVelocity(x, y, z);
      const res = this.findCriticalPoint(realVelFn, basePos, options);

      if (res.converged) {
        const dx = res.root[0] - basePos[0];
        const dy = res.root[1] - basePos[1];
        const dz = res.root[2] - basePos[2];
        const distSq = dx * dx + dy * dy + dz * dz;

        if (distSq <= searchRadiusSq) {
          recoveredRoots.push(res.root);
          const realJ = computeJacobian3DHighOrder(realVelFn, res.root[0], res.root[1], res.root[2]);
          const profile = analyzeVelocityTensor(realJ);
          morseCounts[profile.type] = (morseCounts[profile.type] || 0) + 1;
        }
      }
    }

    const nRec = recoveredRoots.length;
    const stabilityProbability = nRec / numRealizations;

    const meanDisplacedRoot = [0.0, 0.0, 0.0];
    if (nRec > 0) {
      for (const r of recoveredRoots) {
        meanDisplacedRoot[0] += r[0];
        meanDisplacedRoot[1] += r[1];
        meanDisplacedRoot[2] += r[2];
      }
      meanDisplacedRoot[0] /= nRec;
      meanDisplacedRoot[1] /= nRec;
      meanDisplacedRoot[2] /= nRec;
    } else {
      meanDisplacedRoot[0] = basePos[0];
      meanDisplacedRoot[1] = basePos[1];
      meanDisplacedRoot[2] = basePos[2];
    }

    const spatialCov = [
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0]
    ];

    if (nRec >= 2) {
      for (const r of recoveredRoots) {
        const d = [r[0] - meanDisplacedRoot[0], r[1] - meanDisplacedRoot[1], r[2] - meanDisplacedRoot[2]];
        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            spatialCov[i][j] += (d[i] * d[j]) / (nRec - 1);
          }
        }
      }
    }

    const { eigenvalues } = jacobiDiagonalize3x3(spatialCov);
    const semiAxes = [
      Math.sqrt(Math.max(0.0, eigenvalues[0])),
      Math.sqrt(Math.max(0.0, eigenvalues[1])),
      Math.sqrt(Math.max(0.0, eigenvalues[2]))
    ];

    const totalDispersion = Math.hypot(semiAxes[0], semiAxes[1], semiAxes[2]);

    const morseFrequencies = {};
    for (const [key, count] of Object.entries(morseCounts)) {
      morseFrequencies[key] = nRec > 0 ? count / nRec : 0.0;
    }

    const baselineTypeFreq = morseFrequencies[baseProfile.type] || 0.0;
    const bifurcationProbability = 1.0 - baselineTypeFreq;

    return {
      baselineRoot: basePos,
      baselineClassification: baseProfile.type,
      stabilityProbability,
      recoveredRoots,
      meanDisplacedRoot,
      spatialDispersionMpc: totalDispersion,
      spatialCovarianceMatrix: spatialCov,
      confidenceEllipsoidSemiAxesMpc: semiAxes,
      morseFrequencies,
      bifurcationProbability
    };
  }
}

// ============================================================================
// 6. COSMIC VARIANCE PROPAGATION ALONG STREAMLINES
// ============================================================================

/**
 * Propagates spatial trajectory covariance Sigma_x(t) and velocity dispersion along streamlines
 * using matrix ODE integration: d\Sigma / dt = J \Sigma + \Sigma J^T + Q(t).
 */
export class CosmicVarianceStreamlinePropagator {
  /**
   * @param {BayesianHMCField} hmcField
   * @param {SpatialCovarianceKernel} [kernel]
   */
  constructor(hmcField, kernel = null) {
    if (!(hmcField instanceof BayesianHMCField)) {
      throw new TypeError('CosmicVarianceStreamlinePropagator: hmcField must be an instance of BayesianHMCField.');
    }
    this.hmcField = hmcField;
    this.kernel = kernel || new SpatialCovarianceKernel();
  }

  /**
   * Integrates the continuous matrix differential equation along a trajectory.
   *
   * @param {Array<[number, number, number]>} trajectoryPoints - Discretized streamline path points
   * @param {Array<Array<number>>} [initialCovariance] - Initial 3x3 seed uncertainty Sigma_x(0)
   * @param {number} [dt=0.1] - Integration step
   * @returns {{
   *   trajectoryPoints: Array<[number, number, number]>,
   *   covariances: Array<Array<Array<number>>>,
   *   tubeRadiiMpc: number[],
   *   velocityDispersionsKms: number[],
   *   lyapunovExponents: number[]
   * }}
   */
  propagateCovarianceAlongPath(trajectoryPoints, initialCovariance = null, dt = 0.1) {
    const N = trajectoryPoints.length;
    if (N === 0) {
      return {
        trajectoryPoints: [],
        covariances: [],
        tubeRadiiMpc: [],
        velocityDispersionsKms: [],
        lyapunovExponents: []
      };
    }

    const initSigma = initialCovariance || [
      [1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0],
      [0.0, 0.0, 1.0]
    ];

    let currentSigma = [
      [initSigma[0][0], initSigma[0][1], initSigma[0][2]],
      [initSigma[1][0], initSigma[1][1], initSigma[1][2]],
      [initSigma[2][0], initSigma[2][1], initSigma[2][2]]
    ];

    const covariances = [];
    const tubeRadii = [];
    const velDispersions = [];
    const lyapunovExponents = [];

    const initTrace = Math.sqrt(Math.max(1e-12, currentSigma[0][0] + currentSigma[1][1] + currentSigma[2][2]));

    for (let step = 0; step < N; step++) {
      const [x, y, z] = trajectoryPoints[step];

      covariances.push(currentSigma.map(r => [...r]));
      const tr = Math.max(0.0, currentSigma[0][0] + currentSigma[1][1] + currentSigma[2][2]);
      const rTube = Math.sqrt(tr);
      tubeRadii.push(rTube);

      const [sx, sy, sz] = this.hmcField.sampleVelocityDispersion(x, y, z);
      velDispersions.push(Math.hypot(sx, sy, sz) / Math.sqrt(3.0));

      const tElapsed = Math.max(1e-6, (step + 1) * dt);
      const lyap = (1.0 / tElapsed) * Math.log(Math.max(1e-12, rTube / initTrace));
      lyapunovExponents.push(lyap);

      if (step === N - 1) break;

      const J = this.hmcField.sampleMeanVelocityJacobian(x, y, z);

      const Q = [
        [(sx * sx * dt) / this.kernel.L_corr, 0.0, 0.0],
        [0.0, (sy * sy * dt) / this.kernel.L_corr, 0.0],
        [0.0, 0.0, (sz * sz * dt) / this.kernel.L_corr]
      ];

      const dSigma = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
      ];

      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          let jSig = 0.0;
          let sigJt = 0.0;
          for (let k = 0; k < 3; k++) {
            jSig += J[i][k] * currentSigma[k][j];
            sigJt += currentSigma[i][k] * J[j][k];
          }
          dSigma[i][j] = jSig + sigJt + Q[i][j];
        }
      }

      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          currentSigma[i][j] += dSigma[i][j] * dt;
        }
      }

      currentSigma = projectMatrixToNearestPSD(currentSigma, 1e-8);
    }

    return {
      trajectoryPoints,
      covariances,
      tubeRadiiMpc: tubeRadii,
      velocityDispersionsKms: velDispersions,
      lyapunovExponents
    };
  }

  /**
   * Propagates a Monte Carlo streamline ensemble bundle across M stochastic field realizations.
   *
   * @param {[number, number, number]} seedPoint - Starting point [x, y, z] in Mpc/h
   * @param {number} [numRealizations=50]
   * @param {Object} [traceOptions={}]
   * @returns {{
   *   meanTrajectory: Array<[number, number, number]>,
   *   ensembleTrajectories: Array<Array<[number, number, number]>>,
   *   empiricalTubeRadiiMpc: number[],
   *   stepCount: number
   * }}
   */
  traceEnsembleStreamlineBundle(seedPoint, numRealizations = 50, traceOptions = {}) {
    const maxSteps = traceOptions.maxSteps ?? 100;
    const dt = traceOptions.dt ?? 0.2;
    const sampler = new HMCRealizationSampler(this.hmcField);

    const trajectories = [];

    for (let k = 0; k < numRealizations; k++) {
      const traj = [];
      let [x, y, z] = seedPoint;
      traj.push([x, y, z]);

      for (let s = 0; s < maxSteps; s++) {
        const [vx, vy, vz] = sampler.sampleVelocityAtPoint(x, y, z);
        const speed = Math.hypot(vx, vy, vz);
        if (speed < 1e-3) break;

        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
        traj.push([x, y, z]);
      }

      trajectories.push(traj);
    }

    const minLength = Math.min(...trajectories.map(t => t.length));
    const meanTrajectory = [];
    const empiricalTubeRadii = [];

    for (let s = 0; s < minLength; s++) {
      let mx = 0.0, my = 0.0, mz = 0.0;
      for (let k = 0; k < numRealizations; k++) {
        mx += trajectories[k][s][0];
        my += trajectories[k][s][1];
        mz += trajectories[k][s][2];
      }
      mx /= numRealizations;
      my /= numRealizations;
      mz /= numRealizations;
      meanTrajectory.push([mx, my, mz]);

      let varSum = 0.0;
      for (let k = 0; k < numRealizations; k++) {
        const dx = trajectories[k][s][0] - mx;
        const dy = trajectories[k][s][1] - my;
        const dz = trajectories[k][s][2] - mz;
        varSum += dx * dx + dy * dy + dz * dz;
      }
      const tubeRadius = Math.sqrt(varSum / Math.max(1, numRealizations - 1));
      empiricalTubeRadii.push(tubeRadius);
    }

    return {
      meanTrajectory,
      ensembleTrajectories: trajectories,
      empiricalTubeRadiiMpc: empiricalTubeRadii,
      stepCount: minLength
    };
  }
}

// ============================================================================
// 7. STATISTICAL MCMC DIAGNOSTICS FOR 10,000-STEP HMC CHAINS
// ============================================================================

/**
 * Computes the Gelman-Rubin \hat{R} potential scale reduction factor across multiple MCMC chains.
 *
 * @param {Array<Array<number>>} chains - Array of M chains, each containing N samples
 * @returns {number} \hat{R} diagnostic (values < 1.05 indicate good convergence)
 */
export function computeGelmanRubinDiagnostic(chains) {
  const M = chains.length;
  if (M < 2) return 1.0;
  const N = chains[0].length;
  if (N < 2) return 1.0;

  const chainMeans = new Float64Array(M);
  const chainVars = new Float64Array(M);

  for (let m = 0; m < M; m++) {
    const c = chains[m];
    let sum = 0.0;
    for (let i = 0; i < N; i++) sum += c[i];
    const mean = sum / N;
    chainMeans[m] = mean;

    let sumSq = 0.0;
    for (let i = 0; i < N; i++) sumSq += (c[i] - mean) ** 2;
    chainVars[m] = sumSq / (N - 1);
  }

  // Grand mean
  let grandMean = 0.0;
  for (let m = 0; m < M; m++) grandMean += chainMeans[m];
  grandMean /= M;

  // Between-chain variance B / N
  let bOverN = 0.0;
  for (let m = 0; m < M; m++) bOverN += (chainMeans[m] - grandMean) ** 2;
  bOverN /= (M - 1);
  const B = bOverN * N;

  // Within-chain variance W
  let W = 0.0;
  for (let m = 0; m < M; m++) W += chainVars[m];
  W /= M;

  if (W <= 1e-15) return 1.0;

  // Marginal posterior variance estimate \hat{Var}(theta)
  const varHat = ((N - 1) / N) * W + (1.0 / N) * B;
  const rHat = Math.sqrt(varHat / W);

  return rHat;
}

/**
 * Computes integrated autocorrelation time \tau_{int} of an MCMC chain using Sokal's adaptive windowing.
 *
 * @param {number[]|Float64Array} chain - Series of scalar samples
 * @param {number} [cWindow=5.0] - Window cutoff multiplier
 * @returns {number} Integrated autocorrelation time \tau_{int}
 */
export function computeIntegratedAutocorrelationTime(chain, cWindow = 5.0) {
  const N = chain.length;
  if (N < 4) return 1.0;

  let sum = 0.0;
  for (let i = 0; i < N; i++) sum += chain[i];
  const mean = sum / N;

  let var0 = 0.0;
  for (let i = 0; i < N; i++) var0 += (chain[i] - mean) ** 2;
  var0 /= N;

  if (var0 <= 1e-15) return 1.0;

  let tau = 0.5;
  const maxLag = Math.floor(N / 2);

  for (let lag = 1; lag < maxLag; lag++) {
    let autoCov = 0.0;
    for (let i = 0; i < N - lag; i++) {
      autoCov += (chain[i] - mean) * (chain[i + lag] - mean);
    }
    autoCov /= N;

    const rho = autoCov / var0;
    if (rho < 0.0) break; // First negative cutoff

    tau += rho;
    if (lag >= cWindow * tau) break;
  }

  return Math.max(0.5, tau);
}

/**
 * Computes Effective Sample Size N_{eff} = N / (2 \tau_{int}).
 *
 * @param {number[]|Float64Array} chain
 * @returns {number} Effective sample size
 */
export function computeEffectiveSampleSize(chain) {
  const N = chain.length;
  const tau = computeIntegratedAutocorrelationTime(chain);
  return Math.max(1.0, Math.min(N, N / (2.0 * tau)));
}

/**
 * Computes the Energy Bayesian Fraction of Missing Information (E-BFMI) for an HMC chain.
 *
 * E-BFMI = \frac{\sum_{i=1}^{N-1} (E_i - E_{i-1})^2}{\sum_{i=1}^N (E_i - \bar{E})^2}
 * Values > 0.3 typically indicate good momentum energy exploration.
 *
 * @param {number[]|Float64Array} energies - Total Hamiltonian energies along chain
 * @returns {number} E-BFMI statistic
 */
export function computeEnergyBayesianFractionMissingInformation(energies) {
  const N = energies.length;
  if (N < 2) return 1.0;

  let numer = 0.0;
  for (let i = 1; i < N; i++) {
    const dE = energies[i] - energies[i - 1];
    numer += dE * dE;
  }
  numer /= (N - 1);

  let sumE = 0.0;
  for (let i = 0; i < N; i++) sumE += energies[i];
  const meanE = sumE / N;

  let denom = 0.0;
  for (let i = 0; i < N; i++) {
    const diff = energies[i] - meanE;
    denom += diff * diff;
  }
  denom /= (N - 1);

  if (denom <= 1e-15) return 1.0;
  return numer / denom;
}
