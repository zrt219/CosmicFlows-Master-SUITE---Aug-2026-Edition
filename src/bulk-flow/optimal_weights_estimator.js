/**
 * @file optimal_weights_estimator.js
 * @description Kaiser-Watkins-Feldman (KWF) Minimum-Variance Optimal Weighted Bulk Flow Estimator.
 *
 * Implements:
 * 1. Watkins-Feldman-Hudson (2009) minimum-variance optimal weights w_{i, n} for galaxy survey geometry.
 * 2. Survey velocity covariance matrix C_{mn} = \hat{r}_{m, i} \hat{r}_{n, j} \xi_{ij}(r_m - r_n) + \delta_{mn} \sigma_n^2.
 * 3. Cross-covariance with ideal top-hat sphere window Q_{i, n} = \langle u_n V_{ideal, i} \rangle.
 * 4. Survey window matrix W^2(k) and cosmic variance cancellation for reduced cosmic variance.
 *
 * @module bulk-flow/optimal_weights_estimator
 */

import { sphericalTopHatWindow } from './cosmic_variance_calculator.js';
import { choleskyDecomposition } from '../uncertainty/covariance_estimator.js';
import { cartesianToSupergalacticAngles } from './bulk_flow_estimator.js';

/**
 * Computes parallel and perpendicular velocity correlation functions xi_par(r) and xi_perp(r).
 * @param {number} r - Distance in Mpc/h
 * @param {function(number): number} powerSpectrumFn - P(k)
 * @param {number} h0f - H0 * f
 * @returns {{ xiPar: number, xiPerp: number }}
 */
export function computeVelocityCorrelationFunctions(r, powerSpectrumFn, h0f) {
  if (r < 1e-4) {
    // Zero separation
    let int0 = 0.0;
    const numK = 500;
    const dLogK = Math.log(10.0 / 1e-4) / (numK - 1);
    for (let i = 0; i < numK; i++) {
      const k = 1e-4 * Math.exp(i * dLogK);
      int0 += powerSpectrumFn(k) * k * dLogK;
    }
    const var0 = (h0f * h0f) / (6.0 * Math.PI * Math.PI) * int0;
    return { xiPar: var0, xiPerp: var0 };
  }

  const numK = 600;
  const kMin = 1e-4;
  const kMax = 10.0;
  const dLogK = (Math.log(kMax) - Math.log(kMin)) / (numK - 1);

  let intPar = 0.0;
  let intPerp = 0.0;

  for (let i = 0; i < numK; i++) {
    const k = kMin * Math.exp(i * dLogK);
    const kr = k * r;
    const Pk = powerSpectrumFn(k);

    // Spherical Bessel functions j0(kr), j1(kr), j2(kr)
    const sinKr = Math.sin(kr);
    const cosKr = Math.cos(kr);
    const j0 = sinKr / kr;
    const j1 = (sinKr - kr * cosKr) / (kr * kr);

    const fPar = j0 - 2.0 * j1 / kr;
    const fPerp = j1 / kr;

    intPar += Pk * fPar * k * dLogK;
    intPerp += Pk * fPerp * k * dLogK;
  }

  const factor = (h0f * h0f) / (2.0 * Math.PI * Math.PI);
  return {
    xiPar: factor * intPar,
    xiPerp: factor * intPerp
  };
}

/**
 * Minimum Variance Optimal Weights Bulk Flow Estimator.
 */
export class OptimalWeightsBulkFlowEstimator {
  /**
   * @param {Array<{ x: number, y: number, z: number, u: number, error: number }>} galaxies
   * @param {Object} [options={}]
   * @param {number} [options.cosmicDispersion=187.0]
   * @param {number} [options.targetRadiusMpc=50.0] - Ideal top-hat sphere radius R
   */
  constructor(galaxies, options = {}) {
    this.galaxies = galaxies;
    this.n = galaxies.length;
    this.cosmicDispersion = options.cosmicDispersion ?? 187.0;
    this.targetR = options.targetRadiusMpc ?? 50.0;
  }

  /**
   * Estimates bulk flow vector using minimum-variance weights.
   * Solves G * w_i = Q_i for each Cartesian component i in {0, 1, 2}.
   *
   * @param {function(number): number} powerSpectrumFn
   * @param {number} h0f
   * @returns {{
   *   bulkFlowVector: [number, number, number],
   *   magnitude: number,
   *   sglDeg: number,
   *   sgbDeg: number,
   *   covarianceMatrix: number[][],
   *   uncertainty1Sigma: [number, number, number],
   *   weights: Array<[number, number, number]>
   * }}
   */
  estimate(powerSpectrumFn, h0f) {
    const N = this.n;
    if (N < 3) throw new Error('Requires at least 3 galaxies for optimal weights estimation.');

    // Precompute positions and unit vectors
    const rMags = new Float64Array(N);
    const nHats = Array.from({ length: N }, () => new Float64Array(3));

    for (let i = 0; i < N; i++) {
      const g = this.galaxies[i];
      const r = Math.hypot(g.x, g.y, g.z);
      rMags[i] = Math.max(1e-4, r);
      nHats[i][0] = g.x / rMags[i];
      nHats[i][1] = g.y / rMags[i];
      nHats[i][2] = g.z / rMags[i];
    }

    // 1. Build Galaxy Covariance Matrix G_{mn} = <u_m u_n> + delta_{mn} (sigma_m^2 + sigma_*^2)
    const G = Array.from({ length: N }, () => new Float64Array(N));
    const sigStarSq = this.cosmicDispersion * this.cosmicDispersion;

    for (let m = 0; m < N; m++) {
      const errM = this.galaxies[m].error ?? 150.0;
      G[m][m] = errM * errM + sigStarSq;

      for (let n = m + 1; n < N; n++) {
        // Separation vector r_mn = r_m - r_n
        const dx = this.galaxies[m].x - this.galaxies[n].x;
        const dy = this.galaxies[m].y - this.galaxies[n].y;
        const dz = this.galaxies[m].z - this.galaxies[n].z;
        const sep = Math.hypot(dx, dy, dz);

        const { xiPar, xiPerp } = computeVelocityCorrelationFunctions(sep, powerSpectrumFn, h0f);

        // Dot products
        const rHatDot = (nHats[m][0] * nHats[n][0] + nHats[m][1] * nHats[n][1] + nHats[m][2] * nHats[n][2]);
        let rMNDotM = 0.0, rMNDotN = 0.0;
        if (sep > 1e-4) {
          rMNDotM = (dx * nHats[m][0] + dy * nHats[m][1] + dz * nHats[m][2]) / sep;
          rMNDotN = (dx * nHats[n][0] + dy * nHats[n][1] + dz * nHats[n][2]) / sep;
        }

        const covVal = xiPerp * rHatDot + (xiPar - xiPerp) * rMNDotM * rMNDotN;
        G[m][n] = covVal;
        G[n][m] = covVal;
      }
    }

    // 2. Build Cross-Covariance vector Q_{i, n} = <u_n V_{ideal, i}>
    const Q = [
      new Float64Array(N), // Q_x
      new Float64Array(N), // Q_y
      new Float64Array(N)  // Q_z
    ];

    const numK = 300;
    const kMin = 1e-4;
    const kMax = 10.0;
    const dLogK = (Math.log(kMax) - Math.log(kMin)) / (numK - 1);

    for (let n = 0; n < N; n++) {
      const rn = rMags[n];
      let integral = 0.0;

      for (let i = 0; i < numK; i++) {
        const k = kMin * Math.exp(i * dLogK);
        const Pk = powerSpectrumFn(k);
        const krn = k * rn;
        const j1_over_krn = krn < 1e-4 ? 1.0 / 3.0 : (Math.sin(krn) - krn * Math.cos(krn)) / (krn * krn * krn);
        const W_tophat = sphericalTopHatWindow(k * this.targetR);

        integral += Pk * j1_over_krn * W_tophat * k * dLogK;
      }

      const factorQ = ((h0f * h0f) / (2.0 * Math.PI * Math.PI)) * integral;

      Q[0][n] = nHats[n][0] * factorQ;
      Q[1][n] = nHats[n][1] * factorQ;
      Q[2][n] = nHats[n][2] * factorQ;
    }

    // 3. Solve G * w_i = Q_i via Cholesky decomposition
    const L = choleskyDecomposition(G.map(r => Array.from(r)), 1e-4);

    const weights = [
      this._solveCholesky(L, Q[0]),
      this._solveCholesky(L, Q[1]),
      this._solveCholesky(L, Q[2])
    ];

    // 4. Compute bulk flow estimate V_i = sum_n w_{i, n} u_n
    const V = [0.0, 0.0, 0.0];
    for (let comp = 0; comp < 3; comp++) {
      let sum = 0.0;
      for (let n = 0; n < N; n++) {
        sum += weights[comp][n] * this.galaxies[n].u;
      }
      V[comp] = sum;
    }

    const mag = Math.hypot(V[0], V[1], V[2]);
    const angles = cartesianToSupergalacticAngles(V[0], V[1], V[2]);

    // 5. Error covariance matrix C_{ij} = sum_{m, n} w_{i, m} G_{mn} w_{j, n}
    const cov = Array.from({ length: 3 }, () => [0, 0, 0]);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let sum = 0.0;
        for (let m = 0; m < N; m++) {
          for (let n = 0; n < N; n++) {
            sum += weights[i][m] * G[m][n] * weights[j][n];
          }
        }
        cov[i][j] = sum;
      }
    }

    const sigmaVx = Math.sqrt(Math.max(0, cov[0][0]));
    const sigmaVy = Math.sqrt(Math.max(0, cov[1][1]));
    const sigmaVz = Math.sqrt(Math.max(0, cov[2][2]));

    // Format weights per galaxy
    const galaxyWeights = [];
    for (let n = 0; n < N; n++) {
      galaxyWeights.push([weights[0][n], weights[1][n], weights[2][n]]);
    }

    return {
      bulkFlowVector: V,
      magnitude: mag,
      sglDeg: angles.sgl,
      sgbDeg: angles.sgb,
      covarianceMatrix: cov,
      uncertainty1Sigma: [sigmaVx, sigmaVy, sigmaVz],
      weights: galaxyWeights
    };
  }

  /**
   * Solves L L^T x = b via forward and backward substitution.
   * @private
   */
  _solveCholesky(L, b) {
    const N = b.length;
    const y = new Float64Array(N);

    // Forward solve L y = b
    for (let i = 0; i < N; i++) {
      let sum = 0.0;
      for (let k = 0; k < i; k++) sum += L[i][k] * y[k];
      y[i] = (b[i] - sum) / L[i][i];
    }

    // Backward solve L^T x = y
    const x = new Float64Array(N);
    for (let i = N - 1; i >= 0; i--) {
      let sum = 0.0;
      for (let k = i + 1; k < N; k++) sum += L[k][i] * x[k];
      x[i] = (y[i] - sum) / L[i][i];
    }

    return x;
  }
}
