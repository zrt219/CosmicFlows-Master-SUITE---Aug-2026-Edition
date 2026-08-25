/**
 * @file mcmc_convergence_diagnostics.js
 * @description Comprehensive Scientific MCMC Convergence Diagnostics and Gelman-Rubin Verification Suite for Cosmological Chains (CF4++).
 *
 * Implements:
 * 1. Univariate, Split, Folded, and Rank-Normalized Gelman-Rubin Potential Scale Reduction Factor (PSRF / R-hat) [Vehtari et al. 2021].
 * 2. Multivariate Potential Scale Reduction Factor (MPSRF) [Brooks & Gelman 1998] with generalized symmetric eigensolver.
 * 3. Effective Sample Size (ESS) Suite: Bulk-ESS, Tail-ESS, Quantile-ESS, Folded-ESS, and Cross-Chain Variogram Autocorrelation [Geyer 1992].
 * 4. Integrated Autocorrelation Time (tau_int) Estimators: Sokal (1997) Adaptive Window, Goodman & Weare (2010), Geyer IPSE / IMSE / ICSE, and Spectral Density Estimators (Tukey-Hanning, Parzen, Bartlett, Quadratic Spectral AR(1) plug-in, Yule-Walker AR(p)).
 * 5. Geweke Z-score stationarity diagnostic with spectral density at zero frequency and sliding window trajectory.
 * 6. Heidelberger & Welch stationarity (Cramer-von Mises) and relative half-width precision test.
 * 7. Raftery & Lewis quantile diagnostic for required chain length, burn-in, and thinning factor.
 * 8. Energy Bayesian Fraction of Missing Information (E-BFMI / Betancourt 2016) and Kolmogorov CUSUM excursion test.
 * 9. Parameter cross-correlation matrix and multi-lag cross-covariance for multi-dimensional chains.
 * 10. Multi-chain synthetic cosmological sampler (VAR(1) / AR(1) Gaussian and Student-t with prescribed covariance, burn-in drift, and multimodality).
 * 11. End-to-End CF4++ 10,000-Step Chain Convergence Verification Dossier Engine & Quality Gate Certifier.
 *
 * @module uncertainty/mcmc_convergence_diagnostics
 */

// ============================================================================
// SECTION 1: SCIENTIFIC CONSTANTS, COSMOLOGY INVARIANTS & STATISTICAL BASICS
// ============================================================================

/**
 * Standard CF4++ cosmological and MCMC convergence threshold parameters.
 * Reference: CosmicFlows-4 (Dupuy & Courtois 2023, Tully et al. 2023).
 */
export const CF4_MCMC_CONSTANTS = Object.freeze({
  H0: 74.6, // km/s/Mpc (Hubble constant)
  OMEGA_M: 0.31, // Matter density parameter
  OMEGA_LAMBDA: 0.69, // Dark energy density parameter
  SIGMA_8: 0.81, // Matter fluctuation amplitude
  GROWTH_RATE_F: 0.52, // Linear growth rate f ~ Omega_m^0.55
  VELOCITY_SCALE: 52.0, // Characteristic velocity dispersion (km/s)
  DEFAULT_CORRELATION_LENGTH_MPC: 25.0, // Spatial correlation length (Mpc/h)
  DEFAULT_CHAIN_LENGTH: 10000, // Standard HMC realization length
  MIN_ACCEPTABLE_CHAINS: 2, // Minimum chains for Gelman-Rubin
  RECOMMENDED_CHAINS: 4, // Recommended chains for production diagnostics
  R_HAT_EXCELLENT: 1.01, // Excellent convergence threshold
  R_HAT_CONVERGED: 1.05, // Standard convergence threshold
  R_HAT_MARGINAL: 1.10, // Marginal threshold
  BULK_ESS_THRESHOLD: 400.0, // Minimum Bulk-ESS for reliable mean/variance
  TAIL_ESS_THRESHOLD: 200.0, // Minimum Tail-ESS for 95% credible intervals
  GEWEKE_ALPHA: 0.05, // Geweke stationarity significance level
  GEWEKE_Z_CRITICAL: 1.959963984540054, // 95% two-tailed standard normal critical value
  EBFMI_THRESHOLD: 0.30, // Betancourt (2016) HMC energy transition threshold
  MCSE_MAX_REL_SD: 0.05, // Maximum tolerable relative MCSE (MCSE / SD <= 5%)
  SOKAL_C_DEFAULT: 5.0, // Sokal adaptive window cutoff constant
  HW_EPSILON_DEFAULT: 0.02, // Heidelberger-Welch relative half-width threshold (2%)
  RAFTERY_Q_DEFAULT: 0.025, // Raftery-Lewis target lower quantile (2.5%)
  RAFTERY_R_DEFAULT: 0.005, // Raftery-Lewis precision half-width (+-0.5%)
  RAFTERY_S_DEFAULT: 0.95 // Raftery-Lewis probability confidence (95%)
});

/**
 * Standard cosmological parameter metadata for CF4++ Bayesian posteriors.
 */
export const CF4_PARAMETER_DEFINITIONS = Object.freeze({
  H0: Object.freeze({ name: 'H0', unit: 'km/s/Mpc', nominal: 74.6, min: 50.0, max: 100.0, description: 'Hubble expansion rate' }),
  OMEGA_M: Object.freeze({ name: 'Omega_M', unit: 'dimensionless', nominal: 0.31, min: 0.10, max: 0.60, description: 'Total matter density' }),
  OMEGA_LAMBDA: Object.freeze({ name: 'Omega_Lambda', unit: 'dimensionless', nominal: 0.69, min: 0.40, max: 0.90, description: 'Dark energy density' }),
  SIGMA_8: Object.freeze({ name: 'sigma_8', unit: 'dimensionless', nominal: 0.81, min: 0.50, max: 1.20, description: 'Matter power spectrum normalization' }),
  V_BULK: Object.freeze({ name: 'V_bulk', unit: 'km/s', nominal: 388.0, min: 100.0, max: 800.0, description: 'Dipole bulk flow velocity at 150 Mpc/h' }),
  L_APEX: Object.freeze({ name: 'l_apex', unit: 'deg', nominal: 310.0, min: 0.0, max: 360.0, description: 'Galactic longitude of bulk flow apex' }),
  B_APEX: Object.freeze({ name: 'b_apex', unit: 'deg', nominal: 22.0, min: -90.0, max: 90.0, description: 'Galactic latitude of bulk flow apex' }),
  DELTA_H: Object.freeze({ name: 'delta_H', unit: 'dimensionless', nominal: 0.0, min: -0.15, max: 0.15, description: 'Hubble bubble monopole perturbation' }),
  ALPHA_FIELD: Object.freeze({ name: 'alpha_field', unit: 'dimensionless', nominal: 1.0, min: 0.5, max: 1.5, description: 'Velocity field amplitude scaling' })
});

/**
 * Error function approximation using Abramowitz and Stegun (formula 7.1.26).
 * Maximum error < 1.5e-7 across entire real line.
 * @param {number} x
 * @returns {number}
 */
export function erf(x) {
  const sign = x < 0 ? -1 : 1;
  const absX = Math.abs(x);
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;

  const t = 1.0 / (1.0 + p * absX);
  const poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
  const y = 1.0 - poly * Math.exp(-absX * absX);
  return sign * y;
}

/**
 * Complementary error function erfc(x) = 1 - erf(x).
 * @param {number} x
 * @returns {number}
 */
export function erfc(x) {
  return 1.0 - erf(x);
}

/**
 * Standard Normal (Gaussian) Cumulative Distribution Function \Phi(x).
 * @param {number} x
 * @returns {number}
 */
export function normalCDF(x) {
  return 0.5 * (1.0 + erf(x / Math.SQRT2));
}

/**
 * Standard Normal Probability Density Function \phi(x).
 * @param {number} x
 * @returns {number}
 */
export function normalPDF(x) {
  return (1.0 / Math.sqrt(2.0 * Math.PI)) * Math.exp(-0.5 * x * x);
}

/**
 * High-precision Inverse Normal CDF (Probit function) using Acklam's rational approximation.
 * Relative error < 1.15e-9 over the full range p in (0, 1).
 * @param {number} p - Probability strictly in (0, 1)
 * @returns {number} - Standard normal quantile z such that \Phi(z) = p
 */
export function probit(p) {
  if (p <= 0.0 || p >= 1.0 || isNaN(p)) {
    if (p === 0.0) return -Infinity;
    if (p === 1.0) return Infinity;
    throw new RangeError(`Probit input p must be strictly in (0, 1), received ${p}`);
  }

  const a = [
    -3.969683028665376e1,
    2.209460984245205e2,
    -2.759285104469687e2,
    1.383577518672690e2,
    -3.066479806614716e1,
    2.506628277459239e0
  ];
  const b = [
    -5.447609879822406e1,
    1.615858368580409e2,
    -1.556989798598866e2,
    6.680131188771972e1,
    -1.328068155288572e1
  ];
  const c = [
    -7.784894002430293e-3,
    -3.223964580411365e-1,
    -2.400758277161838e0,
    -2.549732539343734e0,
    4.374664141464968e0,
    2.938163982698783e0
  ];
  const d = [
    7.784695709041462e-3,
    3.224671290700398e-1,
    2.445134137142996e0,
    3.754408661907416e0
  ];

  const pLow = 0.02425;
  const pHigh = 1.0 - pLow;

  if (p < pLow) {
    const q = Math.sqrt(-2.0 * Math.log(p));
    return (
      (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    );
  } else if (p <= pHigh) {
    const q = p - 0.5;
    const r = q * q;
    return (
      (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    );
  } else {
    const q = Math.sqrt(-2.0 * Math.log(1.0 - p));
    return -(
      (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    );
  }
}

/**
 * Two-tailed p-value from standard normal z-score.
 * @param {number} z
 * @returns {number}
 */
export function normalTwoTailedPValue(z) {
  const absZ = Math.abs(z);
  return 2.0 * (1.0 - normalCDF(absZ));
}

/**
 * Log-Gamma function ln \Gamma(x) via Lanczos approximation.
 * @param {number} x
 * @returns {number}
 */
export function logGamma(x) {
  if (x <= 0.0) return NaN;
  const c = [
    57.1562356658629235,
    -59.5979603554754912,
    14.1360979747417471,
    -0.491913816097620199,
    0.339946499848118887e-4,
    0.465236289270485756e-4,
    -0.983744753048795646e-4,
    0.158088703224912494e-3,
    -0.210264441724104883e-3,
    0.217439618115212643e-3,
    -0.164318106536763890e-3,
    0.844182239838527433e-4,
    -0.261908384015814087e-4,
    0.368991826595316234e-5
  ];
  let sum = 0.99999999999999709182;
  for (let i = 0; i < c.length; i++) {
    sum += c[i] / (x + i + 1.0);
  }
  const t = x + 4.7421875;
  return (x + 0.5) * Math.log(t) - t + Math.log(Math.SQRT2 * Math.sqrt(Math.PI) * sum / x);
}

/**
 * Student-t PDF with nu degrees of freedom.
 * @param {number} t
 * @param {number} nu
 * @returns {number}
 */
export function studentTPDF(t, nu) {
  const num = Math.exp(logGamma((nu + 1) / 2) - logGamma(nu / 2));
  const den = Math.sqrt(nu * Math.PI) * Math.pow(1.0 + (t * t) / nu, (nu + 1) / 2);
  return num / den;
}

// ============================================================================
// SECTION 2: VECTOR / MATRIX ALGEBRA & EIGENVALUES FOR MPSRF
// ============================================================================

/**
 * Computes sample mean of a 1D array.
 * @param {ArrayLike<number>} x
 * @returns {number}
 */
export function computeMean(x) {
  const n = x.length;
  if (n === 0) return 0.0;
  let sum = 0.0;
  for (let i = 0; i < n; i++) sum += x[i];
  return sum / n;
}

/**
 * Computes sample variance with Bessel correction (N - 1).
 * @param {ArrayLike<number>} x
 * @param {number} [mean] - Optional precomputed mean
 * @returns {number}
 */
export function computeVariance(x, mean = undefined) {
  const n = x.length;
  if (n <= 1) return 0.0;
  const m = mean !== undefined ? mean : computeMean(x);
  let sumSq = 0.0;
  for (let i = 0; i < n; i++) {
    const d = x[i] - m;
    sumSq += d * d;
  }
  return sumSq / (n - 1);
}

/**
 * Computes sample standard deviation.
 * @param {ArrayLike<number>} x
 * @param {number} [mean]
 * @returns {number}
 */
export function computeStd(x, mean = undefined) {
  return Math.sqrt(Math.max(0.0, computeVariance(x, mean)));
}

/**
 * Computes sample skewness (Fisher-Pearson standardized 3rd moment).
 * @param {ArrayLike<number>} x
 * @returns {number}
 */
export function computeSkewness(x) {
  const n = x.length;
  if (n < 3) return 0.0;
  const m = computeMean(x);
  const s = computeStd(x, m);
  if (s < 1e-15) return 0.0;
  let sumCube = 0.0;
  for (let i = 0; i < n; i++) {
    const d = (x[i] - m) / s;
    sumCube += d * d * d;
  }
  return (n / ((n - 1) * (n - 2))) * sumCube;
}

/**
 * Computes sample excess kurtosis (4th standardized moment - 3).
 * @param {ArrayLike<number>} x
 * @returns {number}
 */
export function computeKurtosis(x) {
  const n = x.length;
  if (n < 4) return 0.0;
  const m = computeMean(x);
  const s = computeStd(x, m);
  if (s < 1e-15) return 0.0;
  let sum4 = 0.0;
  for (let i = 0; i < n; i++) {
    const d = (x[i] - m) / s;
    sum4 += d * d * d * d;
  }
  const factor1 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3));
  const factor2 = (3.0 * (n - 1) * (n - 1)) / ((n - 2) * (n - 3));
  return factor1 * sum4 - factor2;
}

/**
 * Computes sample quantiles using Hyndman & Fan Type 7 linear interpolation.
 * @param {ArrayLike<number>} x
 * @param {number} q - Quantile in [0, 1]
 * @returns {number}
 */
export function computeQuantile(x, q) {
  const n = x.length;
  if (n === 0) return NaN;
  if (n === 1) return x[0];
  if (q <= 0.0) {
    let minVal = x[0];
    for (let i = 1; i < n; i++) if (x[i] < minVal) minVal = x[i];
    return minVal;
  }
  if (q >= 1.0) {
    let maxVal = x[0];
    for (let i = 1; i < n; i++) if (x[i] > maxVal) maxVal = x[i];
    return maxVal;
  }

  const sorted = Array.from(x).sort((a, b) => a - b);
  const index = (n - 1) * q;
  const lo = Math.floor(index);
  const hi = Math.ceil(index);
  const frac = index - lo;
  return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
}

/**
 * Computes Median Absolute Deviation (MAD) scaled to match normal standard deviation (factor 1.4826).
 * @param {ArrayLike<number>} x
 * @returns {number}
 */
export function computeMAD(x) {
  const n = x.length;
  if (n === 0) return 0.0;
  const median = computeQuantile(x, 0.5);
  const devs = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    devs[i] = Math.abs(x[i] - median);
  }
  const madRaw = computeQuantile(devs, 0.5);
  return madRaw * 1.482602218505602;
}

/**
 * Computes average ranks (1-based) with exact tie handling.
 * @param {ArrayLike<number>} arr
 * @returns {Float64Array}
 */
export function computeRanks(arr) {
  const n = arr.length;
  const indices = new Int32Array(n);
  for (let i = 0; i < n; i++) indices[i] = i;

  indices.sort((a, b) => arr[a] - arr[b]);

  const ranks = new Float64Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    while (j + 1 < n && arr[indices[j + 1]] === arr[indices[i]]) {
      j++;
    }
    const avgRank = (i + 1 + (j + 1)) / 2.0;
    for (let k = i; k <= j; k++) {
      ranks[indices[k]] = avgRank;
    }
    i = j + 1;
  }
  return ranks;
}

/**
 * Creates a D x D matrix flattened into a 1D Float64Array.
 * @param {number} d
 * @returns {Float64Array}
 */
export function createMatrix(d) {
  return new Float64Array(d * d);
}

/**
 * Creates a D x D identity matrix.
 * @param {number} d
 * @returns {Float64Array}
 */
export function eye(d) {
  const mat = new Float64Array(d * d);
  for (let i = 0; i < d; i++) mat[i * d + i] = 1.0;
  return mat;
}

/**
 * Matrix multiplication C = A * B for D x D matrices.
 * @param {Float64Array} A
 * @param {Float64Array} B
 * @param {number} d
 * @returns {Float64Array}
 */
export function matMul(A, B, d) {
  const C = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    const rowOffset = i * d;
    for (let k = 0; k < d; k++) {
      const aVal = A[rowOffset + k];
      if (Math.abs(aVal) > 1e-20) {
        const bOffset = k * d;
        for (let j = 0; j < d; j++) {
          C[rowOffset + j] += aVal * B[bOffset + j];
        }
      }
    }
  }
  return C;
}

/**
 * Computes matrix transpose A^T (D x D).
 * @param {Float64Array} A
 * @param {number} d
 * @returns {Float64Array}
 */
export function matTranspose(A, d) {
  const AT = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      AT[j * d + i] = A[i * d + j];
    }
  }
  return AT;
}

/**
 * Computes matrix trace Tr(A) = \sum A_{ii}.
 * @param {Float64Array} A
 * @param {number} d
 * @returns {number}
 */
export function matTrace(A, d) {
  let tr = 0.0;
  for (let i = 0; i < d; i++) tr += A[i * d + i];
  return tr;
}

/**
 * Computes Frobenius norm of matrix A.
 * @param {Float64Array} A
 * @param {number} d
 * @returns {number}
 */
export function matFrobeniusNorm(A, d) {
  let sumSq = 0.0;
  const n = d * d;
  for (let i = 0; i < n; i++) sumSq += A[i] * A[i];
  return Math.sqrt(sumSq);
}

/**
 * Performs Cholesky decomposition of symmetric positive-definite matrix A (D x D).
 * Returns lower triangular matrix L such that A = L * L^T.
 * @param {Float64Array} A
 * @param {number} d
 * @param {number} [jitter=1e-10]
 * @returns {Float64Array}
 */
export function choleskyDecompose(A, d, jitter = 1e-10) {
  const L = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0.0;
      for (let k = 0; k < j; k++) {
        sum += L[i * d + k] * L[j * d + k];
      }
      if (i === j) {
        let val = A[i * d + i] - sum;
        if (val <= 0.0) {
          val = jitter;
        }
        L[i * d + j] = Math.sqrt(val);
      } else {
        const denom = L[j * d + j];
        L[i * d + j] = denom > 1e-15 ? (A[i * d + j] - sum) / denom : 0.0;
      }
    }
  }
  return L;
}

/**
 * Inverts a lower triangular matrix L (D x D).
 * @param {Float64Array} L
 * @param {number} d
 * @returns {Float64Array}
 */
export function invertLowerTriangular(L, d) {
  const invL = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    const diag = L[i * d + i];
    if (Math.abs(diag) < 1e-15) {
      throw new Error(`Singular triangular matrix on diagonal ${i}`);
    }
    invL[i * d + i] = 1.0 / diag;
    for (let j = 0; j < i; j++) {
      let sum = 0.0;
      for (let k = j; k < i; k++) {
        sum += L[i * d + k] * invL[k * d + j];
      }
      invL[i * d + j] = -sum / diag;
    }
  }
  return invL;
}

/**
 * Inverts a symmetric positive-definite matrix A using Cholesky factorization.
 * @param {Float64Array} A
 * @param {number} d
 * @returns {Float64Array}
 */
export function invertSymmetricPositiveDefinite(A, d) {
  const L = choleskyDecompose(A, d);
  const invL = invertLowerTriangular(L, d);
  const invA = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      let sum = 0.0;
      for (let k = 0; k < d; k++) {
        sum += invL[k * d + i] * invL[k * d + j];
      }
      invA[i * d + j] = sum;
    }
  }
  return invA;
}

/**
 * Symmetric Jacobi Eigenvalue Decomposition for arbitrary D x D symmetric matrix A.
 * Computes eigenvalues and orthogonal eigenvectors such that A = V * diag(eigenvalues) * V^T.
 * @param {Float64Array} A
 * @param {number} d
 * @param {number} [maxIter=100]
 * @param {number} [tol=1e-12]
 * @returns {{ eigenvalues: Float64Array, eigenvectors: Float64Array }}
 */
export function jacobiEigenvalueDecomposition(A, d, maxIter = 100, tol = 1e-12) {
  const S = new Float64Array(A);
  const V = eye(d);

  for (let iter = 0; iter < maxIter; iter++) {
    let maxOffDiag = 0.0;
    let p = 0;
    let q = 1;

    for (let i = 0; i < d; i++) {
      for (let j = i + 1; j < d; j++) {
        const absVal = Math.abs(S[i * d + j]);
        if (absVal > maxOffDiag) {
          maxOffDiag = absVal;
          p = i;
          q = j;
        }
      }
    }

    if (maxOffDiag < tol) break;

    const spq = S[p * d + q];
    const spp = S[p * d + p];
    const sqq = S[q * d + q];
    const theta = 0.5 * (sqq - spp) / spq;
    let t;
    if (theta >= 0.0) {
      t = 1.0 / (theta + Math.sqrt(1.0 + theta * theta));
    } else {
      t = -1.0 / (-theta + Math.sqrt(1.0 + theta * theta));
    }
    const c = 1.0 / Math.sqrt(1.0 + t * t);
    const s = t * c;
    const tau = s / (1.0 + c);

    S[p * d + q] = 0.0;
    S[q * d + p] = 0.0;
    S[p * d + p] = spp - t * spq;
    S[q * d + q] = sqq + t * spq;

    for (let r = 0; r < d; r++) {
      if (r !== p && r !== q) {
        const spr = S[p * d + r];
        const sqr = S[q * d + r];
        S[p * d + r] = spr - s * (sqr + tau * spr);
        S[r * d + p] = S[p * d + r];
        S[q * d + r] = sqr + s * (spr - tau * sqr);
        S[r * d + q] = S[q * d + r];
      }
    }

    for (let r = 0; r < d; r++) {
      const vpr = V[r * d + p];
      const vqr = V[r * d + q];
      V[r * d + p] = vpr - s * (vqr + tau * vpr);
      V[r * d + q] = vqr + s * (vpr - tau * vqr);
    }
  }

  const eigenvalues = new Float64Array(d);
  for (let i = 0; i < d; i++) eigenvalues[i] = S[i * d + i];

  const indices = new Int32Array(d);
  for (let i = 0; i < d; i++) indices[i] = i;
  indices.sort((a, b) => eigenvalues[b] - eigenvalues[a]);

  const sortedEigVals = new Float64Array(d);
  const sortedEigVecs = new Float64Array(d * d);

  for (let i = 0; i < d; i++) {
    const origIdx = indices[i];
    sortedEigVals[i] = eigenvalues[origIdx];
    for (let r = 0; r < d; r++) {
      sortedEigVecs[r * d + i] = V[r * d + origIdx];
    }
  }

  return { eigenvalues: sortedEigVals, eigenvectors: sortedEigVecs };
}

/**
 * Solves generalized symmetric eigenvalue problem A * v = \lambda * B * v for symmetric positive-definite B.
 * Returns sorted generalized eigenvalues \lambda_1 \ge \dots \ge \lambda_d.
 * @param {Float64Array} A - Symmetric matrix (D x D)
 * @param {Float64Array} B - Symmetric positive-definite matrix (D x D)
 * @param {number} d
 * @returns {Float64Array}
 */
export function solveGeneralizedSymmetricEigenvalues(A, B, d) {
  const L = choleskyDecompose(B, d);
  const invL = invertLowerTriangular(L, d);
  const T = matMul(invL, A, d);

  const C = new Float64Array(d * d);
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      let sum = 0.0;
      for (let k = 0; k < d; k++) {
        sum += T[i * d + k] * invL[j * d + k];
      }
      C[i * d + j] = sum;
    }
  }

  for (let i = 0; i < d; i++) {
    for (let j = i + 1; j < d; j++) {
      const avg = 0.5 * (C[i * d + j] + C[j * d + i]);
      C[i * d + j] = avg;
      C[j * d + i] = avg;
    }
  }

  const { eigenvalues } = jacobiEigenvalueDecomposition(C, d);
  return eigenvalues;
}

// ============================================================================
// SECTION 3: AUTOCOVARIANCE, SPECTRAL DENSITY & SOKAL/GEYER INTEGRATED AUTOCORRELATION TIME
// ============================================================================

/**
 * Computes sample autocovariance sequence \gamma_k for a 1D chain up to maxLag.
 * @param {ArrayLike<number>} chain
 * @param {number} [maxLag]
 * @param {number} [mean]
 * @returns {Float64Array}
 */
export function computeSampleAutocovariance(chain, maxLag = undefined, mean = undefined) {
  const n = chain.length;
  if (n === 0) return new Float64Array(0);
  const m = mean !== undefined ? mean : computeMean(chain);
  const K = maxLag !== undefined ? Math.min(maxLag, n - 1) : n - 1;
  const gamma = new Float64Array(K + 1);

  for (let k = 0; k <= K; k++) {
    let sum = 0.0;
    for (let t = 0; t < n - k; t++) {
      sum += (chain[t] - m) * (chain[t + k] - m);
    }
    gamma[k] = sum / n;
  }
  return gamma;
}

/**
 * Computes sample autocorrelation sequence \rho_k = \gamma_k / \gamma_0 for a 1D chain.
 * @param {ArrayLike<number>} chain
 * @param {number} [maxLag]
 * @returns {Float64Array}
 */
export function computeSampleAutocorrelation(chain, maxLag = undefined) {
  const gamma = computeSampleAutocovariance(chain, maxLag);
  if (gamma.length === 0 || Math.abs(gamma[0]) < 1e-20) {
    const rho = new Float64Array(gamma.length);
    if (rho.length > 0) rho[0] = 1.0;
    return rho;
  }
  const gamma0 = gamma[0];
  const rho = new Float64Array(gamma.length);
  for (let k = 0; k < gamma.length; k++) {
    rho[k] = gamma[k] / gamma0;
  }
  return rho;
}

/**
 * Computes multi-chain pooled autocorrelation function using cross-chain variograms [Vehtari et al. 2021].
 * Formula: \rho_k = 1 - V_k / (2 \hat{V}^+)
 * where V_k = \frac{1}{M(N - k)} \sum_{m=1}^M \sum_{n=k+1}^N (\theta_{m,n} - \theta_{m,n-k})^2.
 *
 * @param {number[][]|Float64Array[]} chains - Array of M chains each of length N
 * @param {number} [maxLag]
 * @returns {{ rho: Float64Array, vHatPlus: number, W: number, B: number }}
 */
export function computePooledAutocorrelation(chains, maxLag = undefined) {
  const M = chains.length;
  if (M < 1) throw new Error('Need at least 1 chain for autocorrelation.');
  const N = chains[0].length;
  if (N < 4) throw new Error('Chain length must be >= 4.');

  const chainMeans = new Float64Array(M);
  const chainVars = new Float64Array(M);
  let grandMean = 0.0;

  for (let m = 0; m < M; m++) {
    const mean = computeMean(chains[m]);
    chainMeans[m] = mean;
    grandMean += mean;
    chainVars[m] = computeVariance(chains[m], mean);
  }
  grandMean /= M;

  let W = 0.0;
  for (let m = 0; m < M; m++) W += chainVars[m];
  W /= M;

  let B = 0.0;
  if (M > 1) {
    for (let m = 0; m < M; m++) {
      const d = chainMeans[m] - grandMean;
      B += d * d;
    }
    B = (N / (M - 1)) * B;
  } else {
    B = 0.0;
  }

  const vHatPlus = M > 1 ? ((N - 1) / N) * W + ((M + 1) / (M * N)) * B : W;

  const K = maxLag !== undefined ? Math.min(maxLag, N - 1) : N - 1;
  const rho = new Float64Array(K + 1);
  rho[0] = 1.0;

  if (vHatPlus <= 1e-15) {
    return { rho, vHatPlus: 0.0, W: 0.0, B: 0.0 };
  }

  for (let k = 1; k <= K; k++) {
    let sumDiffSq = 0.0;
    for (let m = 0; m < M; m++) {
      const chain = chains[m];
      for (let n = k; n < N; n++) {
        const diff = chain[n] - chain[n - k];
        sumDiffSq += diff * diff;
      }
    }
    const Vk = sumDiffSq / (M * (N - k));
    rho[k] = 1.0 - Vk / (2.0 * vHatPlus);
  }

  return { rho, vHatPlus, W, B };
}

/**
 * Computes Sokal (1997) adaptive window cutoff integrated autocorrelation time \tau_{int}.
 * Window condition: find smallest K such that K >= c * \tau(K).
 *
 * @param {ArrayLike<number>} rho - Autocorrelation sequence \rho_0 = 1, \rho_1, ...
 * @param {number} [c=5.0] - Window scale constant (typically 4 to 6)
 * @returns {{ tauInt: number, windowCutoff: number, converged: boolean }}
 */
export function computeSokalAutocorrelationTime(rho, c = 5.0) {
  const maxK = rho.length - 1;
  let runningSum = 0.0;

  for (let k = 1; k <= maxK; k++) {
    runningSum += rho[k];
    const tau = 1.0 + 2.0 * runningSum;
    if (tau > 0 && k >= c * tau) {
      return {
        tauInt: Math.max(1.0, tau),
        windowCutoff: k,
        converged: true
      };
    }
  }

  const fallbackK = Math.floor(maxK / 2);
  let fallbackSum = 0.0;
  for (let k = 1; k <= fallbackK; k++) {
    fallbackSum += rho[k];
  }
  const fallbackTau = Math.max(1.0, 1.0 + 2.0 * fallbackSum);
  return {
    tauInt: fallbackTau,
    windowCutoff: fallbackK,
    converged: false
  };
}

/**
 * Computes Goodman & Weare (2010) affine-invariant ensemble MCMC autocorrelation time estimator.
 * @param {ArrayLike<number>} rho
 * @param {number} [c=5.0]
 * @returns {{ tauInt: number, windowCutoff: number }}
 */
export function computeGoodmanWeareAutocorrelationTime(rho, c = 5.0) {
  return computeSokalAutocorrelationTime(rho, c);
}

/**
 * Computes Geyer (1992) Initial Sequence Estimators for integrated autocorrelation time:
 * 1. Initial Positive Sequence Estimator (IPSE)
 * 2. Initial Monotone Sequence Estimator (IMSE)
 * 3. Initial Convex Sequence Estimator (ICSE)
 *
 * Uses paired autocorrelation sums \hat{P}_m = \hat{\rho}_{2m} + \hat{\rho}_{2m+1}.
 *
 * @param {ArrayLike<number>} rho - Autocorrelation sequence starting at lag 0
 * @returns {{
 *   tauIPSE: number,
 *   tauIMSE: number,
 *   tauICSE: number,
 *   cutoffLag: number,
 *   pairedSums: Float64Array
 * }}
 */
export function computeGeyerSequenceEstimators(rho) {
  const K = rho.length;
  const maxPairs = Math.floor((K - 1) / 2);
  if (maxPairs < 1) {
    return {
      tauIPSE: 1.0,
      tauIMSE: 1.0,
      tauICSE: 1.0,
      cutoffLag: 0,
      pairedSums: new Float64Array(0)
    };
  }

  const P = new Float64Array(maxPairs);
  let maxPosIndex = -1;

  for (let m = 0; m < maxPairs; m++) {
    const lag1 = 2 * m;
    const lag2 = 2 * m + 1;
    P[m] = rho[lag1] + rho[lag2];
    if (P[m] > 0 && maxPosIndex === m - 1) {
      maxPosIndex = m;
    }
  }

  if (maxPosIndex < 0) {
    return {
      tauIPSE: 1.0,
      tauIMSE: 1.0,
      tauICSE: 1.0,
      cutoffLag: 0,
      pairedSums: P
    };
  }

  // 1. Initial Positive Sequence (IPSE)
  let sumIPSE = -rho[0];
  for (let m = 0; m <= maxPosIndex; m++) {
    sumIPSE += 2.0 * P[m];
  }
  const tauIPSE = Math.max(1.0, sumIPSE);

  // 2. Initial Monotone Sequence (IMSE)
  const P_mono = new Float64Array(maxPosIndex + 1);
  P_mono[0] = P[0];
  for (let m = 1; m <= maxPosIndex; m++) {
    P_mono[m] = Math.min(P[m], P_mono[m - 1]);
  }

  let sumIMSE = -rho[0];
  for (let m = 0; m <= maxPosIndex; m++) {
    sumIMSE += 2.0 * P_mono[m];
  }
  const tauIMSE = Math.max(1.0, sumIMSE);

  // 3. Initial Convex Sequence (ICSE)
  const P_conv = new Float64Array(P_mono);
  for (let m = 1; m < maxPosIndex; m++) {
    const target = 0.5 * (P_conv[m - 1] + P_conv[m + 1]);
    if (P_conv[m] > target) {
      P_conv[m] = target;
    }
  }

  let sumICSE = -rho[0];
  for (let m = 0; m <= maxPosIndex; m++) {
    sumICSE += 2.0 * P_conv[m];
  }
  const tauICSE = Math.max(1.0, sumICSE);

  return {
    tauIPSE,
    tauIMSE,
    tauICSE,
    cutoffLag: 2 * maxPosIndex + 1,
    pairedSums: P
  };
}

/**
 * Spectral density at zero frequency S(0) using Tukey-Hanning cosine bell lag window.
 * @param {ArrayLike<number>} chain
 * @param {number} [maxLag]
 * @returns {number}
 */
export function computeSpectralDensityTukeyHanning(chain, maxLag = undefined) {
  const n = chain.length;
  if (n <= 1) return 0.0;
  const K = maxLag !== undefined ? maxLag : Math.min(Math.floor(Math.sqrt(n) * 3), n - 1);
  const gamma = computeSampleAutocovariance(chain, K);

  let S0 = gamma[0];
  for (let k = 1; k <= K; k++) {
    const weight = 0.5 * (1.0 + Math.cos((Math.PI * k) / (K + 1)));
    S0 += 2.0 * weight * gamma[k];
  }
  return Math.max(1e-15, S0);
}

/**
 * Spectral density at zero frequency S(0) using Parzen cubic spline lag window.
 * @param {ArrayLike<number>} chain
 * @param {number} [maxLag]
 * @returns {number}
 */
export function computeSpectralDensityParzen(chain, maxLag = undefined) {
  const n = chain.length;
  if (n <= 1) return 0.0;
  const K = maxLag !== undefined ? maxLag : Math.min(Math.floor(Math.sqrt(n) * 4), n - 1);
  const gamma = computeSampleAutocovariance(chain, K);

  let S0 = gamma[0];
  const halfM = K / 2.0;
  for (let k = 1; k <= K; k++) {
    const z = k / (K + 1.0);
    let weight = 0.0;
    if (k <= halfM) {
      weight = 1.0 - 6.0 * z * z + 6.0 * z * z * z;
    } else {
      weight = 2.0 * Math.pow(1.0 - z, 3.0);
    }
    S0 += 2.0 * weight * gamma[k];
  }
  return Math.max(1e-15, S0);
}

/**
 * Spectral density at zero frequency S(0) using Bartlett triangular lag window.
 * @param {ArrayLike<number>} chain
 * @param {number} [maxLag]
 * @returns {number}
 */
export function computeSpectralDensityBartlett(chain, maxLag = undefined) {
  const n = chain.length;
  if (n <= 1) return 0.0;
  const K = maxLag !== undefined ? maxLag : Math.min(Math.floor(Math.sqrt(n) * 3), n - 1);
  const gamma = computeSampleAutocovariance(chain, K);

  let S0 = gamma[0];
  for (let k = 1; k <= K; k++) {
    const weight = 1.0 - k / (K + 1.0);
    S0 += 2.0 * weight * gamma[k];
  }
  return Math.max(1e-15, S0);
}

/**
 * Spectral density at zero frequency S(0) using Quadratic Spectral (QS) kernel with Andrews (1991) AR(1) optimal bandwidth.
 * @param {ArrayLike<number>} chain
 * @returns {number}
 */
export function computeSpectralDensityAndrewsQS(chain) {
  const n = chain.length;
  if (n <= 2) return computeVariance(chain);

  const m = computeMean(chain);
  let num = 0.0;
  let den = 0.0;
  for (let t = 1; t < n; t++) {
    num += (chain[t] - m) * (chain[t - 1] - m);
    den += (chain[t - 1] - m) * (chain[t - 1] - m);
  }
  const rhoAR = den > 1e-15 ? Math.max(-0.99, Math.min(0.99, num / den)) : 0.0;

  const alpha1 = (4.0 * rhoAR * rhoAR) / Math.pow(1.0 - rhoAR, 4.0);
  const bandwidth = Math.max(1.0, 1.3221 * Math.pow(alpha1 * n, 0.2));

  const K = Math.min(Math.floor(bandwidth * 3.0), n - 1);
  const gamma = computeSampleAutocovariance(chain, K, m);

  let S0 = gamma[0];
  for (let k = 1; k <= K; k++) {
    const x = (6.0 * Math.PI * (k / bandwidth)) / 5.0;
    const sinX = Math.sin(x);
    const cosX = Math.cos(x);
    const kernelVal = (25.0 / (12.0 * Math.PI * Math.PI * (k / bandwidth) * (k / bandwidth))) * (sinX / x - cosX);
    S0 += 2.0 * kernelVal * gamma[k];
  }
  return Math.max(1e-15, S0);
}

/**
 * Computes Yule-Walker Autoregressive Spectral Density S(0) with AIC model order selection (p in 1..maxOrder).
 * @param {ArrayLike<number>} chain
 * @param {number} [maxOrder=10]
 * @returns {{ s0: number, selectedOrder: number, aic: Float64Array }}
 */
export function computeSpectralDensityYuleWalker(chain, maxOrder = 10) {
  const n = chain.length;
  if (n <= 4) return { s0: computeVariance(chain), selectedOrder: 0, aic: new Float64Array([0]) };

  const pMax = Math.min(maxOrder, Math.floor(n / 4));
  const gamma = computeSampleAutocovariance(chain, pMax);
  const gamma0 = gamma[0];

  if (gamma0 < 1e-15) {
    return { s0: 0.0, selectedOrder: 0, aic: new Float64Array([0]) };
  }

  const aic = new Float64Array(pMax + 1);
  aic[0] = n * Math.log(gamma0) + 2 * 1;

  let bestOrder = 0;
  let minAIC = aic[0];
  let bestCoeffs = new Float64Array(0);
  let bestSigmaSq = gamma0;

  let phi = new Float64Array(pMax + 1);
  let prevPhi = new Float64Array(pMax + 1);
  let sigmaSq = gamma0;

  for (let p = 1; p <= pMax; p++) {
    let num = gamma[p];
    for (let j = 1; j < p; j++) {
      num -= prevPhi[j] * gamma[p - j];
    }
    const kCoeff = num / sigmaSq;
    phi[p] = kCoeff;
    for (let j = 1; j < p; j++) {
      phi[j] = prevPhi[j] - kCoeff * prevPhi[p - j];
    }
    sigmaSq *= 1.0 - kCoeff * kCoeff;
    if (sigmaSq <= 1e-15) sigmaSq = 1e-15;

    for (let j = 1; j <= p; j++) prevPhi[j] = phi[j];

    const currentAIC = n * Math.log(sigmaSq) + 2 * (p + 1);
    aic[p] = currentAIC;

    if (currentAIC < minAIC) {
      minAIC = currentAIC;
      bestOrder = p;
      bestSigmaSq = sigmaSq;
      bestCoeffs = phi.slice(1, p + 1);
    }
  }

  let sumPhi = 0.0;
  for (let j = 0; j < bestCoeffs.length; j++) sumPhi += bestCoeffs[j];
  const denom = Math.max(1e-6, 1.0 - sumPhi);
  const s0 = bestSigmaSq / (denom * denom);

  return { s0, selectedOrder: bestOrder, aic };
}

// ============================================================================
// SECTION 4: EFFECTIVE SAMPLE SIZE (ESS) SUITE (BULK, TAIL, QUANTILE, MCSE)
// ============================================================================

/**
 * Computes univariate Effective Sample Size (ESS) using Geyer's Initial Monotone Sequence Estimator (IMSE).
 * Formula: N_{eff} = M * N / \hat{\tau}_{int}.
 *
 * @param {number[][]|Float64Array[]} chains - Array of M chains each of length N
 * @returns {{
 *   ess: number,
 *   tauInt: number,
 *   rho: Float64Array,
 *   cutoffLag: number
 * }}
 */
export function computeEffectiveSampleSize(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const { rho } = computePooledAutocorrelation(chains);
  const { tauIMSE, cutoffLag } = computeGeyerSequenceEstimators(rho);

  const tauInt = Math.max(1.0, tauIMSE);
  const ess = (M * N) / tauInt;

  return {
    ess,
    tauInt,
    rho,
    cutoffLag
  };
}

/**
 * Computes Rank-Normalized Bulk Effective Sample Size (Bulk-ESS) [Vehtari et al. 2021].
 * Normalizes values across all chains via rank transformation and probit score, then computes ESS.
 *
 * @param {number[][]|Float64Array[]} chains
 * @returns {{
 *   bulkESS: number,
 *   tauInt: number,
 *   transformedChains: Float64Array[]
 * }}
 */
export function computeBulkESS(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const totalSamples = M * N;

  const pooled = new Float64Array(totalSamples);
  let idx = 0;
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      pooled[idx++] = chains[m][n];
    }
  }

  const ranks = computeRanks(pooled);

  const transformedChains = [];
  let rankIdx = 0;
  for (let m = 0; m < M; m++) {
    const chainTrans = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const r = ranks[rankIdx++];
      const p = (r - 0.375) / (totalSamples + 0.25);
      chainTrans[n] = probit(p);
    }
    transformedChains.push(chainTrans);
  }

  const { ess, tauInt } = computeEffectiveSampleSize(transformedChains);
  return {
    bulkESS: ess,
    tauInt,
    transformedChains
  };
}

/**
 * Computes Tail Effective Sample Size (Tail-ESS) [Vehtari et al. 2021].
 * Computes indicator ESS at the 5% and 95% sample quantiles: I(\theta \le q_{0.05}) and I(\theta \le q_{0.95}).
 * Returns the minimum of the two quantile ESS values.
 *
 * @param {number[][]|Float64Array[]} chains
 * @param {number} [alpha=0.05]
 * @returns {{
 *   tailESS: number,
 *   essLower: number,
 *   essUpper: number,
 *   qLower: number,
 *   qUpper: number
 * }}
 */
export function computeTailESS(chains, alpha = 0.05) {
  const M = chains.length;
  const N = chains[0].length;
  const totalSamples = M * N;

  const pooled = new Float64Array(totalSamples);
  let idx = 0;
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      pooled[idx++] = chains[m][n];
    }
  }

  const qLower = computeQuantile(pooled, alpha);
  const qUpper = computeQuantile(pooled, 1.0 - alpha);

  const chainsLower = [];
  const chainsUpper = [];

  for (let m = 0; m < M; m++) {
    const indL = new Float64Array(N);
    const indU = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      indL[n] = chains[m][n] <= qLower ? 1.0 : 0.0;
      indU[n] = chains[m][n] <= qUpper ? 1.0 : 0.0;
    }
    chainsLower.push(indL);
    chainsUpper.push(indU);
  }

  const resLower = computeEffectiveSampleSize(chainsLower);
  const resUpper = computeEffectiveSampleSize(chainsUpper);

  const tailESS = Math.min(resLower.ess, resUpper.ess);

  return {
    tailESS,
    essLower: resLower.ess,
    essUpper: resUpper.ess,
    qLower,
    qUpper
  };
}

/**
 * Computes Quantile-specific ESS for an arbitrary probability quantile p \in (0, 1).
 * @param {number[][]|Float64Array[]} chains
 * @param {number} p
 * @returns {{ ess: number, tauInt: number, quantileValue: number }}
 */
export function computeQuantileESS(chains, p = 0.5) {
  const M = chains.length;
  const N = chains[0].length;
  const totalSamples = M * N;

  const pooled = new Float64Array(totalSamples);
  let idx = 0;
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      pooled[idx++] = chains[m][n];
    }
  }

  const qVal = computeQuantile(pooled, p);
  const indicatorChains = [];
  for (let m = 0; m < M; m++) {
    const ind = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      ind[n] = chains[m][n] <= qVal ? 1.0 : 0.0;
    }
    indicatorChains.push(ind);
  }

  const { ess, tauInt } = computeEffectiveSampleSize(indicatorChains);
  return {
    ess,
    tauInt,
    quantileValue: qVal
  };
}

/**
 * Computes Monte Carlo Standard Error (MCSE) for the posterior mean, standard deviation, and quantiles.
 * Formulas:
 *   MCSE(\mu) = \sigma / \sqrt{N_{eff,bulk}}
 *   MCSE(\sigma) = \sigma / \sqrt{2 * (N_{eff,bulk} - 1)}
 *
 * @param {number[][]|Float64Array[]} chains
 * @returns {{
 *   mcseMean: number,
 *   mcseStd: number,
 *   relMcseMean: number,
 *   mean: number,
 *   sd: number,
 *   bulkESS: number,
 *   tailESS: number
 * }}
 */
export function computeMCSEDiagnostics(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const total = M * N;

  const pooled = new Float64Array(total);
  let idx = 0;
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      pooled[idx++] = chains[m][n];
    }
  }

  const mean = computeMean(pooled);
  const sd = computeStd(pooled, mean);

  const { bulkESS } = computeBulkESS(chains);
  const { tailESS } = computeTailESS(chains);

  const safeBulkESS = Math.max(2.0, bulkESS);
  const mcseMean = sd / Math.sqrt(safeBulkESS);
  const mcseStd = sd / Math.sqrt(2.0 * (safeBulkESS - 1.0));
  const relMcseMean = sd > 1e-15 ? mcseMean / sd : 0.0;

  return {
    mcseMean,
    mcseStd,
    relMcseMean,
    mean,
    sd,
    bulkESS,
    tailESS
  };
}

// ============================================================================
// SECTION 5: GELMAN-RUBIN \hat{R} SUITE (UNIVARIATE, SPLIT, RANK, MULTIVARIATE)
// ============================================================================

/**
 * Computes standard univariate Gelman-Rubin potential scale reduction factor \hat{R}.
 *
 * @param {number[][]|Float64Array[]} chains - M chains of length N
 * @returns {{
 *   rHat: number,
 *   w: number,
 *   b: number,
 *   vHat: number,
 *   df: number,
 *   isConverged: boolean
 * }}
 */
export function computeStandardGelmanRubin(chains) {
  const M = chains.length;
  if (M < 2) throw new Error('Gelman-Rubin requires at least M >= 2 independent chains.');
  const N = chains[0].length;
  if (N < 4) throw new Error('Chain length N must be at least 4.');

  const chainMeans = new Float64Array(M);
  const chainVars = new Float64Array(M);
  let grandMean = 0.0;

  for (let m = 0; m < M; m++) {
    const mean = computeMean(chains[m]);
    chainMeans[m] = mean;
    grandMean += mean;
    chainVars[m] = computeVariance(chains[m], mean);
  }
  grandMean /= M;

  let W = 0.0;
  for (let m = 0; m < M; m++) W += chainVars[m];
  W /= M;

  let B = 0.0;
  for (let m = 0; m < M; m++) {
    const d = chainMeans[m] - grandMean;
    B += d * d;
  }
  B = (N / (M - 1)) * B;

  if (W <= 1e-15) {
    return { rHat: 1.0, w: W, b: B, vHat: 0.0, df: Infinity, isConverged: true };
  }

  const vHat = ((N - 1) / N) * W + ((M + 1) / (M * N)) * B;

  let varVHat = 0.0;
  for (let m = 0; m < M; m++) {
    const dW = chainVars[m] - W;
    varVHat += dW * dW;
  }
  varVHat = (varVHat / (M - 1)) * Math.pow((N - 1) / N, 2.0);

  const df = varVHat > 1e-15 ? (2.0 * vHat * vHat) / varVHat : Infinity;
  const dfCorrection = df < Infinity && df > 2 ? (df + 3.0) / (df + 1.0) : 1.0;

  const rHat = Math.sqrt(Math.max(1.0, (vHat / W) * dfCorrection));
  const isConverged = rHat < CF4_MCMC_CONSTANTS.R_HAT_CONVERGED;

  return { rHat, w: W, b: B, vHat, df, isConverged };
}

/**
 * Splits each of M chains in half, creating 2M chains of length N/2 to compute Split-\hat{R}.
 * Catches non-stationarity within individual chains.
 *
 * @param {number[][]|Float64Array[]} chains
 * @returns {{
 *   splitRHat: number,
 *   numSplitChains: number,
 *   splitChainLength: number,
 *   w: number,
 *   b: number,
 *   vHat: number,
 *   isConverged: boolean
 * }}
 */
export function computeSplitGelmanRubin(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const halfN = Math.floor(N / 2);

  if (halfN < 4) {
    throw new Error(`Chain length N=${N} too short for split R-hat (half-chain=${halfN}).`);
  }

  const splitChains = [];
  for (let m = 0; m < M; m++) {
    const c = chains[m];
    splitChains.push(c.slice(0, halfN));
    splitChains.push(c.slice(halfN, halfN * 2));
  }

  const res = computeStandardGelmanRubin(splitChains);
  return {
    splitRHat: res.rHat,
    numSplitChains: splitChains.length,
    splitChainLength: halfN,
    w: res.w,
    b: res.b,
    vHat: res.vHat,
    isConverged: res.isConverged
  };
}

/**
 * Computes Rank-Normalized Split-\hat{R} [Vehtari et al. 2021].
 * Robust against heavy tails, skewness, and multi-modality.
 *
 * @param {number[][]|Float64Array[]} chains
 * @returns {{
 *   rankSplitRHat: number,
 *   bulkRHat: number,
 *   foldedRHat: number,
 *   isConverged: boolean
 * }}
 */
export function computeRankNormalizedGelmanRubin(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const halfN = Math.floor(N / 2);

  const splitChains = [];
  for (let m = 0; m < M; m++) {
    const c = chains[m];
    splitChains.push(c.slice(0, halfN));
    splitChains.push(c.slice(halfN, halfN * 2));
  }

  const S = splitChains.length * halfN;
  const pooled = new Float64Array(S);
  let idx = 0;
  for (let m = 0; m < splitChains.length; m++) {
    for (let n = 0; n < halfN; n++) {
      pooled[idx++] = splitChains[m][n];
    }
  }

  // 1. Bulk R-hat on probit ranks
  const ranks = computeRanks(pooled);
  const zChains = [];
  let rIdx = 0;
  for (let m = 0; m < splitChains.length; m++) {
    const z = new Float64Array(halfN);
    for (let n = 0; n < halfN; n++) {
      const r = ranks[rIdx++];
      const p = (r - 0.375) / (S + 0.25);
      z[n] = probit(p);
    }
    zChains.push(z);
  }
  const bulkRes = computeStandardGelmanRubin(zChains);
  const bulkRHat = bulkRes.rHat;

  // 2. Folded transformation for Scale / Variance convergence: \zeta = |\theta - median(\theta)|
  const medianVal = computeQuantile(pooled, 0.5);
  const foldedPooled = new Float64Array(S);
  for (let i = 0; i < S; i++) {
    foldedPooled[i] = Math.abs(pooled[i] - medianVal);
  }

  const foldedRanks = computeRanks(foldedPooled);
  const foldedZChains = [];
  let fIdx = 0;
  for (let m = 0; m < splitChains.length; m++) {
    const z = new Float64Array(halfN);
    for (let n = 0; n < halfN; n++) {
      const r = foldedRanks[fIdx++];
      const p = (r - 0.375) / (S + 0.25);
      z[n] = probit(p);
    }
    foldedZChains.push(z);
  }
  const foldedRes = computeStandardGelmanRubin(foldedZChains);
  const foldedRHat = foldedRes.rHat;

  const rankSplitRHat = Math.max(bulkRHat, foldedRHat);
  const isConverged = rankSplitRHat < CF4_MCMC_CONSTANTS.R_HAT_CONVERGED;

  return {
    rankSplitRHat,
    bulkRHat,
    foldedRHat,
    isConverged
  };
}

/**
 * Computes Multivariate Potential Scale Reduction Factor (MPSRF) [Brooks & Gelman 1998].
 * Evaluates joint convergence across D parameters simultaneously.
 *
 * @param {Array<number[][]|Float64Array[]>} multiDimChains - Array of D parameters, each having M chains of length N
 * @returns {{
 *   mpsrf: number,
 *   maxEigenvalue: number,
 *   eigenvalues: Float64Array,
 *   wMatrix: Float64Array,
 *   bMatrix: Float64Array,
 *   vHatMatrix: Float64Array,
 *   isConverged: boolean
 * }}
 */
export function computeMultivariateGelmanRubin(multiDimChains) {
  const D = multiDimChains.length;
  if (D < 1) throw new Error('Multi-dimensional chains array must have D >= 1 parameters.');
  const M = multiDimChains[0].length;
  if (M < 2) throw new Error('Multivariate Gelman-Rubin requires M >= 2 chains.');
  const N = multiDimChains[0][0].length;
  if (N < 4) throw new Error('Chain length N must be at least 4.');

  // 1. Compute chain mean vectors (M vectors of length D) and grand mean vector (length D)
  const chainMeans = [];
  const grandMean = new Float64Array(D);

  for (let m = 0; m < M; m++) {
    const meanVec = new Float64Array(D);
    for (let d = 0; d < D; d++) {
      meanVec[d] = computeMean(multiDimChains[d][m]);
      grandMean[d] += meanVec[d];
    }
    chainMeans.push(meanVec);
  }
  for (let d = 0; d < D; d++) grandMean[d] /= M;

  // 2. Compute within-chain covariance matrix W (D x D)
  const W = createMatrix(D);
  for (let m = 0; m < M; m++) {
    const meanVec = chainMeans[m];
    for (let n = 0; n < N; n++) {
      for (let i = 0; i < D; i++) {
        const di = multiDimChains[i][m][n] - meanVec[i];
        for (let j = 0; j < D; j++) {
          const dj = multiDimChains[j][m][n] - meanVec[j];
          W[i * D + j] += di * dj;
        }
      }
    }
  }
  const normW = 1.0 / (M * (N - 1));
  for (let idx = 0; idx < D * D; idx++) W[idx] *= normW;

  // 3. Compute between-chain covariance matrix B / N (D x D)
  const BoverN = createMatrix(D);
  for (let m = 0; m < M; m++) {
    const meanVec = chainMeans[m];
    for (let i = 0; i < D; i++) {
      const di = meanVec[i] - grandMean[i];
      for (let j = 0; j < D; j++) {
        const dj = meanVec[j] - grandMean[j];
        BoverN[i * D + j] += di * dj;
      }
    }
  }
  const normB = 1.0 / (M - 1);
  for (let idx = 0; idx < D * D; idx++) BoverN[idx] *= normB;

  // 4. Compute estimated pooled posterior covariance matrix \hat{V}
  const vHatMatrix = createMatrix(D);
  for (let idx = 0; idx < D * D; idx++) {
    vHatMatrix[idx] = ((N - 1) / N) * W[idx] + ((M + 1) / M) * BoverN[idx];
  }

  // 5. Solve generalized symmetric eigenvalue problem (B/N) v = \lambda W v
  const Wreg = new Float64Array(W);
  for (let i = 0; i < D; i++) Wreg[i * D + i] += 1e-10;

  const eigenvalues = solveGeneralizedSymmetricEigenvalues(BoverN, Wreg, D);
  const maxEigenvalue = Math.max(0.0, eigenvalues[0]);

  const mpsrf = Math.sqrt(Math.max(1.0, (N - 1) / N + ((M + 1) / M) * maxEigenvalue));
  const isConverged = mpsrf < CF4_MCMC_CONSTANTS.R_HAT_CONVERGED;

  return {
    mpsrf,
    maxEigenvalue,
    eigenvalues,
    wMatrix: W,
    bMatrix: BoverN,
    vHatMatrix,
    isConverged
  };
}

// ============================================================================
// SECTION 6: CHAIN HEALTH DIAGNOSTICS (GEWEKE, HEIDELBERGER-WELCH, RAFTERY-LEWIS, E-BFMI, CUSUM)
// ============================================================================

/**
 * Computes Geweke Diagnostic Z-score for stationarity.
 * Tests H0: mean(first fracA) == mean(last fracB) using spectral density estimates.
 *
 * @param {ArrayLike<number>} chain - 1D MCMC chain
 * @param {number} [fracA=0.10] - First window fraction (default 10%)
 * @param {number} [fracB=0.50] - Second window fraction (default 50%)
 * @returns {{
 *   zScore: number,
 *   pValue: number,
 *   meanA: number,
 *   meanB: number,
 *   s0A: number,
 *   s0B: number,
 *   isStationary: boolean
 * }}
 */
export function computeGewekeDiagnostic(chain, fracA = 0.10, fracB = 0.50) {
  const n = chain.length;
  if (n < 50) throw new Error(`Geweke diagnostic requires n >= 50 samples, received ${n}`);

  const nA = Math.floor(fracA * n);
  const nB = Math.floor(fracB * n);
  const startB = n - nB;

  const segA = chain.slice(0, nA);
  const segB = chain.slice(startB);

  const meanA = computeMean(segA);
  const meanB = computeMean(segB);

  const s0A = computeSpectralDensityAndrewsQS(segA);
  const s0B = computeSpectralDensityAndrewsQS(segB);

  const seA = Math.sqrt(s0A / nA);
  const seB = Math.sqrt(s0B / nB);
  const denom = Math.sqrt(seA * seA + seB * seB);

  const zScore = denom > 1e-12 ? (meanA - meanB) / denom : 0.0;
  const pValue = normalTwoTailedPValue(zScore);
  const isStationary = Math.abs(zScore) < CF4_MCMC_CONSTANTS.GEWEKE_Z_CRITICAL;

  return {
    zScore,
    pValue,
    meanA,
    meanB,
    s0A,
    s0B,
    isStationary
  };
}

/**
 * Computes sliding Geweke path diagnostic over progressive chain discard fractions.
 * @param {ArrayLike<number>} chain
 * @param {number} [numIntervals=20]
 * @returns {Array<{ startFraction: number, startIndex: number, zScore: number, pValue: number, isStationary: boolean }>}
 */
export function computeSlidingGewekePath(chain, numIntervals = 20) {
  const n = chain.length;
  const results = [];
  const maxDiscardFraction = 0.5;

  for (let i = 0; i < numIntervals; i++) {
    const discardFrac = (i / numIntervals) * maxDiscardFraction;
    const startIdx = Math.floor(discardFrac * n);
    const subChain = chain.slice(startIdx);
    if (subChain.length < 50) break;

    const geweke = computeGewekeDiagnostic(subChain, 0.10, 0.50);
    results.push({
      startFraction: discardFrac,
      startIndex: startIdx,
      zScore: geweke.zScore,
      pValue: geweke.pValue,
      isStationary: geweke.isStationary
    });
  }
  return results;
}

/**
 * Heidelberger and Welch Stationarity & Relative Half-Width Precision Test (1983).
 * Uses Cramer-von Mises functional on the partial sum Brownian bridge.
 *
 * @param {ArrayLike<number>} chain
 * @param {number} [pAlpha=0.05]
 * @param {number} [targetEpsilon=0.02] - Relative half-width error threshold
 * @returns {{
 *   stationarityPassed: boolean,
 *   discardFraction: number,
 *   cvmStatistic: number,
 *   cvmPValue: number,
 *   halfWidthPassed: boolean,
 *   halfWidth: number,
 *   relativeHalfWidth: number,
 *   retainedMean: number
 * }}
 */
export function computeHeidelbergerWelchDiagnostic(chain, pAlpha = 0.05, targetEpsilon = 0.02) {
  const n = chain.length;
  if (n < 100) throw new Error('Heidelberger-Welch requires at least n >= 100 samples.');

  const discardSteps = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
  const cvmCritical = 0.461;

  for (let s = 0; s < discardSteps.length; s++) {
    const discard = discardSteps[s];
    const startIdx = Math.floor(discard * n);
    const subChain = chain.slice(startIdx);
    const N_sub = subChain.length;

    const mean = computeMean(subChain);
    const s0 = computeSpectralDensityTukeyHanning(subChain);

    let sumDev = 0.0;
    let cvmSum = 0.0;
    for (let t = 0; t < N_sub; t++) {
      sumDev += subChain[t] - mean;
      const bridge = (sumDev - ((t + 1) / N_sub) * 0.0) / Math.sqrt(N_sub * s0);
      cvmSum += bridge * bridge;
    }
    const cvmStat = cvmSum / N_sub;

    if (cvmStat < cvmCritical || s === discardSteps.length - 1) {
      const stationarityPassed = cvmStat < cvmCritical;
      const cvmPValue = Math.exp(-1.2337 * Math.max(0.0, cvmStat - 0.1));

      const se = Math.sqrt(s0 / N_sub);
      const halfWidth = 1.95996 * se;
      const relHalfWidth = Math.abs(mean) > 1e-12 ? halfWidth / Math.abs(mean) : halfWidth;
      const halfWidthPassed = relHalfWidth <= targetEpsilon;

      return {
        stationarityPassed,
        discardFraction: discard,
        cvmStatistic: cvmStat,
        cvmPValue,
        halfWidthPassed,
        halfWidth,
        relativeHalfWidth: relHalfWidth,
        retainedMean: mean
      };
    }
  }

  return {
    stationarityPassed: false,
    discardFraction: 0.5,
    cvmStatistic: 1.0,
    cvmPValue: 0.0,
    halfWidthPassed: false,
    halfWidth: 0.0,
    relativeHalfWidth: 1.0,
    retainedMean: computeMean(chain)
  };
}

/**
 * Raftery & Lewis (1992) Quantile Run-Length Diagnostic.
 * Computes required sample size N, burn-in M, and thinning factor I to estimate
 * quantile q with precision r at confidence level s.
 *
 * @param {ArrayLike<number>} chain
 * @param {number} [q=0.025] - Target quantile
 * @param {number} [r=0.005] - Precision tolerance
 * @param {number} [s=0.95] - Probability confidence
 * @returns {{
 *   requiredN: number,
 *   burnInM: number,
 *   minN: number,
 *   thinningFactorI: number,
 *   quantileEstimate: number
 * }}
 */
export function computeRafteryLewisDiagnostic(chain, q = 0.025, r = 0.005, s = 0.95) {
  const n = chain.length;
  if (n < 100) throw new Error('Raftery-Lewis requires at least n >= 100 samples.');

  const qVal = computeQuantile(chain, q);
  const phiZ = probit(0.5 * (1.0 + s));

  const binary = new Uint8Array(n);
  for (let i = 0; i < n; i++) binary[i] = chain[i] <= qVal ? 1 : 0;

  let n00 = 0, n01 = 0, n10 = 0, n11 = 0;
  for (let t = 0; t < n - 1; t++) {
    if (binary[t] === 0 && binary[t + 1] === 0) n00++;
    else if (binary[t] === 0 && binary[t + 1] === 1) n01++;
    else if (binary[t] === 1 && binary[t + 1] === 0) n10++;
    else if (binary[t] === 1 && binary[t + 1] === 1) n11++;
  }

  const alpha = (n01 + 1.0) / (n00 + n01 + 2.0);
  const beta = (n10 + 1.0) / (n10 + n11 + 2.0);
  const lambda = 1.0 - alpha - beta;

  const minN = Math.ceil((phiZ * phiZ * q * (1.0 - q)) / (r * r));

  const numerator = (2.0 - alpha - beta) * alpha * beta;
  const denominator = Math.pow(alpha + beta, 3.0);
  const factor = denominator > 1e-15 ? numerator / denominator : 1.0;
  const requiredN = Math.ceil((phiZ * phiZ * factor) / (r * r));

  const eps = 0.001;
  const absLambda = Math.abs(lambda);
  const burnInM = absLambda > 1e-6 && absLambda < 1.0 ? Math.ceil(Math.log(eps * (alpha + beta) / Math.max(alpha, beta)) / Math.log(absLambda)) : 10;

  const thinningFactorI = minN > 0 ? requiredN / minN : 1.0;

  return {
    requiredN,
    burnInM: Math.max(0, burnInM),
    minN,
    thinningFactorI,
    quantileEstimate: qVal
  };
}

/**
 * Computes Energy Bayesian Fraction of Missing Information (E-BFMI / Betancourt 2016).
 * Diagnoses momentum distribution exploration in HMC/NUTS algorithms.
 *
 * @param {ArrayLike<number>} energies - Hamiltonian energy values H(q, p) at each step
 * @returns {{
 *   ebfmi: number,
 *   energyVariance: number,
 *   energyDiffVariance: number,
 *   isAdequate: boolean
 * }}
 */
export function computeEBFMIDiagnostic(energies) {
  const n = energies.length;
  if (n < 4) throw new Error('E-BFMI requires at least 4 energy samples.');

  let diffSumSq = 0.0;
  for (let i = 1; i < n; i++) {
    const diff = energies[i] - energies[i - 1];
    diffSumSq += diff * diff;
  }
  const energyDiffVariance = diffSumSq / (n - 1);

  const meanE = computeMean(energies);
  let eSumSq = 0.0;
  for (let i = 0; i < n; i++) {
    const d = energies[i] - meanE;
    eSumSq += d * d;
  }
  const energyVariance = eSumSq / (n - 1);

  const ebfmi = energyVariance > 1e-15 ? energyDiffVariance / energyVariance : 0.0;
  const isAdequate = ebfmi >= CF4_MCMC_CONSTANTS.EBFMI_THRESHOLD;

  return {
    ebfmi,
    energyVariance,
    energyDiffVariance,
    isAdequate
  };
}

/**
 * Computes Cumulative Sum (CUSUM) Trajectory and Kolmogorov Maximum Path Excursion.
 * Detects persistent trends or drift in the MCMC chain.
 *
 * @param {ArrayLike<number>} chain
 * @returns {{
 *   cusumPath: Float64Array,
 *   maxExcursion: number,
 *   maxExcursionIndex: number,
 *   isStationary: boolean
 * }}
 */
export function computeCUSUMDiagnostic(chain) {
  const n = chain.length;
  if (n === 0) throw new Error('Empty chain for CUSUM.');

  const mean = computeMean(chain);
  const sd = computeStd(chain, mean);
  const denom = sd > 1e-15 ? sd * Math.sqrt(n) : 1.0;

  const cusumPath = new Float64Array(n);
  let running = 0.0;
  let maxExc = 0.0;
  let maxIdx = 0;

  for (let t = 0; t < n; t++) {
    running += chain[t] - mean;
    const val = running / denom;
    cusumPath[t] = val;
    const absVal = Math.abs(val);
    if (absVal > maxExc) {
      maxExc = absVal;
      maxIdx = t;
    }
  }

  const isStationary = maxExc < 1.358;

  return {
    cusumPath,
    maxExcursion: maxExc,
    maxExcursionIndex: maxIdx,
    isStationary
  };
}

/**
 * Computes empirical Pearson cross-correlation matrix across D parameters for multi-dimensional chains.
 * @param {Array<number[][]|Float64Array[]>} multiDimChains - [d][m][n]
 * @returns {Float64Array} - D x D correlation matrix
 */
export function computeCrossCorrelationMatrix(multiDimChains) {
  const D = multiDimChains.length;
  const M = multiDimChains[0].length;
  const N = multiDimChains[0][0].length;
  const total = M * N;

  const means = new Float64Array(D);
  const sds = new Float64Array(D);
  const pooled = [];

  for (let d = 0; d < D; d++) {
    const flat = new Float64Array(total);
    let idx = 0;
    for (let m = 0; m < M; m++) {
      for (let n = 0; n < N; n++) {
        flat[idx++] = multiDimChains[d][m][n];
      }
    }
    means[d] = computeMean(flat);
    sds[d] = computeStd(flat, means[d]);
    pooled.push(flat);
  }

  const corrMat = createMatrix(D);
  for (let i = 0; i < D; i++) {
    corrMat[i * D + i] = 1.0;
    for (let j = i + 1; j < D; j++) {
      let sumProd = 0.0;
      const flatI = pooled[i];
      const flatJ = pooled[j];
      const mi = means[i];
      const mj = means[j];

      for (let k = 0; k < total; k++) {
        sumProd += (flatI[k] - mi) * (flatJ[k] - mj);
      }
      const cov = sumProd / (total - 1);
      const denom = sds[i] * sds[j];
      const r = denom > 1e-15 ? cov / denom : 0.0;
      corrMat[i * D + j] = r;
      corrMat[j * D + i] = r;
    }
  }

  return corrMat;
}

// ============================================================================
// SECTION 7: SYNTHETIC COSMOLOGICAL MCMC GENERATOR FOR VERIFICATION
// ============================================================================

/**
 * Deterministic PRNG with Box-Muller normal and Student-t sampling.
 */
export class DeterministicPRNG {
  constructor(seed = 123456789) {
    this.seed = seed;
    this.hasSpare = false;
    this.spare = 0.0;
  }

  /** Uniform random number in [0, 1) */
  random() {
    this.seed = (this.seed * 1664525 + 1013904223) % 4294967296;
    return this.seed / 4294967296.0;
  }

  /** Standard normal N(0, 1) sample via Box-Muller */
  standardNormal() {
    if (this.hasSpare) {
      this.hasSpare = false;
      return this.spare;
    }
    let u, v, s;
    do {
      u = this.random() * 2.0 - 1.0;
      v = this.random() * 2.0 - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s === 0.0);

    const mul = Math.sqrt(-2.0 * Math.log(s) / s);
    this.spare = v * mul;
    this.hasSpare = true;
    return u * mul;
  }

  /** Student-t sample with nu degrees of freedom */
  studentT(nu = 4) {
    const z = this.standardNormal();
    let chi2 = 0.0;
    for (let i = 0; i < nu; i++) {
      const g = this.standardNormal();
      chi2 += g * g;
    }
    return z / Math.sqrt(chi2 / nu);
  }
}

/**
 * Generates synthetic multi-chain MCMC cosmological parameters with prescribed autocorrelation,
 * covariance structure, and optional non-stationary burn-in drift.
 *
 * @param {Object} options
 * @param {number} [options.numChains=4] - Number of independent chains M
 * @param {number} [options.chainLength=10000] - Number of steps N
 * @param {number[]} [options.means=[74.6, 0.31, 0.69, 0.81, 388.0, 310.0, 22.0]] - Means: [H0, Om, OL, s8, Vbulk, l_apex, b_apex]
 * @param {number[]} [options.sigmas=[1.2, 0.015, 0.015, 0.02, 28.0, 4.5, 3.2]] - Marginal standard deviations
 * @param {number} [options.ar1Persistence=0.70] - AR(1) persistence \phi \in [0, 0.99]
 * @param {number} [options.burnInTransientSteps=0] - Number of initial non-stationary drift steps
 * @param {number} [options.seed=42]
 * @returns {Array<Float64Array[]>} - Structure: [paramIndex][chainIndex][sampleIndex]
 */
export function generateSyntheticCosmologicalChains(options = {}) {
  const numChains = options.numChains || 4;
  const chainLength = options.chainLength || 10000;
  const means = options.means || [74.6, 0.31, 0.69, 0.81, 388.0, 310.0, 22.0];
  const sigmas = options.sigmas || [1.2, 0.015, 0.015, 0.02, 28.0, 4.5, 3.2];
  const phi = options.ar1Persistence !== undefined ? options.ar1Persistence : 0.70;
  const burnSteps = options.burnInTransientSteps || 0;
  const seed = options.seed || 42;

  const numParams = means.length;
  const prng = new DeterministicPRNG(seed);

  const innovScale = Math.sqrt(Math.max(1e-6, 1.0 - phi * phi));

  const result = [];
  for (let d = 0; d < numParams; d++) {
    const chainList = [];
    for (let m = 0; m < numChains; m++) {
      chainList.push(new Float64Array(chainLength));
    }
    result.push(chainList);
  }

  for (let m = 0; m < numChains; m++) {
    const currentVals = new Float64Array(numParams);
    for (let d = 0; d < numParams; d++) {
      currentVals[d] = means[d] + prng.standardNormal() * sigmas[d] * 2.0;
    }

    for (let t = 0; t < chainLength; t++) {
      let driftFactor = 1.0;
      if (burnSteps > 0 && t < burnSteps) {
        driftFactor = 1.0 + 3.0 * (1.0 - t / burnSteps);
      }

      for (let d = 0; d < numParams; d++) {
        const innov = prng.standardNormal() * sigmas[d] * innovScale;
        currentVals[d] = means[d] + phi * (currentVals[d] - means[d]) + innov;

        let sampleVal = currentVals[d];
        if (burnSteps > 0 && t < burnSteps) {
          sampleVal += (prng.random() - 0.5) * sigmas[d] * driftFactor;
        }

        result[d][m][t] = sampleVal;
      }
    }
  }

  return result;
}

// ============================================================================
// SECTION 8: CF4++ 10,000-STEP CONVERGENCE VERIFICATION DOSSIER & REPORT ENGINE
// ============================================================================

/**
 * Complete parameter-level convergence diagnostic record.
 * @typedef {Object} ParameterConvergenceRecord
 * @property {string} name
 * @property {number} mean
 * @property {number} sd
 * @property {number} median
 * @property {number[]} credibleInterval68
 * @property {number[]} credibleInterval95
 * @property {number} rHatStandard
 * @property {number} rHatSplit
 * @property {number} rHatRankSplit
 * @property {number} rHatFolded
 * @property {number} bulkESS
 * @property {number} tailESS
 * @property {number} tauInt
 * @property {number} gewekeZ
 * @property {number} gewekePValue
 * @property {boolean} gewekeStationary
 * @property {boolean} hwStationary
 * @property {boolean} hwHalfWidthPassed
 * @property {number} rafteryThinningI
 * @property {number} mcseMean
 * @property {number} relMcseMean
 * @property {boolean} isCertified
 */

/**
 * Evaluates full convergence diagnostic battery for a single multi-chain cosmological parameter.
 *
 * @param {number[][]|Float64Array[]} chains - M chains of length N
 * @param {string} paramName - Parameter name (e.g. 'H0', 'Omega_m')
 * @returns {ParameterConvergenceRecord}
 */
export function evaluateParameterConvergence(chains, paramName = 'param') {
  const M = chains.length;
  const N = chains[0].length;
  const total = M * N;

  const pooled = new Float64Array(total);
  let idx = 0;
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      pooled[idx++] = chains[m][n];
    }
  }

  const mean = computeMean(pooled);
  const sd = computeStd(pooled, mean);
  const median = computeQuantile(pooled, 0.5);

  const ci68 = [computeQuantile(pooled, 0.15865), computeQuantile(pooled, 0.84135)];
  const ci95 = [computeQuantile(pooled, 0.02275), computeQuantile(pooled, 0.97725)];

  const stdGelman = computeStandardGelmanRubin(chains);
  const splitGelman = computeSplitGelmanRubin(chains);
  const rankGelman = computeRankNormalizedGelmanRubin(chains);

  const { bulkESS } = computeBulkESS(chains);
  const { tailESS } = computeTailESS(chains);
  const { tauInt } = computeEffectiveSampleSize(chains);

  const geweke = computeGewekeDiagnostic(chains[0]);
  const hw = computeHeidelbergerWelchDiagnostic(chains[0]);
  const raftery = computeRafteryLewisDiagnostic(chains[0]);
  const mcse = computeMCSEDiagnostics(chains);

  const isCertified = (
    rankGelman.rankSplitRHat <= CF4_MCMC_CONSTANTS.R_HAT_CONVERGED &&
    bulkESS >= CF4_MCMC_CONSTANTS.BULK_ESS_THRESHOLD &&
    tailESS >= CF4_MCMC_CONSTANTS.TAIL_ESS_THRESHOLD &&
    geweke.isStationary &&
    mcse.relMcseMean <= CF4_MCMC_CONSTANTS.MCSE_MAX_REL_SD
  );

  return {
    name: paramName,
    mean,
    sd,
    median,
    credibleInterval68: ci68,
    credibleInterval95: ci95,
    rHatStandard: stdGelman.rHat,
    rHatSplit: splitGelman.splitRHat,
    rHatRankSplit: rankGelman.rankSplitRHat,
    rHatFolded: rankGelman.foldedRHat,
    bulkESS,
    tailESS,
    tauInt,
    gewekeZ: geweke.zScore,
    gewekePValue: geweke.pValue,
    gewekeStationary: geweke.isStationary,
    hwStationary: hw.stationarityPassed,
    hwHalfWidthPassed: hw.halfWidthPassed,
    rafteryThinningI: raftery.thinningFactorI,
    mcseMean: mcse.mcseMean,
    relMcseMean: mcse.relMcseMean,
    isCertified
  };
}

/**
 * Master Verification Dossier Engine for multi-chain CF4++ MCMC realizations.
 * Runs comprehensive univariate, split, rank-normalized, and multivariate diagnostics,
 * computes overall certification score (0 - 100), and outputs Markdown and JSON dossiers.
 *
 * @param {Array<number[][]|Float64Array[]>} multiDimChains - [paramIndex][chainIndex][sampleIndex]
 * @param {Object} [options]
 * @param {string[]} [options.paramNames] - Array of parameter names
 * @param {ArrayLike<number>} [options.energyTrace] - Optional HMC energy trace
 * @param {string} [options.datasetCitation='Dupuy & Courtois (2023)']
 * @returns {{
 *   dossierId: string,
 *   timestamp: string,
 *   citation: string,
 *   numParameters: number,
 *   numChains: number,
 *   chainLength: number,
 *   mpsrf: number,
 *   ebfmi: number|null,
 *   parameters: ParameterConvergenceRecord[],
 *   overallConvergenceScore: number,
 *   certificationStatus: 'CERTIFIED_CONVERGED' | 'MARGINAL_CONVERGENCE' | 'FAILED_NON_CONVERGED',
 *   failedGates: string[],
 *   markdownReport: string,
 *   jsonSummary: Object
 * }}
 */
export function generateCF4ConvergenceDossier(multiDimChains, options = {}) {
  const D = multiDimChains.length;
  if (D < 1) throw new Error('Multi-dimensional chains must have at least 1 parameter.');
  const M = multiDimChains[0].length;
  const N = multiDimChains[0][0].length;

  const defaultNames = ['H0', 'Omega_M', 'Omega_Lambda', 'sigma_8', 'V_bulk', 'l_apex', 'b_apex'];
  const paramNames = options.paramNames || defaultNames.slice(0, D);
  const citation = options.datasetCitation || 'CosmicFlows-4++ (Dupuy & Courtois 2023)';
  const energyTrace = options.energyTrace || null;

  const paramRecords = [];
  const failedGates = [];

  for (let d = 0; d < D; d++) {
    const name = d < paramNames.length ? paramNames[d] : `param_${d}`;
    const record = evaluateParameterConvergence(multiDimChains[d], name);
    paramRecords.push(record);

    if (record.rHatRankSplit > CF4_MCMC_CONSTANTS.R_HAT_CONVERGED) {
      failedGates.push(`Gate 1 (Rank-Split R-hat): ${name} R-hat = ${record.rHatRankSplit.toFixed(4)} > ${CF4_MCMC_CONSTANTS.R_HAT_CONVERGED}`);
    }
    if (record.bulkESS < CF4_MCMC_CONSTANTS.BULK_ESS_THRESHOLD) {
      failedGates.push(`Gate 2 (Bulk-ESS): ${name} Bulk-ESS = ${record.bulkESS.toFixed(1)} < ${CF4_MCMC_CONSTANTS.BULK_ESS_THRESHOLD}`);
    }
    if (record.tailESS < CF4_MCMC_CONSTANTS.TAIL_ESS_THRESHOLD) {
      failedGates.push(`Gate 3 (Tail-ESS): ${name} Tail-ESS = ${record.tailESS.toFixed(1)} < ${CF4_MCMC_CONSTANTS.TAIL_ESS_THRESHOLD}`);
    }
    if (!record.gewekeStationary) {
      failedGates.push(`Gate 4 (Geweke Stationarity): ${name} |Z| = ${Math.abs(record.gewekeZ).toFixed(2)} > 1.96`);
    }
    if (record.relMcseMean > CF4_MCMC_CONSTANTS.MCSE_MAX_REL_SD) {
      failedGates.push(`Gate 5 (MCSE Precision): ${name} MCSE/SD = ${(record.relMcseMean * 100).toFixed(2)}% > 5%`);
    }
  }

  let mpsrfVal = 1.0;
  if (D > 1) {
    const multiGelman = computeMultivariateGelmanRubin(multiDimChains);
    mpsrfVal = multiGelman.mpsrf;
    if (mpsrfVal > CF4_MCMC_CONSTANTS.R_HAT_CONVERGED) {
      failedGates.push(`Gate 6 (Multivariate MPSRF): MPSRF = ${mpsrfVal.toFixed(4)} > ${CF4_MCMC_CONSTANTS.R_HAT_CONVERGED}`);
    }
  }

  let ebfmiVal = null;
  if (energyTrace && energyTrace.length >= 4) {
    const ebfmiDiag = computeEBFMIDiagnostic(energyTrace);
    ebfmiVal = ebfmiDiag.ebfmi;
    if (!ebfmiDiag.isAdequate) {
      failedGates.push(`Gate 7 (E-BFMI Energy): E-BFMI = ${ebfmiVal.toFixed(3)} < ${CF4_MCMC_CONSTANTS.EBFMI_THRESHOLD}`);
    }
  }

  let score = 100.0;
  let maxRHat = 1.0;
  let minBulkESS = Infinity;
  for (let d = 0; d < D; d++) {
    if (paramRecords[d].rHatRankSplit > maxRHat) maxRHat = paramRecords[d].rHatRankSplit;
    if (paramRecords[d].bulkESS < minBulkESS) minBulkESS = paramRecords[d].bulkESS;
  }
  if (maxRHat > 1.01) {
    score -= Math.min(40.0, (maxRHat - 1.01) * 400.0);
  }
  if (minBulkESS < 1000.0) {
    score -= Math.min(30.0, (1000.0 - minBulkESS) / 25.0);
  }
  if (mpsrfVal > 1.01) {
    score -= Math.min(20.0, (mpsrfVal - 1.01) * 200.0);
  }
  if (failedGates.length > 0) {
    score -= failedGates.length * 5.0;
  }
  score = Math.max(0.0, Math.min(100.0, score));

  let certificationStatus = 'CERTIFIED_CONVERGED';
  if (failedGates.length > 0 || score < 75.0) {
    certificationStatus = score >= 50.0 ? 'MARGINAL_CONVERGENCE' : 'FAILED_NON_CONVERGED';
  }

  const now = new Date().toISOString();
  const dossierId = `CF4-MCMC-DIAG-${Date.now()}`;

  let md = `# CF4++ MCMC Convergence & Gelman-Rubin Verification Dossier\n\n`;
  md += `**Dossier ID:** \`${dossierId}\`  \n`;
  md += `**Timestamp:** ${now}  \n`;
  md += `**Citation Reference:** ${citation}  \n`;
  md += `**Chains Analyzed:** $M = ${M}$, Chain Length $N = ${N}$ ($S_{total} = ${M * N}$ samples)  \n`;
  md += `**Multivariate PSRF (MPSRF):** \`${mpsrfVal.toFixed(4)}\`  \n`;
  if (ebfmiVal !== null) {
    md += `**HMC Energy E-BFMI:** \`${ebfmiVal.toFixed(4)}\`  \n`;
  }
  md += `**Overall Convergence Score:** \`${score.toFixed(1)} / 100\`  \n`;
  md += `**Certification Gate Status:** **${certificationStatus}**  \n\n`;

  md += `## 1. Cosmological Parameter Posterior Convergence Table\n\n`;
  md += `| Parameter | Mean ± SD | 68.3% CI | 95.4% CI | $\\hat{R}_{\\text{rank}}$ | Bulk ESS | Tail ESS | $\\tau_{\\text{int}}$ | Geweke $Z$ | MCSE/SD | Status |\n`;
  md += `|---|---|---|---|---|---|---|---|---|---|---|\n`;

  for (let d = 0; d < D; d++) {
    const p = paramRecords[d];
    const meanSd = `${p.mean.toFixed(2)} ± ${p.sd.toFixed(2)}`;
    const ci68Str = `[${p.credibleInterval68[0].toFixed(2)}, ${p.credibleInterval68[1].toFixed(2)}]`;
    const ci95Str = `[${p.credibleInterval95[0].toFixed(2)}, ${p.credibleInterval95[1].toFixed(2)}]`;
    const statusStr = p.isCertified ? '✅ PASS' : '❌ FAIL';
    const relMcseStr = `${(p.relMcseMean * 100).toFixed(1)}%`;

    md += `| **${p.name}** | ${meanSd} | ${ci68Str} | ${ci95Str} | ${p.rHatRankSplit.toFixed(4)} | ${Math.round(p.bulkESS)} | ${Math.round(p.tailESS)} | ${p.tauInt.toFixed(1)} | ${p.gewekeZ.toFixed(2)} | ${relMcseStr} | ${statusStr} |\n`;
  }

  md += `\n## 2. Quality Verification Gates\n\n`;
  if (failedGates.length === 0) {
    md += `> [!NOTE]\n> All CF4++ scientific convergence gates PASSED with high confidence ($\\hat{R} < 1.05$, $\\text{ESS} > 400$, Geweke $|Z| < 1.96$, $\\text{MPSRF} < 1.05$).\n\n`;
  } else {
    md += `> [!WARNING]\n> ${failedGates.length} convergence gate(s) failed or exhibited marginal tension:\n`;
    for (let f = 0; f < failedGates.length; f++) {
      md += `> - ${failedGates[f]}\n`;
    }
    md += `\n`;
  }

  const jsonSummary = {
    dossierId,
    timestamp: now,
    citation,
    numParameters: D,
    numChains: M,
    chainLength: N,
    totalSamples: M * N,
    mpsrf: mpsrfVal,
    ebfmi: ebfmiVal,
    score,
    certificationStatus,
    failedGates,
    parameterTable: paramRecords
  };

  return {
    dossierId,
    timestamp: now,
    citation,
    numParameters: D,
    numChains: M,
    chainLength: N,
    mpsrf: mpsrfVal,
    ebfmi: ebfmiVal,
    parameters: paramRecords,
    overallConvergenceScore: score,
    certificationStatus,
    failedGates,
    markdownReport: md,
    jsonSummary
  };
}
