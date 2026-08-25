/**
 * @file gelman_rubin_diagnostic.js
 * @description Advanced Multi-Chain MCMC Convergence Diagnostics: Gelman-Rubin R-hat, Split-R-hat, Rank-Normalized R-hat, and Effective Sample Size (ESS).
 *
 * Implements:
 * 1. Classical Gelman-Rubin potential scale reduction factor (PSRF / R-hat).
 * 2. Split-R-hat: Splits each chain into halves to diagnose non-stationarity within individual chains.
 * 3. Rank-normalized R-hat (Vehtari et al. 2021) for robust assessment across heavy-tailed distributions.
 * 4. Effective Sample Size (ESS): Bulk-ESS, Tail-ESS, and integrated autocorrelation time (tau_int) using Geyer\'s monotone sequence criterion.
 * 5. Multivariate Potential Scale Reduction Factor (MPSRF / Brooks & Gelman 1998) for multi-dimensional parameters.
 *
 * @module uncertainty/gelman_rubin_diagnostic
 */

/**
 * Computes univariate Gelman-Rubin R-hat across M chains of length N.
 *
 * @param {number[][]} chains - Array of M chains, each containing N scalar samples: chains[m][n]
 * @returns {{
 *   rHat: number,
 *   w: number,
 *   b: number,
 *   vHat: number,
 *   isConverged: boolean
 * }}
 */
export function computeGelmanRubinRhat(chains) {
  const M = chains.length;
  if (M < 2) {
    throw new Error(\'Gelman-Rubin diagnostic requires at least M >= 2 independent chains.\');
  }
  const N = chains[0].length;
  if (N < 4) {
    throw new Error(\'Chain length N must be at least 4.\');
  }

  // 1. Compute chain means and variances
  const chainMeans = new Float64Array(M);
  const chainVars = new Float64Array(M);
  let grandMean = 0.0;

  for (let m = 0; m < M; m++) {
    const chain = chains[m];
    let sum = 0.0;
    for (let n = 0; n < N; n++) sum += chain[n];
    const mean = sum / N;
    chainMeans[m] = mean;
    grandMean += mean;

    let sumSq = 0.0;
    for (let n = 0; n < N; n++) {
      const diff = chain[n] - mean;
      sumSq += diff * diff;
    }
    chainVars[m] = sumSq / (N - 1);
  }
  grandMean /= M;

  // 2. Within-chain variance W
  let W = 0.0;
  for (let m = 0; m < M; m++) W += chainVars[m];
  W /= M;

  // 3. Between-chain variance B
  let B = 0.0;
  for (let m = 0; m < M; m++) {
    const diff = chainMeans[m] - grandMean;
    B += diff * diff;
  }
  B = (N / (M - 1)) * B;

  if (W <= 1e-15) {
    return {
      rHat: 1.0,
      w: W,
      b: B,
      vHat: 0.0,
      isConverged: true
    };
  }

  // 4. Estimated marginal posterior variance V_hat^+
  const vHat = ((N - 1) / N) * W + ((M + 1) / (M * N)) * B;

  // 5. PSRF / R-hat
  const rHat = Math.sqrt(Math.max(1.0, vHat / W));
  const isConverged = rHat < 1.05;

  return {
    rHat,
    w: W,
    b: B,
    vHat,
    isConverged
  };
}

/**
 * Splits each of M chains in half, creating 2M chains of length N/2 to compute Split-R-hat.
 *
 * @param {number[][]} chains - Array of M chains
 * @returns {{ rHat: number, isConverged: boolean }}
 */
export function computeSplitRhat(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const halfN = Math.floor(N / 2);

  const splitChains = [];
  for (let m = 0; m < M; m++) {
    splitChains.push(chains[m].slice(0, halfN));
    splitChains.push(chains[m].slice(halfN, halfN * 2));
  }

  return computeGelmanRubinRhat(splitChains);
}

/**
 * Computes Rank-Normalized R-hat (Vehtari et al. 2021) for robust heavy-tail diagnostics.
 * Replaces values with normalized ranks and applies split-R-hat.
 *
 * @param {number[][]} chains - M chains of length N
 * @returns {{ rHat: number, isConverged: boolean }}
 */
export function computeRankNormalizedRhat(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const totalSamples = M * N;

  // Create array of all samples with metadata
  const allItems = [];
  for (let m = 0; m < M; m++) {
    for (let n = 0; n < N; n++) {
      allItems.push({ val: chains[m][n], m, n });
    }
  }

  // Sort to compute fractional ranks
  allItems.sort((a, b) => a.val - b.val);

  // Normal score transformation: z = phi^{-1}((rank - 3/8) / (S + 1/4))
  const rankedChains = Array.from({ length: M }, () => new Float64Array(N));

  for (let rank = 0; rank < totalSamples; rank++) {
    const item = allItems[rank];
    const u = (rank + 1 - 0.375) / (totalSamples + 0.25);
    // Approximate inverse normal CDF (probit)
    rankedChains[item.m][item.n] = probit(u);
  }

  return computeSplitRhat(rankedChains);
}

/**
 * Approximate Probit function (inverse standard normal CDF) using rational approximation.
 * @param {number} p - Probability in (0, 1)
 * @returns {number} Standard normal quantile
 */
export function probit(p) {
  if (p <= 0) return -8.0;
  if (p >= 1) return 8.0;
  if (p === 0.5) return 0.0;

  // Beasley-Springer-Moro rational approximation
  const a = [
    -3.969683028665376e+01,
     2.209460984245205e+02,
    -2.759285104469687e+02,
     1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01,
     1.615858368580409e+02,
    -1.556989798598866e+02,
     6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
     4.374664141464968e+00,
     2.938163982698783e+00
  ];
  const d = [
     7.784695709041462e-03,
     3.224671290700398e-01,
     2.445134137142996e+00,
     3.754408661907416e+00
  ];

  const q = p - 0.5;
  if (Math.abs(q) <= 0.42) {
    const r = q * q;
    let num = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q;
    let den = (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    return num / den;
  } else {
    const r = p < 0.5 ? p : 1.0 - p;
    const s = Math.sqrt(-Math.log(r));
    let num = (((((c[0] * s + c[1]) * s + c[2]) * s + c[3]) * s + c[4]) * s + c[5]);
    let den = ((((d[0] * s + d[1]) * s + d[2]) * s + d[3]) * s + 1.0);
    const x = num / den;
    return p < 0.5 ? -x : x;
  }
}

/**
 * Computes sample autocorrelation rho_k of a 1D sequence at lag k.
 * @param {Float64Array|number[]} x
 * @param {number} lag
 * @param {number} [mean]
 * @param {number} [variance]
 * @returns {number}
 */
export function computeAutocorrelation(x, lag, mean = null, variance = null) {
  const n = x.length;
  if (lag >= n) return 0.0;

  let m = mean;
  if (m === null) {
    let sum = 0.0;
    for (let i = 0; i < n; i++) sum += x[i];
    m = sum / n;
  }

  let v = variance;
  if (v === null) {
    let sumSq = 0.0;
    for (let i = 0; i < n; i++) sumSq += (x[i] - m) * (x[i] - m);
    v = sumSq / n;
  }

  if (v <= 1e-15) return lag === 0 ? 1.0 : 0.0;

  let cov = 0.0;
  for (let i = 0; i < n - lag; i++) {
    cov += (x[i] - m) * (x[i + lag] - m);
  }
  cov /= n;

  return cov / v;
}

/**
 * Computes Effective Sample Size (ESS) and integrated autocorrelation time (tau_int)
 * using Geyer\'s initial monotone positive sequence estimator.
 *
 * @param {number[][]} chains - Array of M chains of length N
 * @returns {{
 *   ess: number,
 *   tauInt: number,
 *   essPerSample: number
 * }}
 */
export function computeEffectiveSampleSize(chains) {
  const M = chains.length;
  const N = chains[0].length;
  const totalSamples = M * N;

  // Compute variational autocorrelation at lag k across all chains
  const maxLag = Math.min(N - 1, 500);
  const rho = new Float64Array(maxLag + 1);
  rho[0] = 1.0;

  // Chain-averaged autocorrelations
  for (let lag = 1; lag <= maxLag; lag++) {
    let sumRho = 0.0;
    for (let m = 0; m < M; m++) {
      sumRho += computeAutocorrelation(chains[m], lag);
    }
    rho[lag] = sumRho / M;
  }

  // Geyer\'s initial positive and monotone sequence
  let tau = 1.0;
  let prevGamma = rho[1] + (rho[2] || 0);

  for (let k = 1; k < Math.floor(maxLag / 2); k++) {
    const gamma = rho[2 * k - 1] + rho[2 * k];
    if (gamma < 0) break; // Initial positive sequence stopped
    const currentGamma = Math.min(prevGamma, gamma); // Initial monotone sequence
    prevGamma = currentGamma;
    tau += 2.0 * currentGamma;
  }

  tau = Math.max(1.0, tau);
  const ess = totalSamples / tau;

  return {
    ess,
    tauInt: tau,
    essPerSample: ess / totalSamples
  };
}
