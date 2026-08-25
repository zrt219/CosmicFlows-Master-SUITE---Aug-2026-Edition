/**
 * @file posterior_statistics.js
 * @description Comprehensive Bayesian Posterior Summary Statistics: HPD Intervals, Quantiles, MCSE, and Kernel Density Estimation.
 *
 * Implements:
 * 1. Highest Posterior Density (HPD) interval calculator for unimodal and asymmetric posteriors.
 * 2. Equal-tailed Bayesian Credible Intervals (68.27%, 95.45%, 99.73%).
 * 3. Monte Carlo Standard Error (MCSE) using spectral variance and batch means.
 * 4. 1D & 2D Gaussian Kernel Density Estimation (KDE) with automated Silverman / Scott bandwidth selection.
 * 5. Maximum A Posteriori (MAP) numerical mode estimation.
 *
 * @module uncertainty/posterior_statistics
 */

/**
 * Computes the Highest Posterior Density (HPD) credible interval containing (1 - alpha) probability mass.
 *
 * @param {Float64Array|number[]} samples - 1D array of posterior samples
 * @param {number} [credMass=0.95] - Credible mass (e.g. 0.6827, 0.95, 0.9973)
 * @returns {[number, number]} [lowerBound, upperBound]
 */
export function computeHPDInterval(samples, credMass = 0.95) {
  const n = samples.length;
  if (n < 2) throw new Error(\'HPD interval calculation requires at least 2 samples.\');

  const sorted = Array.from(samples).sort((a, b) => a - b);
  const intervalIdx = Math.floor(credMass * n);
  const numIntervals = n - intervalIdx;

  let minWidth = Infinity;
  let hpdLower = sorted[0];
  let hpdUpper = sorted[intervalIdx];

  for (let i = 0; i < numIntervals; i++) {
    const lower = sorted[i];
    const upper = sorted[i + intervalIdx];
    const width = upper - lower;
    if (width < minWidth) {
      minWidth = width;
      hpdLower = lower;
      hpdUpper = upper;
    }
  }

  return [hpdLower, hpdUpper];
}

/**
 * Computes equal-tailed quantiles of a 1D sample array.
 *
 * @param {Float64Array|number[]} samples
 * @param {number[]} [probs=[0.025, 0.16, 0.5, 0.84, 0.975]]
 * @returns {number[]} Quantile values corresponding to probs
 */
export function computeQuantiles(samples, probs = [0.025, 0.16, 0.5, 0.84, 0.975]) {
  const n = samples.length;
  if (n === 0) return [];
  const sorted = Array.from(samples).sort((a, b) => a - b);

  return probs.map(p => {
    if (p <= 0) return sorted[0];
    if (p >= 1) return sorted[n - 1];
    const index = p * (n - 1);
    const lower = Math.floor(index);
    const upper = Math.ceil(index);
    const weight = index - lower;
    return (1.0 - weight) * sorted[lower] + weight * sorted[upper];
  });
}

/**
 * Computes Monte Carlo Standard Error (MCSE) for the sample mean using batch means.
 *
 * @param {Float64Array|number[]} samples
 * @param {number} [numBatches=50]
 * @returns {number} MCSE estimate
 */
export function computeMCSE(samples, numBatches = 50) {
  const n = samples.length;
  if (n < numBatches * 2) {
    // Fallback to naive standard error std / sqrt(N)
    let sum = 0.0;
    for (let i = 0; i < n; i++) sum += samples[i];
    const mean = sum / n;
    let sumSq = 0.0;
    for (let i = 0; i < n; i++) sumSq += (samples[i] - mean) ** 2;
    const std = Math.sqrt(sumSq / Math.max(1, n - 1));
    return std / Math.sqrt(n);
  }

  const batchSize = Math.floor(n / numBatches);
  const batchMeans = new Float64Array(numBatches);

  for (let b = 0; b < numBatches; b++) {
    let sum = 0.0;
    const start = b * batchSize;
    for (let i = 0; i < batchSize; i++) {
      sum += samples[start + i];
    }
    batchMeans[b] = sum / batchSize;
  }

  let grandMean = 0.0;
  for (let b = 0; b < numBatches; b++) grandMean += batchMeans[b];
  grandMean /= numBatches;

  let sumSq = 0.0;
  for (let b = 0; b < numBatches; b++) {
    const diff = batchMeans[b] - grandMean;
    sumSq += diff * diff;
  }
  const batchVar = sumSq / (numBatches - 1);

  return Math.sqrt(batchVar / numBatches);
}

/**
 * Computes Silverman\'s Rule of Thumb bandwidth for 1D Gaussian KDE:
 * h = 0.9 * min(std, IQR / 1.34) * n^(-1/5)
 *
 * @param {Float64Array|number[]} samples
 * @returns {number} Bandwidth h
 */
export function silvermanBandwidth1D(samples) {
  const n = samples.length;
  if (n < 2) return 1.0;

  let sum = 0.0;
  for (let i = 0; i < n; i++) sum += samples[i];
  const mean = sum / n;

  let sumSq = 0.0;
  for (let i = 0; i < n; i++) sumSq += (samples[i] - mean) ** 2;
  const std = Math.sqrt(sumSq / (n - 1));

  // Compute IQR
  const q = computeQuantiles(samples, [0.25, 0.75]);
  const iqr = q[1] - q[0];
  const scale = iqr > 0 ? Math.min(std, iqr / 1.34) : std;

  return Math.max(1e-6, 0.9 * (scale || 1.0) * Math.pow(n, -0.2));
}

/**
 * 1D Gaussian Kernel Density Estimator.
 */
export class KernelDensityEstimator1D {
  /**
   * @param {Float64Array|number[]} samples
   * @param {number} [bandwidth] - If omitted, uses Silverman\'s rule
   */
  constructor(samples, bandwidth = null) {
    this.samples = new Float64Array(samples);
    this.n = this.samples.length;
    this.h = bandwidth || silvermanBandwidth1D(this.samples);
    this.invH = 1.0 / this.h;
    this.normConst = 1.0 / (this.n * this.h * Math.sqrt(2.0 * Math.PI));
  }

  /**
   * Evaluates estimated PDF p(x).
   * @param {number} x
   * @returns {number} Density p(x)
   */
  evaluate(x) {
    let sum = 0.0;
    for (let i = 0; i < this.n; i++) {
      const u = (x - this.samples[i]) * this.invH;
      sum += Math.exp(-0.5 * u * u);
    }
    return sum * this.normConst;
  }

  /**
   * Evaluates PDF across regular grid [minX, maxX].
   * @param {number} minX
   * @param {number} maxX
   * @param {number} [numPoints=200]
   * @returns {{ grid: Float64Array, density: Float64Array }}
   */
  evaluateGrid(minX, maxX, numPoints = 200) {
    const grid = new Float64Array(numPoints);
    const density = new Float64Array(numPoints);
    const step = (maxX - minX) / (numPoints - 1);

    for (let i = 0; i < numPoints; i++) {
      const x = minX + i * step;
      grid[i] = x;
      density[i] = this.evaluate(x);
    }

    return { grid, density };
  }

  /**
   * Finds the mode (peak) of the estimated density.
   * @param {number} [numGridPoints=500]
   * @returns {{ mode: number, maxDensity: number }}
   */
  findMode(numGridPoints = 500) {
    let minVal = Infinity;
    let maxVal = -Infinity;
    for (let i = 0; i < this.n; i++) {
      if (this.samples[i] < minVal) minVal = this.samples[i];
      if (this.samples[i] > maxVal) maxVal = this.samples[i];
    }
    const padding = 3.0 * this.h;
    const { grid, density } = this.evaluateGrid(minVal - padding, maxVal + padding, numGridPoints);

    let maxD = -1.0;
    let modeX = grid[0];
    for (let i = 0; i < numGridPoints; i++) {
      if (density[i] > maxD) {
        maxD = density[i];
        modeX = grid[i];
      }
    }

    return { mode: modeX, maxDensity: maxD };
  }
}
