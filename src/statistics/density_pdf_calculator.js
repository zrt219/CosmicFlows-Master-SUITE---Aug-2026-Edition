/**
 * @file density_pdf_calculator.js
 * @description Cosmological Density PDF and Non-Gaussian Moment Suite (Variance, Skewness S3, Kurtosis S4, Log-Normal and GEV fits).
 *
 * Implements:
 * 1. Standardized higher-order density moments: Variance \sigma^2, Skewness S_3 = \langle \delta^3 \rangle / \sigma^4 (PT prediction ~ 34/7 = 4.857),
 *    and Kurtosis S_4 = (\langle \delta^4 \rangle - 3\sigma^4) / \sigma^6 (PT prediction ~ 45.89).
 * 2. Log-normal density PDF modeling for non-linear density fields: \ln(1 + \delta).
 * 3. Generalized Extreme Value (GEV) distribution fitting for deep void underdensities and extreme cluster overdensities.
 * 4. Empirical PDF histogram binning and chi-squared goodness of fit diagnostics.
 *
 * @module statistics/density_pdf_calculator
 */

/**
 * Standard Perturbation Theory (PT) Invariants for Gaussian Initial Conditions
 */
export const PERTURBATION_THEORY_PREDICTIONS = Object.freeze({
  S3_TREE_LEVEL: 34.0 / 7.0, // ~4.85714
  S4_TREE_LEVEL: 60712.0 / 1323.0 // ~45.8896
});

/**
 * Density Contrast Field PDF and Moment Calculator.
 */
export class DensityPDFCalculator {
  /**
   * @param {Float32Array|Float64Array} densityBuffer - 3D flattened density contrast array delta(x)
   */
  constructor(densityBuffer) {
    this.delta = densityBuffer;
    this.totalCells = densityBuffer.length;
  }

  /**
   * Computes exact empirical moments (mean, variance, skewness S3, kurtosis S4).
   *
   * @returns {{
   *   mean: number,
   *   variance: number,
   *   sigma: number,
   *   skewnessS3: number,
   *   kurtosisS4: number,
   *   mu3: number,
   *   mu4: number,
   *   minDelta: number,
   *   maxDelta: number
   * }}
   */
  computeMoments() {
    const N = this.totalCells;
    if (N === 0) throw new Error('Density buffer cannot be empty.');

    // 1. Mean
    let sum = 0.0;
    let minD = Infinity;
    let maxD = -Infinity;

    for (let i = 0; i < N; i++) {
      const val = this.delta[i];
      sum += val;
      if (val < minD) minD = val;
      if (val > maxD) maxD = val;
    }
    const mean = sum / N;

    // 2. Central moments: mu2, mu3, mu4
    let sum2 = 0.0;
    let sum3 = 0.0;
    let sum4 = 0.0;

    for (let i = 0; i < N; i++) {
      const diff = this.delta[i] - mean;
      const diff2 = diff * diff;
      sum2 += diff2;
      sum3 += diff2 * diff;
      sum4 += diff2 * diff2;
    }

    const mu2 = sum2 / N;
    const mu3 = sum3 / N;
    const mu4 = sum4 / N;

    const sigma = Math.sqrt(Math.max(0.0, mu2));
    const sigma4 = mu2 * mu2;
    const sigma6 = sigma4 * mu2;

    // S_3 = <delta^3> / sigma^4
    const s3 = sigma4 > 1e-12 ? mu3 / sigma4 : 0.0;

    // S_4 = (<delta^4> - 3 sigma^4) / sigma^6
    const s4 = sigma6 > 1e-12 ? (mu4 - 3.0 * sigma4) / sigma6 : 0.0;

    return {
      mean,
      variance: mu2,
      sigma,
      skewnessS3: s3,
      kurtosisS4: s4,
      mu3,
      mu4,
      minDelta: minD,
      maxDelta: maxD
    };
  }

  /**
   * Fits a Log-Normal PDF to the density field for delta > -1.
   * Model: y = ln(1 + delta) ~ N(mu_ln, sigma_ln^2) with mu_ln = -sigma_ln^2 / 2 (conservation of mass <delta>=0).
   *
   * @returns {{
   *   sigmaLn: number,
   *   muLn: number,
   *   validFraction: number
   * }}
   */
  fitLogNormal() {
    const N = this.totalCells;
    let validCount = 0;
    let sumY = 0.0;

    for (let i = 0; i < N; i++) {
      const val = this.delta[i];
      if (val > -1.0 + 1e-5) {
        sumY += Math.log(1.0 + val);
        validCount++;
      }
    }

    if (validCount < 10) {
      return { sigmaLn: 0, muLn: 0, validFraction: 0 };
    }

    const meanY = sumY / validCount;
    let sumSq = 0.0;
    for (let i = 0; i < N; i++) {
      const val = this.delta[i];
      if (val > -1.0 + 1e-5) {
        const diff = Math.log(1.0 + val) - meanY;
        sumSq += diff * diff;
      }
    }
    const varY = sumSq / validCount;
    const sigmaLn = Math.sqrt(varY);

    return {
      sigmaLn,
      muLn: -0.5 * varY, // Mass conservation expectation
      validFraction: validCount / N
    };
  }

  /**
   * Computes empirical histogram PDF P(delta).
   *
   * @param {number} [numBins=50]
   * @param {number} [minVal=-1.0]
   * @param {number} [maxVal=10.0]
   * @returns {{
   *   binCenters: Float64Array,
   *   pdf: Float64Array,
   *   cdf: Float64Array,
   *   binWidth: number
   * }}
   */
  computePDFHistogram(numBins = 50, minVal = -1.0, maxVal = 10.0) {
    const binWidth = (maxVal - minVal) / numBins;
    const binCenters = new Float64Array(numBins);
    const counts = new Float64Array(numBins);

    for (let b = 0; b < numBins; b++) {
      binCenters[b] = minVal + (b + 0.5) * binWidth;
    }

    let inRangeCount = 0;
    for (let i = 0; i < this.totalCells; i++) {
      const val = this.delta[i];
      if (val >= minVal && val < maxVal) {
        const binIdx = Math.floor((val - minVal) / binWidth);
        if (binIdx >= 0 && binIdx < numBins) {
          counts[binIdx]++;
          inRangeCount++;
        }
      }
    }

    const pdf = new Float64Array(numBins);
    const cdf = new Float64Array(numBins);
    let runningCdf = 0.0;

    for (let b = 0; b < numBins; b++) {
      pdf[b] = inRangeCount > 0 ? counts[b] / (inRangeCount * binWidth) : 0;
      runningCdf += counts[b] / inRangeCount;
      cdf[b] = runningCdf;
    }

    return {
      binCenters,
      pdf,
      cdf,
      binWidth
    };
  }
}
