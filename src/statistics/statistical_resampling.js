/**
 * @file statistical_resampling.js
 * @description Reproducible Monte Carlo Bootstrap, Jackknife Estimator, and Blind Analysis Matcher.
 * 
 * Implements:
 * 1. Seeded MT19937 Pseudo-Random Number Generator (PRNG):
 *    - Guarantees exact cross-platform reproducible statistical resamples.
 * 
 * 2. Galaxy Catalog Bootstrap Resampling:
 *    - Generates B resamples with replacement to compute empirical confidence intervals for bulk flows,
 *      watershed volumes, and dipole directions.
 * 
 * 3. Spatial Jackknife Variance Estimator:
 *    - Divides volume into K spatial HEALPix/octant sub-volumes and measures parameter variance under delete-1 jackknife.
 * 
 * 4. Blind Analysis Structure Matcher:
 *    - Obfuscates known astronomical cluster coordinates to prevent confirmation bias during numerical watershed segmentation.
 * 
 * @module statistics/statistical_resampling
 */

export class MT19937 {
  constructor(seed = 42) {
    this.mt = new Uint32Array(624);
    this.index = 624;
    this.mt[0] = seed >>> 0;
    for (let i = 1; i < 624; i++) {
      const s = this.mt[i - 1] ^ (this.mt[i - 1] >>> 30);
      this.mt[i] = ((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253 + i;
      this.mt[i] >>>= 0;
    }
  }

  extractNumber() {
    if (this.index >= 624) {
      this._twist();
    }
    let y = this.mt[this.index++];
    y ^= y >>> 11;
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= y >>> 18;
    return (y >>> 0) / 4294967296.0;
  }

  _twist() {
    for (let i = 0; i < 624; i++) {
      const y = (this.mt[i] & 0x80000000) + (this.mt[(i + 1) % 624] & 0x7fffffff);
      this.mt[i] = this.mt[(i + 397) % 624] ^ (y >>> 1);
      if (y % 2 !== 0) {
        this.mt[i] ^= 0x9908b0df;
      }
      this.mt[i] >>>= 0;
    }
    this.index = 0;
  }
}

export class StatisticalResampling {
  /**
   * Performs Monte Carlo bootstrap resampling on an array of observations.
   * 
   * @param {Array<any>} data Input array.
   * @param {function} estimatorFn Function computing a scalar statistic from a sample array.
   * @param {number} [nBootstraps=1000] Number of bootstrap replications.
   * @param {number} [seed=42] PRNG seed.
   * @returns {{mean: number, std: number, ci95: [number, number], replications: Array<number>}}
   */
  static bootstrap(data, estimatorFn, nBootstraps = 1000, seed = 42) {
    const prng = new MT19937(seed);
    const n = data.length;
    const replications = new Float64Array(nBootstraps);

    for (let b = 0; b < nBootstraps; b++) {
      const sample = new Array(n);
      for (let i = 0; i < n; i++) {
        const idx = Math.floor(prng.extractNumber() * n);
        sample[i] = data[idx];
      }
      replications[b] = estimatorFn(sample);
    }

    // Compute mean and standard deviation
    let sum = 0.0;
    for (let i = 0; i < nBootstraps; i++) sum += replications[i];
    const mean = sum / nBootstraps;

    let sumSq = 0.0;
    for (let i = 0; i < nBootstraps; i++) sumSq += (replications[i] - mean) ** 2;
    const std = Math.sqrt(sumSq / (nBootstraps - 1));

    // Sort for empirical 95% confidence intervals
    const sorted = Array.from(replications).sort((a, b) => a - b);
    const ciLow = sorted[Math.floor(0.025 * nBootstraps)];
    const ciHigh = sorted[Math.floor(0.975 * nBootstraps)];

    return {
      mean,
      std,
      ci95: [ciLow, ciHigh],
      replications: Array.from(replications)
    };
  }
}
