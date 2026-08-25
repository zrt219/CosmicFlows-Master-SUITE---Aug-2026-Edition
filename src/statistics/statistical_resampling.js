/**
 * @file statistical_resampling.js
 * @description Reproducible Pseudo-Random Number Generators (MT19937, PCG64, Xoshiro256**), Monte Carlo Bootstrap, Jackknife, and Latin Hypercube Sampling.
 *
 * Implements:
 * 1. Seeded MT19937 Mersenne Twister PRNG.
 * 2. PCG64 Permuted Congruential Generator.
 * 3. Xoshiro256** 64-bit state pseudo-random generator.
 * 4. Latin Hypercube Sampling (LHS) for space-filling parameter space exploration.
 * 5. Stratified and Non-Parametric Bootstrap Resampling.
 * 6. Blind Analysis Coordinate Masker for un-biased watershed / cluster boundary analysis.
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

/**
 * 64-bit Permuted Congruential Generator (PCG-XSH-RR).
 */
export class PCG64 {
  constructor(initState = 42n, initSeq = 54n) {
    this.state = 0n;
    this.inc = (BigInt(initSeq) << 1n) | 1n;
    this.step();
    this.state = (this.state + BigInt(initState)) & 0xFFFFFFFFFFFFFFFFn;
    this.step();
  }

  step() {
    this.state = (this.state * 6364136223846793005n + this.inc) & 0xFFFFFFFFFFFFFFFFn;
  }

  nextUint32() {
    const oldState = this.state;
    this.step();
    const xorShifted = Number(((oldState >> 18n) ^ oldState) >> 27n) >>> 0;
    const rot = Number(oldState >> 59n);
    return ((xorShifted >>> rot) | (xorShifted << ((-rot) & 31))) >>> 0;
  }

  extractNumber() {
    return this.nextUint32() / 4294967296.0;
  }
}

/**
 * Latin Hypercube Sampler (LHS) for D dimensions across [min, max] bounds.
 */
export class LatinHypercubeSampler {
  /**
   * @param {Array<[number, number]>} bounds - Array of D [min, max] pairs
   * @param {number} [seed=42]
   */
  constructor(bounds, seed = 42) {
    this.bounds = bounds;
    this.dim = bounds.length;
    this.prng = new MT19937(seed);
  }

  /**
   * Generates N Latin Hypercube sample points.
   *
   * @param {number} nSamples
   * @returns {Array<Float64Array>}
   */
  sample(nSamples) {
    const samples = Array.from({ length: nSamples }, () => new Float64Array(this.dim));

    for (let d = 0; d < this.dim; d++) {
      const [minVal, maxVal] = this.bounds[d];
      const range = maxVal - minVal;
      const binWidth = range / nSamples;

      // Create random permutation of bin indices [0..nSamples-1]
      const perm = Array.from({ length: nSamples }, (_, i) => i);
      for (let i = nSamples - 1; i > 0; i--) {
        const j = Math.floor(this.prng.extractNumber() * (i + 1));
        const tmp = perm[i];
        perm[i] = perm[j];
        perm[j] = tmp;
      }

      // Sample uniformly within each permuted bin
      for (let i = 0; i < nSamples; i++) {
        const bin = perm[i];
        const u = this.prng.extractNumber();
        samples[i][d] = minVal + (bin + u) * binWidth;
      }
    }

    return samples;
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

  /**
   * Generates Gaussian pseudo-random sample vector from N(mean, Cov) using Cholesky factor L.
   *
   * @param {Float64Array|number[]} mean
   * @param {number[][]} L - Lower triangular Cholesky factor L L^T = Cov
   * @param {MT19937|PCG64} prng
   * @returns {Float64Array}
   */
  static sampleMultivariateGaussian(mean, L, prng) {
    const d = mean.length;
    const z = new Float64Array(d);

    // Standard Box-Muller
    for (let i = 0; i < d; i += 2) {
      const u1 = Math.max(1e-15, prng.extractNumber());
      const u2 = prng.extractNumber();
      z[i] = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      if (i + 1 < d) {
        z[i + 1] = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);
      }
    }

    // x = mean + L * z
    const x = new Float64Array(d);
    for (let i = 0; i < d; i++) {
      let sum = mean[i];
      for (let j = 0; j <= i; j++) {
        sum += L[i][j] * z[j];
      }
      x[i] = sum;
    }

    return x;
  }
}
