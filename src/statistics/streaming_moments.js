/**
 * @file streaming_moments.js
 * @description Cosmological Pairwise Velocity Statistics: Mean Infall v_12(r), Velocity Dispersion Tensor \sigma_12^2(r), and Streaming Models.
 *
 * Implements:
 * 1. Mean pairwise streaming velocity v_12(r) = \langle (v_1 - v_2) \cdot \hat{r}_{12} \rangle.
 * 2. Pairwise velocity dispersion tensor: Parallel \sigma_\parallel^2(r) and Perpendicular \sigma_\perp^2(r).
 * 3. Davis-Peebles linear infall model: v_12(r) = -2/3 * f * H0 * r * \xi(r) / (1 + \xi(r)).
 * 4. Pair-weighted streaming velocity profiles across radial bins in Mpc/h.
 *
 * @module statistics/streaming_moments
 */

/**
 * Pairwise Velocity Streaming Moments Calculator.
 */
export class StreamingMomentsCalculator {
  /**
   * @param {Array<{ x: number, y: number, z: number, vx: number, vy: number, vz: number }>} particles
   * @param {Float64Array|number[]} binEdges - Radial separation bin edges in Mpc/h
   */
  constructor(particles, binEdges) {
    this.particles = particles;
    this.n = particles.length;
    this.binEdges = new Float64Array(binEdges);
    this.numBins = this.binEdges.length - 1;
    this.rCenters = new Float64Array(this.numBins);
    for (let b = 0; b < this.numBins; b++) {
      this.rCenters[b] = 0.5 * (this.binEdges[b] + this.binEdges[b + 1]);
    }
  }

  /**
   * Computes empirical pairwise moments: mean infall v_12(r), parallel dispersion sigma_par(r),
   * and perpendicular dispersion sigma_perp(r).
   *
   * @returns {{
   *   r: Float64Array,
   *   meanInfallV12: Float64Array,
   *   sigmaParallel: Float64Array,
   *   sigmaPerpendicular: Float64Array,
   *   pairCounts: Float64Array
   * }}
   */
  computePairwiseMoments() {
    const N = this.n;
    const numBins = this.numBins;
    const minR = this.binEdges[0];
    const maxR = this.binEdges[numBins];
    const minRSq = minR * minR;
    const maxRSq = maxR * maxR;

    const pairCounts = new Float64Array(numBins);
    const sumVPar = new Float64Array(numBins);
    const sumVParSq = new Float64Array(numBins);
    const sumVPerpSq = new Float64Array(numBins);

    for (let i = 0; i < N; i++) {
      const p1 = this.particles[i];
      for (let j = i + 1; j < N; j++) {
        const p2 = this.particles[j];

        // Separation vector r_12 = r_1 - r_2
        const dx = p1.x - p2.x;
        const dy = p1.y - p2.y;
        const dz = p1.z - p2.z;
        const distSq = dx * dx + dy * dy + dz * dz;

        if (distSq >= minRSq && distSq < maxRSq) {
          const dist = Math.sqrt(distSq);
          const nHat = [dx / dist, dy / dist, dz / dist];

          // Relative velocity v_12 = v_1 - v_2
          const dvx = p1.vx - p2.vx;
          const dvy = p1.vy - p2.vy;
          const dvz = p1.vz - p2.vz;

          // Parallel relative velocity v_par = v_12 . nHat
          const vPar = dvx * nHat[0] + dvy * nHat[1] + dvz * nHat[2];

          // Perpendicular relative velocity squared = |v_12|^2 - v_par^2
          const vTotSq = dvx * dvx + dvy * dvy + dvz * dvz;
          const vPerpSq = Math.max(0.0, vTotSq - vPar * vPar);

          // Find bin
          for (let b = 0; b < numBins; b++) {
            if (dist >= this.binEdges[b] && dist < this.binEdges[b + 1]) {
              pairCounts[b] += 1.0;
              sumVPar[b] += vPar;
              sumVParSq[b] += vPar * vPar;
              sumVPerpSq[b] += vPerpSq;
              break;
            }
          }
        }
      }
    }

    const meanInfall = new Float64Array(numBins);
    const sigmaPar = new Float64Array(numBins);
    const sigmaPerp = new Float64Array(numBins);

    for (let b = 0; b < numBins; b++) {
      const count = pairCounts[b];
      if (count > 1) {
        const meanP = sumVPar[b] / count;
        meanInfall[b] = meanP;

        const varPar = (sumVParSq[b] - count * meanP * meanP) / (count - 1);
        sigmaPar[b] = Math.sqrt(Math.max(0.0, varPar));

        // 2 perpendicular degrees of freedom => divide by 2
        const varPerp = sumVPerpSq[b] / (2.0 * count);
        sigmaPerp[b] = Math.sqrt(Math.max(0.0, varPerp));
      }
    }

    return {
      r: this.rCenters,
      meanInfallV12: meanInfall,
      sigmaParallel: sigmaPar,
      sigmaPerpendicular: sigmaPerp,
      pairCounts
    };
  }

  /**
   * Computes Davis-Peebles linear theoretical infall profile:
   * v_12(r) = - (2/3) * H0 * f * r * xi(r) / (1 + xi(r)).
   *
   * @param {Float64Array|number[]} rValues - Separation distances in Mpc/h
   * @param {Float64Array|number[]} xiValues - Correlation function values xi(r)
   * @param {number} [h0f=47.18] - H0 * f (e.g. 74.6 * 0.315^0.55)
   * @returns {Float64Array} Theoretical v_12(r) in km/s
   */
  static davisPeeblesModel(rValues, xiValues, h0f = 47.18) {
    const n = rValues.length;
    const v12 = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      const r = rValues[i];
      const xi = xiValues[i];
      const factor = xi / (1.0 + Math.max(0.0, xi));
      v12[i] = -(2.0 / 3.0) * h0f * r * factor;
    }

    return v12;
  }
}
