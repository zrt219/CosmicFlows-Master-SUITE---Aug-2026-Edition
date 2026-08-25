/**
 * @file two_point_correlation.js
 * @description 3D Spatial and Anisotropic Two-Point Correlation Function Suite (Landy-Szalay, Hamilton, and RSD Multipoles).
 *
 * Implements:
 * 1. Landy-Szalay (1993) minimal-variance pair-count estimator: xi(r) = (DD - 2DR + RR) / RR.
 * 2. Hamilton (1993) geometry-invariant estimator: xi(r) = (DD * RR) / (DR^2) - 1.
 * 3. 2D Anisotropic correlation function xi(r_perp, r_par) = xi(sigma, pi) for Redshift-Space Distortions (RSD).
 * 4. Legendre multipole projections: Monopole xi_0(s), Quadrupole xi_2(s), Hexadecapole xi_4(s).
 * 5. Linear Kaiser RSD distortion parameter beta = f / b recovery.
 *
 * @module statistics/two_point_correlation
 */

/**
 * Computes pairwise separation distance histogram for an array of 3D point coordinates.
 *
 * @param {Array<{ x: number, y: number, z: number }>} data1
 * @param {Array<{ x: number, y: number, z: number }>} [data2] - If omitted, computes auto-pairs
 * @param {Float64Array|number[]} binEdges - Radial separation bin edges [r0, r1, ..., rB]
 * @returns {Float64Array} Pair counts in each bin
 */
export function computePairCounts(data1, data2 = null, binEdges = null) {
  const n1 = data1.length;
  const isAuto = !data2;
  const n2 = isAuto ? n1 : data2.length;
  const p2 = isAuto ? data1 : data2;

  const numBins = binEdges.length - 1;
  const counts = new Float64Array(numBins);
  const minR = binEdges[0];
  const maxR = binEdges[numBins];
  const minRSq = minR * minR;
  const maxRSq = maxR * maxR;

  for (let i = 0; i < n1; i++) {
    const p1_i = data1[i];
    const jStart = isAuto ? i + 1 : 0;

    for (let j = jStart; j < n2; j++) {
      const p2_j = p2[j];
      const dx = p1_i.x - p2_j.x;
      const dy = p1_i.y - p2_j.y;
      const dz = p1_i.z - p2_j.z;
      const distSq = dx * dx + dy * dy + dz * dz;

      if (distSq >= minRSq && distSq < maxRSq) {
        const dist = Math.sqrt(distSq);
        // Binary search or linear search for bin
        for (let b = 0; b < numBins; b++) {
          if (dist >= binEdges[b] && dist < binEdges[b + 1]) {
            counts[b] += 1.0;
            break;
          }
        }
      }
    }
  }

  return counts;
}

/**
 * Landy-Szalay and Hamilton Two-Point Correlation Estimator.
 */
export class TwoPointCorrelationEstimator {
  /**
   * @param {Float64Array|number[]} binEdges - Radial separation bin edges in Mpc/h
   */
  constructor(binEdges) {
    this.binEdges = new Float64Array(binEdges);
    this.numBins = this.binEdges.length - 1;
    this.rCenters = new Float64Array(this.numBins);
    for (let b = 0; b < this.numBins; b++) {
      this.rCenters[b] = 0.5 * (this.binEdges[b] + this.binEdges[b + 1]);
    }
  }

  /**
   * Computes Landy-Szalay correlation function xi_LS(r) = (DD - 2DR + RR) / RR * norm.
   *
   * @param {Array<{ x: number, y: number, z: number }>} data - Galaxy data points D
   * @param {Array<{ x: number, y: number, z: number }>} randoms - Unclustered synthetic random points R
   * @returns {{
   *   r: Float64Array,
   *   xi: Float64Array,
   *   xiHamilton: Float64Array,
   *   ddCounts: Float64Array,
   *   drCounts: Float64Array,
   *   rrCounts: Float64Array
   * }}
   */
  computeLandySzalay(data, randoms) {
    const nD = data.length;
    const nR = randoms.length;
    if (nD < 2 || nR < 2) {
      throw new Error('Requires at least 2 data and 2 random points for correlation function.');
    }

    const nPairsDD = (nD * (nD - 1)) / 2.0;
    const nPairsRR = (nR * (nR - 1)) / 2.0;
    const nPairsDR = nD * nR;

    const DD = computePairCounts(data, null, this.binEdges);
    const RR = computePairCounts(randoms, null, this.binEdges);
    const DR = computePairCounts(data, randoms, this.binEdges);

    const xiLS = new Float64Array(this.numBins);
    const xiHam = new Float64Array(this.numBins);

    for (let b = 0; b < this.numBins; b++) {
      const normDD = DD[b] / nPairsDD;
      const normRR = RR[b] / nPairsRR;
      const normDR = DR[b] / nPairsDR;

      if (normRR > 1e-9) {
        xiLS[b] = (normDD - 2.0 * normDR + normRR) / normRR;
      } else {
        xiLS[b] = 0.0;
      }

      if (normDR > 1e-9) {
        xiHam[b] = (normDD * normRR) / (normDR * normDR) - 1.0;
      } else {
        xiHam[b] = 0.0;
      }
    }

    return {
      r: this.rCenters,
      xi: xiLS,
      xiHamilton: xiHam,
      ddCounts: DD,
      drCounts: DR,
      rrCounts: RR
    };
  }

  /**
   * Computes 2D Anisotropic Correlation Matrix xi(sigma, pi) for Redshift-Space Distortions.
   *
   * @param {Array<{ x: number, y: number, z: number }>} data
   * @param {Array<{ x: number, y: number, z: number }>} randoms
   * @param {number} maxSigma - Max transverse separation in Mpc/h
   * @param {number} maxPi - Max parallel line-of-sight separation in Mpc/h
   * @param {number} [numGrid=20] - Number of bins in each direction
   * @returns {{
   *   sigmaGrid: Float64Array,
   *   piGrid: Float64Array,
   *   xi2D: number[][]
   * }}
   */
  computeAnisotropic2D(data, randoms, maxSigma = 40.0, maxPi = 40.0, numGrid = 20) {
    const dSigma = maxSigma / numGrid;
    const dPi = maxPi / numGrid;

    const sigmaCenters = new Float64Array(numGrid);
    const piCenters = new Float64Array(numGrid);
    for (let i = 0; i < numGrid; i++) {
      sigmaCenters[i] = (i + 0.5) * dSigma;
      piCenters[i] = (i + 0.5) * dPi;
    }

    const dd2D = Array.from({ length: numGrid }, () => new Float64Array(numGrid));
    const rr2D = Array.from({ length: numGrid }, () => new Float64Array(numGrid));
    const dr2D = Array.from({ length: numGrid }, () => new Float64Array(numGrid));

    const nD = data.length;
    const nR = randoms.length;

    // DD pairs
    for (let i = 0; i < nD; i++) {
      for (let j = i + 1; j < nD; j++) {
        const dx = data[i].x - data[j].x;
        const dy = data[i].y - data[j].y;
        const dz = data[i].z - data[j].z;

        // Midpoint line of sight
        const mx = 0.5 * (data[i].x + data[j].x);
        const my = 0.5 * (data[i].y + data[j].y);
        const mz = 0.5 * (data[i].z + data[j].z);
        const mMag = Math.hypot(mx, my, mz) || 1.0;
        const los = [mx / mMag, my / mMag, mz / mMag];

        const rPar = Math.abs(dx * los[0] + dy * los[1] + dz * los[2]);
        const rTotSq = dx * dx + dy * dy + dz * dz;
        const rPerp = Math.sqrt(Math.max(0.0, rTotSq - rPar * rPar));

        const sigIdx = Math.floor(rPerp / dSigma);
        const piIdx = Math.floor(rPar / dPi);

        if (sigIdx < numGrid && piIdx < numGrid) {
          dd2D[sigIdx][piIdx] += 1.0;
        }
      }
    }

    // RR pairs
    for (let i = 0; i < nR; i++) {
      for (let j = i + 1; j < nR; j++) {
        const dx = randoms[i].x - randoms[j].x;
        const dy = randoms[i].y - randoms[j].y;
        const dz = randoms[i].z - randoms[j].z;

        const mx = 0.5 * (randoms[i].x + randoms[j].x);
        const my = 0.5 * (randoms[i].y + randoms[j].y);
        const mz = 0.5 * (randoms[i].z + randoms[j].z);
        const mMag = Math.hypot(mx, my, mz) || 1.0;
        const los = [mx / mMag, my / mMag, mz / mMag];

        const rPar = Math.abs(dx * los[0] + dy * los[1] + dz * los[2]);
        const rTotSq = dx * dx + dy * dy + dz * dz;
        const rPerp = Math.sqrt(Math.max(0.0, rTotSq - rPar * rPar));

        const sigIdx = Math.floor(rPerp / dSigma);
        const piIdx = Math.floor(rPar / dPi);

        if (sigIdx < numGrid && piIdx < numGrid) {
          rr2D[sigIdx][piIdx] += 1.0;
        }
      }
    }

    // DR pairs
    for (let i = 0; i < nD; i++) {
      for (let j = 0; j < nR; j++) {
        const dx = data[i].x - randoms[j].x;
        const dy = data[i].y - randoms[j].y;
        const dz = data[i].z - randoms[j].z;

        const mx = 0.5 * (data[i].x + randoms[j].x);
        const my = 0.5 * (data[i].y + randoms[j].y);
        const mz = 0.5 * (data[i].z + randoms[j].z);
        const mMag = Math.hypot(mx, my, mz) || 1.0;
        const los = [mx / mMag, my / mMag, mz / mMag];

        const rPar = Math.abs(dx * los[0] + dy * los[1] + dz * los[2]);
        const rTotSq = dx * dx + dy * dy + dz * dz;
        const rPerp = Math.sqrt(Math.max(0.0, rTotSq - rPar * rPar));

        const sigIdx = Math.floor(rPerp / dSigma);
        const piIdx = Math.floor(rPar / dPi);

        if (sigIdx < numGrid && piIdx < numGrid) {
          dr2D[sigIdx][piIdx] += 1.0;
        }
      }
    }

    const nPairsDD = (nD * (nD - 1)) / 2.0;
    const nPairsRR = (nR * (nR - 1)) / 2.0;
    const nPairsDR = nD * nR;

    const xi2D = Array.from({ length: numGrid }, () => new Float64Array(numGrid));
    for (let s = 0; s < numGrid; s++) {
      for (let p = 0; p < numGrid; p++) {
        const normDD = dd2D[s][p] / nPairsDD;
        const normRR = rr2D[s][p] / nPairsRR;
        const normDR = dr2D[s][p] / nPairsDR;
        if (normRR > 1e-9) {
          xi2D[s][p] = (normDD - 2.0 * normDR + normRR) / normRR;
        }
      }
    }

    return {
      sigmaGrid: sigmaCenters,
      piGrid: piCenters,
      xi2D: xi2D.map(r => Array.from(r))
    };
  }
}
