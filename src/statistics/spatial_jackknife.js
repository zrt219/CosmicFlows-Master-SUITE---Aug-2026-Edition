/**
 * @file spatial_jackknife.js
 * @description General Spatial Jackknife Partitioning and Covariance Estimator for 3D Cosmological Surveys.
 *
 * Implements:
 * 1. 3D Cartesian block partitioning (Nx x Ny x Nz sub-cubes).
 * 2. Spatial Delete-1 Jackknife with automated bias estimation and bias-corrected parameter values.
 * 3. Jackknife pseudo-values J_k = K * theta_nom - (K - 1) * theta_{(k)}.
 * 4. Multi-parameter Jackknife error covariance matrix C_ij and correlation matrix.
 *
 * @module statistics/spatial_jackknife
 */

import { computeSampleCovariance } from '../uncertainty/covariance_estimator.js';

/**
 * Assigns a 3D position [x, y, z] to a 3D Cartesian grid partition of size [Nx, Ny, Nz].
 *
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number[]} origin - [x0, y0, z0]
 * @param {number[]} boxSize - [Lx, Ly, Lz]
 * @param {number[]} partitions - [Nx, Ny, Nz] (e.g. [2, 2, 2] = 8 blocks)
 * @returns {number} Sub-volume index k in [0, Nx*Ny*Nz - 1]
 */
export function getCartesianPartition(x, y, z, origin, boxSize, partitions) {
  const [x0, y0, z0] = origin;
  const [Lx, Ly, Lz] = boxSize;
  const [Nx, Ny, Nz] = partitions;

  let ix = Math.floor(((x - x0) / Lx) * Nx);
  let iy = Math.floor(((y - y0) / Ly) * Ny);
  let iz = Math.floor(((z - z0) / Lz) * Nz);

  ix = Math.max(0, Math.min(Nx - 1, ix));
  iy = Math.max(0, Math.min(Ny - 1, iy));
  iz = Math.max(0, Math.min(Nz - 1, iz));

  return (ix * Ny + iy) * Nz + iz;
}

/**
 * Spatial Jackknife Analysis Engine.
 */
export class SpatialJackknife {
  /**
   * @param {Array<{ x: number, y: number, z: number, [key: string]: any }>} data
   * @param {Object} [partitionOptions={}]
   * @param {number[]} [partitionOptions.origin=[-200, -200, -200]]
   * @param {number[]} [partitionOptions.boxSize=[400, 400, 400]]
   * @param {number[]} [partitionOptions.partitions=[2, 2, 2]] - Default 8 blocks
   */
  constructor(data, partitionOptions = {}) {
    this.data = data;
    this.n = data.length;

    this.origin = partitionOptions.origin ?? [-200, -200, -200];
    this.boxSize = partitionOptions.boxSize ?? [400, 400, 400];
    this.partitions = partitionOptions.partitions ?? [2, 2, 2];
    this.numRegions = this.partitions[0] * this.partitions[1] * this.partitions[2];

    // Partition data into spatial sub-buckets
    this.buckets = Array.from({ length: this.numRegions }, () => []);
    for (let i = 0; i < this.n; i++) {
      const item = data[i];
      const pIdx = getCartesianPartition(item.x, item.y, item.z, this.origin, this.boxSize, this.partitions);
      this.buckets[pIdx].push(item);
    }
  }

  /**
   * Executes Delete-1 Jackknife estimation for a vector-valued or scalar estimator.
   *
   * @param {function(Array): Float64Array|number[]} estimatorFn - Function mapping sub-sample to parameter vector
   * @returns {{
   *   nominal: number[],
   *   jackknifeMean: number[],
   *   bias: number[],
   *   biasCorrected: number[],
   *   covarianceMatrix: number[][],
   *   standardErrors: number[],
   *   pseudoValues: number[][],
   *   activeRegions: number
   * }}
   */
  evaluate(estimatorFn) {
    const nominal = Array.from(estimatorFn(this.data));
    const dim = nominal.length;

    const subEstimates = [];

    for (let k = 0; k < this.numRegions; k++) {
      // Create sample with region k excluded
      const subSample = [];
      for (let r = 0; r < this.numRegions; r++) {
        if (r === k) continue;
        subSample.push(...this.buckets[r]);
      }

      if (subSample.length >= Math.max(3, dim)) {
        const est = Array.from(estimatorFn(subSample));
        subEstimates.push(est);
      }
    }

    const K = subEstimates.length;
    if (K < 2) {
      throw new Error('Not enough non-empty spatial partitions for jackknife estimation.');
    }

    // Jackknife mean theta_bar
    const jackMean = new Float64Array(dim);
    for (let k = 0; k < K; k++) {
      for (let d = 0; d < dim; d++) {
        jackMean[d] += subEstimates[k][d];
      }
    }
    for (let d = 0; d < dim; d++) jackMean[d] /= K;

    // Bias = (K - 1) * (jackMean - nominal)
    const bias = new Float64Array(dim);
    const biasCorrected = new Float64Array(dim);
    for (let d = 0; d < dim; d++) {
      bias[d] = (K - 1) * (jackMean[d] - nominal[d]);
      biasCorrected[d] = nominal[d] - bias[d];
    }

    // Jackknife covariance: C = ((K - 1) / K) * sum_k (theta_k - jackMean)(theta_k - jackMean)^T
    const cov = Array.from({ length: dim }, () => new Float64Array(dim));
    const factor = (K - 1) / K;

    for (let i = 0; i < dim; i++) {
      for (let j = i; j < dim; j++) {
        let sum = 0.0;
        for (let k = 0; k < K; k++) {
          sum += (subEstimates[k][i] - jackMean[i]) * (subEstimates[k][j] - jackMean[j]);
        }
        const val = factor * sum;
        cov[i][j] = val;
        cov[j][i] = val;
      }
    }

    const standardErrors = new Float64Array(dim);
    for (let d = 0; d < dim; d++) {
      standardErrors[d] = Math.sqrt(Math.max(0.0, cov[d][d]));
    }

    // Pseudo-values: J_k = K * nominal - (K - 1) * theta_k
    const pseudoValues = [];
    for (let k = 0; k < K; k++) {
      const pv = new Float64Array(dim);
      for (let d = 0; d < dim; d++) {
        pv[d] = K * nominal[d] - (K - 1) * subEstimates[k][d];
      }
      pseudoValues.push(Array.from(pv));
    }

    return {
      nominal,
      jackknifeMean: Array.from(jackMean),
      bias: Array.from(bias),
      biasCorrected: Array.from(biasCorrected),
      covarianceMatrix: cov.map(r => Array.from(r)),
      standardErrors: Array.from(standardErrors),
      pseudoValues,
      activeRegions: K
    };
  }
}
