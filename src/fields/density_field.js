/**
 * @file density_field.js
 * @description Cosmological density contrast delta(x) evaluation, logarithmic density A = ln(1 + delta),
 * exact Cloud-in-Cell (CIC) mass assignment from galaxy/tracer catalogues, and density gradients.
 * 
 * Astrophysical Formulations:
 * 1. Density Contrast:
 *    delta(x) = (rho(x) - <rho>) / <rho> >= -1.0
 * 
 * 2. Logarithmic Density:
 *    A(x) = ln(1 + delta(x)) (for lognormal cosmological fields, bounded with delta_floor > -1.0)
 * 
 * 3. Cloud-in-Cell (CIC) Assignment:
 *    For particle at r_p = (x_p, y_p, z_p) with mass w_p:
 *    Distributes weight across 8 surrounding grid vertices via trilinear partition:
 *      W(dx, dy, dz) = (1 - |dx|/hx) * (1 - |dy|/hy) * (1 - |dz|/hz)
 *    Preserves total mass: sum_cells(rho_i * V_cell) = sum_particles(w_p)
 * 
 * @module fields/density_field
 */

import { GridIndexer } from './grid_indexer.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { FiniteDifference } from '../interpolation/finite_difference.js';

/**
 * DensityField represents a continuous 3D cosmological density field.
 */
export class DensityField {
  /**
   * Constructs a DensityField.
   * 
   * @param {GridIndexer} gridIndexer Associated grid indexer.
   * @param {Float32Array|Float64Array} deltaBuffer 1D flat buffer of overdensity values delta(x).
   * @param {Object} [options] Configuration options.
   * @param {number} [options.deltaFloor=-0.999] Floor value for logarithmic density evaluation.
   */
  constructor(gridIndexer, deltaBuffer, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('DensityField: gridIndexer must be an instance of GridIndexer.');
    }

    const total = gridIndexer.totalCells;
    if (!deltaBuffer || deltaBuffer.length !== total) {
      throw new Error(`DensityField: deltaBuffer size (${deltaBuffer?.length}) must equal grid total cells (${total}).`);
    }

    this.grid = gridIndexer;
    this.delta = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      this.delta[i] = deltaBuffer[i];
    }

    this.deltaFloor = typeof options.deltaFloor === 'number' ? options.deltaFloor : -0.999;
    this.interpolator = new TrilinearInterpolator(this.grid);
    this.fd = new FiniteDifference(this.grid, this.interpolator);
  }

  /**
   * Factory method to create a DensityField from discrete galaxy/tracer particles using CIC mass assignment.
   * 
   * @param {GridIndexer} gridIndexer Target grid indexer.
   * @param {ArrayLike<number>|Array<[number, number, number]>} particles Array of 3D positions [x0, y0, z0, x1, ...] or [[x0,y0,z0], ...].
   * @param {ArrayLike<number>} [weights] Particle masses/weights. If omitted, all particles have weight 1.0.
   * @returns {DensityField} Resulting overdensity field.
   */
  static fromParticlesCIC(gridIndexer, particles, weights = null) {
    const totalCells = gridIndexer.totalCells;
    const massGrid = new Float64Array(totalCells);

    let numParticles = 0;
    const isFlatArray = typeof particles[0] === 'number';

    if (isFlatArray) {
      numParticles = Math.floor(particles.length / 3);
    } else {
      numParticles = particles.length;
    }

    let totalParticleMass = 0.0;

    for (let p = 0; p < numParticles; p++) {
      const px = isFlatArray ? particles[p * 3] : particles[p][0];
      const py = isFlatArray ? particles[p * 3 + 1] : particles[p][1];
      const pz = isFlatArray ? particles[p * 3 + 2] : particles[p][2];
      const pw = weights ? weights[p] : 1.0;

      totalParticleMass += pw;

      // Convert to grid index space
      const [gx, gy, gz] = gridIndexer.coordToGridIndex(px, py, pz);

      const i0 = Math.floor(gx);
      const j0 = Math.floor(gy);
      const k0 = Math.floor(gz);

      const tx = gx - i0;
      const ty = gy - j0;
      const tz = gz - k0;

      // Weights for 8 vertices
      const w000 = (1.0 - tx) * (1.0 - ty) * (1.0 - tz);
      const w100 = tx * (1.0 - ty) * (1.0 - tz);
      const w010 = (1.0 - tx) * ty * (1.0 - tz);
      const w110 = tx * ty * (1.0 - tz);
      const w001 = (1.0 - tx) * (1.0 - ty) * tz;
      const w101 = tx * (1.0 - ty) * tz;
      const w011 = (1.0 - tx) * ty * tz;
      const w111 = tx * ty * tz;

      const deposit = (ix, iy, iz, weight) => {
        const idx = gridIndexer.handleBoundaryIndex(ix, iy, iz);
        if (idx !== null) {
          const linearIdx = gridIndexer.getLinearIndex(idx[0], idx[1], idx[2]);
          massGrid[linearIdx] += pw * weight;
        }
      };

      deposit(i0,     j0,     k0,     w000);
      deposit(i0 + 1, j0,     k0,     w100);
      deposit(i0,     j0 + 1, k0,     w010);
      deposit(i0 + 1, j0 + 1, k0,     w110);
      deposit(i0,     j0,     k0 + 1, w001);
      deposit(i0 + 1, j0,     k0 + 1, w101);
      deposit(i0,     j0 + 1, k0 + 1, w011);
      deposit(i0 + 1, j0 + 1, k0 + 1, w111);
    }

    // Compute mean mass per cell
    let gridTotalMass = 0.0;
    for (let i = 0; i < totalCells; i++) {
      gridTotalMass += massGrid[i];
    }

    const meanMassPerCell = gridTotalMass / totalCells;
    const deltaBuffer = new Float64Array(totalCells);

    if (meanMassPerCell > 0.0) {
      for (let i = 0; i < totalCells; i++) {
        deltaBuffer[i] = (massGrid[i] - meanMassPerCell) / meanMassPerCell;
      }
    }

    return new DensityField(gridIndexer, deltaBuffer);
  }

  /**
   * Factory method to create DensityField from an analytic function f(x, y, z).
   * 
   * @param {GridIndexer} gridIndexer Grid indexer.
   * @param {Function} analyticFn Analytic scalar function (x, y, z) => delta.
   * @returns {DensityField} Density field.
   */
  static fromAnalyticFunction(gridIndexer, analyticFn) {
    const total = gridIndexer.totalCells;
    const delta = new Float64Array(total);

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const idx = gridIndexer.getLinearIndex(ix, iy, iz);
          const [x, y, z] = gridIndexer.getNodeCoord(ix, iy, iz);
          delta[idx] = analyticFn(x, y, z);
        }
      }
    }

    return new DensityField(gridIndexer, delta);
  }

  /**
   * Samples the continuous density contrast delta(x) at physical coordinates (x, y, z).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {number} Interpolated density contrast delta.
   */
  sampleDensity(x, y, z) {
    return this.interpolator.interpolateScalar(this.delta, x, y, z);
  }

  /**
   * Samples the logarithmic density A(x) = ln(1 + delta(x)) at physical coordinates (x, y, z).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @param {number} [floorOverride] Optional custom delta floor.
   * @returns {number} Logarithmic density ln(1 + delta).
   */
  sampleLogDensity(x, y, z, floorOverride = null) {
    const floor = typeof floorOverride === 'number' ? floorOverride : this.deltaFloor;
    const d = Math.max(floor, this.sampleDensity(x, y, z));
    return Math.log(1.0 + d);
  }

  /**
   * Computes the exact spatial gradient of the density contrast grad(delta) = [d delta/dx, d delta/dy, d delta/dz].
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {[number, number, number]} Density gradient vector.
   */
  sampleGradient(x, y, z) {
    return this.interpolator.analyticGradient(this.delta, x, y, z);
  }

  /**
   * Computes the spatial gradient of the logarithmic density grad(A) = grad(delta) / (1 + delta).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {[number, number, number]} Log-density gradient vector.
   */
  sampleLogGradient(x, y, z) {
    const [gx, gy, gz] = this.sampleGradient(x, y, z);
    const d = Math.max(this.deltaFloor, this.sampleDensity(x, y, z));
    const denom = 1.0 + d;
    return [gx / denom, gy / denom, gz / denom];
  }

  /**
   * Computes statistical moments of the density contrast field across all grid voxels.
   * 
   * @returns {{min: number, max: number, mean: number, variance: number, stdDev: number, skewness: number}} Summary metrics.
   */
  computeStatistics() {
    const total = this.delta.length;
    let min = Infinity;
    let max = -Infinity;
    let sum = 0.0;

    for (let i = 0; i < total; i++) {
      const v = this.delta[i];
      if (v < min) min = v;
      if (v > max) max = v;
      sum += v;
    }

    const mean = sum / total;
    let sumVar = 0.0;
    let sumSkew = 0.0;

    for (let i = 0; i < total; i++) {
      const diff = this.delta[i] - mean;
      sumVar += diff * diff;
      sumSkew += diff * diff * diff;
    }

    const variance = sumVar / total;
    const stdDev = Math.sqrt(variance);
    const skewness = stdDev > 0 ? (sumSkew / total) / Math.pow(stdDev, 3) : 0.0;

    return {
      min,
      max,
      mean,
      variance,
      stdDev,
      skewness
    };
  }
}
