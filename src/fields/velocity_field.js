/**
 * @file velocity_field.js
 * @description Cosmic peculiar velocity field representation (vx, vy, vz),
 * continuous 3D evaluation, x52.0 scale factor enforcement (Wiener Filter displacement
 * to km/s conversion), boundary gradient tensors, and kinematic strain/vorticity decomposition.
 * 
 * Astrophysical Context:
 * In Cosmicflows Wiener Filter / Constrained Realizations reconstructions,
 * dimensionless displacement fields Psi are converted into peculiar velocities v (km/s):
 *   v = H0 * f(Omega_m) * Psi ~= 52.0 * Psi (for H0=74.6, Omega_m=0.31 -> H0*f ~= 52.0 km/s / (h^-1 Mpc))
 * 
 * @module fields/velocity_field
 */

import { GridIndexer } from './grid_indexer.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { FiniteDifference } from '../interpolation/finite_difference.js';

/**
 * Standard Cosmicflows Wiener Filter displacement-to-velocity scale factor (km/s / unit).
 */
export const DEFAULT_VELOCITY_SCALE_FACTOR = 52.0;

/**
 * VelocityField represents a 3D cosmological peculiar velocity vector field.
 */
export class VelocityField {
  /**
   * Constructs a VelocityField.
   * 
   * @param {GridIndexer} gridIndexer Associated grid indexer.
   * @param {Float32Array|Float64Array} vxBuffer X-component peculiar velocity array.
   * @param {Float32Array|Float64Array} vyBuffer Y-component peculiar velocity array.
   * @param {Float32Array|Float64Array} vzBuffer Z-component peculiar velocity array.
   * @param {Object} [options] Configuration options.
   * @param {number} [options.scaleFactor=52.0] Multiplicative velocity scale factor.
   * @param {boolean} [options.enforceScaleFactor=true] Whether to apply the scale factor.
   */
  constructor(gridIndexer, vxBuffer, vyBuffer, vzBuffer, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('VelocityField: gridIndexer must be an instance of GridIndexer.');
    }

    const totalCells = gridIndexer.totalCells;
    if (!vxBuffer || vxBuffer.length !== totalCells ||
        !vyBuffer || vyBuffer.length !== totalCells ||
        !vzBuffer || vzBuffer.length !== totalCells) {
      throw new Error(`VelocityField: Buffer sizes must match grid total cells (${totalCells}).`);
    }

    this.grid = gridIndexer;
    this.scaleFactor = typeof options.scaleFactor === 'number' ? options.scaleFactor : DEFAULT_VELOCITY_SCALE_FACTOR;
    this.isScaleFactorEnforced = options.enforceScaleFactor !== false;

    // Allocate internal Float64Array storage and apply scale factor if enforced
    const scale = this.isScaleFactorEnforced ? this.scaleFactor : 1.0;
    this.vx = new Float64Array(totalCells);
    this.vy = new Float64Array(totalCells);
    this.vz = new Float64Array(totalCells);

    for (let i = 0; i < totalCells; i++) {
      this.vx[i] = vxBuffer[i] * scale;
      this.vy[i] = vyBuffer[i] * scale;
      this.vz[i] = vzBuffer[i] * scale;
    }

    this.interpolator = new TrilinearInterpolator(this.grid);
    this.fd = new FiniteDifference(this.grid, this.interpolator);
  }

  /**
   * Factory method to create VelocityField from analytic functions.
   * 
   * @param {GridIndexer} gridIndexer Grid indexer.
   * @param {Function} vectorFn Analytic function (x, y, z) => [vx, vy, vz].
   * @param {Object} [options] Options including scaleFactor.
   * @returns {VelocityField} Created VelocityField.
   */
  static fromAnalyticFunction(gridIndexer, vectorFn, options = {}) {
    const total = gridIndexer.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    for (let iz = 0; iz < gridIndexer.nz; iz++) {
      for (let iy = 0; iy < gridIndexer.ny; iy++) {
        for (let ix = 0; ix < gridIndexer.nx; ix++) {
          const idx = gridIndexer.getLinearIndex(ix, iy, iz);
          const [x, y, z] = gridIndexer.getNodeCoord(ix, iy, iz);
          const [v_x, v_y, v_z] = vectorFn(x, y, z);
          vx[idx] = v_x;
          vy[idx] = v_y;
          vz[idx] = v_z;
        }
      }
    }

    // Analytic function is already in true physical units, so scaleFactor can default to 1.0 or as given
    const opts = { scaleFactor: 1.0, enforceScaleFactor: false, ...options };
    return new VelocityField(gridIndexer, vx, vy, vz, opts);
  }

  /**
   * Factory method to create VelocityField from interleaved [vx0, vy0, vz0, vx1, ...] buffer.
   * 
   * @param {GridIndexer} gridIndexer Grid indexer.
   * @param {ArrayLike<number>} interleavedBuffer Interleaved buffer (length = 3 * totalCells).
   * @param {Object} [options] Options.
   * @returns {VelocityField} Created VelocityField.
   */
  static fromInterleavedArray(gridIndexer, interleavedBuffer, options = {}) {
    const total = gridIndexer.totalCells;
    if (interleavedBuffer.length !== total * 3) {
      throw new Error(`VelocityField: Interleaved buffer length must be ${total * 3}.`);
    }

    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      vx[i] = interleavedBuffer[i * 3];
      vy[i] = interleavedBuffer[i * 3 + 1];
      vz[i] = interleavedBuffer[i * 3 + 2];
    }

    return new VelocityField(gridIndexer, vx, vy, vz, options);
  }

  /**
   * Evaluates continuous peculiar velocity vector [vx, vy, vz] at physical coordinates (x, y, z).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {[number, number, number]} Interpolated velocity vector [vx, vy, vz] in km/s.
   */
  sampleVelocity(x, y, z) {
    return this.interpolator.interpolateVector(this.vx, this.vy, this.vz, x, y, z);
  }

  /**
   * Computes peculiar velocity magnitude |v| = sqrt(vx^2 + vy^2 + vz^2) at (x, y, z).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {number} Velocity speed |v| in km/s.
   */
  sampleSpeed(x, y, z) {
    const [vx, vy, vz] = this.sampleVelocity(x, y, z);
    return Math.sqrt(vx * vx + vy * vy + vz * vz);
  }

  /**
   * Evaluates the continuous exact velocity gradient Jacobian tensor J_ij = dv_i / dx_j at (x, y, z).
   * 
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @param {boolean} [useAnalytic=true] True for exact analytic piecewise derivative, false for 14-pt finite difference.
   * @returns {Array<Array<number>>} 3x3 Jacobian matrix.
   */
  sampleGradientTensor(x, y, z, useAnalytic = true) {
    if (useAnalytic) {
      return this.interpolator.analyticVectorJacobian(this.vx, this.vy, this.vz, x, y, z);
    }
    return this.fd.velocityJacobian(this.vx, this.vy, this.vz, x, y, z, '14point');
  }

  /**
   * Computes velocity divergence div(v) = Tr(J) at (x, y, z).
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {number} Velocity divergence (km/s / (h^-1 Mpc)).
   */
  sampleDivergence(x, y, z) {
    const J = this.sampleGradientTensor(x, y, z, true);
    return J[0][0] + J[1][1] + J[2][2];
  }

  /**
   * Computes vorticity curl(v) at (x, y, z).
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {[number, number, number]} Vorticity vector [wx, wy, wz].
   */
  sampleVorticity(x, y, z) {
    const J = this.sampleGradientTensor(x, y, z, true);
    return [
      J[2][1] - J[1][2],
      J[0][2] - J[2][0],
      J[1][0] - J[0][1]
    ];
  }

  /**
   * Computes the 3x3 symmetric rate-of-strain tensor S_ij = 0.5 * (J_ij + J_ji) at (x, y, z).
   * 
   * @param {number} x SGX coordinate.
   * @param {number} y SGY coordinate.
   * @param {number} z SGZ coordinate.
   * @returns {Array<Array<number>>} 3x3 Rate-of-strain tensor.
   */
  sampleRateOfStrain(x, y, z) {
    const J = this.sampleGradientTensor(x, y, z, true);
    const S = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        S[i][j] = 0.5 * (J[i][j] + J[j][i]);
      }
    }
    return S;
  }

  /**
   * Diagonalizes a 3x3 symmetric matrix (such as the Rate-of-Strain tensor S_ij)
   * using the exact cyclic Jacobi eigenvalue algorithm.
   * 
   * @param {Array<Array<number>>} matrix 3x3 symmetric matrix.
   * @param {number} [maxIter=50] Maximum iterations.
   * @returns {{eigenvalues: [number, number, number], eigenvectors: Array<Array<number>>}} Sorted eigenvalues (descending) & eigenvectors.
   */
  static diagonalizeSymmetric3x3(matrix, maxIter = 50) {
    // Copy input matrix
    const A = [
      [matrix[0][0], matrix[0][1], matrix[0][2]],
      [matrix[1][0], matrix[1][1], matrix[1][2]],
      [matrix[2][0], matrix[2][1], matrix[2][2]]
    ];

    // Identity matrix for eigenvectors
    const V = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];

    for (let iter = 0; iter < maxIter; iter++) {
      let maxOff = 0.0;
      let p = 0, q = 1;

      for (let i = 0; i < 3; i++) {
        for (let j = i + 1; j < 3; j++) {
          const absVal = Math.abs(A[i][j]);
          if (absVal > maxOff) {
            maxOff = absVal;
            p = i;
            q = j;
          }
        }
      }

      if (maxOff < 1e-15) break;

      const app = A[p][p];
      const aqq = A[q][q];
      const apq = A[p][q];

      const phi = 0.5 * Math.atan2(2.0 * apq, aqq - app);
      const c = Math.cos(phi);
      const s = Math.sin(phi);

      // Perform Jacobi rotation on A
      A[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
      A[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
      A[p][q] = 0.0;
      A[q][p] = 0.0;

      const r = 3 - p - q; // the third index
      const apr = A[p][r];
      const aqr = A[q][r];
      A[p][r] = c * apr - s * aqr;
      A[r][p] = A[p][r];
      A[q][r] = s * apr + c * aqr;
      A[r][q] = A[q][r];

      // Update eigenvectors
      for (let i = 0; i < 3; i++) {
        const vip = V[i][p];
        const viq = V[i][q];
        V[i][p] = c * vip - s * viq;
        V[i][q] = s * vip + c * viq;
      }
    }

    // Sort eigenvalues descending: lambda1 >= lambda2 >= lambda3
    const pairs = [
      { val: A[0][0], vec: [V[0][0], V[1][0], V[2][0]] },
      { val: A[1][1], vec: [V[0][1], V[1][1], V[2][1]] },
      { val: A[2][2], vec: [V[0][2], V[1][2], V[2][2]] }
    ];
    pairs.sort((a, b) => b.val - a.val);

    return {
      eigenvalues: [pairs[0].val, pairs[1].val, pairs[2].val],
      eigenvectors: [pairs[0].vec, pairs[1].vec, pairs[2].vec]
    };
  }

  /**
   * Computes the mean bulk flow vector V_bulk in a spherical volume centered at [cx, cy, cz] with radius R.
   * 
   * @param {Array<number>} center [cx, cy, cz] coordinates.
   * @param {number} radius Spherical radius R.
   * @param {number} [nSamples=20] Number of sampling steps per axis.
   * @returns {{bulkFlow: [number, number, number], speed: number, count: number}} Bulk flow velocity vector and magnitude.
   */
  computeBulkFlow(center, radius, nSamples = 20) {
    const [cx, cy, cz] = center;
    const r2 = radius * radius;
    let sumVx = 0.0;
    let sumVy = 0.0;
    let sumVz = 0.0;
    let count = 0;

    const step = (2.0 * radius) / nSamples;

    for (let ix = 0; ix <= nSamples; ix++) {
      const x = cx - radius + ix * step;
      for (let iy = 0; iy <= nSamples; iy++) {
        const y = cy - radius + iy * step;
        for (let iz = 0; iz <= nSamples; iz++) {
          const z = cz - radius + iz * step;
          const dist2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
          if (dist2 <= r2) {
            const [vx, vy, vz] = this.sampleVelocity(x, y, z);
            sumVx += vx;
            sumVy += vy;
            sumVz += vz;
            count++;
          }
        }
      }
    }

    if (count === 0) {
      const [vx, vy, vz] = this.sampleVelocity(cx, cy, cz);
      return { bulkFlow: [vx, vy, vz], speed: Math.sqrt(vx * vx + vy * vy + vz * vz), count: 1 };
    }

    const bVx = sumVx / count;
    const bVy = sumVy / count;
    const bVz = sumVz / count;
    const speed = Math.sqrt(bVx * bVx + bVy * bVy + bVz * bVz);

    return {
      bulkFlow: [bVx, bVy, bVz],
      speed,
      count
    };
  }
}
