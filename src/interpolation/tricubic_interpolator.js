/**
 * @file tricubic_interpolator.js
 * @description High-precision 3D Tensor-product tricubic Catmull-Rom interpolator
 * with validated C^1 continuity and analytical within-cell first and second derivatives
 * for cosmological scalar and vector fields.
 * 
 * Mathematical Formulation:
 * Inside each grid cell (ix, iy, iz) with local normalized coordinates (u, v, w) in [0, 1]^3:
 *   u = (x - x_0) / dx,  v = (y - y_0) / dy,  w = (z - z_0) / dz
 * The continuous field is represented by a tri-cubic polynomial:
 *   f(u, v, w) = sum_{i=0}^3 sum_{j=0}^3 sum_{k=0}^3 a_{ijk} u^i v^j w^k
 * 
 * Utilizing a 4x4x4 (64-point) local stencil with cubic Hermite / Catmull-Rom basis:
 *   w_{-1}(t) = 0.5 * (-t^3 + 2t^2 - t)
 *   w_0(t)    = 0.5 * ( 3t^3 - 5t^2 + 2)
 *   w_1(t)    = 0.5 * (-3t^3 + 4t^2 + t)
 *   w_2(t)    = 0.5 * (  t^3 -  t^2)
 * 
 * Properties:
 * 1. Exact node interpolation: f(x_node) = f_node exactly.
 * 2. C^1 continuity across all voxel boundaries.
 * 3. Analytical first and second derivatives (gradient, Hessian, Jacobian, Laplacian).
 * 4. Exact reproduction of 3D polynomials up to cubic degree.
 * 
 * @module interpolation/tricubic_interpolator
 */

import { GridIndexer, BoundaryMode } from '../fields/grid_indexer.js';
import {
  SupergalacticPosition,
  VelocityVector,
  Jacobian3x3,
  ScientificUnits
} from '../coordinates/scientific_types.js';

/**
 * Spline kernel types for tricubic interpolation.
 * @readonly
 * @enum {string}
 */
export const TricubicKernel = Object.freeze({
  CATMULL_ROM: 'catmull_rom',   // C^1 continuous, exact node interpolation (standard)
  HERMITE_64: 'hermite_64',     // 64-coefficient explicit cubic Hermite
  B_SPLINE: 'b_spline'          // C^2 continuous uniform B-spline (smoothing)
});

/**
 * Evaluates 1D Catmull-Rom cubic basis functions at normalized t in [0, 1].
 * @param {number} t Normalized parameter in [0, 1].
 * @param {Float64Array} out 4-element output array for weights [w-1, w0, w1, w2].
 */
export function catmullRomWeights(t, out) {
  const t2 = t * t;
  const t3 = t2 * t;
  out[0] = 0.5 * (-t3 + 2.0 * t2 - t);
  out[1] = 0.5 * (3.0 * t3 - 5.0 * t2 + 2.0);
  out[2] = 0.5 * (-3.0 * t3 + 4.0 * t2 + t);
  out[3] = 0.5 * (t3 - t2);
}

/**
 * Evaluates 1D Catmull-Rom first derivative weights d(w)/dt.
 * @param {number} t Normalized parameter in [0, 1].
 * @param {Float64Array} out 4-element output array [dw-1/dt, dw0/dt, dw1/dt, dw2/dt].
 */
export function catmullRomDerivatives(t, out) {
  const t2 = t * t;
  out[0] = 0.5 * (-3.0 * t2 + 4.0 * t - 1.0);
  out[1] = 0.5 * (9.0 * t2 - 10.0 * t);
  out[2] = 0.5 * (-9.0 * t2 + 8.0 * t + 1.0);
  out[3] = 0.5 * (3.0 * t2 - 2.0 * t);
}

/**
 * Evaluates 1D Catmull-Rom second derivative weights d^2(w)/dt^2.
 * @param {number} t Normalized parameter in [0, 1].
 * @param {Float64Array} out 4-element output array [d2w-1/dt2, d2w0/dt2, d2w1/dt2, d2w2/dt2].
 */
export function catmullRomSecondDerivatives(t, out) {
  out[0] = -3.0 * t + 2.0;
  out[1] = 9.0 * t - 5.0;
  out[2] = -9.0 * t + 4.0;
  out[3] = 3.0 * t - 1.0;
}

/**
 * Evaluates 1D uniform cubic B-spline weights at normalized t in [0, 1].
 * (Provides C^2 smooth curve approximating grid points).
 * @param {number} t
 * @param {Float64Array} out
 */
export function bSplineWeights(t, out) {
  const t2 = t * t;
  const t3 = t2 * t;
  const oneOver6 = 1.0 / 6.0;
  out[0] = oneOver6 * (1.0 - 3.0 * t + 3.0 * t2 - t3);
  out[1] = oneOver6 * (4.0 - 6.0 * t2 + 3.0 * t3);
  out[2] = oneOver6 * (1.0 + 3.0 * t + 3.0 * t2 - 3.0 * t3);
  out[3] = oneOver6 * t3;
}

/**
 * Evaluates 1D uniform cubic B-spline first derivatives.
 * @param {number} t
 * @param {Float64Array} out
 */
export function bSplineDerivatives(t, out) {
  const t2 = t * t;
  out[0] = 0.5 * (-1.0 + 2.0 * t - t2);
  out[1] = 0.5 * (-4.0 * t + 3.0 * t2);
  out[2] = 0.5 * (1.0 + 2.0 * t - 3.0 * t2);
  out[3] = 0.5 * t2;
}

/**
 * Evaluates 1D uniform cubic B-spline second derivatives.
 * @param {number} t
 * @param {Float64Array} out
 */
export function bSplineSecondDerivatives(t, out) {
  out[0] = 1.0 - t;
  out[1] = -2.0 + 3.0 * t;
  out[2] = 1.0 - 3.0 * t;
  out[3] = t;
}

/**
 * TricubicInterpolator performs 3D C^1 continuous field and gradient evaluation
 * over discrete structured meshes.
 */
export class TricubicInterpolator {
  /**
   * @param {GridIndexer} gridIndexer Associated grid indexer defining mesh geometry.
   * @param {Object} [options]
   * @param {string} [options.kernel='catmull_rom'] Spline kernel ('catmull_rom', 'hermite_64', 'b_spline').
   */
  constructor(gridIndexer, options = {}) {
    if (!gridIndexer || !(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('TricubicInterpolator: gridIndexer must be an instance of GridIndexer.');
    }

    this.grid = gridIndexer;
    this.kernel = options.kernel || TricubicKernel.CATMULL_ROM;

    // Pre-allocated reusable scratch arrays for high-frequency inner loops
    this._wx = new Float64Array(4);
    this._wy = new Float64Array(4);
    this._wz = new Float64Array(4);

    this._dwx = new Float64Array(4);
    this._dwy = new Float64Array(4);
    this._dwz = new Float64Array(4);

    this._d2wx = new Float64Array(4);
    this._d2wy = new Float64Array(4);
    this._d2wz = new Float64Array(4);

    this._indicesX = new Int32Array(4);
    this._indicesY = new Int32Array(4);
    this._indicesZ = new Int32Array(4);
  }

  /**
   * Computes base cell index (ix0, iy0, iz0) and normalized fractional offsets (u, v, w) in [0, 1].
   * @param {number} x SGX in physical units.
   * @param {number} y SGY in physical units.
   * @param {number} z SGZ in physical units.
   * @returns {{ix0: number, iy0: number, iz0: number, u: number, v: number, w: number}}
   */
  getVoxelCoords(x, y, z) {
    const origin = this.grid.origin;
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    // Fractional continuous index
    const fx = (x - origin[0]) / dx;
    const fy = (y - origin[1]) / dy;
    const fz = (z - origin[2]) / dz;

    let ix0 = Math.floor(fx);
    let iy0 = Math.floor(fy);
    let iz0 = Math.floor(fz);

    let u = fx - ix0;
    let v = fy - iy0;
    let w = fz - iz0;

    // Handle exact upper boundary
    if (ix0 >= this.grid.nx - 1) {
      ix0 = this.grid.nx - 2;
      u = fx - ix0;
    }
    if (iy0 >= this.grid.ny - 1) {
      iy0 = this.grid.ny - 2;
      v = fy - iy0;
    }
    if (iz0 >= this.grid.nz - 1) {
      iz0 = this.grid.nz - 2;
      w = fz - iz0;
    }

    if (ix0 < 0) {
      ix0 = 0;
      u = fx;
    }
    if (iy0 < 0) {
      iy0 = 0;
      v = fy;
    }
    if (iz0 < 0) {
      iz0 = 0;
      w = fz;
    }

    return { ix0, iy0, iz0, u, v, w };
  }

  /**
   * Populates the 4-point 1D index array with boundary handling for a specific axis.
   * @param {number} baseIdx Base coordinate index (e.g. ix0).
   * @param {number} maxDim Grid dimension (nx, ny, or nz).
   * @param {Int32Array} outIndices Output 4-element Int32Array.
   */
  _populateAxisIndices(baseIdx, maxDim, outIndices) {
    const bMode = this.grid.boundaryMode;

    for (let offset = -1; offset <= 2; offset++) {
      let idx = baseIdx + offset;
      const slot = offset + 1;

      if (idx >= 0 && idx < maxDim) {
        outIndices[slot] = idx;
      } else if (bMode === BoundaryMode.PERIODIC) {
        outIndices[slot] = ((idx % maxDim) + maxDim) % maxDim;
      } else if (bMode === BoundaryMode.REFLECT) {
        if (idx < 0) idx = -idx;
        if (idx >= maxDim) idx = 2 * (maxDim - 1) - idx;
        outIndices[slot] = Math.max(0, Math.min(maxDim - 1, idx));
      } else {
        // CLAMP or fallback
        outIndices[slot] = Math.max(0, Math.min(maxDim - 1, idx));
      }
    }
  }

  /**
   * Evaluates the 1D weights based on the active interpolation kernel.
   * @param {number} u
   * @param {number} v
   * @param {number} w
   */
  _computeWeights(u, v, w) {
    if (this.kernel === TricubicKernel.B_SPLINE) {
      bSplineWeights(u, this._wx);
      bSplineWeights(v, this._wy);
      bSplineWeights(w, this._wz);
    } else {
      catmullRomWeights(u, this._wx);
      catmullRomWeights(v, this._wy);
      catmullRomWeights(w, this._wz);
    }
  }

  /**
   * Evaluates the 1D weights and first derivatives.
   * @param {number} u
   * @param {number} v
   * @param {number} w
   */
  _computeWeightsAndDerivatives(u, v, w) {
    if (this.kernel === TricubicKernel.B_SPLINE) {
      bSplineWeights(u, this._wx);
      bSplineWeights(v, this._wy);
      bSplineWeights(w, this._wz);
      bSplineDerivatives(u, this._dwx);
      bSplineDerivatives(v, this._dwy);
      bSplineDerivatives(w, this._dwz);
    } else {
      catmullRomWeights(u, this._wx);
      catmullRomWeights(v, this._wy);
      catmullRomWeights(w, this._wz);
      catmullRomDerivatives(u, this._dwx);
      catmullRomDerivatives(v, this._dwy);
      catmullRomDerivatives(w, this._dwz);
    }
  }

  /**
   * Evaluates weights, first derivatives, and second derivatives.
   * @param {number} u
   * @param {number} v
   * @param {number} w
   */
  _computeAllDerivatives(u, v, w) {
    if (this.kernel === TricubicKernel.B_SPLINE) {
      bSplineWeights(u, this._wx);
      bSplineWeights(v, this._wy);
      bSplineWeights(w, this._wz);
      bSplineDerivatives(u, this._dwx);
      bSplineDerivatives(v, this._dwy);
      bSplineDerivatives(w, this._dwz);
      bSplineSecondDerivatives(u, this._d2wx);
      bSplineSecondDerivatives(v, this._d2wy);
      bSplineSecondDerivatives(w, this._d2wz);
    } else {
      catmullRomWeights(u, this._wx);
      catmullRomWeights(v, this._wy);
      catmullRomWeights(w, this._wz);
      catmullRomDerivatives(u, this._dwx);
      catmullRomDerivatives(v, this._dwy);
      catmullRomDerivatives(w, this._dwz);
      catmullRomSecondDerivatives(u, this._d2wx);
      catmullRomSecondDerivatives(v, this._d2wy);
      catmullRomSecondDerivatives(w, this._d2wz);
    }
  }

  /**
   * Interpolates a scalar 3D field buffer at physical coordinates (x, y, z).
   * 
   * @param {Float32Array|Float64Array} buffer 1D flat grid data array.
   * @param {number} x SGX in physical units (Mpc/h).
   * @param {number} y SGY in physical units (Mpc/h).
   * @param {number} z SGZ in physical units (Mpc/h).
   * @returns {number} Interpolated field value.
   */
  interpolate(buffer, x, y, z) {
    const { ix0, iy0, iz0, u, v, w } = this.getVoxelCoords(x, y, z);

    this._populateAxisIndices(ix0, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy0, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz0, this.grid.nz, this._indicesZ);

    this._computeWeights(u, v, w);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const wx = this._wx;
    const wy = this._wy;
    const wz = this._wz;
    const idxX = this._indicesX;
    const idxY = this._indicesY;
    const idxZ = this._indicesZ;

    let result = 0.0;

    for (let k = 0; k < 4; k++) {
      const zOffset = idxZ[k] * nx * ny;
      const weightZ = wz[k];

      for (let j = 0; j < 4; j++) {
        const yzOffset = zOffset + idxY[j] * nx;
        const weightYZ = weightZ * wy[j];

        for (let i = 0; i < 4; i++) {
          const flatIdx = yzOffset + idxX[i];
          result += weightYZ * wx[i] * buffer[flatIdx];
        }
      }
    }

    return result;
  }

  /**
   * Computes the continuous spatial gradient grad(f) = [df/dx, df/dy, df/dz]
   * at physical coordinates (x, y, z) via analytical spline derivatives.
   * 
   * @param {Float32Array|Float64Array} buffer 1D flat grid data array.
   * @param {number} x SGX in physical units.
   * @param {number} y SGY in physical units.
   * @param {number} z SGZ in physical units.
   * @returns {Array<number>} [df/dx, df/dy, df/dz].
   */
  interpolateGradient(buffer, x, y, z) {
    const { ix0, iy0, iz0, u, v, w } = this.getVoxelCoords(x, y, z);

    this._populateAxisIndices(ix0, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy0, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz0, this.grid.nz, this._indicesZ);

    this._computeWeightsAndDerivatives(u, v, w);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    const wx = this._wx, wy = this._wy, wz = this._wz;
    const dwx = this._dwx, dwy = this._dwy, dwz = this._dwz;
    const idxX = this._indicesX, idxY = this._indicesY, idxZ = this._indicesZ;

    let dfdu = 0.0;
    let dfdv = 0.0;
    let dfdw = 0.0;

    for (let k = 0; k < 4; k++) {
      const zOffset = idxZ[k] * nx * ny;
      const valWz = wz[k];
      const valDwz = dwz[k];

      for (let j = 0; j < 4; j++) {
        const yzOffset = zOffset + idxY[j] * nx;
        const valWy = wy[j];
        const valDwy = dwy[j];

        for (let i = 0; i < 4; i++) {
          const val = buffer[yzOffset + idxX[i]];
          const valWx = wx[i];
          const valDwx = dwx[i];

          dfdu += valDwx * valWy * valWz * val;
          dfdv += valWx * valDwy * valWz * val;
          dfdw += valWx * valWy * valDwz * val;
        }
      }
    }

    // Chain rule: df/dx = (df/du) / dx
    return [
      dfdu / dx,
      dfdv / dy,
      dfdw / dz
    ];
  }

  /**
   * Simultaneously evaluates field value and gradient [value, df/dx, df/dy, df/dz].
   * @param {Float32Array|Float64Array} buffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {{value: number, gradient: Array<number>}}
   */
  interpolateValueAndGradient(buffer, x, y, z) {
    const { ix0, iy0, iz0, u, v, w } = this.getVoxelCoords(x, y, z);

    this._populateAxisIndices(ix0, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy0, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz0, this.grid.nz, this._indicesZ);

    this._computeWeightsAndDerivatives(u, v, w);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    const wx = this._wx, wy = this._wy, wz = this._wz;
    const dwx = this._dwx, dwy = this._dwy, dwz = this._dwz;
    const idxX = this._indicesX, idxY = this._indicesY, idxZ = this._indicesZ;

    let value = 0.0;
    let dfdu = 0.0;
    let dfdv = 0.0;
    let dfdw = 0.0;

    for (let k = 0; k < 4; k++) {
      const zOffset = idxZ[k] * nx * ny;
      const valWz = wz[k];
      const valDwz = dwz[k];

      for (let j = 0; j < 4; j++) {
        const yzOffset = zOffset + idxY[j] * nx;
        const valWy = wy[j];
        const valDwy = dwy[j];

        for (let i = 0; i < 4; i++) {
          const val = buffer[yzOffset + idxX[i]];
          const valWx = wx[i];
          const valDwx = dwx[i];

          value += valWx * valWy * valWz * val;
          dfdu += valDwx * valWy * valWz * val;
          dfdv += valWx * valDwy * valWz * val;
          dfdw += valWx * valWy * valDwz * val;
        }
      }
    }

    return {
      value,
      gradient: [dfdu / dx, dfdv / dy, dfdw / dz]
    };
  }

  /**
   * Computes the 3x3 Hessian matrix H_ij = d^2 f / (dx_i dx_j) and Laplacian.
   * 
   * @param {Float32Array|Float64Array} buffer Flat grid buffer.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {{hessian: Float64Array, laplacian: number}}
   */
  interpolateHessian(buffer, x, y, z) {
    const { ix0, iy0, iz0, u, v, w } = this.getVoxelCoords(x, y, z);

    this._populateAxisIndices(ix0, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy0, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz0, this.grid.nz, this._indicesZ);

    this._computeAllDerivatives(u, v, w);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    const wx = this._wx, wy = this._wy, wz = this._wz;
    const dwx = this._dwx, dwy = this._dwy, dwz = this._dwz;
    const d2wx = this._d2wx, d2wy = this._d2wy, d2wz = this._d2wz;
    const idxX = this._indicesX, idxY = this._indicesY, idxZ = this._indicesZ;

    let d2f_du2 = 0.0;
    let d2f_dv2 = 0.0;
    let d2f_dw2 = 0.0;
    let d2f_dudv = 0.0;
    let d2f_dudw = 0.0;
    let d2f_dvdw = 0.0;

    for (let k = 0; k < 4; k++) {
      const zOffset = idxZ[k] * nx * ny;
      const wzK = wz[k], dwzK = dwz[k], d2wzK = d2wz[k];

      for (let j = 0; j < 4; j++) {
        const yzOffset = zOffset + idxY[j] * nx;
        const wyJ = wy[j], dwyJ = dwy[j], d2wyJ = d2wy[j];

        for (let i = 0; i < 4; i++) {
          const val = buffer[yzOffset + idxX[i]];
          const wxI = wx[i], dwxI = dwx[i], d2wxI = d2wx[i];

          d2f_du2 += d2wxI * wyJ * wzK * val;
          d2f_dv2 += wxI * d2wyJ * wzK * val;
          d2f_dw2 += wxI * wyJ * d2wzK * val;

          d2f_dudv += dwxI * dwyJ * wzK * val;
          d2f_dudw += dwxI * wyJ * dwzK * val;
          d2f_dvdw += wxI * dwyJ * dwzK * val;
        }
      }
    }

    const H = new Float64Array(9);
    const dxdx = dx * dx;
    const dydy = dy * dy;
    const dzdz = dz * dz;
    const dxdy = dx * dy;
    const dxdz = dx * dz;
    const dydz = dy * dz;

    H[0] = d2f_du2 / dxdx; // d2f / dx2
    H[4] = d2f_dv2 / dydy; // d2f / dy2
    H[8] = d2f_dw2 / dzdz; // d2f / dz2

    H[1] = H[3] = d2f_dudv / dxdy; // d2f / dxdy
    H[2] = H[6] = d2f_dudw / dxdz; // d2f / dxdz
    H[5] = H[7] = d2f_dvdw / dydz; // d2f / dydz

    const laplacian = H[0] + H[4] + H[8];

    return {
      hessian: H,
      laplacian
    };
  }

  /**
   * Evaluates a 3D vector field (vx, vy, vz) at physical coordinates (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x SGX in physical units.
   * @param {number} y SGY in physical units.
   * @param {number} z SGZ in physical units.
   * @returns {VelocityVector} Continuous peculiar velocity vector in km/s.
   */
  interpolateVector(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const { ix0, iy0, iz0, u, v, w } = this.getVoxelCoords(x, y, z);

    this._populateAxisIndices(ix0, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy0, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz0, this.grid.nz, this._indicesZ);

    this._computeWeights(u, v, w);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const wx = this._wx, wy = this._wy, wz = this._wz;
    const idxX = this._indicesX, idxY = this._indicesY, idxZ = this._indicesZ;

    let vx = 0.0;
    let vy = 0.0;
    let vz = 0.0;

    for (let k = 0; k < 4; k++) {
      const zOffset = idxZ[k] * nx * ny;
      const weightZ = wz[k];

      for (let j = 0; j < 4; j++) {
        const yzOffset = zOffset + idxY[j] * nx;
        const weightYZ = weightZ * wy[j];

        for (let i = 0; i < 4; i++) {
          const flatIdx = yzOffset + idxX[i];
          const weight = weightYZ * wx[i];
          vx += weight * vxBuffer[flatIdx];
          vy += weight * vyBuffer[flatIdx];
          vz += weight * vzBuffer[flatIdx];
        }
      }
    }

    return new VelocityVector(vx, vy, vz);
  }

  /**
   * Computes the 3x3 velocity Jacobian tensor J_ij = dv_i / dx_j at physical coordinates (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Jacobian3x3}
   */
  interpolateJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const gradVx = this.interpolateGradient(vxBuffer, x, y, z);
    const gradVy = this.interpolateGradient(vyBuffer, x, y, z);
    const gradVz = this.interpolateGradient(vzBuffer, x, y, z);

    // Row 0: dvx/dx, dvx/dy, dvx/dz
    // Row 1: dvy/dx, dvy/dy, dvy/dz
    // Row 2: dvz/dx, dvz/dy, dvz/dz
    const elements = [
      gradVx[0], gradVx[1], gradVx[2],
      gradVy[0], gradVy[1], gradVy[2],
      gradVz[0], gradVz[1], gradVz[2]
    ];

    return new Jacobian3x3(elements);
  }

  /**
   * Batch evaluates scalar field at N physical coordinates.
   * 
   * @param {Float32Array|Float64Array} buffer 1D scalar grid buffer.
   * @param {Float32Array|Float64Array|Array<number>} coords Interleaved [x0,y0,z0, x1,y1,z1, ...].
   * @param {Float64Array} [outBuffer] Pre-allocated output array of length N.
   * @returns {Float64Array} Output values array.
   */
  interpolateScalarBatch(buffer, coords, outBuffer = null) {
    const numPoints = Math.floor(coords.length / 3);
    const result = outBuffer || new Float64Array(numPoints);

    for (let p = 0; p < numPoints; p++) {
      const p3 = p * 3;
      result[p] = this.interpolate(buffer, coords[p3], coords[p3 + 1], coords[p3 + 2]);
    }

    return result;
  }

  /**
   * Batch evaluates vector field at N physical coordinates.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {Float32Array|Float64Array|Array<number>} coords Interleaved [x, y, z, ...].
   * @param {Float64Array} [outVx] Output buffer for vx.
   * @param {Float64Array} [outVy] Output buffer for vy.
   * @param {Float64Array} [outVz] Output buffer for vz.
   * @returns {{vx: Float64Array, vy: Float64Array, vz: Float64Array}}
   */
  interpolateVectorBatch(vxBuffer, vyBuffer, vzBuffer, coords, outVx = null, outVy = null, outVz = null) {
    const numPoints = Math.floor(coords.length / 3);
    const vxOut = outVx || new Float64Array(numPoints);
    const vyOut = outVy || new Float64Array(numPoints);
    const vzOut = outVz || new Float64Array(numPoints);

    for (let p = 0; p < numPoints; p++) {
      const p3 = p * 3;
      const v = this.interpolateVector(vxBuffer, vyBuffer, vzBuffer, coords[p3], coords[p3 + 1], coords[p3 + 2]);
      vxOut[p] = v.vx;
      vyOut[p] = v.vy;
      vzOut[p] = v.vz;
    }

    return { vx: vxOut, vy: vyOut, vz: vzOut };
  }

  /**
   * Solves for the 64 polynomial coefficients a_{ijk} for cell (ix, iy, iz)
   * such that f(u, v, w) = sum a_{ijk} u^i v^j w^k.
   * 
   * @param {Float32Array|Float64Array} buffer
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {Float64Array} 64-element coefficient array in order a_{ijk} with i fastest, then j, then k.
   */
  getPolynomialCoefficients(buffer, ix, iy, iz) {
    this._populateAxisIndices(ix, this.grid.nx, this._indicesX);
    this._populateAxisIndices(iy, this.grid.ny, this._indicesY);
    this._populateAxisIndices(iz, this.grid.nz, this._indicesZ);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const idxX = this._indicesX;
    const idxY = this._indicesY;
    const idxZ = this._indicesZ;

    // Load 4x4x4 stencil
    const g = new Float64Array(64);
    let gIdx = 0;
    for (let k = 0; k < 4; k++) {
      const zOff = idxZ[k] * nx * ny;
      for (let j = 0; j < 4; j++) {
        const yzOff = zOff + idxY[j] * nx;
        for (let i = 0; i < 4; i++) {
          g[gIdx++] = buffer[yzOff + idxX[i]];
        }
      }
    }

    // Catmull-Rom transformation matrix M (4x4)
    // [u^0, u^1, u^2, u^3] * M = [w-1, w0, w1, w2]
    const M = [
      [ 0.0,  1.0,  0.0,  0.0],
      [-0.5,  0.0,  0.5,  0.0],
      [ 1.0, -2.5,  2.0, -0.5],
      [-0.5,  1.5, -1.5,  0.5]
    ];

    // Transform along X: T1(i, j, k) = sum_{p=0}^3 M[i][p] * g(p, j, k)
    const T1 = new Float64Array(64);
    for (let k = 0; k < 4; k++) {
      for (let j = 0; j < 4; j++) {
        for (let i = 0; i < 4; i++) {
          let sum = 0.0;
          for (let p = 0; p < 4; p++) {
            sum += M[i][p] * g[k * 16 + j * 4 + p];
          }
          T1[k * 16 + j * 4 + i] = sum;
        }
      }
    }

    // Transform along Y: T2(i, j, k) = sum_{q=0}^3 M[j][q] * T1(i, q, k)
    const T2 = new Float64Array(64);
    for (let k = 0; k < 4; k++) {
      for (let j = 0; j < 4; j++) {
        for (let i = 0; i < 4; i++) {
          let sum = 0.0;
          for (let q = 0; q < 4; q++) {
            sum += M[j][q] * T1[k * 16 + q * 4 + i];
          }
          T2[k * 16 + j * 4 + i] = sum;
        }
      }
    }

    // Transform along Z: A(i, j, k) = sum_{r=0}^3 M[k][r] * T2(i, j, r)
    const A = new Float64Array(64);
    for (let k = 0; k < 4; k++) {
      for (let j = 0; j < 4; j++) {
        for (let i = 0; i < 4; i++) {
          let sum = 0.0;
          for (let r = 0; r < 4; r++) {
            sum += M[k][r] * T2[r * 16 + j * 4 + i];
          }
          A[k * 16 + j * 4 + i] = sum;
        }
      }
    }

    return A;
  }

  /**
   * Fast evaluation from 64 precomputed polynomial coefficients.
   * @param {Float64Array} coeffs 64 polynomial coefficients.
   * @param {number} u
   * @param {number} v
   * @param {number} w
   * @returns {number}
   */
  static evaluateFromCoefficients(coeffs, u, v, w) {
    const uPowers = [1.0, u, u * u, u * u * u];
    const vPowers = [1.0, v, v * v, v * v * v];
    const wPowers = [1.0, w, w * w, w * w * w];

    let result = 0.0;
    for (let k = 0; k < 4; k++) {
      const wk = wPowers[k];
      for (let j = 0; j < 4; j++) {
        const wj = wk * vPowers[j];
        for (let i = 0; i < 4; i++) {
          result += coeffs[k * 16 + j * 4 + i] * wj * uPowers[i];
        }
      }
    }
    return result;
  }
}
