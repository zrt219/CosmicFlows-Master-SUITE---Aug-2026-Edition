/**
 * @file trilinear_interpolator.js
 * @description Continuous trilinear interpolation of 3D scalar and vector cosmological fields,
 * with exact analytic spatial derivatives (gradient vector, velocity Jacobian matrix),
 * and second derivatives (scalar Hessian matrix, vector Hessian tensor).
 * 
 * Mathematical Formulation:
 * Inside each grid voxel with normalized fractional coordinates (u, v, w) in [0, 1]^3:
 *   f(u, v, w) = a0 + a1*u + a2*v + a3*w + a4*u*v + a5*u*w + a6*v*w + a7*u*v*w
 * 
 * Exact First Spatial Derivatives:
 *   df/dx = (1/dx) * (a1 + a4*v + a5*w + a7*v*w)
 *   df/dy = (1/dy) * (a2 + a4*u + a6*w + a7*u*w)
 *   df/dz = (1/dz) * (a3 + a5*u + a6*v + a7*u*v)
 * 
 * Exact Second Spatial Derivatives (Hessian):
 *   d2f/dx2 = 0, d2f/dy2 = 0, d2f/dz2 = 0 (trilinear piecewise linear along axes)
 *   d2f/dxdy = (1 / (dx*dy)) * (a4 + a7*w)
 *   d2f/dxdz = (1 / (dx*dz)) * (a5 + a7*v)
 *   d2f/dydz = (1 / (dy*dz)) * (a6 + a7*u)
 * 
 * @module interpolation/trilinear_interpolator
 */

import { GridIndexer, BoundaryMode } from '../fields/grid_indexer.js';

/**
 * TrilinearInterpolator evaluates continuous scalar and vector fields,
 * their analytic gradients/Jacobians, and Hessian tensors on 3D structured grids.
 */
export class TrilinearInterpolator {
  /**
   * Constructs a TrilinearInterpolator.
   * 
   * @param {GridIndexer} gridIndexer Associated grid indexer.
   */
  constructor(gridIndexer) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('TrilinearInterpolator: gridIndexer must be an instance of GridIndexer.');
    }
    this.grid = gridIndexer;
  }

  /**
   * Computes the 8 polynomial coefficients a0..a7 from the 8 corner voxel values.
   * 
   * @param {number} c000 Value at (i,   j,   k)
   * @param {number} c100 Value at (i+1, j,   k)
   * @param {number} c010 Value at (i,   j+1, k)
   * @param {number} c110 Value at (i+1, j+1, k)
   * @param {number} c001 Value at (i,   j,   k+1)
   * @param {number} c101 Value at (i+1, j,   k+1)
   * @param {number} c011 Value at (i,   j+1, k+1)
   * @param {number} c111 Value at (i+1, j+1, k+1)
   * @returns {Float64Array} 8 polynomial coefficients [a0, a1, a2, a3, a4, a5, a6, a7].
   */
  static computePolynomialCoefficients(c000, c100, c010, c110, c001, c101, c011, c111) {
    const a0 = c000;
    const a1 = c100 - c000;
    const a2 = c010 - c000;
    const a3 = c001 - c000;
    const a4 = c110 - c100 - c010 + c000;
    const a5 = c101 - c100 - c001 + c000;
    const a6 = c011 - c010 - c001 + c000;
    const a7 = c111 - c110 - c101 - c011 + c100 + c010 + c001 - c000;

    return Float64Array.of(a0, a1, a2, a3, a4, a5, a6, a7);
  }

  /**
   * Internal helper to extract the 8 corner values of a scalar field buffer.
   * 
   * @private
   * @param {ArrayLike<number>} buffer 1D flat buffer of grid values.
   * @param {number} gx Grid index along SGX.
   * @param {number} gy Grid index along SGY.
   * @param {number} gz Grid index along SGZ.
   * @returns {{coeffs: Float64Array, u: number, v: number, w: number}} Cell interpolation context.
   */
  _getCellContext(buffer, gx, gy, gz) {
    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;

    // Floor indices for base corner
    let i0 = Math.floor(gx);
    let j0 = Math.floor(gy);
    let k0 = Math.floor(gz);

    // Compute fractional cell coordinates in [0, 1]
    let u = gx - i0;
    let v = gy - j0;
    let w = gz - k0;

    // Handle edge boundaries for vertex grid
    if (i0 >= nx - 1) {
      i0 = nx - 2;
      u = 1.0;
    } else if (i0 < 0) {
      i0 = 0;
      u = 0.0;
    }

    if (j0 >= ny - 1) {
      j0 = ny - 2;
      v = 1.0;
    } else if (j0 < 0) {
      j0 = 0;
      v = 0.0;
    }

    if (k0 >= nz - 1) {
      k0 = nz - 2;
      w = 1.0;
    } else if (k0 < 0) {
      k0 = 0;
      w = 0.0;
    }

    const getVal = (ix, iy, iz) => {
      const idx = this.grid.handleBoundaryIndex(ix, iy, iz);
      if (idx === null) return 0.0;
      return buffer[this.grid.getLinearIndex(idx[0], idx[1], idx[2])];
    };

    const c000 = getVal(i0,     j0,     k0);
    const c100 = getVal(i0 + 1, j0,     k0);
    const c010 = getVal(i0,     j0 + 1, k0);
    const c110 = getVal(i0 + 1, j0 + 1, k0);
    const c001 = getVal(i0,     j0,     k0 + 1);
    const c101 = getVal(i0 + 1, j0,     k0 + 1);
    const c011 = getVal(i0,     j0 + 1, k0 + 1);
    const c111 = getVal(i0 + 1, j0 + 1, k0 + 1);

    const coeffs = TrilinearInterpolator.computePolynomialCoefficients(
      c000, c100, c010, c110, c001, c101, c011, c111
    );

    return { coeffs, u, v, w };
  }

  /**
   * Interpolates a continuous scalar field value at physical coordinates (x, y, z).
   * 
   * @param {ArrayLike<number>} buffer 1D flat buffer of scalar grid values.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {number} Interpolated scalar value.
   */
  interpolateScalar(buffer, x, y, z) {
    const [gx, gy, gz] = this.grid.coordToGridIndex(x, y, z);
    const { coeffs, u, v, w } = this._getCellContext(buffer, gx, gy, gz);
    const [a0, a1, a2, a3, a4, a5, a6, a7] = coeffs;

    return a0 + a1 * u + a2 * v + a3 * w +
           a4 * u * v + a5 * u * w + a6 * v * w +
           a7 * u * v * w;
  }

  /**
   * Computes the exact analytic spatial gradient of a scalar field [df/dx, df/dy, df/dz] at (x, y, z).
   * 
   * @param {ArrayLike<number>} buffer 1D flat buffer of scalar grid values.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {[number, number, number]} Analytic gradient vector [df/dx, df/dy, df/dz].
   */
  analyticGradient(buffer, x, y, z) {
    const [gx, gy, gz] = this.grid.coordToGridIndex(x, y, z);
    const { coeffs, u, v, w } = this._getCellContext(buffer, gx, gy, gz);
    const [, a1, a2, a3, a4, a5, a6, a7] = coeffs;

    const df_du = a1 + a4 * v + a5 * w + a7 * v * w;
    const df_dv = a2 + a4 * u + a6 * w + a7 * u * w;
    const df_dw = a3 + a5 * u + a6 * v + a7 * u * v;

    const df_dx = df_du * this.grid.invDx;
    const df_dy = df_dv * this.grid.invDy;
    const df_dz = df_dw * this.grid.invDz;

    return [df_dx, df_dy, df_dz];
  }

  /**
   * Computes the exact analytic Hessian matrix (3x3 symmetric matrix of 2nd spatial derivatives)
   * of a scalar field at (x, y, z).
   * 
   * [ [d2f/dx2,   d2f/dxdy,  d2f/dxdz],
   *   [d2f/dydx,  d2f/dy2,   d2f/dydz],
   *   [d2f/dzdx,  d2f/dzdy,  d2f/dz2] ]
   * 
   * @param {ArrayLike<number>} buffer 1D flat buffer of scalar grid values.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {Array<Array<number>>} 3x3 Hessian matrix.
   */
  analyticHessian(buffer, x, y, z) {
    const [gx, gy, gz] = this.grid.coordToGridIndex(x, y, z);
    const { coeffs, u, v, w } = this._getCellContext(buffer, gx, gy, gz);
    const [, , , , a4, a5, a6, a7] = coeffs;

    const invDxDy = this.grid.invDx * this.grid.invDy;
    const invDxDz = this.grid.invDx * this.grid.invDz;
    const invDyDz = this.grid.invDy * this.grid.invDz;

    const d2f_dxdy = (a4 + a7 * w) * invDxDy;
    const d2f_dxdz = (a5 + a7 * v) * invDxDz;
    const d2f_dydz = (a6 + a7 * u) * invDyDz;

    return [
      [0.0, d2f_dxdy, d2f_dxdz],
      [d2f_dxdy, 0.0, d2f_dydz],
      [d2f_dxdz, d2f_dydz, 0.0]
    ];
  }

  /**
   * Interpolates continuous 3D vector field components [vx, vy, vz] at physical coordinates (x, y, z).
   * 
   * @param {ArrayLike<number>} vxBuffer Flat buffer for X component.
   * @param {ArrayLike<number>} vyBuffer Flat buffer for Y component.
   * @param {ArrayLike<number>} vzBuffer Flat buffer for Z component.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {[number, number, number]} Interpolated vector [vx, vy, vz].
   */
  interpolateVector(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const vx = this.interpolateScalar(vxBuffer, x, y, z);
    const vy = this.interpolateScalar(vyBuffer, x, y, z);
    const vz = this.interpolateScalar(vzBuffer, x, y, z);
    return [vx, vy, vz];
  }

  /**
   * Computes the exact analytic 3x3 Jacobian matrix J_ij = dv_i / dx_j for a vector field at (x, y, z).
   * 
   * J = [ [dvx/dx, dvx/dy, dvx/dz],
   *       [dvy/dx, dvy/dy, dvy/dz],
   *       [dvz/dx, dvz/dy, dvz/dz] ]
   * 
   * @param {ArrayLike<number>} vxBuffer Flat buffer for X component.
   * @param {ArrayLike<number>} vyBuffer Flat buffer for Y component.
   * @param {ArrayLike<number>} vzBuffer Flat buffer for Z component.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {Array<Array<number>>} 3x3 Jacobian matrix.
   */
  analyticVectorJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const gradVx = this.analyticGradient(vxBuffer, x, y, z);
    const gradVy = this.analyticGradient(vyBuffer, x, y, z);
    const gradVz = this.analyticGradient(vzBuffer, x, y, z);

    return [
      gradVx,
      gradVy,
      gradVz
    ];
  }

  /**
   * Computes the exact 3x3x3 Hessian tensor H_ijk = d^2 v_i / (dx_j dx_k) for a vector field at (x, y, z).
   * 
   * @param {ArrayLike<number>} vxBuffer Flat buffer for X component.
   * @param {ArrayLike<number>} vyBuffer Flat buffer for Y component.
   * @param {ArrayLike<number>} vzBuffer Flat buffer for Z component.
   * @param {number} x SGX physical coordinate.
   * @param {number} y SGY physical coordinate.
   * @param {number} z SGZ physical coordinate.
   * @returns {Array<Array<Array<number>>>} 3x3x3 Hessian tensor [H_x, H_y, H_z].
   */
  analyticVectorHessian(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    return [
      this.analyticHessian(vxBuffer, x, y, z),
      this.analyticHessian(vyBuffer, x, y, z),
      this.analyticHessian(vzBuffer, x, y, z)
    ];
  }
}
