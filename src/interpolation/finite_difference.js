/**
 * @file finite_difference.js
 * @description 6-point and 14-point isotropic central finite differencing
 * for cosmological 3D scalar and vector velocity fields, computing exact gradients,
 * velocity Jacobian tensors, divergence, vorticity, strain rate, and convergence checks.
 * 
 * Mathematical Formulation:
 * 
 * 1. Standard 6-Point Central Stencil (2nd order O(h^2)):
 *    df/dx = (f(x + hx, y, z) - f(x - hx, y, z)) / (2 * hx)
 *    df/dy = (f(x, y + hy, z) - f(x, y - hy, z)) / (2 * hy)
 *    df/dz = (f(x, y, z + hz) - f(x, y, z - hz)) / (2 * hz)
 * 
 * 2. 14-Point Isotropic Central Stencil (Eliminates lattice anisotropy up to O(h^4)):
 *    df/dx = (1 / (2 * hx)) * [ w_face * (f(+hx, 0, 0) - f(-hx, 0, 0)) +
 *            w_corner * sum_{sy, sz = +-1} (f(+hx, sy*hy, sz*hz) - f(-hx, sy*hy, sz*hz)) ]
 *    where w_face = 2/3, w_corner = 1/12, normalized so w_face + 4*w_corner = 1.
 * 
 * 3. Kinematic Quantities:
 *    - Jacobian: J_ij = dv_i / dx_j
 *    - Divergence: div(v) = Tr(J) = J_00 + J_11 + J_22
 *    - Vorticity: omega = curl(v) = [J_21 - J_12, J_02 - J_20, J_10 - J_01]
 *    - Rate of Strain Tensor: S_ij = 0.5 * (J_ij + J_ji)
 *    - Spin Tensor: Omega_ij = 0.5 * (J_ij - J_ji)
 * 
 * @module interpolation/finite_difference
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { TrilinearInterpolator } from './trilinear_interpolator.js';

/**
 * FiniteDifference provides 6-point and 14-point difference operators,
 * velocity field decomposition, and convergence verification.
 */
export class FiniteDifference {
  /**
   * Constructs a FiniteDifference operator.
   * 
   * @param {GridIndexer} gridIndexer Associated grid indexer.
   * @param {TrilinearInterpolator} [interpolator] Optional shared interpolator.
   */
  constructor(gridIndexer, interpolator = null) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('FiniteDifference: gridIndexer must be an instance of GridIndexer.');
    }
    this.grid = gridIndexer;
    this.interpolator = interpolator || new TrilinearInterpolator(gridIndexer);
  }

  /**
   * Computes the 3D scalar gradient using the standard 6-point central difference stencil.
   * 
   * @param {ArrayLike<number>|Function} field Scalar buffer or evaluation function f(x, y, z).
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {Array<number>} [stepSize] Custom step sizes [hx, hy, hz]. Defaults to grid [dx, dy, dz].
   * @returns {[number, number, number]} Gradient vector [df/dx, df/dy, df/dz].
   */
  gradient6Point(field, x, y, z, stepSize = null) {
    const hx = stepSize ? stepSize[0] : this.grid.dx;
    const hy = stepSize ? stepSize[1] : this.grid.dy;
    const hz = stepSize ? stepSize[2] : this.grid.dz;

    const evalFn = typeof field === 'function'
      ? field
      : (px, py, pz) => this.interpolator.interpolateScalar(field, px, py, pz);

    const df_dx = (evalFn(x + hx, y, z) - evalFn(x - hx, y, z)) / (2.0 * hx);
    const df_dy = (evalFn(x, y + hy, z) - evalFn(x, y - hy, z)) / (2.0 * hy);
    const df_dz = (evalFn(x, y, z + hz) - evalFn(x, y, z - hz)) / (2.0 * hz);

    return [df_dx, df_dy, df_dz];
  }

  /**
   * Computes the 3D scalar gradient using the 14-point isotropic central difference stencil.
   * 
   * @param {ArrayLike<number>|Function} field Scalar buffer or evaluation function f(x, y, z).
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {Array<number>} [stepSize] Custom step sizes [hx, hy, hz]. Defaults to grid [dx, dy, dz].
   * @returns {[number, number, number]} Gradient vector [df/dx, df/dy, df/dz].
   */
  gradient14Point(field, x, y, z, stepSize = null) {
    const hx = stepSize ? stepSize[0] : this.grid.dx;
    const hy = stepSize ? stepSize[1] : this.grid.dy;
    const hz = stepSize ? stepSize[2] : this.grid.dz;

    const evalFn = typeof field === 'function'
      ? field
      : (px, py, pz) => this.interpolator.interpolateScalar(field, px, py, pz);

    const wFace = 2.0 / 3.0;
    const wCorner = 1.0 / 12.0;

    // Gradient df/dx
    const faceX = evalFn(x + hx, y, z) - evalFn(x - hx, y, z);
    const cornerX =
      (evalFn(x + hx, y + hy, z + hz) - evalFn(x - hx, y + hy, z + hz)) +
      (evalFn(x + hx, y + hy, z - hz) - evalFn(x - hx, y + hy, z - hz)) +
      (evalFn(x + hx, y - hy, z + hz) - evalFn(x - hx, y - hy, z + hz)) +
      (evalFn(x + hx, y - hy, z - hz) - evalFn(x - hx, y - hy, z - hz));
    const df_dx = (wFace * faceX + wCorner * cornerX) / (2.0 * hx);

    // Gradient df/dy
    const faceY = evalFn(x, y + hy, z) - evalFn(x, y - hy, z);
    const cornerY =
      (evalFn(x + hx, y + hy, z + hz) - evalFn(x + hx, y - hy, z + hz)) +
      (evalFn(x + hx, y + hy, z - hz) - evalFn(x + hx, y - hy, z - hz)) +
      (evalFn(x - hx, y + hy, z + hz) - evalFn(x - hx, y - hy, z + hz)) +
      (evalFn(x - hx, y + hy, z - hz) - evalFn(x - hx, y - hy, z - hz));
    const df_dy = (wFace * faceY + wCorner * cornerY) / (2.0 * hy);

    // Gradient df/dz
    const faceZ = evalFn(x, y, z + hz) - evalFn(x, y, z - hz);
    const cornerZ =
      (evalFn(x + hx, y + hy, z + hz) - evalFn(x + hx, y + hy, z - hz)) +
      (evalFn(x + hx, y - hy, z + hz) - evalFn(x + hx, y - hy, z - hz)) +
      (evalFn(x - hx, y + hy, z + hz) - evalFn(x - hx, y + hy, z - hz)) +
      (evalFn(x - hx, y - hy, z + hz) - evalFn(x - hx, y - hy, z - hz));
    const df_dz = (wFace * faceZ + wCorner * cornerZ) / (2.0 * hz);

    return [df_dx, df_dy, df_dz];
  }

  /**
   * Computes the 3x3 velocity gradient Jacobian matrix J_ij = dv_i / dx_j.
   * 
   * @param {ArrayLike<number>} vxBuffer Flat X velocity buffer.
   * @param {ArrayLike<number>} vyBuffer Flat Y velocity buffer.
   * @param {ArrayLike<number>} vzBuffer Flat Z velocity buffer.
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {string} [stencil='14point'] '6point' or '14point'.
   * @param {Array<number>} [stepSize] Custom step size.
   * @returns {Array<Array<number>>} 3x3 Jacobian tensor.
   */
  velocityJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil = '14point', stepSize = null) {
    const gradMethod = stencil === '6point'
      ? (buf) => this.gradient6Point(buf, x, y, z, stepSize)
      : (buf) => this.gradient14Point(buf, x, y, z, stepSize);

    const gradVx = gradMethod(vxBuffer);
    const gradVy = gradMethod(vyBuffer);
    const gradVz = gradMethod(vzBuffer);

    return [
      gradVx,
      gradVy,
      gradVz
    ];
  }

  /**
   * Computes velocity divergence div(v) = dvx/dx + dvy/dy + dvz/dz.
   * 
   * @param {ArrayLike<number>} vxBuffer Flat X velocity buffer.
   * @param {ArrayLike<number>} vyBuffer Flat Y velocity buffer.
   * @param {ArrayLike<number>} vzBuffer Flat Z velocity buffer.
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {string} [stencil='14point'] Stencil type ('6point' or '14point').
   * @returns {number} Velocity divergence.
   */
  divergence(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil = '14point') {
    const J = this.velocityJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil);
    return J[0][0] + J[1][1] + J[2][2];
  }

  /**
   * Computes vorticity vector omega = curl(v) = [dvz/dy - dvy/dz, dvx/dz - dvz/dx, dvy/dx - dvx/dy].
   * 
   * @param {ArrayLike<number>} vxBuffer Flat X velocity buffer.
   * @param {ArrayLike<number>} vyBuffer Flat Y velocity buffer.
   * @param {ArrayLike<number>} vzBuffer Flat Z velocity buffer.
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {string} [stencil='14point'] Stencil type ('6point' or '14point').
   * @returns {[number, number, number]} Vorticity vector [wx, wy, wz].
   */
  vorticity(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil = '14point') {
    const J = this.velocityJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil);
    const wx = J[2][1] - J[1][2];
    const wy = J[0][2] - J[2][0];
    const wz = J[1][0] - J[0][1];
    return [wx, wy, wz];
  }

  /**
   * Computes the symmetric rate-of-strain tensor S_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i).
   * 
   * @param {ArrayLike<number>} vxBuffer Flat X velocity buffer.
   * @param {ArrayLike<number>} vyBuffer Flat Y velocity buffer.
   * @param {ArrayLike<number>} vzBuffer Flat Z velocity buffer.
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {string} [stencil='14point'] Stencil type.
   * @returns {Array<Array<number>>} 3x3 symmetric strain rate tensor.
   */
  rateOfStrain(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil = '14point') {
    const J = this.velocityJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z, stencil);
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
   * Performs automated Richardson extrapolation and convergence error check
   * for a scalar derivative function at (x, y, z).
   * 
   * @param {Function} evalFn Evaluation function f(x, y, z).
   * @param {number} x Physical coordinate X.
   * @param {number} y Physical coordinate Y.
   * @param {number} z Physical coordinate Z.
   * @param {number} [baseH=1.0] Base step size h.
   * @param {string} [stencil='14point'] '6point' or '14point'.
   * @returns {Object} Convergence telemetry { d_h, d_half_h, extrapolated, errorEstimate, orderObserved }.
   */
  checkConvergence(evalFn, x, y, z, baseH = 1.0, stencil = '14point') {
    const gradFn = stencil === '6point'
      ? (h) => this.gradient6Point(evalFn, x, y, z, [h, h, h])
      : (h) => this.gradient14Point(evalFn, x, y, z, [h, h, h]);

    const g1 = gradFn(baseH);
    const g2 = gradFn(baseH * 0.5);
    const g3 = gradFn(baseH * 0.25);

    // Richardson extrapolation for 2nd order (p=2) -> extrapolated = (4*g2 - g1)/3
    const extrapolatedX = (4.0 * g2[0] - g1[0]) / 3.0;
    const errorEstimateX = Math.abs(g2[0] - g1[0]) / 3.0;

    // Calculate observed convergence rate p = log2( |g1 - g2| / |g2 - g3| )
    const diff1 = Math.abs(g1[0] - g2[0]);
    const diff2 = Math.abs(g2[0] - g3[0]);
    const orderObserved = (diff1 > 1e-14 && diff2 > 1e-14)
      ? Math.log2(diff1 / diff2)
      : 2.0;

    return {
      grad_h: g1,
      grad_half_h: g2,
      grad_quarter_h: g3,
      extrapolated: [
        (4.0 * g2[0] - g1[0]) / 3.0,
        (4.0 * g2[1] - g1[1]) / 3.0,
        (4.0 * g2[2] - g1[2]) / 3.0
      ],
      errorEstimate: [
        Math.abs(g2[0] - g1[0]) / 3.0,
        Math.abs(g2[1] - g1[1]) / 3.0,
        Math.abs(g2[2] - g1[2]) / 3.0
      ],
      orderObserved
    };
  }
}
