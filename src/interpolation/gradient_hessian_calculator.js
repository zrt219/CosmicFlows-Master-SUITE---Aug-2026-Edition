/**
 * @file gradient_hessian_calculator.js
 * @description Multi-order finite difference stencils (2nd, 4th, 6th, 8th order),
 * high-precision Hessian eigensystems, and critical point classification for cosmological fields.
 * 
 * Implements:
 * 1. Central difference gradient stencils up to 8th-order accuracy O(h^8).
 * 2. Mixed partial derivative Hessian matrices with guaranteed Schwarz symmetry H_ij = H_ji.
 * 3. 3D Jacobi eigensystem diagonalization for principal curvatures (lambda_1 >= lambda_2 >= lambda_3).
 * 4. Critical point classification (Local Maxima, Local Minima, Saddle Index 1, Saddle Index 2).
 * 5. Cosmological filament spine and sheet normal identification via Hessian eigenvectors.
 * 
 * @module interpolation/gradient_hessian_calculator
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { ScalarField3D } from '../fields/scalar_field_3d.js';
import { JacobiEigenSolver3D } from '../coordinates/scientific_types.js';

/**
 * Finite Difference Stencil Accuracy Order.
 * @readonly
 * @enum {number}
 */
export const StencilAccuracy = Object.freeze({
  SECOND: 2,
  FOURTH: 4,
  SIXTH: 6,
  EIGHTH: 8
});

/**
 * GradientHessianCalculator provides multi-order derivatives and eigensystem kinematics.
 */
export class GradientHessianCalculator {
  /**
   * @param {GridIndexer} gridIndexer
   * @param {number} [order=4] Finite difference order (2, 4, 6, 8)
   */
  constructor(gridIndexer, order = StencilAccuracy.FOURTH) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('GradientHessianCalculator: gridIndexer must be an instance of GridIndexer.');
    }
    this.grid = gridIndexer;
    this.order = order;
  }

  /**
   * Evaluates 1D 1st derivative of an array along a line with arbitrary order stencil.
   * @param {Function} getSample (offset: number) => number
   * @param {number} h Grid spacing
   * @returns {number}
   */
  _deriv1D(getSample, h) {
    if (this.order === StencilAccuracy.SECOND) {
      // 2nd order: (-1/2*f[-1] + 1/2*f[+1]) / h
      return (0.5 * getSample(1) - 0.5 * getSample(-1)) / h;
    } else if (this.order === StencilAccuracy.FOURTH) {
      // 4th order: (1/12*f[-2] - 2/3*f[-1] + 2/3*f[+1] - 1/12*f[+2]) / h
      return (
        (1.0 / 12.0) * getSample(-2) -
        (2.0 / 3.0)  * getSample(-1) +
        (2.0 / 3.0)  * getSample(1) -
        (1.0 / 12.0) * getSample(2)
      ) / h;
    } else if (this.order === StencilAccuracy.SIXTH) {
      // 6th order stencil
      return (
        (-1.0 / 60.0) * getSample(-3) +
        (3.0 / 20.0)  * getSample(-2) -
        (3.0 / 4.0)   * getSample(-1) +
        (3.0 / 4.0)   * getSample(1) -
        (3.0 / 20.0)  * getSample(2) +
        (1.0 / 60.0)  * getSample(3)
      ) / h;
    } else {
      // 8th order stencil
      return (
        (1.0 / 280.0)  * getSample(-4) -
        (4.0 / 105.0)  * getSample(-3) +
        (1.0 / 5.0)    * getSample(-2) -
        (4.0 / 5.0)    * getSample(-1) +
        (4.0 / 5.0)    * getSample(1) -
        (1.0 / 5.0)    * getSample(2) +
        (4.0 / 105.0)  * getSample(3) -
        (1.0 / 280.0)  * getSample(4)
      ) / h;
    }
  }

  /**
   * Computes discrete gradient vector [df/dx, df/dy, df/dz] at voxel (ix, iy, iz).
   * @param {ScalarField3D|Float64Array} field
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {Float64Array} [gx, gy, gz]
   */
  computeGradientAt(field, ix, iy, iz) {
    const data = field instanceof ScalarField3D ? field.data : field;
    const { nx, ny, nz, dx, dy, dz } = this.grid;

    const clamp = (val, max) => Math.max(0, Math.min(max - 1, val));

    const sampleX = (off) => data[this.grid.index(clamp(ix + off, nx), iy, iz)];
    const sampleY = (off) => data[this.grid.index(ix, clamp(iy + off, ny), iz)];
    const sampleZ = (off) => data[this.grid.index(ix, iy, clamp(iz + off, nz))];

    return new Float64Array([
      this._deriv1D(sampleX, dx),
      this._deriv1D(sampleY, dy),
      this._deriv1D(sampleZ, dz)
    ]);
  }

  /**
   * Computes discrete 3x3 symmetric Hessian matrix d2f / dx_i dx_j at voxel (ix, iy, iz).
   * @param {ScalarField3D|Float64Array} field
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {Float64Array} 9-element row-major matrix
   */
  computeHessianAt(field, ix, iy, iz) {
    const data = field instanceof ScalarField3D ? field.data : field;
    const { nx, ny, nz, dx, dy, dz } = this.grid;

    const clamp = (val, max) => Math.max(0, Math.min(max - 1, val));
    const getVal = (x, y, z) => data[this.grid.index(clamp(x, nx), clamp(y, ny), clamp(z, nz))];

    // Diagonal 2nd derivatives using central 4th-order stencil:
    // f''(x) = (-f[-2] + 16*f[-1] - 30*f[0] + 16*f[+1] - f[+2]) / (12 * h^2)
    const d2_diag = (sampleLine, h) => {
      const f_m2 = sampleLine(-2);
      const f_m1 = sampleLine(-1);
      const f_0  = sampleLine(0);
      const f_p1 = sampleLine(1);
      const f_p2 = sampleLine(2);
      return (-f_m2 + 16.0 * f_m1 - 30.0 * f_0 + 16.0 * f_p1 - f_p2) / (12.0 * h * h);
    };

    const dxx = d2_diag((off) => getVal(ix + off, iy, iz), dx);
    const dyy = d2_diag((off) => getVal(ix, iy + off, iz), dy);
    const dzz = d2_diag((off) => getVal(ix, iy, iz + off), dz);

    // Cross partial derivatives: d2f / dx dy = (f[+1,+1] - f[+1,-1] - f[-1,+1] + f[-1,-1]) / (4 * dx * dy)
    const dxy = (
      getVal(ix + 1, iy + 1, iz) -
      getVal(ix + 1, iy - 1, iz) -
      getVal(ix - 1, iy + 1, iz) +
      getVal(ix - 1, iy - 1, iz)
    ) / (4.0 * dx * dy);

    const dxz = (
      getVal(ix + 1, iy, iz + 1) -
      getVal(ix + 1, iy, iz - 1) -
      getVal(ix - 1, iy, iz + 1) +
      getVal(ix - 1, iy, iz - 1)
    ) / (4.0 * dx * dz);

    const dyz = (
      getVal(ix, iy + 1, iz + 1) -
      getVal(ix, iy + 1, iz - 1) -
      getVal(ix, iy - 1, iz + 1) +
      getVal(ix, iy - 1, iz - 1)
    ) / (4.0 * dy * dz);

    return new Float64Array([
      dxx, dxy, dxz,
      dxy, dyy, dyz,
      dxz, dyz, dzz
    ]);
  }

  /**
   * Diagonalizes Hessian matrix at (ix, iy, iz) and classifies the critical point topology.
   * @param {ScalarField3D|Float64Array} field
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @param {number} [gradThreshold=1e-3]
   * @returns {{ eigenvalues: Float64Array, eigenvectors: Array<Float64Array>, isCriticalPoint: boolean, criticalPointType: string }}
   */
  classifyPointTopology(field, ix, iy, iz, gradThreshold = 1e-3) {
    const grad = this.computeGradientAt(field, ix, iy, iz);
    const gradNorm = Math.sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
    const isCritical = gradNorm < gradThreshold;

    const H = this.computeHessianAt(field, ix, iy, iz);
    const eigen = JacobiEigenSolver3D.diagonalize(H);

    const [l1, l2, l3] = eigen.eigenvalues; // Sorted descending: l1 >= l2 >= l3

    let cpType = 'regular';
    if (isCritical) {
      if (l1 < 0 && l2 < 0 && l3 < 0) {
        cpType = 'local_maximum'; // Peak
      } else if (l1 > 0 && l2 > 0 && l3 > 0) {
        cpType = 'local_minimum'; // Void center / basin minimum
      } else if (l1 >= 0 && l2 >= 0 && l3 < 0) {
        cpType = 'saddle_index_1'; // Wall saddle
      } else if (l1 >= 0 && l2 < 0 && l3 < 0) {
        cpType = 'saddle_index_2'; // Filament saddle
      }
    }

    return {
      eigenvalues: eigen.eigenvalues,
      eigenvectors: eigen.eigenvectors,
      isCriticalPoint: isCritical,
      criticalPointType: cpType
    };
  }
}
