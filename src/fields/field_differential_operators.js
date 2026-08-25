/**
 * @file field_differential_operators.js
 * @description Comprehensive cosmological differential operator library and kinetic sanity
 * diagnostics suite for the ZRT CosmicFlows Computational Research Workbench.
 * 
 * Computes:
 * - Continuous and discrete gradient: grad(f) = [df/dx, df/dy, df/dz]
 * - Velocity divergence: div(v) = dvx/dx + dvy/dy + dvz/dz (theta)
 * - Vorticity vector: curl(v) = [dvz/dy - dvy/dz, dvx/dz - dvz/dx, dvy/dx - dvx/dy] (omega)
 * - Scalar and vector Laplacian: div(grad(f)) = d2f/dx2 + d2f/dy2 + d2f/dz2
 * - Symmetric strain rate tensor: S_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i)
 * - Antisymmetric vorticity tensor: Omega_ij = 0.5 * (dv_i/dx_j - dv_j/dx_i)
 * - Traceless kinematic shear tensor: sigma_ij = S_ij - (1/3) * theta * delta_ij
 * - Principal strain eigenvalues (lambda_1 >= lambda_2 >= lambda_3) and eigenvectors
 * - Okubo-Weiss parameter: Q = 0.25 * (||omega||^2 - 2 * ||S||^2)
 * - Kinetic energy density: E_k = 0.5 * rho_bar * (1 + delta) * ||v||^2
 * - Enstrophy density: E = 0.5 * ||omega||^2
 * - Kinetic helicity density: h = v . omega
 * - Cosmic web classification (voids, sheets, filaments, knots) via deformation eigenvalue counting
 * - Full diagnostic kinetic sanity maps and linear continuity audit reports
 * 
 * @module fields/field_differential_operators
 */

import { GridIndexer, BoundaryMode } from './grid_indexer.js';
import { TricubicInterpolator, TricubicKernel } from '../interpolation/tricubic_interpolator.js';
import {
  SupergalacticPosition,
  VelocityVector,
  DensityContrast,
  Jacobian3x3,
  EigenSystem3D,
  RunManifest,
  ScientificUnits,
  CF4_PUBLIC_VELOCITY_SCALE,
  CosmologicalParameters,
  sha256Hex
} from '../coordinates/scientific_types.js';

/**
 * Finite difference stencil order.
 * @readonly
 * @enum {number}
 */
export const StencilOrder = Object.freeze({
  SECOND: 2,   // 2nd-order central difference (3-point stencil)
  FOURTH: 4,   // 4th-order central difference (5-point stencil)
  SIXTH: 6     // 6th-order central difference (7-point stencil)
});

/**
 * Statistical summary container for field arrays.
 */
export class FieldStatistics {
  /**
   * @param {Float32Array|Float64Array} buffer
   */
  constructor(buffer) {
    if (!buffer || buffer.length === 0) {
      throw new TypeError('FieldStatistics: buffer must be non-empty TypedArray.');
    }

    const n = buffer.length;
    let min = Infinity;
    let max = -Infinity;
    let sum = 0.0;
    let sumSq = 0.0;

    for (let i = 0; i < n; i++) {
      const v = buffer[i];
      if (v < min) min = v;
      if (v > max) max = v;
      sum += v;
      sumSq += v * v;
    }

    const mean = sum / n;
    const variance = Math.max(0.0, (sumSq / n) - (mean * mean));
    const stdDev = Math.sqrt(variance);
    const rms = Math.sqrt(sumSq / n);

    // Compute skewness and kurtosis
    let sumCube = 0.0;
    let sumQuad = 0.0;
    for (let i = 0; i < n; i++) {
      const diff = buffer[i] - mean;
      const diff2 = diff * diff;
      sumCube += diff2 * diff;
      sumQuad += diff2 * diff2;
    }

    const skewness = stdDev > 1e-12 ? (sumCube / n) / Math.pow(stdDev, 3) : 0.0;
    const kurtosis = stdDev > 1e-12 ? (sumQuad / n) / Math.pow(stdDev, 4) - 3.0 : 0.0;

    this.count = n;
    this.min = min;
    this.max = max;
    this.mean = mean;
    this.stdDev = stdDev;
    this.rms = rms;
    this.variance = variance;
    this.skewness = skewness;
    this.kurtosis = kurtosis;
    Object.freeze(this);
  }

  toJSON() {
    return {
      count: this.count,
      min: this.min,
      max: this.max,
      mean: this.mean,
      stdDev: this.stdDev,
      rms: this.rms,
      skewness: this.skewness,
      kurtosis: this.kurtosis
    };
  }
}

/**
 * FieldDifferentialOperators provides high-order continuous and discrete differential
 * operators for cosmological velocity, density, and potential fields.
 */
export class FieldDifferentialOperators {
  /**
   * @param {GridIndexer} gridIndexer Mesh geometry and boundary conditions.
   * @param {Object} [options]
   * @param {number} [options.stencilOrder=4] Finite difference order (2, 4, 6).
   * @param {TricubicInterpolator} [options.interpolator] Optional custom interpolator.
   */
  constructor(gridIndexer, options = {}) {
    if (!gridIndexer || !(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('FieldDifferentialOperators: gridIndexer must be an instance of GridIndexer.');
    }

    this.grid = gridIndexer;
    this.stencilOrder = options.stencilOrder || StencilOrder.FOURTH;
    this.interpolator = options.interpolator || new TricubicInterpolator(this.grid);
  }

  // ==========================================================================
  // Continuous Point-Wise Operators (via Tricubic Spline Derivatives)
  // ==========================================================================

  /**
   * Continuous gradient grad(f) = [df/dx, df/dy, df/dz] at arbitrary (x, y, z).
   * 
   * @param {Float32Array|Float64Array} buffer 1D flat scalar grid.
   * @param {number} x SGX in physical units (Mpc/h).
   * @param {number} y SGY in physical units (Mpc/h).
   * @param {number} z SGZ in physical units (Mpc/h).
   * @returns {Array<number>} [df/dx, df/dy, df/dz].
   */
  gradient(buffer, x, y, z) {
    return this.interpolator.interpolateGradient(buffer, x, y, z);
  }

  /**
   * Continuous velocity divergence div(v) = dvx/dx + dvy/dy + dvz/dz at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number} Divergence in (km/s)/(Mpc/h).
   */
  divergence(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const dvx = this.interpolator.interpolateGradient(vxBuffer, x, y, z);
    const dvy = this.interpolator.interpolateGradient(vyBuffer, x, y, z);
    const dvz = this.interpolator.interpolateGradient(vzBuffer, x, y, z);
    return dvx[0] + dvy[1] + dvz[2];
  }

  /**
   * Continuous vorticity vector curl(v) = [dvz/dy - dvy/dz, dvx/dz - dvz/dx, dvy/dx - dvx/dy] at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Array<number>} [omega_x, omega_y, omega_z] in (km/s)/(Mpc/h).
   */
  curl(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const dvx = this.interpolator.interpolateGradient(vxBuffer, x, y, z);
    const dvy = this.interpolator.interpolateGradient(vyBuffer, x, y, z);
    const dvz = this.interpolator.interpolateGradient(vzBuffer, x, y, z);

    return [
      dvz[1] - dvy[2], // dvz/dy - dvy/dz
      dvx[2] - dvz[0], // dvx/dz - dvz/dx
      dvy[0] - dvx[1]  // dvy/dx - dvx/dy
    ];
  }

  /**
   * Continuous scalar Laplacian lap(f) = d2f/dx2 + d2f/dy2 + d2f/dz2 at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} buffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number} Laplacian value.
   */
  laplacian(buffer, x, y, z) {
    const { laplacian } = this.interpolator.interpolateHessian(buffer, x, y, z);
    return laplacian;
  }

  /**
   * Continuous 3x3 velocity Jacobian tensor J_ij = dv_i / dx_j at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Jacobian3x3}
   */
  jacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    return this.interpolator.interpolateJacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z);
  }

  /**
   * Continuous symmetric strain rate tensor S_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i) at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Float64Array} 9-element symmetric matrix.
   */
  strainRateTensor(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const J = this.jacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z);
    return J.symmetricStrainRate();
  }

  /**
   * Continuous traceless shear tensor sigma_ij = S_ij - (1/3)*div(v)*delta_ij at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Float64Array} 9-element symmetric traceless matrix.
   */
  shearTensor(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const J = this.jacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z);
    return J.shearTensor();
  }

  /**
   * Diagonalizes the continuous strain rate tensor into principal strain eigenvalues
   * and eigenvectors at (x, y, z).
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {EigenSystem3D}
   */
  strainEigenSystem(vxBuffer, vyBuffer, vzBuffer, x, y, z) {
    const J = this.jacobian(vxBuffer, vyBuffer, vzBuffer, x, y, z);
    return J.diagonalizeStrain();
  }

  // ==========================================================================
  // High-Order Discrete Grid Differential Operators
  // ==========================================================================

  /**
   * Discrete partial derivative df/dx along SGX, SGY, or SGZ using high-order finite differences.
   * 
   * @param {Float32Array|Float64Array} buffer 1D flat grid array.
   * @param {number} axis 0 for SGX, 1 for SGY, 2 for SGZ.
   * @param {Float64Array} [outBuffer] Pre-allocated output buffer.
   * @returns {Float64Array} Output derivative grid.
   */
  discreteDerivative(buffer, axis, outBuffer = null) {
    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const total = this.grid.totalCells;
    const out = outBuffer || new Float64Array(total);
    const h = axis === 0 ? this.grid.dx : (axis === 1 ? this.grid.dy : this.grid.dz);

    const order = this.stencilOrder;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const idx = this.grid.getLinearIndex(ix, iy, iz);

          if (order === StencilOrder.SIXTH) {
            // 6th-order central difference: (-f_-3 + 9 f_-2 - 45 f_-1 + 45 f_1 - 9 f_2 + f_3) / (60 h)
            const vM3 = this._sampleOffset(buffer, ix, iy, iz, axis, -3);
            const vM2 = this._sampleOffset(buffer, ix, iy, iz, axis, -2);
            const vM1 = this._sampleOffset(buffer, ix, iy, iz, axis, -1);
            const vP1 = this._sampleOffset(buffer, ix, iy, iz, axis, 1);
            const vP2 = this._sampleOffset(buffer, ix, iy, iz, axis, 2);
            const vP3 = this._sampleOffset(buffer, ix, iy, iz, axis, 3);
            out[idx] = (-vM3 + 9.0 * vM2 - 45.0 * vM1 + 45.0 * vP1 - 9.0 * vP2 + vP3) / (60.0 * h);
          } else if (order === StencilOrder.FOURTH) {
            // 4th-order central difference: (-f_2 + 8 f_1 - 8 f_-1 + f_-2) / (12 h)
            const vM2 = this._sampleOffset(buffer, ix, iy, iz, axis, -2);
            const vM1 = this._sampleOffset(buffer, ix, iy, iz, axis, -1);
            const vP1 = this._sampleOffset(buffer, ix, iy, iz, axis, 1);
            const vP2 = this._sampleOffset(buffer, ix, iy, iz, axis, 2);
            out[idx] = (-vP2 + 8.0 * vP1 - 8.0 * vM1 + vM2) / (12.0 * h);
          } else {
            // 2nd-order central difference: (f_1 - f_-1) / (2 h)
            const vM1 = this._sampleOffset(buffer, ix, iy, iz, axis, -1);
            const vP1 = this._sampleOffset(buffer, ix, iy, iz, axis, 1);
            out[idx] = (vP1 - vM1) / (2.0 * h);
          }
        }
      }
    }

    return out;
  }

  /**
   * Samples a grid value with index offset along an axis respecting boundary mode.
   * @param {Float32Array|Float64Array} buffer
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @param {number} axis
   * @param {number} offset
   * @returns {number}
   */
  _sampleOffset(buffer, ix, iy, iz, axis, offset) {
    let ox = ix, oy = iy, oz = iz;
    if (axis === 0) ox += offset;
    else if (axis === 1) oy += offset;
    else if (axis === 2) oz += offset;

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const bMode = this.grid.boundaryMode;

    if (ox >= 0 && ox < nx && oy >= 0 && oy < ny && oz >= 0 && oz < nz) {
      return buffer[this.grid.getLinearIndex(ox, oy, oz)];
    }

    if (bMode === BoundaryMode.PERIODIC) {
      ox = ((ox % nx) + nx) % nx;
      oy = ((oy % ny) + ny) % ny;
      oz = ((oz % nz) + nz) % nz;
      return buffer[this.grid.getLinearIndex(ox, oy, oz)];
    }

    if (bMode === BoundaryMode.REFLECT) {
      if (ox < 0) ox = -ox;
      if (ox >= nx) ox = 2 * (nx - 1) - ox;
      if (oy < 0) oy = -oy;
      if (oy >= ny) oy = 2 * (ny - 1) - oy;
      if (oz < 0) oz = -oz;
      if (oz >= nz) oz = 2 * (nz - 1) - oz;
      ox = Math.max(0, Math.min(nx - 1, ox));
      oy = Math.max(0, Math.min(ny - 1, oy));
      oz = Math.max(0, Math.min(nz - 1, oz));
      return buffer[this.grid.getLinearIndex(ox, oy, oz)];
    }

    // CLAMP or fallback
    ox = Math.max(0, Math.min(nx - 1, ox));
    oy = Math.max(0, Math.min(ny - 1, oy));
    oz = Math.max(0, Math.min(nz - 1, oz));
    return buffer[this.grid.getLinearIndex(ox, oy, oz)];
  }

  /**
   * Computes discrete 3D divergence grid theta = div(v) = dvx/dx + dvy/dy + dvz/dz.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {Float64Array} [outBuffer]
   * @returns {Float64Array} Divergence field grid.
   */
  computeDivergenceGrid(vxBuffer, vyBuffer, vzBuffer, outBuffer = null) {
    const total = this.grid.totalCells;
    const div = outBuffer || new Float64Array(total);

    const dvx_dx = this.discreteDerivative(vxBuffer, 0);
    const dvy_dy = this.discreteDerivative(vyBuffer, 1);
    const dvz_dz = this.discreteDerivative(vzBuffer, 2);

    for (let i = 0; i < total; i++) {
      div[i] = dvx_dx[i] + dvy_dy[i] + dvz_dz[i];
    }

    return div;
  }

  /**
   * Computes discrete 3D vorticity vector field grids: (omega_x, omega_y, omega_z) and magnitude ||omega||.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @returns {{wx: Float64Array, wy: Float64Array, wz: Float64Array, magnitude: Float64Array}}
   */
  computeVorticityGrid(vxBuffer, vyBuffer, vzBuffer) {
    const total = this.grid.totalCells;

    const dvx_dy = this.discreteDerivative(vxBuffer, 1);
    const dvx_dz = this.discreteDerivative(vxBuffer, 2);

    const dvy_dx = this.discreteDerivative(vyBuffer, 0);
    const dvy_dz = this.discreteDerivative(vyBuffer, 2);

    const dvz_dx = this.discreteDerivative(vzBuffer, 0);
    const dvz_dy = this.discreteDerivative(vzBuffer, 1);

    const wx = new Float64Array(total);
    const wy = new Float64Array(total);
    const wz = new Float64Array(total);
    const magnitude = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const ox = dvz_dy[i] - dvy_dz[i];
      const oy = dvx_dz[i] - dvz_dx[i];
      const oz = dvy_dx[i] - dvx_dy[i];

      wx[i] = ox;
      wy[i] = oy;
      wz[i] = oz;
      magnitude[i] = Math.sqrt(ox * ox + oy * oy + oz * oz);
    }

    return { wx, wy, wz, magnitude };
  }

  /**
   * Computes discrete specific kinetic energy grid: e_k = 0.5 * ||v||^2 in (km/s)^2.
   * When delta is provided without background mass density rho_bar, computes density-weighted proxy:
   *   e_k_weighted = 0.5 * (1 + delta) * ||v||^2.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {Float32Array|Float64Array} [densityBuffer] Optional density contrast delta.
   * @returns {Float64Array} Specific kinetic energy grid in (km/s)^2.
   */
  computeSpecificKineticEnergyGrid(vxBuffer, vyBuffer, vzBuffer, densityBuffer = null) {
    const total = this.grid.totalCells;
    const energy = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const vx = vxBuffer[i];
      const vy = vyBuffer[i];
      const vz = vzBuffer[i];
      const vSq = vx * vx + vy * vy + vz * vz;
      const weight = densityBuffer ? (1.0 + Math.max(-1.0, densityBuffer[i])) : 1.0;
      energy[i] = 0.5 * weight * vSq;
    }

    return energy;
  }

  /**
   * Computes physical kinetic energy density grid: E_k = 0.5 * rho_bar * (1 + delta) * ||v||^2
   * in (M_sun / (Mpc/h)^3) * (km/s)^2 requiring explicit background cosmology.
   */
  computePhysicalKineticEnergyDensityGrid(vxBuffer, vyBuffer, vzBuffer, densityBuffer = null, omegaM = 0.315) {
    const total = this.grid.totalCells;
    const RHO_CRIT_0 = 2.77536627e11; // M_sun / (Mpc/h)^3
    const rhoBar = omegaM * RHO_CRIT_0;
    const specificKE = this.computeSpecificKineticEnergyGrid(vxBuffer, vyBuffer, vzBuffer, densityBuffer);
    const energyDensity = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      energyDensity[i] = rhoBar * specificKE[i];
    }

    return energyDensity;
  }

  /**
   * Alias for backward compatibility.
   */
  computeKineticEnergyGrid(vxBuffer, vyBuffer, vzBuffer, densityBuffer = null, omegaM = 0.315) {
    return this.computePhysicalKineticEnergyDensityGrid(vxBuffer, vyBuffer, vzBuffer, densityBuffer, omegaM);
  }

  /**
   * Computes discrete enstrophy grid: E = 0.5 * ||omega||^2.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @returns {Float64Array} Enstrophy grid in ((km/s)/(Mpc/h))^2.
   */
  computeEnstrophyGrid(vxBuffer, vyBuffer, vzBuffer) {
    const { magnitude } = this.computeVorticityGrid(vxBuffer, vyBuffer, vzBuffer);
    const total = this.grid.totalCells;
    const enstrophy = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const m = magnitude[i];
      enstrophy[i] = 0.5 * m * m;
    }

    return enstrophy;
  }

  /**
   * Computes discrete kinetic helicity density grid: h = v . omega.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @returns {Float64Array} Helicity grid in (km/s)^2 / (Mpc/h).
   */
  computeHelicityGrid(vxBuffer, vyBuffer, vzBuffer) {
    const { wx, wy, wz } = this.computeVorticityGrid(vxBuffer, vyBuffer, vzBuffer);
    const total = this.grid.totalCells;
    const helicity = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      helicity[i] = vxBuffer[i] * wx[i] + vyBuffer[i] * wy[i] + vzBuffer[i] * wz[i];
    }

    return helicity;
  }

  /**
   * Computes discrete 3D Q-criterion grid: Q = 0.5 * (||Omega||^2 - ||S||^2)
   * where Omega is the antisymmetric vorticity tensor and S is the symmetric strain tensor.
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @returns {Float64Array} 3D Q-criterion parameter grid.
   */
  computeQCriterionGrid(vxBuffer, vyBuffer, vzBuffer) {
    const total = this.grid.totalCells;

    const dvx_dx = this.discreteDerivative(vxBuffer, 0);
    const dvx_dy = this.discreteDerivative(vxBuffer, 1);
    const dvx_dz = this.discreteDerivative(vxBuffer, 2);

    const dvy_dx = this.discreteDerivative(vyBuffer, 0);
    const dvy_dy = this.discreteDerivative(vyBuffer, 1);
    const dvy_dz = this.discreteDerivative(vyBuffer, 2);

    const dvz_dx = this.discreteDerivative(vzBuffer, 0);
    const dvz_dy = this.discreteDerivative(vzBuffer, 1);
    const dvz_dz = this.discreteDerivative(vzBuffer, 2);

    const Q = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      // Antisymmetric vorticity tensor components: Omega_ij = 0.5 * (dv_i/dx_j - dv_j/dx_i)
      const w01 = 0.5 * (dvx_dy[i] - dvy_dx[i]);
      const w02 = 0.5 * (dvx_dz[i] - dvz_dx[i]);
      const w12 = 0.5 * (dvy_dz[i] - dvz_dy[i]);
      const normOmegaSq = 2.0 * (w01 * w01 + w02 * w02 + w12 * w12);

      // Symmetric Strain tensor components: S_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i)
      const s00 = dvx_dx[i];
      const s11 = dvy_dy[i];
      const s22 = dvz_dz[i];
      const s01 = 0.5 * (dvx_dy[i] + dvy_dx[i]);
      const s02 = 0.5 * (dvx_dz[i] + dvz_dx[i]);
      const s12 = 0.5 * (dvy_dz[i] + dvz_dy[i]);
      const normStrainSq = s00 * s00 + s11 * s11 + s22 * s22 + 2.0 * (s01 * s01 + s02 * s02 + s12 * s12);

      Q[i] = 0.5 * (normOmegaSq - normStrainSq);
    }

    return Q;
  }

  /**
   * Alias for computeQCriterionGrid for 3D flows.
   */
  computeOkuboWeissGrid(vxBuffer, vyBuffer, vzBuffer) {
    return this.computeQCriterionGrid(vxBuffer, vyBuffer, vzBuffer);
  }

  /**
   * Cosmic Web Classification Map:
   * Computes the 3 eigenvalues of the velocity shear tensor or deformation tensor at every grid voxel,
   * and classifies each voxel into:
   *   - Void (0 eigenvalues > threshold) [Code: 0]
   *   - Sheet (1 eigenvalue > threshold) [Code: 1]
   *   - Filament (2 eigenvalues > threshold) [Code: 2]
   *   - Knot / Cluster (3 eigenvalues > threshold) [Code: 3]
   * 
   * @param {Float32Array|Float64Array} vxBuffer
   * @param {Float32Array|Float64Array} vyBuffer
   * @param {Float32Array|Float64Array} vzBuffer
   * @param {number} [threshold=0.0] Eigenvalue threshold lambda_th.
   * @returns {{
   *   classificationGrid: Uint8Array,
   *   lambda1Grid: Float64Array,
   *   lambda2Grid: Float64Array,
   *   lambda3Grid: Float64Array,
   *   volumeFractions: {voidPct: number, sheetPct: number, filamentPct: number, knotPct: number}
   * }}
   */
  classifyCosmicWebGrid(vxBuffer, vyBuffer, vzBuffer, threshold = 0.0) {
    const total = this.grid.totalCells;

    const dvx_dx = this.discreteDerivative(vxBuffer, 0);
    const dvx_dy = this.discreteDerivative(vxBuffer, 1);
    const dvx_dz = this.discreteDerivative(vxBuffer, 2);

    const dvy_dx = this.discreteDerivative(vyBuffer, 0);
    const dvy_dy = this.discreteDerivative(vyBuffer, 1);
    const dvy_dz = this.discreteDerivative(vyBuffer, 2);

    const dvz_dx = this.discreteDerivative(vzBuffer, 0);
    const dvz_dy = this.discreteDerivative(vzBuffer, 1);
    const dvz_dz = this.discreteDerivative(vzBuffer, 2);

    const classGrid = new Uint8Array(total);
    const l1Grid = new Float64Array(total);
    const l2Grid = new Float64Array(total);
    const l3Grid = new Float64Array(total);

    const counts = [0, 0, 0, 0];
    const S = new Float64Array(9);

    for (let i = 0; i < total; i++) {
      // In cosmological conventions, deformation / tidal tensor is proportional to -0.5 * (dv_i/dx_j + dv_j/dx_i) / (H0*f)
      // Positive eigenvalue corresponds to gravitational collapse/compression.
      // S_ij = -0.5 * (dv_i/dx_j + dv_j/dx_i)
      S[0] = -dvx_dx[i];
      S[4] = -dvy_dy[i];
      S[8] = -dvz_dz[i];
      S[1] = S[3] = -0.5 * (dvx_dy[i] + dvy_dx[i]);
      S[2] = S[6] = -0.5 * (dvx_dz[i] + dvz_dx[i]);
      S[5] = S[7] = -0.5 * (dvy_dz[i] + dvz_dy[i]);

      const eigen = EigenSystem3D.fromSymmetricMatrix(S);
      const l1 = eigen.lambda1;
      const l2 = eigen.lambda2;
      const l3 = eigen.lambda3;

      l1Grid[i] = l1;
      l2Grid[i] = l2;
      l3Grid[i] = l3;

      const code = eigen.classifyWebStructure(threshold).code;
      classGrid[i] = code;
      counts[code]++;
    }

    const volumeFractions = {
      voidPct: (counts[0] / total) * 100.0,
      sheetPct: (counts[1] / total) * 100.0,
      filamentPct: (counts[2] / total) * 100.0,
      knotPct: (counts[3] / total) * 100.0
    };

    return {
      classificationGrid: classGrid,
      lambda1Grid: l1Grid,
      lambda2Grid: l2Grid,
      lambda3Grid: l3Grid,
      volumeFractions
    };
  }

  // ==========================================================================
  // Diagnostic Kinetic Sanity Maps & Linear Continuity Audit Suite
  // ==========================================================================

  /**
   * Generates a complete suite of diagnostic kinetic sanity maps, statistical summaries,
   * linear perturbation theory continuity verification (theta ~= -H0*f*delta),
   * irrotationality audit, and W3C PROV-O reproducible manifest.
   * 
   * @param {Object} fields Input fields container.
   * @param {Float32Array|Float64Array} fields.vx Velocity X buffer (km/s).
   * @param {Float32Array|Float64Array} fields.vy Velocity Y buffer (km/s).
   * @param {Float32Array|Float64Array} fields.vz Velocity Z buffer (km/s).
   * @param {Float32Array|Float64Array} [fields.density] Density contrast delta buffer.
   * @param {Object} [options]
   * @param {number} [options.h0f=52.0] Linear growth rate product H0 * f(Omega_m) in (km/s)/(Mpc/h).
   * @param {number} [options.webThreshold=0.0] Threshold for cosmic web eigenvalue classification.
   * @returns {Object} Full kinetic sanity dossier.
   */
  generateKineticSanityMaps(fields, options = {}) {
    if (!fields || !fields.vx || !fields.vy || !fields.vz) {
      throw new TypeError('generateKineticSanityMaps: fields must contain vx, vy, and vz buffers.');
    }

    const startTime = Date.now();
    const vx = fields.vx;
    const vy = fields.vy;
    const vz = fields.vz;
    const density = fields.density || null;
    const cosmo = options.cosmology instanceof CosmologicalParameters 
      ? options.cosmology 
      : new CosmologicalParameters(typeof options.h0f === 'number' ? { H0: options.h0f / 0.524 } : {});
    const continuityCoeff = cosmo.continuityCoefficient;
    const threshold = typeof options.webThreshold === 'number' ? options.webThreshold : 0.0;

    const total = this.grid.totalCells;

    // 1. Divergence Map
    const divergenceGrid = this.computeDivergenceGrid(vx, vy, vz);
    const divStats = new FieldStatistics(divergenceGrid);

    // 2. Vorticity Map
    const vorticity = this.computeVorticityGrid(vx, vy, vz);
    const vortStats = new FieldStatistics(vorticity.magnitude);

    // 3. Kinetic Energy Map
    const kineticEnergyGrid = this.computeKineticEnergyGrid(vx, vy, vz, density);
    const keStats = new FieldStatistics(kineticEnergyGrid);

    // 4. Enstrophy Map
    const enstrophyGrid = this.computeEnstrophyGrid(vx, vy, vz);
    const enstrophyStats = new FieldStatistics(enstrophyGrid);

    // 5. Helicity Map
    const helicityGrid = this.computeHelicityGrid(vx, vy, vz);
    const helicityStats = new FieldStatistics(helicityGrid);

    // 6. Okubo-Weiss Map
    const okuboWeissGrid = this.computeOkuboWeissGrid(vx, vy, vz);
    const owStats = new FieldStatistics(okuboWeissGrid);

    // 7. Cosmic Web Classification
    const web = this.classifyCosmicWebGrid(vx, vy, vz, threshold);

    // 8. Linear Perturbation Theory Continuity Residual:
    // In linear theory, theta = div(v) = -H0 * f * delta  =>  delta_linear = -div(v) / (H0 * f)
    // Residual = div(v) + H0 * f * delta
    let continuityReport = null;
    if (density) {
      const residualGrid = new Float64Array(total);
      let resSqSum = 0.0;
      let divSqSum = 0.0;

      for (let i = 0; i < total; i++) {
        const expectedDiv = -continuityCoeff * density[i];
        const actualDiv = divergenceGrid[i];
        const res = actualDiv - expectedDiv;
        residualGrid[i] = res;
        resSqSum += res * res;
        divSqSum += actualDiv * actualDiv;
      }

      const l2Residual = Math.sqrt(resSqSum / total);
      const l2Div = Math.sqrt(divSqSum / total);
      const relativeResidual = l2Div > 1e-12 ? l2Residual / l2Div : 0.0;
      const resStats = new FieldStatistics(residualGrid);

      continuityReport = {
        continuityCoefficient: continuityCoeff,
        cosmology: { H0: cosmo.H0, omegaM: cosmo.omegaM, a: cosmo.scaleFactorA, f: cosmo.growthRateF },
        l2Residual,
        l2Divergence: l2Div,
        relativeResidual,
        residualStats: resStats.toJSON(),
        isLinearConsistent: relativeResidual < 0.25 // < 25% residual on discrete 16^3 finite difference grid
      };
    }

    // 9. Irrotationality Metric (Potential Flow Assessment)
    // Ratio of RMS vorticity to RMS divergence: in purely irrotational potential flow, ratio -> 0
    const vortToDivRatio = divStats.rms > 1e-12 ? vortStats.rms / divStats.rms : 0.0;
    const isPotentialFlowDominant = vortToDivRatio < 0.25;

    const durationMs = Date.now() - startTime;

    // Cryptographic SHA-256 fingerprint of diagnostics
    const hashPayload = JSON.stringify({
      nx: this.grid.nx,
      ny: this.grid.ny,
      nz: this.grid.nz,
      divMean: divStats.mean,
      vortRms: vortStats.rms,
      keMean: keStats.mean,
      fractions: web.volumeFractions
    });
    const diagnosticsHash = sha256Hex(hashPayload);

    const manifest = new RunManifest({
      runId: `kinetic-sanity-${Date.now()}`,
      algorithmName: 'FieldDifferentialOperators.generateKineticSanityMaps',
      parameters: {
        continuityCoefficient: continuityCoeff,
        H0: cosmo.H0,
        omegaM: cosmo.omegaM,
        stencilOrder: this.stencilOrder,
        boundaryMode: this.grid.boundaryMode,
        webThreshold: threshold
      },
      inputDatasetHash: sha256Hex(`grid-${total}`),
      outputFieldHash: diagnosticsHash,
      executionDurationMs: durationMs
    });

    return {
      grid: this.grid,
      maps: {
        divergence: divergenceGrid,
        vorticityX: vorticity.wx,
        vorticityY: vorticity.wy,
        vorticityZ: vorticity.wz,
        vorticityMagnitude: vorticity.magnitude,
        kineticEnergy: kineticEnergyGrid,
        enstrophy: enstrophyGrid,
        helicity: helicityGrid,
        okuboWeiss: okuboWeissGrid,
        cosmicWeb: web.classificationGrid,
        lambda1: web.lambda1Grid,
        lambda2: web.lambda2Grid,
        lambda3: web.lambda3Grid
      },
      statistics: {
        divergence: divStats.toJSON(),
        vorticity: vortStats.toJSON(),
        kineticEnergy: keStats.toJSON(),
        enstrophy: enstrophyStats.toJSON(),
        helicity: helicityStats.toJSON(),
        okuboWeiss: owStats.toJSON()
      },
      cosmicWeb: {
        volumeFractions: web.volumeFractions,
        threshold
      },
      continuityAudit: continuityReport,
      flowIrrotationality: {
        vorticityToDivergenceRatio: vortToDivRatio,
        isPotentialFlowDominant,
        assessment: isPotentialFlowDominant
          ? 'Flow is potential/irrotational dominant (compatible with linear cosmological perturbation theory)'
          : 'High vorticity detected (turbulent or non-linear multi-streaming regime)'
      },
      manifest: manifest.toJSON()
    };
  }
}
