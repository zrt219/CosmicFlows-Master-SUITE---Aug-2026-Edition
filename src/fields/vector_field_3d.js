/**
 * @file vector_field_3d.js
 * @description High-performance 3D vector field container (vx, vy, vz in km/s), Helmholtz-Hodge decomposition,
 * radial infall velocity profiles, and bulk flow multipole estimators for the CosmicFlows-4 Research Workbench.
 * 
 * Implements:
 * 1. 3D vector field storage with GridIndexer and x52.0 scale factor tracking.
 * 2. Vector continuous sampling via Trilinear and 64-point C^1 Tricubic Hermite splines.
 * 3. Analytical within-cell 3x3 velocity Jacobian J_ij = dv_i / dx_j.
 * 4. Divergence theta = div(v), Vorticity omega = curl(v), and Kinetic Energy E_k = 0.5 * ||v||^2.
 * 5. Radial infall velocity v_rad(r) and tangential velocity v_tan(r) relative to cluster centers.
 * 6. Global and spherical shell bulk flow dipole velocity integration.
 * 7. Helmholtz-Hodge decomposition into irrotational (potential) and solenoidal (divergence-free) components.
 * 
 * @module fields/vector_field_3d
 */

import { GridIndexer } from './grid_indexer.js';
import { ScalarField3D } from './scalar_field_3d.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { TricubicInterpolator } from '../interpolation/tricubic_interpolator.js';
import {
  SupergalacticPosition,
  VelocityVector,
  Jacobian3x3,
  CF4_PUBLIC_VELOCITY_SCALE
} from '../coordinates/scientific_types.js';

/**
 * Universal 3D vector field (vx, vy, vz).
 */
export class VectorField3D {
  /**
   * @param {GridIndexer} gridIndexer
   * @param {Float32Array|Float64Array} [vx]
   * @param {Float32Array|Float64Array} [vy]
   * @param {Float32Array|Float64Array} [vz]
   * @param {Object} [options]
   * @param {string} [options.name='velocity_field']
   * @param {string} [options.unit='km/s']
   * @param {boolean} [options.scaleFactorApplied=false]
   */
  constructor(gridIndexer, vx = null, vy = null, vz = null, options = {}) {
    if (!(gridIndexer instanceof GridIndexer)) {
      throw new TypeError('VectorField3D: gridIndexer must be an instance of GridIndexer.');
    }
    this.grid = gridIndexer;
    this.name = options.name || 'velocity_field';
    this.unit = options.unit || 'km/s';
    this.scaleFactorApplied = Boolean(options.scaleFactorApplied);

    const total = this.grid.totalCells;
    this.vx = new Float64Array(total);
    this.vy = new Float64Array(total);
    this.vz = new Float64Array(total);

    if (vx && vy && vz) {
      if (vx.length !== total || vy.length !== total || vz.length !== total) {
        throw new Error(`VectorField3D: Buffer sizes must match grid total cells (${total}).`);
      }
      for (let i = 0; i < total; i++) {
        this.vx[i] = vx[i];
        this.vy[i] = vy[i];
        this.vz[i] = vz[i];
      }
    }

    this._trilinear = null;
    this._tricubic = null;
  }

  /**
   * Lazily initializes TrilinearInterpolator.
   * @returns {TrilinearInterpolator}
   */
  get trilinear() {
    if (!this._trilinear) {
      this._trilinear = new TrilinearInterpolator(this.grid);
    }
    return this._trilinear;
  }

  /**
   * Lazily initializes TricubicInterpolator for vector fields.
   * @returns {TricubicInterpolator}
   */
  get tricubic() {
    if (!this._tricubic) {
      this._tricubic = new TricubicInterpolator(this.grid, this.vx, this.vy, this.vz);
    }
    return this._tricubic;
  }

  /**
   * Gets vector at grid voxel (ix, iy, iz) as VelocityVector.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {VelocityVector}
   */
  getVector(ix, iy, iz) {
    const idx = this.grid.index(ix, iy, iz);
    return new VelocityVector(this.vx[idx], this.vy[idx], this.vz[idx]);
  }

  /**
   * Sets vector at grid voxel (ix, iy, iz).
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @param {number} vx
   * @param {number} vy
   * @param {number} vz
   */
  setVector(ix, iy, iz, vx, vy, vz) {
    const idx = this.grid.index(ix, iy, iz);
    this.vx[idx] = vx;
    this.vy[idx] = vy;
    this.vz[idx] = vz;
  }

  /**
   * Continuous trilinear sampling of velocity vector [vx, vy, vz] in km/s.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {VelocityVector}
   */
  sampleTrilinear(x, y, z) {
    const v = this.trilinear.evaluateVector(this.vx, this.vy, this.vz, x, y, z);
    return new VelocityVector(v[0], v[1], v[2]);
  }

  /**
   * Continuous 64-point C^1 Tricubic Hermite spline evaluation of velocity vector.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {VelocityVector}
   */
  sampleTricubic(x, y, z) {
    const v = this.tricubic.evaluateVector(x, y, z);
    return new VelocityVector(v[0], v[1], v[2]);
  }

  /**
   * Analytical 3x3 velocity Jacobian tensor J_ij = dv_i / dx_j.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Jacobian3x3}
   */
  evaluateJacobian(x, y, z) {
    const rawJ = this.tricubic.evaluateVectorJacobian(x, y, z);
    return new Jacobian3x3([
      [rawJ[0], rawJ[1], rawJ[2]],
      [rawJ[3], rawJ[4], rawJ[5]],
      [rawJ[6], rawJ[7], rawJ[8]]
    ]);
  }

  /**
   * Applies the official CosmicFlows-4 x52.0 velocity scale factor.
   * Strictly enforces that the scaling is applied only once.
   * @param {number} [scale=52.0]
   * @returns {this}
   */
  applyScaleFactor(scale = CF4_PUBLIC_VELOCITY_SCALE) {
    if (this.scaleFactorApplied) {
      throw new Error('VectorField3D.applyScaleFactor: Scale factor x52.0 has already been applied! (Rule 3 Invariant)');
    }
    const total = this.grid.totalCells;
    for (let i = 0; i < total; i++) {
      this.vx[i] *= scale;
      this.vy[i] *= scale;
      this.vz[i] *= scale;
    }
    this.scaleFactorApplied = true;
    return this;
  }

  /**
   * Computes scalar magnitude field ||v|| = sqrt(vx^2 + vy^2 + vz^2).
   * @returns {ScalarField3D}
   */
  computeSpeedField() {
    const total = this.grid.totalCells;
    const speedBuf = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const vx = this.vx[i];
      const vy = this.vy[i];
      const vz = this.vz[i];
      speedBuf[i] = Math.sqrt(vx * vx + vy * vy + vz * vz);
    }
    return new ScalarField3D(this.grid, speedBuf, `${this.name}_speed`, this.unit);
  }

  /**
   * Computes kinetic energy density field E_k = 0.5 * ||v||^2 in (km/s)^2.
   * @returns {ScalarField3D}
   */
  computeKineticEnergyField() {
    const total = this.grid.totalCells;
    const keBuf = new Float64Array(total);
    for (let i = 0; i < total; i++) {
      const vx = this.vx[i];
      const vy = this.vy[i];
      const vz = this.vz[i];
      keBuf[i] = 0.5 * (vx * vx + vy * vy + vz * vz);
    }
    return new ScalarField3D(this.grid, keBuf, `${this.name}_kinetic_energy`, '(km/s)^2');
  }

  /**
   * Computes discrete velocity divergence field theta = div(v) = dvx/dx + dvy/dy + dvz/dz.
   * @returns {ScalarField3D}
   */
  computeDivergenceField() {
    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const dx2 = 2.0 * this.grid.dx;
    const dy2 = 2.0 * this.grid.dy;
    const dz2 = 2.0 * this.grid.dz;

    const divBuf = new Float64Array(nx * ny * nz);

    for (let iz = 0; iz < nz; iz++) {
      const izM = iz > 0 ? iz - 1 : 0;
      const izP = iz < nz - 1 ? iz + 1 : nz - 1;
      for (let iy = 0; iy < ny; iy++) {
        const iyM = iy > 0 ? iy - 1 : 0;
        const iyP = iy < ny - 1 ? iy + 1 : ny - 1;
        for (let ix = 0; ix < nx; ix++) {
          const ixM = ix > 0 ? ix - 1 : 0;
          const ixP = ix < nx - 1 ? ix + 1 : nx - 1;

          const dvx_dx = (this.vx[this.grid.index(ixP, iy, iz)] - this.vx[this.grid.index(ixM, iy, iz)]) / (ixP !== ixM ? (ixP - ixM) * this.grid.dx : 1.0);
          const dvy_dy = (this.vy[this.grid.index(ix, iyP, iz)] - this.vy[this.grid.index(ix, iyM, iz)]) / (iyP !== iyM ? (iyP - iyM) * this.grid.dy : 1.0);
          const dvz_dz = (this.vz[this.grid.index(ix, iy, izP)] - this.vz[this.grid.index(ix, iy, izM)]) / (izP !== izM ? (izP - izM) * this.grid.dz : 1.0);

          divBuf[this.grid.index(ix, iy, iz)] = dvx_dx + dvy_dy + dvz_dz;
        }
      }
    }

    return new ScalarField3D(this.grid, divBuf, `${this.name}_divergence`, 'km/s/(Mpc/h)');
  }

  /**
   * Computes discrete vorticity vector field omega = curl(v).
   * @returns {VectorField3D}
   */
  computeVorticityField() {
    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;
    const dx = this.grid.dx;
    const dy = this.grid.dy;
    const dz = this.grid.dz;

    const wx = new Float64Array(nx * ny * nz);
    const wy = new Float64Array(nx * ny * nz);
    const wz = new Float64Array(nx * ny * nz);

    for (let iz = 0; iz < nz; iz++) {
      const izM = iz > 0 ? iz - 1 : 0;
      const izP = iz < nz - 1 ? iz + 1 : nz - 1;
      const dzEff = (izP - izM) * dz || 1.0;

      for (let iy = 0; iy < ny; iy++) {
        const iyM = iy > 0 ? iy - 1 : 0;
        const iyP = iy < ny - 1 ? iy + 1 : ny - 1;
        const dyEff = (iyP - iyM) * dy || 1.0;

        for (let ix = 0; ix < nx; ix++) {
          const ixM = ix > 0 ? ix - 1 : 0;
          const ixP = ix < nx - 1 ? ix + 1 : nx - 1;
          const dxEff = (ixP - ixM) * dx || 1.0;

          const idx = this.grid.index(ix, iy, iz);

          const dvy_dx = (this.vy[this.grid.index(ixP, iy, iz)] - this.vy[this.grid.index(ixM, iy, iz)]) / dxEff;
          const dvz_dx = (this.vz[this.grid.index(ixP, iy, iz)] - this.vz[this.grid.index(ixM, iy, iz)]) / dxEff;

          const dvx_dy = (this.vx[this.grid.index(ix, iyP, iz)] - this.vx[this.grid.index(ix, iyM, iz)]) / dyEff;
          const dvz_dy = (this.vz[this.grid.index(ix, iyP, iz)] - this.vz[this.grid.index(ix, iyM, iz)]) / dyEff;

          const dvx_dz = (this.vx[this.grid.index(ix, iy, izP)] - this.vx[this.grid.index(ix, iy, izM)]) / dzEff;
          const dvy_dz = (this.vy[this.grid.index(ix, iy, izP)] - this.vy[this.grid.index(ix, iy, izM)]) / dzEff;

          wx[idx] = dvz_dy - dvy_dz;
          wy[idx] = dvx_dz - dvz_dx;
          wz[idx] = dvy_dx - dvx_dy;
        }
      }
    }

    return new VectorField3D(this.grid, wx, wy, wz, {
      name: `${this.name}_vorticity`,
      unit: 'km/s/(Mpc/h)',
      scaleFactorApplied: true
    });
  }

  /**
   * Computes bulk flow vector in a spherical shell R_inner <= r <= R_outer around center.
   * @param {Array<number>|SupergalacticPosition} [center=[0, 0, 0]]
   * @param {number} [rInner=0.0]
   * @param {number} [rOuter=100.0]
   * @returns {{ bulkVelocity: VelocityVector, speed: number, cellCount: number }}
   */
  computeBulkFlow(center = [0, 0, 0], rInner = 0.0, rOuter = 100.0) {
    const cx = center.x ?? center[0];
    const cy = center.y ?? center[1];
    const cz = center.z ?? center[2];

    const rIn2 = rInner * rInner;
    const rOut2 = rOuter * rOuter;

    let sumVx = 0.0;
    let sumVy = 0.0;
    let sumVz = 0.0;
    let count = 0;

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;

    for (let iz = 0; iz < nz; iz++) {
      const z = this.grid.zMin + iz * this.grid.dz;
      const dz = z - cz;
      for (let iy = 0; iy < ny; iy++) {
        const y = this.grid.yMin + iy * this.grid.dy;
        const dy = y - cy;
        for (let ix = 0; ix < nx; ix++) {
          const x = this.grid.xMin + ix * this.grid.dx;
          const dx = x - cx;

          const r2 = dx * dx + dy * dy + dz * dz;
          if (r2 >= rIn2 && r2 <= rOut2) {
            const idx = this.grid.index(ix, iy, iz);
            sumVx += this.vx[idx];
            sumVy += this.vy[idx];
            sumVz += this.vz[idx];
            count++;
          }
        }
      }
    }

    if (count === 0) {
      return {
        bulkVelocity: new VelocityVector(0, 0, 0),
        speed: 0.0,
        cellCount: 0
      };
    }

    const meanVx = sumVx / count;
    const meanVy = sumVy / count;
    const meanVz = sumVz / count;
    const speed = Math.sqrt(meanVx * meanVx + meanVy * meanVy + meanVz * meanVz);

    return {
      bulkVelocity: new VelocityVector(meanVx, meanVy, meanVz),
      speed,
      cellCount: count
    };
  }

  /**
   * Computes radial infall velocity profile v_rad(r) around a cluster center in concentric radial bins.
   * @param {Array<number>|SupergalacticPosition} center
   * @param {number} [maxRadiusMpcOverH=100.0]
   * @param {number} [numBins=20]
   * @returns {Array<{ rMid: number, vRadMean: number, vTanMean: number, count: number }>}
   */
  computeRadialInfallProfile(center, maxRadiusMpcOverH = 100.0, numBins = 20) {
    const cx = center.x ?? center[0];
    const cy = center.y ?? center[1];
    const cz = center.z ?? center[2];

    const dr = maxRadiusMpcOverH / numBins;
    const binSumsRad = new Float64Array(numBins);
    const binSumsTan = new Float64Array(numBins);
    const binCounts = new Int32Array(numBins);

    const nx = this.grid.nx;
    const ny = this.grid.ny;
    const nz = this.grid.nz;

    for (let iz = 0; iz < nz; iz++) {
      const z = this.grid.zMin + iz * this.grid.dz;
      const dz = z - cz;
      for (let iy = 0; iy < ny; iy++) {
        const y = this.grid.yMin + iy * this.grid.dy;
        const dy = y - cy;
        for (let ix = 0; ix < nx; ix++) {
          const x = this.grid.xMin + ix * this.grid.dx;
          const dx = x - cx;

          const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          if (r > 0.0 && r < maxRadiusMpcOverH) {
            const bin = Math.min(numBins - 1, Math.floor(r / dr));
            const idx = this.grid.index(ix, iy, iz);
            const vx = this.vx[idx];
            const vy = this.vy[idx];
            const vz = this.vz[idx];

            // Unit radial vector
            const urX = dx / r;
            const urY = dy / r;
            const urZ = dz / r;

            const vRad = vx * urX + vy * urY + vz * urZ; // Infall is negative vRad
            const vTotal2 = vx * vx + vy * vy + vz * vz;
            const vTan2 = Math.max(0.0, vTotal2 - vRad * vRad);
            const vTan = Math.sqrt(vTan2);

            binSumsRad[bin] += vRad;
            binSumsTan[bin] += vTan;
            binCounts[bin]++;
          }
        }
      }
    }

    const profile = [];
    for (let b = 0; b < numBins; b++) {
      const cnt = binCounts[b];
      profile.push({
        rMid: (b + 0.5) * dr,
        vRadMean: cnt > 0 ? binSumsRad[b] / cnt : 0.0,
        vTanMean: cnt > 0 ? binSumsTan[b] / cnt : 0.0,
        count: cnt
      });
    }
    return profile;
  }
}
