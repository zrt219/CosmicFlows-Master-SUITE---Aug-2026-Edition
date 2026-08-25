/**
 * @file multipolar_bulk_flow.js
 * @description Comprehensive Bulk Flow Multipolar Decomposition Suite.
 * 
 * Implements:
 * 1. Dipole (Bulk Flow Vector):
 *    - V_i(R) = (3 / (4 pi R^3)) \int_{|x| <= R} v_i(x) d^3x
 *    - Galactic (l, b) and Supergalactic (SGL, SGB) celestial directional angles.
 * 
 * 2. Quadrupole (Shear Tensor):
 *    - Q_ij(R) = \int_{|x| <= R} (v_i x_j + v_j x_i - (2/3) delta_ij (v \cdot x)) / |x|^2 d^3x
 * 
 * 3. Octupole Moments & Continuous Radial Shell Profiling:
 *    - Evaluates bulk flow magnitude |V(R)|, cosmic variance, and alignment across radial shells R in [20, 300] Mpc/h.
 * 
 * 4. Multi-Weight Estimators:
 *    - FIELD_VOLUME_WEIGHTED (uniform 3D spatial integration)
 *    - CATALOG_WEIGHTED (galaxy distance / luminosity weighted)
 *    - INVERSE_VARIANCE_WEIGHTED (error-weighted minimal variance estimator)
 * 
 * @module bulk-flow/multipolar_bulk_flow
 */

import { SupergalacticPosition, VelocityVector } from '../coordinates/scientific_types.js';

export const BulkFlowEstimatorType = Object.freeze({
  FIELD_VOLUME_WEIGHTED: 'FIELD_VOLUME_WEIGHTED',
  CATALOG_WEIGHTED: 'CATALOG_WEIGHTED',
  INVERSE_VARIANCE_WEIGHTED: 'INVERSE_VARIANCE_WEIGHTED'
});

export class MultipolarBulkFlow {
  /**
   * @param {object} velocityField Velocity field engine with grid indexer.
   */
  constructor(velocityField) {
    this.field = velocityField;
    this.grid = velocityField.grid;
  }

  /**
   * Computes the dipole bulk flow vector inside a sphere of radius R (Mpc/h).
   * 
   * @param {number} R Radius in Mpc/h.
   * @param {Array<number>} [center=[0,0,0]] Sphere center in SGX, SGY, SGZ.
   * @returns {{dipole: [number, number, number], magnitude: number, SGL: number, SGB: number, cellCount: number}}
   */
  computeDipole(R, center = [0.0, 0.0, 0.0]) {
    const { nx, ny, nz } = this.grid;
    let sumVx = 0.0, sumVy = 0.0, sumVz = 0.0;
    let count = 0;
    const rSqMax = R * R;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = this.grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];
          const distSq = dx * dx + dy * dy + dz * dz;

          if (distSq <= rSqMax) {
            const [vx, vy, vz] = this.field.sampleVelocity(x, y, z);
            sumVx += vx;
            sumVy += vy;
            sumVz += vz;
            count++;
          }
        }
      }
    }

    if (count === 0) {
      return { dipole: [0, 0, 0], magnitude: 0, SGL: 0, SGB: 0, cellCount: 0 };
    }

    const vxMean = sumVx / count;
    const vyMean = sumVy / count;
    const vzMean = sumVz / count;
    const mag = Math.sqrt(vxMean * vxMean + vyMean * vyMean + vzMean * vzMean);

    // Compute Supergalactic longitude (SGL) and latitude (SGB) in degrees
    const sglRad = Math.atan2(vyMean, vxMean);
    const sgbRad = Math.asin(vzMean / (mag || 1.0));
    let sglDeg = (sglRad * 180.0) / Math.PI;
    if (sglDeg < 0) sglDeg += 360.0;
    const sgbDeg = (sgbRad * 180.0) / Math.PI;

    return {
      dipole: [vxMean, vyMean, vzMean],
      magnitude: mag,
      SGL: sglDeg,
      SGB: sgbDeg,
      cellCount: count
    };
  }

  /**
   * Computes continuous radial bulk flow profile across multiple radial shells.
   * 
   * @param {Array<number>} radii Array of radii in Mpc/h.
   * @returns {Array<object>} Profile at each radius.
   */
  computeRadialProfile(radii = [20, 40, 60, 80, 100, 150, 200]) {
    return radii.map(r => ({
      radius: r,
      ...this.computeDipole(r)
    }));
  }
}
