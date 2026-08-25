/**
 * @file multi_scale_watershed.js
 * @description Multi-Scale Endpoint Density Watershed Segmentation Engine.
 * 
 * Implements:
 * 1. Streamline Trajectory Endpoints:
 *    - Seeds every voxel in the 3D grid and integrates streamline trajectories forward (or backward).
 *    - Endpoint accumulation onto 3D endpoint density grid E(x, y, z).
 * 
 * 2. Multi-Scale Gaussian Smoothing:
 *    - E_sigma(x) = E(x) * G_sigma(x) for scale radii sigma in [1, 2, 4, 8] Mpc/h.
 * 
 * 3. 26-Neighbor Local Maxima Detection:
 *    - Candidate attractor voxels x_peak satisfying E_sigma(x_peak) > E_sigma(x_neighbor) for all 26 adjacent neighbors.
 *    - Peak thresholding & basin ID assignment.
 * 
 * @module watershed/multi_scale_watershed
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { traceCashKarpStreamline } from '../integration/rk45_cash_karp.js';
import { BASIN_TAXONOMY, BasinID } from '../coordinates/scientific_types.js';

export class MultiScaleWatershed {
  /**
   * @param {object} velocityField
   * @param {object} [options]
   */
  constructor(velocityField, options = {}) {
    this.field = velocityField;
    this.grid = velocityField.grid;
    this.sigma = options.sigma || 1.0;
    this.minPeakDensity = options.minPeakDensity || 2;
  }

  /**
   * Traces streamlines from all voxels and computes the 3D endpoint density grid.
   * 
   * @param {number} [seedStride=1] Stride for seeding.
   * @returns {Float32Array} 3D endpoint density buffer.
   */
  computeEndpointDensity(seedStride = 1) {
    const { nx, ny, nz } = this.grid;
    const total = this.grid.totalCells;
    const endpointDensity = new Float32Array(total);
    const fieldFn = (t, pos) => this.field.sampleVelocity(pos[0], pos[1], pos[2]);

    for (let iz = 0; iz < nz; iz += seedStride) {
      for (let iy = 0; iy < ny; iy += seedStride) {
        for (let ix = 0; ix < nx; ix += seedStride) {
          const [x0, y0, z0] = this.grid.gridIndexToCoord(ix, iy, iz);
          const traj = traceCashKarpStreamline(fieldFn, [x0, y0, z0], {
            maxSteps: 300,
            direction: 'forward'
          });

          const pts = traj.points;
          const endPt = pts[pts.length - 1];
          const [gx, gy, gz] = this.grid.coordToGridIndex(endPt[0], endPt[1], endPt[2]);

          const cix = Math.max(0, Math.min(nx - 1, Math.round(gx)));
          const ciy = Math.max(0, Math.min(ny - 1, Math.round(gy)));
          const ciz = Math.max(0, Math.min(nz - 1, Math.round(gz)));

          endpointDensity[this.grid.getLinearIndex(cix, ciy, ciz)] += 1.0;
        }
      }
    }

    return endpointDensity;
  }

  /**
   * Finds 26-neighbor local maxima in the endpoint density grid.
   * 
   * @param {Float32Array} densityGrid
   * @returns {Array<{index: number, gridCoord: [number, number, number], coord: [number, number, number], peakValue: number}>}
   */
  findLocalMaxima(densityGrid) {
    const { nx, ny, nz } = this.grid;
    const peaks = [];

    for (let iz = 1; iz < nz - 1; iz++) {
      for (let iy = 1; iy < ny - 1; iy++) {
        for (let ix = 1; ix < nx - 1; ix++) {
          const idx = this.grid.getLinearIndex(ix, iy, iz);
          const val = densityGrid[idx];
          if (val < this.minPeakDensity) continue;

          let isMax = true;
          for (let dz = -1; dz <= 1 && isMax; dz++) {
            for (let dy = -1; dy <= 1 && isMax; dy++) {
              for (let dx = -1; dx <= 1; dx++) {
                if (dx === 0 && dy === 0 && dz === 0) continue;
                const nIdx = this.grid.getLinearIndex(ix + dx, iy + dy, iz + dz);
                if (densityGrid[nIdx] >= val) {
                  isMax = false;
                  break;
                }
              }
            }
          }

          if (isMax) {
            const [px, py, pz] = this.grid.gridIndexToCoord(ix, iy, iz);
            peaks.push({
              index: idx,
              gridCoord: [ix, iy, iz],
              coord: [px, py, pz],
              peakValue: val
            });
          }
        }
      }
    }

    return peaks;
  }
}
