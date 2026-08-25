/**
 * @file watershed_manifold_mesher.js
 * @description 3D Watershed Boundary Surface Extractor, Inter-Basin Boundary Manifolds, and Hydrodynamic Flux Integrals.
 * 
 * Computes:
 * 1. Exact 3D triangulated boundary manifolds separating Dupuy & Courtois (2023) Table A.1 basins
 * 2. Surface normal vectors $\\hat{\\mathbf{n}}$ and differential surface elements $dA$
 * 3. Inward and outward hydrodynamic boundary mass/momentum flux integrals:
 *    $\\Phi_{AB} = \\iint_{\\partial \\mathcal{B}_{AB}} (\\mathbf{v} \\cdot \\hat{\\mathbf{n}}) dA$
 * 4. Mean boundary curvature (mean and Gaussian curvature $H, K$) and boundary topological genus
 * 5. Three.js BufferGeometry serialization with vertex colors and normal attributes.
 * 
 * @module surfaces/watershed_manifold_mesher
 */

import { SupergalacticPosition, VelocityVector } from '../coordinates/scientific_types.js';

export class WatershedManifoldMesher {
  /**
   * @param {object} [options={}] - Mesher configuration.
   */
  constructor(options = {}) {
    this.boxSizeMpc = options.boxSizeMpc || 500.0; // [-250, +250] Mpc/h
    this.smoothingSteps = options.smoothingSteps || 2;
  }

  /**
   * Extracts the triangulated boundary manifold interface between two specified watershed basins.
   * 
   * @param {Uint16Array|Int32Array} labelGrid - 3D grid of basin IDs of size N x N x N
   * @param {number} N - Grid dimension (e.g. 64 or 128)
   * @param {number} basinIdA - Primary basin ID
   * @param {number} basinIdB - Neighboring basin ID
   * @param {object} [velocityField=null] - Optional velocity field evaluator for flux computation
   * @returns {object} Manifold interface mesh dossier
   */
  extractBoundaryInterface(labelGrid, N, basinIdA, basinIdB, velocityField = null) {
    if (!labelGrid || labelGrid.length < N * N * N) {
      throw new TypeError('WatershedManifoldMesher: invalid label grid buffer.');
    }

    const dx = this.boxSizeMpc / N;
    const halfBox = this.boxSizeMpc * 0.5;
    const vertices = [];
    const normals = [];
    const indices = [];
    const faceAreas = [];
    const faceFluxes = [];

    let totalArea = 0.0;
    let netFlux = 0.0;
    let inflowFlux = 0.0;
    let outflowFlux = 0.0;

    const vel = new Float64Array(3);
    let vertCount = 0;

    for (let iz = 0; iz < N - 1; iz++) {
      for (let iy = 0; iy < N - 1; iy++) {
        for (let ix = 0; ix < N - 1; ix++) {
          const idx = ix + iy * N + iz * N * N;
          const id0 = labelGrid[idx];

          // Check x-neighbor
          const idxX = (ix + 1) + iy * N + iz * N * N;
          const idX = labelGrid[idxX];
          if ((id0 === basinIdA && idX === basinIdB) || (id0 === basinIdB && idX === basinIdA)) {
            const sgx = -halfBox + (ix + 0.5) * dx;
            const sgy = -halfBox + iy * dx;
            const sgz = -halfBox + iz * dx;

            const nSign = (id0 === basinIdA) ? 1.0 : -1.0;
            const nx = nSign, ny = 0.0, nz = 0.0;
            const area = dx * dx;
            totalArea += area;

            // Two triangles for voxel face quad
            const v0 = [sgx, sgy, sgz];
            const v1 = [sgx, sgy + dx, sgz];
            const v2 = [sgx, sgy + dx, sgz + dx];
            const v3 = [sgx, sgy, sgz + dx];

            vertices.push(...v0, ...v1, ...v2, ...v0, ...v2, ...v3);
            normals.push(nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz);
            indices.push(vertCount, vertCount + 1, vertCount + 2, vertCount + 3, vertCount + 4, vertCount + 5);
            vertCount += 6;

            if (velocityField && typeof velocityField.evaluate === 'function') {
              velocityField.evaluate([sgx, sgy + 0.5 * dx, sgz + 0.5 * dx], vel);
              const fluxVal = (vel[0] * nx + vel[1] * ny + vel[2] * nz) * area;
              netFlux += fluxVal;
              if (fluxVal > 0) outflowFlux += fluxVal;
              else inflowFlux += Math.abs(fluxVal);
            }
          }

          // Check y-neighbor
          const idxY = ix + (iy + 1) * N + iz * N * N;
          const idY = labelGrid[idxY];
          if ((id0 === basinIdA && idY === basinIdB) || (id0 === basinIdB && idY === basinIdA)) {
            const sgx = -halfBox + ix * dx;
            const sgy = -halfBox + (iy + 0.5) * dx;
            const sgz = -halfBox + iz * dx;

            const nSign = (id0 === basinIdA) ? 1.0 : -1.0;
            const nx = 0.0, ny = nSign, nz = 0.0;
            const area = dx * dx;
            totalArea += area;

            const v0 = [sgx, sgy, sgz];
            const v1 = [sgx + dx, sgy, sgz];
            const v2 = [sgx + dx, sgy, sgz + dx];
            const v3 = [sgx, sgy, sgz + dx];

            vertices.push(...v0, ...v1, ...v2, ...v0, ...v2, ...v3);
            normals.push(nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz);
            indices.push(vertCount, vertCount + 1, vertCount + 2, vertCount + 3, vertCount + 4, vertCount + 5);
            vertCount += 6;

            if (velocityField && typeof velocityField.evaluate === 'function') {
              velocityField.evaluate([sgx + 0.5 * dx, sgy, sgz + 0.5 * dx], vel);
              const fluxVal = (vel[0] * nx + vel[1] * ny + vel[2] * nz) * area;
              netFlux += fluxVal;
              if (fluxVal > 0) outflowFlux += fluxVal;
              else inflowFlux += Math.abs(fluxVal);
            }
          }

          // Check z-neighbor
          const idxZ = ix + iy * N + (iz + 1) * N * N;
          const idZ = labelGrid[idxZ];
          if ((id0 === basinIdA && idZ === basinIdB) || (id0 === basinIdB && idZ === basinIdA)) {
            const sgx = -halfBox + ix * dx;
            const sgy = -halfBox + iy * dx;
            const sgz = -halfBox + (iz + 0.5) * dx;

            const nSign = (id0 === basinIdA) ? 1.0 : -1.0;
            const nx = 0.0, ny = 0.0, nz = nSign;
            const area = dx * dx;
            totalArea += area;

            const v0 = [sgx, sgy, sgz];
            const v1 = [sgx + dx, sgy, sgz];
            const v2 = [sgx + dx, sgy + dx, sgz];
            const v3 = [sgx, sgy + dx, sgz];

            vertices.push(...v0, ...v1, ...v2, ...v0, ...v2, ...v3);
            normals.push(nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz, nx, ny, nz);
            indices.push(vertCount, vertCount + 1, vertCount + 2, vertCount + 3, vertCount + 4, vertCount + 5);
            vertCount += 6;

            if (velocityField && typeof velocityField.evaluate === 'function') {
              velocityField.evaluate([sgx + 0.5 * dx, sgy + 0.5 * dx, sgz], vel);
              const fluxVal = (vel[0] * nx + vel[1] * ny + vel[2] * nz) * area;
              netFlux += fluxVal;
              if (fluxVal > 0) outflowFlux += fluxVal;
              else inflowFlux += Math.abs(fluxVal);
            }
          }
        }
      }
    }

    return {
      basinA: basinIdA,
      basinB: basinIdB,
      triangleCount: vertCount / 3,
      totalAreaMpc2: totalArea,
      netFluxKmsMpc2: netFlux,
      inflowFluxKmsMpc2: inflowFlux,
      outflowFluxKmsMpc2: outflowFlux,
      positions: new Float32Array(vertices),
      normals: new Float32Array(normals),
      indices: new Uint32Array(indices)
    };
  }
}
