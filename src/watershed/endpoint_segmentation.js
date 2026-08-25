/**
 * @file endpoint_segmentation.js
 * @description Full implementation of the Cosmic Velocity Watershed Segmentation Methodology.
 * Follows Courtois, Dupuy, Hoffman, Tully et al.:
 * 1. Forward streamline integration from all voxel centers through peculiar velocity field v(x)
 * 2. 3D endpoint density grid accumulation
 * 3. 26-neighbor 3D local maxima detection for attractor sinks
 * 4. Watershed basin labeling and attribution mapping
 */

import { VoxelGrid, createGridMetadata } from './watershed_classifier.js';

/**
 * @typedef {Object} SegmentationOptions
 * @property {number} [maxSteps=1000] - Maximum RK4 integration steps per streamline
 * @property {number} [dt=0.05] - Base integration time step (Gyr or normalized units)
 * @property {boolean} [adaptiveStep=true] - Enable adaptive step-size control
 * @property {number} [velocityTol=1e-3] - Convergence threshold for velocity magnitude (km/s)
 * @property {number} [stagnationTol=1e-4] - Displacement stagnation threshold (Mpc/h)
 * @property {number} [densityThreshold=5.0] - Minimum endpoint density to consider as an attractor peak
 * @property {number} [mergeRadiusMpc=20.0] - Spatial radius to cluster adjacent attractor peaks (Mpc/h)
 * @property {function} [progressCallback] - Optional callback (percent, status)
 */

/**
 * Trilinear velocity sampler for a continuous coordinate [x, y, z] inside the box.
 * Coordinates are in Mpc/h, bounded within [-halfExtent, +halfExtent].
 *
 * @param {Float32Array|Float64Array} velocityGrid - Flattened [nx * ny * nz * 3] velocity field
 * @param {number} nx - Dimension X
 * @param {number} ny - Dimension Y
 * @param {number} nz - Dimension Z
 * @param {number} halfExtent - Box half-extent in Mpc/h
 * @param {number} x - Coordinate X
 * @param {number} y - Coordinate Y
 * @param {number} z - Coordinate Z
 * @returns {[number, number, number]} Interpolated velocity [vx, vy, vz]
 */
export function sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, x, y, z) {
  const boxSize = 2.0 * halfExtent;
  const gx = ((x + halfExtent) / boxSize) * (nx - 1);
  const gy = ((y + halfExtent) / boxSize) * (ny - 1);
  const gz = ((z + halfExtent) / boxSize) * (nz - 1);

  if (gx < 0 || gx > nx - 1 || gy < 0 || gy > ny - 1 || gz < 0 || gz > nz - 1) {
    return [0, 0, 0];
  }

  const ix0 = Math.floor(gx);
  const iy0 = Math.floor(gy);
  const iz0 = Math.floor(gz);

  const ix1 = Math.min(ix0 + 1, nx - 1);
  const iy1 = Math.min(iy0 + 1, ny - 1);
  const iz1 = Math.min(iz0 + 1, nz - 1);

  const tx = gx - ix0;
  const ty = gy - iy0;
  const tz = gz - iz0;

  const invTx = 1.0 - tx;
  const invTy = 1.0 - ty;
  const invTz = 1.0 - tz;

  const w000 = invTx * invTy * invTz;
  const w100 = tx * invTy * invTz;
  const w010 = invTx * ty * invTz;
  const w110 = tx * ty * invTz;
  const w001 = invTx * invTy * tz;
  const w101 = tx * invTy * tz;
  const w011 = invTx * ty * tz;
  const w111 = tx * ty * tz;

  const getIdx = (i, j, k) => ((i * ny + j) * nz + k) * 3;

  const idx000 = getIdx(ix0, iy0, iz0);
  const idx100 = getIdx(ix1, iy0, iz0);
  const idx010 = getIdx(ix0, iy1, iz0);
  const idx110 = getIdx(ix1, iy1, iz0);
  const idx001 = getIdx(ix0, iy0, iz1);
  const idx101 = getIdx(ix1, iy0, iz1);
  const idx011 = getIdx(ix0, iy1, iz1);
  const idx111 = getIdx(ix1, iy1, iz1);

  const vx = w000 * velocityGrid[idx000] + w100 * velocityGrid[idx100] +
             w010 * velocityGrid[idx010] + w110 * velocityGrid[idx110] +
             w001 * velocityGrid[idx001] + w101 * velocityGrid[idx101] +
             w011 * velocityGrid[idx011] + w111 * velocityGrid[idx111];

  const vy = w000 * velocityGrid[idx000 + 1] + w100 * velocityGrid[idx100 + 1] +
             w010 * velocityGrid[idx010 + 1] + w110 * velocityGrid[idx110 + 1] +
             w001 * velocityGrid[idx001 + 1] + w101 * velocityGrid[idx101 + 1] +
             w011 * velocityGrid[idx011 + 1] + w111 * velocityGrid[idx111 + 1];

  const vz = w000 * velocityGrid[idx000 + 2] + w100 * velocityGrid[idx100 + 2] +
             w010 * velocityGrid[idx010 + 2] + w110 * velocityGrid[idx110 + 2] +
             w001 * velocityGrid[idx001 + 2] + w101 * velocityGrid[idx101 + 2] +
             w011 * velocityGrid[idx011 + 2] + w111 * velocityGrid[idx111 + 2];

  return [vx, vy, vz];
}

/**
 * Integrates a single streamline from starting position x0 using 4th-Order Runge-Kutta (RK4).
 *
 * @param {Float32Array} velocityGrid
 * @param {number} nx
 * @param {number} ny
 * @param {number} nz
 * @param {number} halfExtent
 * @param {number} x0 - Initial X (Mpc/h)
 * @param {number} y0 - Initial Y (Mpc/h)
 * @param {number} z0 - Initial Z (Mpc/h)
 * @param {SegmentationOptions} options
 * @returns {{ endPoint: [number, number, number], steps: number, status: 'converged'|'boundary'|'stagnant'|'max_steps' }}
 */
export function traceStreamlineRK4(velocityGrid, nx, ny, nz, halfExtent, x0, y0, z0, options = {}) {
  const maxSteps = options.maxSteps ?? 1000;
  const dt = options.dt ?? 0.05;
  const vTol = options.velocityTol ?? 1e-3;
  const stagTol = options.stagnationTol ?? 1e-4;

  let x = x0, y = y0, z = z0;
  let status = 'max_steps';
  let steps = 0;

  for (let s = 0; s < maxSteps; s++) {
    steps++;

    // Boundary check
    if (Math.abs(x) > halfExtent || Math.abs(y) > halfExtent || Math.abs(z) > halfExtent) {
      status = 'boundary';
      break;
    }

    // k1 = v(x)
    const [k1x, k1y, k1z] = sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, x, y, z);
    const speed1 = Math.hypot(k1x, k1y, k1z);
    if (speed1 < vTol) {
      status = 'converged';
      break;
    }

    // k2 = v(x + dt/2 * k1)
    const xMid1 = x + 0.5 * dt * k1x;
    const yMid1 = y + 0.5 * dt * k1y;
    const zMid1 = z + 0.5 * dt * k1z;
    const [k2x, k2y, k2z] = sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, xMid1, yMid1, zMid1);

    // k3 = v(x + dt/2 * k2)
    const xMid2 = x + 0.5 * dt * k2x;
    const yMid2 = y + 0.5 * dt * k2y;
    const zMid2 = z + 0.5 * dt * k2z;
    const [k3x, k3y, k3z] = sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, xMid2, yMid2, zMid2);

    // k4 = v(x + dt * k3)
    const xEnd = x + dt * k3x;
    const yEnd = y + dt * k3y;
    const zEnd = z + dt * k3z;
    const [k4x, k4y, k4z] = sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, xEnd, yEnd, zEnd);

    // RK4 step
    const dx = (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    const dy = (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    const dz = (dt / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);

    const stepDist = Math.hypot(dx, dy, dz);
    if (stepDist < stagTol) {
      status = 'stagnant';
      break;
    }

    x += dx;
    y += dy;
    z += dz;
  }

  return {
    endPoint: [x, y, z],
    steps,
    status
  };
}

/**
 * Full Watershed Segmentation Pipeline from Velocity Field.
 *
 * @param {Float32Array} velocityGrid - Velocity field buffer [nx * ny * nz * 3] in km/s or Mpc/h
 * @param {number} nx - Dimension X (e.g. 128 or 64)
 * @param {number} ny - Dimension Y
 * @param {number} nz - Dimension Z
 * @param {number} halfExtent - Box half-extent in Mpc/h (e.g. 500.0)
 * @param {SegmentationOptions} [options={}] - Pipeline options
 * @returns {{
 *   basinGrid: VoxelGrid,
 *   endpointGrid: Float32Array,
 *   densityGrid: Float32Array,
 *   attractors: Array<{ id: number, name: string, positionMpc: [number, number, number], gridCoords: [number, number, number], peakDensity: number, capturedVoxels: number }>,
 *   telemetry: { totalVoxels: number, convergedCount: number, boundaryCount: number, executionTimeMs: number }
 * }}
 */
export function computeWatershedSegmentation(velocityGrid, nx, ny, nz, halfExtent, options = {}) {
  const startTime = performance.now();
  const totalVoxels = nx * ny * nz;
  const metadata = createGridMetadata({ nx, ny, nz, halfExtentMpc: halfExtent });
  const basinVoxelGrid = new VoxelGrid(nx, ny, nz, null, metadata, 'Uint16Array');

  const endpointArray = new Float32Array(totalVoxels * 3);
  const densityGrid = new Float32Array(totalVoxels);
  const statusArray = new Uint8Array(totalVoxels); // 0=other, 1=converged, 2=boundary, 3=stagnant

  let convergedCount = 0;
  let boundaryCount = 0;

  // Step 1: Trace streamlines from every voxel centroid
  for (let idx = 0; idx < totalVoxels; idx++) {
    const [ix, iy, iz] = basinVoxelGrid.getGridCoords(idx);
    const { x, y, z } = basinVoxelGrid.gridToWorld(ix, iy, iz, 'Mpc');

    const trace = traceStreamlineRK4(velocityGrid, nx, ny, nz, halfExtent, x, y, z, options);
    const endX = trace.endPoint[0];
    const endY = trace.endPoint[1];
    const endZ = trace.endPoint[2];

    endpointArray[idx * 3] = endX;
    endpointArray[idx * 3 + 1] = endY;
    endpointArray[idx * 3 + 2] = endZ;

    if (trace.status === 'converged' || trace.status === 'stagnant') {
      convergedCount++;
      statusArray[idx] = 1;

      // Accumulate into density grid using nearest grid voxel
      const gPos = basinVoxelGrid.worldToGrid(endX, endY, endZ, 'Mpc');
      if (gPos.inside) {
        const destIdx = basinVoxelGrid.getIndex(gPos.ix, gPos.iy, gPos.iz);
        densityGrid[destIdx] += 1.0;
      }
    } else if (trace.status === 'boundary') {
      boundaryCount++;
      statusArray[idx] = 2;
    }
  }

  // Step 2: 26-neighbor Local Maxima Detection on Endpoint Density Grid
  const candidateAttractors = [];
  const densityThresh = options.densityThreshold ?? 5.0;

  for (let ix = 1; ix < nx - 1; ix++) {
    for (let iy = 1; iy < ny - 1; iy++) {
      for (let iz = 1; iz < nz - 1; iz++) {
        const centerIdx = basinVoxelGrid.getIndex(ix, iy, iz);
        const centerVal = densityGrid[centerIdx];
        if (centerVal < densityThresh) continue;

        let isMax = true;
        for (let dx = -1; dx <= 1 && isMax; dx++) {
          for (let dy = -1; dy <= 1 && isMax; dy++) {
            for (let dz = -1; dz <= 1 && dz <= 1; dz++) {
              if (dx === 0 && dy === 0 && dz === 0) continue;
              const nIdx = basinVoxelGrid.getIndex(ix + dx, iy + dy, iz + dz);
              if (densityGrid[nIdx] > centerVal) {
                isMax = false;
                break;
              }
            }
          }
        }

        if (isMax) {
          const worldPos = basinVoxelGrid.gridToWorld(ix, iy, iz, 'Mpc');
          candidateAttractors.push({
            ix,
            iy,
            iz,
            worldPos: [worldPos.x, worldPos.y, worldPos.z],
            peakDensity: centerVal
          });
        }
      }
    }
  }

  // Step 3: Cluster/Merge close local maxima within mergeRadius
  const mergeRadius = options.mergeRadiusMpc ?? 20.0;
  const attractors = [];
  candidateAttractors.sort((a, b) => b.peakDensity - a.peakDensity);

  for (const cand of candidateAttractors) {
    let merged = false;
    for (const attr of attractors) {
      const dist = Math.hypot(
        cand.worldPos[0] - attr.positionMpc[0],
        cand.worldPos[1] - attr.positionMpc[1],
        cand.worldPos[2] - attr.positionMpc[2]
      );
      if (dist < mergeRadius) {
        // Merge into stronger peak (weighted centroid)
        const totalDensity = attr.peakDensity + cand.peakDensity;
        attr.positionMpc[0] = (attr.positionMpc[0] * attr.peakDensity + cand.worldPos[0] * cand.peakDensity) / totalDensity;
        attr.positionMpc[1] = (attr.positionMpc[1] * attr.peakDensity + cand.worldPos[1] * cand.peakDensity) / totalDensity;
        attr.positionMpc[2] = (attr.positionMpc[2] * attr.peakDensity + cand.worldPos[2] * cand.peakDensity) / totalDensity;
        attr.peakDensity = totalDensity;
        merged = true;
        break;
      }
    }

    if (!merged) {
      const basinId = attractors.length + 1;
      const gCoord = basinVoxelGrid.worldToGrid(cand.worldPos[0], cand.worldPos[1], cand.worldPos[2], 'Mpc');
      attractors.push({
        id: basinId,
        name: `Basin-${basinId}`,
        positionMpc: [...cand.worldPos],
        gridCoords: [gCoord.ix, gCoord.iy, gCoord.iz],
        peakDensity: cand.peakDensity,
        capturedVoxels: 0
      });
    }
  }

  // Step 4: Basin Attribution (assign each voxel to closest attractor basin if endpoint is within domain)
  for (let idx = 0; idx < totalVoxels; idx++) {
    if (statusArray[idx] !== 1 || attractors.length === 0) {
      basinVoxelGrid.data[idx] = 0; // Unsegmented / Boundary
      continue;
    }

    const endX = endpointArray[idx * 3];
    const endY = endpointArray[idx * 3 + 1];
    const endZ = endpointArray[idx * 3 + 2];

    let bestDistSq = Infinity;
    let bestBasinId = 0;

    for (const attr of attractors) {
      const dx = endX - attr.positionMpc[0];
      const dy = endY - attr.positionMpc[1];
      const dz = endZ - attr.positionMpc[2];
      const distSq = dx * dx + dy * dy + dz * dz;
      if (distSq < bestDistSq) {
        bestDistSq = distSq;
        bestBasinId = attr.id;
      }
    }

    basinVoxelGrid.data[idx] = bestBasinId;
    if (bestBasinId > 0) {
      attractors[bestBasinId - 1].capturedVoxels++;
    }
  }

  const executionTimeMs = performance.now() - startTime;

  return {
    basinGrid: basinVoxelGrid,
    endpointGrid: endpointArray,
    densityGrid,
    attractors,
    telemetry: {
      totalVoxels,
      convergedCount,
      boundaryCount,
      executionTimeMs
    }
  };
}
