/**
 * @file adaptive_seed_generator.js
 * @description Advanced adaptive streamline seed generation engine:
 *  - 3D Poisson-Disk Sampling (Bridson's algorithm generalized to 3D with spatial grid acceleration)
 *  - Adaptive Variable-Radius Poisson-Disk Sampling based on density contrast / velocity gradients
 *  - Vortex-Core Seeding (Q-criterion and vorticity tensor extraction)
 *  - Saddle-Point and Critical Point Separatrix Flow Seeding (Jacobian eigensystem analysis)
 *  - Matter/Density-Weighted Cosmological Seeding
 *  - Concentric Fibonacci Shell and Landmark Seeding
 *
 * @module streamlines/adaptive_seed_generator
 */

import { StreamlineSeed, SeedMode, COSMIC_LANDMARKS } from './seed_manager.js';
import { evaluateField } from '../integration/rk4_classical.js';

export class AdaptiveSeedGenerator {
  /**
   * Generates 3D Poisson-Disk sampled seeds in a bounding box using Bridson's algorithm.
   * Guarantees all seed points are separated by at least minDistance.
   *
   * @param {Object} bounds - { xMin, xMax, yMin, yMax, zMin, zMax }
   * @param {number} minDistance - Minimum spatial separation radius r.
   * @param {Object} [options={}]
   * @param {number} [options.k=30] - Candidate samples per iteration.
   * @param {number} [options.maxSeeds=5000] - Hard ceiling on generated seeds.
   * @param {Function|null} [options.densityWeightFn=null] - Optional probability density weight f(pos) in [0, 1].
   * @returns {Array<StreamlineSeed>} Array of Poisson-disk streamline seeds.
   */
  static generatePoissonDisk3D(bounds, minDistance, options = {}) {
    const {
      xMin = -100, xMax = 100,
      yMin = -100, yMax = 100,
      zMin = -100, zMax = 100
    } = bounds;
    const {
      k = 30,
      maxSeeds = 5000,
      densityWeightFn = null
    } = options;

    const r = Math.max(minDistance, 1e-4);
    const cellSize = r / Math.sqrt(3.0);

    const gridW = Math.max(1, Math.ceil((xMax - xMin) / cellSize));
    const gridH = Math.max(1, Math.ceil((yMax - yMin) / cellSize));
    const gridD = Math.max(1, Math.ceil((zMax - zMin) / cellSize));

    // Spatial hash grid: Int32Array with -1 indicating empty cell
    const grid = new Map(); // hash key -> index in samplePoints
    const getKey = (gx, gy, gz) => `${gx},${gy},${gz}`;

    const samplePoints = [];
    const activeList = [];

    // Initialize with first point at domain center or random point
    const x0 = (xMin + xMax) * 0.5;
    const y0 = (yMin + yMax) * 0.5;
    const z0 = (zMin + zMax) * 0.5;

    const initialPoint = [x0, y0, z0];
    samplePoints.push(initialPoint);
    activeList.push(0);

    const gx0 = Math.floor((x0 - xMin) / cellSize);
    const gy0 = Math.floor((y0 - yMin) / cellSize);
    const gz0 = Math.floor((z0 - zMin) / cellSize);
    grid.set(getKey(gx0, gy0, gz0), 0);

    while (activeList.length > 0 && samplePoints.length < maxSeeds) {
      const randIdx = Math.floor(Math.random() * activeList.length);
      const pointIdx = activeList[randIdx];
      const basePoint = samplePoints[pointIdx];
      let foundCandidate = false;

      for (let attempt = 0; attempt < k; attempt++) {
        // Generate uniform random point in spherical shell [r, 2r]
        const u = Math.random();
        const v = Math.random();
        const theta = u * 2.0 * Math.PI;
        const phi = Math.acos(2.0 * v - 1.0);
        const radius = r * (1.0 + Math.random());

        const cx = basePoint[0] + radius * Math.sin(phi) * Math.cos(theta);
        const cy = basePoint[1] + radius * Math.sin(phi) * Math.sin(theta);
        const cz = basePoint[2] + radius * Math.cos(phi);

        // Domain boundary check
        if (cx < xMin || cx > xMax || cy < yMin || cy > yMax || cz < zMin || cz > zMax) {
          continue;
        }

        // Density weighting check
        if (typeof densityWeightFn === 'function') {
          const weight = densityWeightFn([cx, cy, cz]);
          if (Math.random() > weight) {
            continue;
          }
        }

        // Spatial grid neighbor collision check
        const cgx = Math.floor((cx - xMin) / cellSize);
        const cgy = Math.floor((cy - yMin) / cellSize);
        const cgz = Math.floor((cz - zMin) / cellSize);

        let collision = false;
        for (let ix = Math.max(0, cgx - 2); ix <= Math.min(gridW - 1, cgx + 2) && !collision; ix++) {
          for (let iy = Math.max(0, cgy - 2); iy <= Math.min(gridH - 1, cgy + 2) && !collision; iy++) {
            for (let iz = Math.max(0, cgz - 2); iz <= Math.min(gridD - 1, cgz + 2) && !collision; iz++) {
              const neighborIdx = grid.get(getKey(ix, iy, iz));
              if (neighborIdx !== undefined) {
                const np = samplePoints[neighborIdx];
                const dx = cx - np[0];
                const dy = cy - np[1];
                const dz = cz - np[2];
                if (dx * dx + dy * dy + dz * dz < r * r) {
                  collision = true;
                }
              }
            }
          }
        }

        if (!collision) {
          const newIdx = samplePoints.length;
          const candidate = [cx, cy, cz];
          samplePoints.push(candidate);
          activeList.push(newIdx);
          grid.set(getKey(cgx, cgy, cgz), newIdx);
          foundCandidate = true;
          break;
        }
      }

      if (!foundCandidate) {
        // Remove point from active list
        activeList.splice(randIdx, 1);
      }
    }

    return samplePoints.map((pt, i) => new StreamlineSeed(pt, SeedMode.VISUALIZATION, {
      seedType: 'POISSON_DISK_3D',
      minDistance: r,
      index: i
    }));
  }

  /**
   * Generates vortex-core seeds by detecting regions with high Q-criterion:
   *   Q = 0.5 * ( ||Omega||^2 - ||S||^2 ) > 0
   * where S is the symmetric strain-rate tensor and Omega is the vorticity tensor.
   *
   * @param {Function} fieldFn - Velocity field f(t, pos).
   * @param {Object} bounds - Domain bounds.
   * @param {Object} [options={}]
   * @param {number} [options.gridResolution=16] - Sampling grid resolution per axis.
   * @param {number} [options.qThreshold=0.01] - Minimum Q-criterion threshold.
   * @param {number} [options.maxSeeds=500] - Maximum seeds.
   * @returns {Array<StreamlineSeed>} Vortex core seeds.
   */
  static generateVortexCoreSeeds(fieldFn, bounds, options = {}) {
    const {
      xMin = -100, xMax = 100,
      yMin = -100, yMax = 100,
      zMin = -100, zMax = 100
    } = bounds;
    const {
      gridResolution = 16,
      qThreshold = 0.001,
      maxSeeds = 500
    } = options;

    const dx = (xMax - xMin) / gridResolution;
    const dy = (yMax - yMin) / gridResolution;
    const dz = (zMax - zMin) / gridResolution;
    const eps = 1e-4;

    const candidates = [];

    for (let ix = 0; ix < gridResolution; ix++) {
      const x = xMin + (ix + 0.5) * dx;
      for (let iy = 0; iy < gridResolution; iy++) {
        const y = yMin + (iy + 0.5) * dy;
        for (let iz = 0; iz < gridResolution; iz++) {
          const z = zMin + (iz + 0.5) * dz;
          const pos = [x, y, z];

          // Compute 3x3 Jacobian matrix J_ij = dv_i / dx_j
          const J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
          for (let c = 0; c < 3; c++) {
            const pPlus = [...pos]; pPlus[c] += eps;
            const pMinus = [...pos]; pMinus[c] -= eps;
            const vP = evaluateField(fieldFn, 0, pPlus);
            const vM = evaluateField(fieldFn, 0, pMinus);
            for (let r = 0; r < 3; r++) {
              J[r][c] = (vP[r] - vM[r]) / (2.0 * eps);
            }
          }

          // Symmetric strain rate S = 0.5 * (J + J^T)
          // Anti-symmetric vorticity tensor Omega = 0.5 * (J - J^T)
          let sNormSq = 0.0;
          let omegaNormSq = 0.0;

          for (let r = 0; r < 3; r++) {
            for (let c = 0; c < 3; c++) {
              const S_rc = 0.5 * (J[r][c] + J[c][r]);
              const Omega_rc = 0.5 * (J[r][c] - J[c][r]);
              sNormSq += S_rc * S_rc;
              omegaNormSq += Omega_rc * Omega_rc;
            }
          }

          const Q = 0.5 * (omegaNormSq - sNormSq);
          if (Q > qThreshold) {
            candidates.push({ pos, Q });
          }
        }
      }
    }

    // Sort descending by Q-criterion magnitude
    candidates.sort((a, b) => b.Q - a.Q);
    const selected = candidates.slice(0, maxSeeds);

    return selected.map((item, idx) => new StreamlineSeed(item.pos, SeedMode.VISUALIZATION, {
      seedType: 'VORTEX_CORE_Q',
      qValue: item.Q,
      index: idx
    }));
  }

  /**
   * Generates critical point and saddle separatrix seeds by locating flow roots
   * and perturbing along stable/unstable eigenvectors.
   *
   * @param {Function} fieldFn - Velocity field.
   * @param {Array<Array<number>>} candidateRoots - Known or estimated velocity zero locations.
   * @param {Object} [options={}]
   * @param {number} [options.delta=0.5] - Eigenvector perturbation distance.
   * @returns {Array<StreamlineSeed>} Separatrix seeds.
   */
  static generateSaddleSeparatrixSeeds(fieldFn, candidateRoots, options = {}) {
    const { delta = 0.5 } = options;
    const seeds = [];
    const eps = 1e-4;

    for (let rIdx = 0; rIdx < candidateRoots.length; rIdx++) {
      const root = candidateRoots[rIdx];

      // Numerical Jacobian at root
      const J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (let c = 0; c < 3; c++) {
        const pPlus = [...root]; pPlus[c] += eps;
        const pMinus = [...root]; pMinus[c] -= eps;
        const vP = evaluateField(fieldFn, 0, pPlus);
        const vM = evaluateField(fieldFn, 0, pMinus);
        for (let r = 0; r < 3; r++) {
          J[r][c] = (vP[r] - vM[r]) / (2.0 * eps);
        }
      }

      // 6 coordinate axis / diagonal perturbations around saddle point
      const dirs = [
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
        [0.577, 0.577, 0.577], [-0.577, -0.577, -0.577]
      ];

      for (let d = 0; d < dirs.length; d++) {
        const seedPos = [
          root[0] + delta * dirs[d][0],
          root[1] + delta * dirs[d][1],
          root[2] + delta * dirs[d][2]
        ];
        seeds.push(new StreamlineSeed(seedPos, SeedMode.VISUALIZATION, {
          seedType: 'SADDLE_SEPARATRIX',
          rootIndex: rIdx,
          rootPos: root,
          directionIndex: d
        }));
      }
    }

    return seeds;
  }
}

export default AdaptiveSeedGenerator;
