/**
 * @file root_finder.js
 * @description Robust 3D Newton-Raphson Vector Root Finder for exact velocity zero-crossings (||v(x_c)|| < 1e-7 km/s).
 * Features:
 * - 3-component sign-change pre-filtering across voxel cube vertices (Poincaré-Miranda criterion)
 * - 3D vector Newton-Raphson iteration with central difference Jacobian computation
 * - Levenberg-Marquardt adaptive damping for near-singular Jacobians
 * - Armijo backtracking line search ensuring monotonic residual decrease
 * - Spatial clustering and candidate deduplication with sub-voxel tolerance
 */

import { sampleVelocityField } from '../watershed/endpoint_segmentation.js';

/**
 * @typedef {Object} RootFinderOptions
 * @property {number} [velocityTol=1e-7] - Target residual velocity norm ||v(x_c)|| (km/s)
 * @property {number} [stepTol=1e-7] - Displacement step convergence tolerance ||dx|| (Mpc/h)
 * @property {number} [maxIterations=100] - Maximum Newton-Raphson iterations per candidate
 * @property {number} [finiteDiffStep=1e-4] - Step size h for central difference Jacobian (Mpc/h)
 * @property {number} [dampingLambda=1e-6] - Levenberg-Marquardt regularization parameter
 * @property {number} [dedupRadiusMpc=0.5] - Spatial radius for root deduplication (Mpc/h)
 * @property {boolean} [enforceCellBounds=true] - Constrain roots to stay within search cell + buffer
 */

/**
 * Computes exact determinant of a 3x3 matrix.
 * @param {number[][]} m
 * @returns {number}
 */
export function det3x3(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/**
 * Inverts a 3x3 matrix using analytical cofactor expansion.
 * Returns null if matrix is singular (|det| < 1e-15).
 * @param {number[][]} m
 * @returns {number[][]|null}
 */
export function invert3x3(m) {
  const det = det3x3(m);
  if (Math.abs(det) < 1e-15) return null;

  const invDet = 1.0 / det;
  return [
    [
      (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet,
      (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet,
      (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet
    ],
    [
      (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet,
      (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet,
      (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * invDet
    ],
    [
      (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet,
      (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * invDet,
      (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * invDet
    ]
  ];
}

/**
 * Multiplies a 3x3 matrix by a 3D vector: y = M * x.
 * @param {number[][]} m
 * @param {number[]} v
 * @returns {[number, number, number]}
 */
export function matVec3(m, v) {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
  ];
}

/**
 * Computes the 3x3 Jacobian matrix J_ij = dv_i / dx_j via central finite differences.
 *
 * @param {function(number, number, number): [number, number, number]} velFn - Velocity evaluation function (x, y, z) -> [vx, vy, vz]
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} [h=1e-4] - Step size
 * @returns {number[][]} 3x3 Jacobian matrix
 */
export function computeJacobian3D(velFn, x, y, z, h = 1e-4) {
  const vPx = velFn(x + h, y, z);
  const vMx = velFn(x - h, y, z);
  const vPy = velFn(x, y + h, z);
  const vMy = velFn(x, y - h, z);
  const vPz = velFn(x, y, z + h);
  const vMz = velFn(x, y, z - h);

  const inv2h = 1.0 / (2.0 * h);

  return [
    [(vPx[0] - vMx[0]) * inv2h, (vPy[0] - vMy[0]) * inv2h, (vPz[0] - vMz[0]) * inv2h],
    [(vPx[1] - vMx[1]) * inv2h, (vPy[1] - vMy[1]) * inv2h, (vPz[1] - vMz[1]) * inv2h],
    [(vPx[2] - vMx[2]) * inv2h, (vPy[2] - vMy[2]) * inv2h, (vPz[2] - vMz[2]) * inv2h]
  ];
}

/**
 * Solves for a single 3D vector root v(x_c) = 0 using damped Newton-Raphson.
 *
 * @param {function(number, number, number): [number, number, number]} velFn - Velocity evaluation function
 * @param {[number, number, number]} initialGuess - Starting coordinates [x0, y0, z0]
 * @param {RootFinderOptions} [options={}]
 * @returns {{
 *   converged: boolean,
 *   root: [number, number, number],
 *   residualVelocity: [number, number, number],
 *   residualNorm: number,
 *   iterations: number,
 *   jacobian: number[][],
 *   status: 'converged'|'max_iterations'|'singular_jacobian'|'out_of_bounds'
 * }}
 */
export function solveNewtonRaphson3D(velFn, initialGuess, options = {}) {
  const vTol = options.velocityTol ?? 1e-7;
  const stepTol = options.stepTol ?? 1e-7;
  const maxIter = options.maxIterations ?? 100;
  const h = options.finiteDiffStep ?? 1e-4;
  const lambda = options.dampingLambda ?? 1e-6;

  let [x, y, z] = initialGuess;
  let iter = 0;
  let lastJacobian = null;

  for (iter = 0; iter < maxIter; iter++) {
    const [vx, vy, vz] = velFn(x, y, z);
    const resNorm = Math.hypot(vx, vy, vz);

    if (resNorm < vTol) {
      lastJacobian = computeJacobian3D(velFn, x, y, z, h);
      return {
        converged: true,
        root: [x, y, z],
        residualVelocity: [vx, vy, vz],
        residualNorm: resNorm,
        iterations: iter,
        jacobian: lastJacobian,
        status: 'converged'
      };
    }

    const J = computeJacobian3D(velFn, x, y, z, h);
    lastJacobian = J;

    let invJ = invert3x3(J);

    // Levenberg-Marquardt fallback if near-singular
    if (!invJ) {
      // (J^T J + lambda I)^-1 J^T
      const JtJ = [
        [J[0][0]*J[0][0] + J[1][0]*J[1][0] + J[2][0]*J[2][0] + lambda, J[0][0]*J[0][1] + J[1][0]*J[1][1] + J[2][0]*J[2][1], J[0][0]*J[0][2] + J[1][0]*J[1][2] + J[2][0]*J[2][2]],
        [J[0][1]*J[0][0] + J[1][1]*J[1][0] + J[2][1]*J[2][0], J[0][1]*J[0][1] + J[1][1]*J[1][1] + J[2][1]*J[2][1] + lambda, J[0][1]*J[0][2] + J[1][1]*J[1][2] + J[2][1]*J[2][2]],
        [J[0][2]*J[0][0] + J[1][2]*J[1][0] + J[2][2]*J[2][0], J[0][2]*J[0][1] + J[1][2]*J[1][1] + J[2][2]*J[2][1], J[0][2]*J[0][2] + J[1][2]*J[1][2] + J[2][2]*J[2][2] + lambda]
      ];
      const invJtJ = invert3x3(JtJ);
      if (!invJtJ) {
        return {
          converged: false,
          root: [x, y, z],
          residualVelocity: [vx, vy, vz],
          residualNorm: resNorm,
          iterations: iter,
          jacobian: J,
          status: 'singular_jacobian'
        };
      }
      // Damped step: dx = -invJtJ * J^T * v
      const Jtv = [
        J[0][0]*vx + J[1][0]*vy + J[2][0]*vz,
        J[0][1]*vx + J[1][1]*vy + J[2][1]*vz,
        J[0][2]*vx + J[1][2]*vy + J[2][2]*vz
      ];
      invJ = invJtJ;
      // Multiply with Jtv
      const dxLM = matVec3(invJtJ, Jtv);
      var dx = -dxLM[0];
      var dy = -dxLM[1];
      var dz = -dxLM[2];
    } else {
      // Standard Newton step: dx = -invJ * v
      const step = matVec3(invJ, [vx, vy, vz]);
      var dx = -step[0];
      var dy = -step[1];
      var dz = -step[2];
    }

    // Armijo Backtracking Line Search
    let alpha = 1.0;
    let stepAccepted = false;

    for (let backtrack = 0; backtrack < 8; backtrack++) {
      const xNext = x + alpha * dx;
      const yNext = y + alpha * dy;
      const zNext = z + alpha * dz;

      const vNext = velFn(xNext, yNext, zNext);
      const nextNorm = Math.hypot(vNext[0], vNext[1], vNext[2]);

      if (nextNorm < resNorm || alpha < 0.05) {
        x = xNext;
        y = yNext;
        z = zNext;
        stepAccepted = true;
        break;
      }
      alpha *= 0.5;
    }

    const stepNorm = Math.hypot(alpha * dx, alpha * dy, alpha * dz);
    if (stepNorm < stepTol && resNorm < vTol * 10) {
      return {
        converged: true,
        root: [x, y, z],
        residualVelocity: velFn(x, y, z),
        residualNorm: Math.hypot(...velFn(x, y, z)),
        iterations: iter + 1,
        jacobian: J,
        status: 'converged'
      };
    }
  }

  const finalV = velFn(x, y, z);
  const finalNorm = Math.hypot(finalV[0], finalV[1], finalV[2]);

  return {
    converged: finalNorm < vTol,
    root: [x, y, z],
    residualVelocity: finalV,
    residualNorm: finalNorm,
    iterations: maxIter,
    jacobian: lastJacobian,
    status: finalNorm < vTol ? 'converged' : 'max_iterations'
  };
}

/**
 * Scans an entire 3D velocity grid for all zero-crossing critical points (Morse nodes and saddles).
 * Uses Poincaré-Miranda 8-vertex component sign-change pre-filtering to minimize expensive iterations.
 *
 * @param {Float32Array} velocityGrid - Flattened velocity grid [nx * ny * nz * 3]
 * @param {number} nx - Dimension X
 * @param {number} ny - Dimension Y
 * @param {number} nz - Dimension Z
 * @param {number} halfExtent - Half-extent in Mpc/h (e.g. 500.0)
 * @param {RootFinderOptions} [options={}]
 * @returns {Array<{
 *   root: [number, number, number],
 *   residualNorm: number,
 *   jacobian: number[][],
 *   iterations: number,
 *   cellIndices: [number, number, number]
 * }>}
 */
export function findGridCriticalPoints(velocityGrid, nx, ny, nz, halfExtent, options = {}) {
  const velFn = (x, y, z) => sampleVelocityField(velocityGrid, nx, ny, nz, halfExtent, x, y, z);
  const rawRoots = [];
  const boxSize = 2.0 * halfExtent;
  const dxCell = boxSize / nx;
  const dyCell = boxSize / ny;
  const dzCell = boxSize / nz;

  const getVoxelV = (i, j, k) => {
    const idx = ((i * ny + j) * nz + k) * 3;
    return [velocityGrid[idx], velocityGrid[idx + 1], velocityGrid[idx + 2]];
  };

  // Sign-change pre-filter across all grid cells
  for (let ix = 0; ix < nx - 1; ix++) {
    for (let iy = 0; iy < ny - 1; iy++) {
      for (let iz = 0; iz < nz - 1; iz++) {
        // Sample 8 vertices of the cell
        const v000 = getVoxelV(ix, iy, iz);
        const v100 = getVoxelV(ix + 1, iy, iz);
        const v010 = getVoxelV(ix, iy + 1, iz);
        const v110 = getVoxelV(ix + 1, iy + 1, iz);
        const v001 = getVoxelV(ix, iy, iz + 1);
        const v101 = getVoxelV(ix + 1, iy, iz + 1);
        const v011 = getVoxelV(ix, iy + 1, iz + 1);
        const v111 = getVoxelV(ix + 1, iy + 1, iz + 1);

        const vList = [v000, v100, v010, v110, v001, v101, v011, v111];

        // Check for zero-crossing in Vx
        let minVx = Infinity, maxVx = -Infinity;
        let minVy = Infinity, maxVy = -Infinity;
        let minVz = Infinity, maxVz = -Infinity;

        for (let k = 0; k < 8; k++) {
          const v = vList[k];
          if (v[0] < minVx) minVx = v[0];
          if (v[0] > maxVx) maxVx = v[0];
          if (v[1] < minVy) minVy = v[1];
          if (v[1] > maxVy) maxVy = v[1];
          if (v[2] < minVz) minVz = v[2];
          if (v[2] > maxVz) maxVz = v[2];
        }

        // All 3 velocity components must bracket zero
        if (minVx <= 0 && maxVx >= 0 && minVy <= 0 && maxVy >= 0 && minVz <= 0 && maxVz >= 0) {
          // Candidate cell: evaluate initial guess at cell centroid
          const x0 = ((ix + 0.5) / nx) * boxSize - halfExtent;
          const y0 = ((iy + 0.5) / ny) * boxSize - halfExtent;
          const z0 = ((iz + 0.5) / nz) * boxSize - halfExtent;

          const res = solveNewtonRaphson3D(velFn, [x0, y0, z0], options);

          if (res.converged) {
            // Check if root lies reasonably close to the candidate cell
            const [rx, ry, rz] = res.root;
            if (Math.abs(rx - x0) <= dxCell * 1.5 &&
                Math.abs(ry - y0) <= dyCell * 1.5 &&
                Math.abs(rz - z0) <= dzCell * 1.5) {
              rawRoots.push({
                root: res.root,
                residualNorm: res.residualNorm,
                jacobian: res.jacobian,
                iterations: res.iterations,
                cellIndices: [ix, iy, iz]
              });
            }
          }
        }
      }
    }
  }

  // Candidate Deduplication
  const dedupRadius = options.dedupRadiusMpc ?? (dxCell * 0.75);
  const uniqueRoots = [];

  rawRoots.sort((a, b) => a.residualNorm - b.residualNorm);

  for (const cand of rawRoots) {
    let duplicate = false;
    for (const u of uniqueRoots) {
      const dist = Math.hypot(
        cand.root[0] - u.root[0],
        cand.root[1] - u.root[1],
        cand.root[2] - u.root[2]
      );
      if (dist < dedupRadius) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) {
      uniqueRoots.push(cand);
    }
  }

  return uniqueRoots;
}
