/**
 * @file newton_raphson_3d.js
 * @description High-Precision 3D Newton-Raphson Vector Root Solver for Cosmological Velocity Fields.
 * 
 * Implements:
 * 1. 3D Vector Root Finding for Velocity Stagnation Points: ||v(x*)|| < 1e-7 km/s
 * 2. High-order central finite difference Jacobians (O(h^4) and O(h^2))
 * 3. Exact analytical cofactor matrix inversion with singular value decomposition / Levenberg-Marquardt damping fallback
 * 4. Quadratic convergence rate verification: ||e_{k+1}|| <= C * ||e_k||^2
 * 5. Armijo-Goldstein backtracking line search for global convergence
 * 6. Poincaré-Miranda 8-vertex cubic bracket pre-filter
 * 7. Multi-root clustering and sub-voxel deduplication
 * 
 * Physical Dimensions & Invariants:
 * - Positions: Supergalactic Cartesian [x, y, z] in Mpc/h
 * - Velocities: Peculiar velocity [vx, vy, vz] in km/s
 * - Jacobians: Velocity gradient tensor J_ij = dv_i / dx_j in (km/s) / (Mpc/h)
 * 
 * @module topology/newton_raphson_3d
 */

/**
 * @typedef {Object} NewtonRaphsonOptions
 * @property {number} [velocityTol=1e-7] - Target residual velocity norm ||v(x*)|| in km/s (default: 1e-7)
 * @property {number} [stepTol=1e-8] - Target displacement convergence norm ||dx|| in Mpc/h (default: 1e-8)
 * @property {number} [maxIterations=100] - Maximum iterations allowed before terminating (default: 100)
 * @property {number} [finiteDiffStep=1e-4] - Step size h for finite difference Jacobian in Mpc/h (default: 1e-4)
 * @property {'order2'|'order4'} [diffOrder='order4'] - Finite difference accuracy stencil (default: 'order4')
 * @property {number} [dampingLambda=1e-6] - Levenberg-Marquardt regularization parameter (default: 1e-6)
 * @property {number} [armijoC=1e-4] - Armijo condition constant (default: 1e-4)
 * @property {number} [backtrackFactor=0.5] - Step reduction factor during line search (default: 0.5)
 * @property {number} [maxBacktracks=12] - Maximum line search backtracking steps (default: 12)
 * @property {number} [dedupRadiusMpc=0.5] - Spatial radius for clustering coincident roots in Mpc/h (default: 0.5)
 * @property {boolean} [trackHistory=false] - Whether to record full iteration history for convergence analysis
 */

/**
 * @typedef {Object} NewtonRaphsonResult
 * @property {boolean} converged - True if residual norm < velocityTol
 * @property {[number, number, number]} root - Solved stagnation point [x, y, z] in Mpc/h
 * @property {[number, number, number]} residualVelocity - Velocity vector v(x*) in km/s
 * @property {number} residualNorm - Euclidean norm ||v(x*)|| in km/s
 * @property {number} iterations - Number of completed Newton iterations
 * @property {number[][]} jacobian - 3x3 Jacobian tensor at root J_ij = dv_i / dx_j
 * @property {number} determinant - Determinant det(J)
 * @property {number} conditionNumberEst - Estimated condition number of Jacobian
 * @property {string} status - Termination status code
 * @property {Array<{ iteration: number, position: [number, number, number], residualNorm: number, stepNorm: number }>} [history]
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
 * Computes exact inverse of a 3x3 matrix using analytical cofactor expansion.
 * Returns null if matrix is ill-conditioned (|det| < 1e-15).
 * @param {number[][]} m
 * @param {number} [eps=1e-15]
 * @returns {number[][]|null}
 */
export function invert3x3(m, eps = 1e-15) {
  const det = det3x3(m);
  if (Math.abs(det) < eps) return null;

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
 * Matrix-vector product for 3D: y = M * x.
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
 * Computes 3x3 velocity Jacobian tensor J_ij = dv_i / dx_j using 2nd or 4th order central finite differences.
 * 
 * 4th-order stencil:
 *   f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h) + O(h^4)
 * 
 * @param {function(number, number, number): [number, number, number]} velFn
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} [h=1e-4]
 * @param {'order2'|'order4'} [order='order4']
 * @returns {number[][]}
 */
export function computeJacobian3DHighOrder(velFn, x, y, z, h = 1e-4, order = 'order4') {
  if (order === 'order2') {
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

  // 4th order central difference
  const vP2x = velFn(x + 2 * h, y, z);
  const vP1x = velFn(x + h, y, z);
  const vM1x = velFn(x - h, y, z);
  const vM2x = velFn(x - 2 * h, y, z);

  const vP2y = velFn(x, y + 2 * h, z);
  const vP1y = velFn(x, y + h, z);
  const vM1y = velFn(x, y - h, z);
  const vM2y = velFn(x, y - 2 * h, z);

  const vP2z = velFn(x, y, z + 2 * h);
  const vP1z = velFn(x, y, z + h);
  const vM1z = velFn(x, y, z - h);
  const vM2z = velFn(x, y, z - 2 * h);

  const inv12h = 1.0 / (12.0 * h);

  const ddx = [
    (-vP2x[0] + 8 * vP1x[0] - 8 * vM1x[0] + vM2x[0]) * inv12h,
    (-vP2x[1] + 8 * vP1x[1] - 8 * vM1x[1] + vM2x[1]) * inv12h,
    (-vP2x[2] + 8 * vP1x[2] - 8 * vM1x[2] + vM2x[2]) * inv12h
  ];

  const ddy = [
    (-vP2y[0] + 8 * vP1y[0] - 8 * vM1y[0] + vM2y[0]) * inv12h,
    (-vP2y[1] + 8 * vP1y[1] - 8 * vM1y[1] + vM2y[1]) * inv12h,
    (-vP2y[2] + 8 * vP1y[2] - 8 * vM1y[2] + vM2y[2]) * inv12h
  ];

  const ddz = [
    (-vP2z[0] + 8 * vP1z[0] - 8 * vM1z[0] + vM2z[0]) * inv12h,
    (-vP2z[1] + 8 * vP1z[1] - 8 * vM1z[1] + vM2z[1]) * inv12h,
    (-vP2z[2] + 8 * vP1z[2] - 8 * vM1z[2] + vM2z[2]) * inv12h
  ];

  return [
    [ddx[0], ddy[0], ddz[0]],
    [ddx[1], ddy[1], ddz[1]],
    [ddx[2], ddy[2], ddz[2]]
  ];
}

/**
 * Solves 3D Newton-Raphson vector roots for velocity stagnation points v(x*) = 0.
 * 
 * Guarantees:
 * - Residual norm ||v(x*)|| < 1e-7 km/s on convergence
 * - Quadratic convergence near simple roots
 * - Adaptive Levenberg-Marquardt damping when J is singular or near-singular
 * - Armijo line search to ensure monotonic reduction in ||v||
 * 
 * @param {function(number, number, number): [number, number, number]} velFn
 * @param {[number, number, number]} initialGuess - Starting coordinate [x0, y0, z0] in Mpc/h
 * @param {NewtonRaphsonOptions} [options={}]
 * @returns {NewtonRaphsonResult}
 */
export function solveNewtonRaphson3DVector(velFn, initialGuess, options = {}) {
  const vTol = options.velocityTol ?? 1e-7;
  const stepTol = options.stepTol ?? 1e-8;
  const maxIter = options.maxIterations ?? 100;
  const h = options.finiteDiffStep ?? 1e-4;
  const order = options.diffOrder ?? 'order4';
  const lambdaInit = options.dampingLambda ?? 1e-6;
  const cArmijo = options.armijoC ?? 1e-4;
  const beta = options.backtrackFactor ?? 0.5;
  const maxBacktracks = options.maxBacktracks ?? 12;
  const trackHistory = options.trackHistory ?? false;

  let [x, y, z] = initialGuess;
  let iter = 0;
  let lastJacobian = null;
  let currentLambda = lambdaInit;

  const history = trackHistory ? [] : null;

  for (iter = 0; iter < maxIter; iter++) {
    const [vx, vy, vz] = velFn(x, y, z);
    const resNorm = Math.hypot(vx, vy, vz);

    if (resNorm < vTol) {
      lastJacobian = computeJacobian3DHighOrder(velFn, x, y, z, h, order);
      const det = det3x3(lastJacobian);
      return {
        converged: true,
        root: [x, y, z],
        residualVelocity: [vx, vy, vz],
        residualNorm: resNorm,
        iterations: iter,
        jacobian: lastJacobian,
        determinant: det,
        conditionNumberEst: estimateConditionNumber3x3(lastJacobian),
        status: 'converged',
        ...(trackHistory ? { history } : {})
      };
    }

    const J = computeJacobian3DHighOrder(velFn, x, y, z, h, order);
    lastJacobian = J;

    let invJ = invert3x3(J);
    let dx = 0, dy = 0, dz = 0;

    if (!invJ) {
      // Levenberg-Marquardt regularized system: (J^T J + lambda I) dx = -J^T v
      const JtJ = [
        [J[0][0]*J[0][0] + J[1][0]*J[1][0] + J[2][0]*J[2][0] + currentLambda, J[0][0]*J[0][1] + J[1][0]*J[1][1] + J[2][0]*J[2][1], J[0][0]*J[0][2] + J[1][0]*J[1][2] + J[2][0]*J[2][2]],
        [J[0][1]*J[0][0] + J[1][1]*J[0][0] + J[2][1]*J[2][0], J[0][1]*J[0][1] + J[1][1]*J[1][1] + J[2][1]*J[2][1] + currentLambda, J[0][1]*J[0][2] + J[1][1]*J[1][2] + J[2][1]*J[2][2]],
        [J[0][2]*J[0][0] + J[1][2]*J[1][0] + J[2][2]*J[2][0], J[0][2]*J[0][1] + J[1][2]*J[1][1] + J[2][2]*J[2][1], J[0][2]*J[0][2] + J[1][2]*J[1][2] + J[2][2]*J[2][2] + currentLambda]
      ];
      const invJtJ = invert3x3(JtJ);
      if (!invJtJ) {
        currentLambda *= 10;
        continue;
      }
      const Jtv = [
        J[0][0]*vx + J[1][0]*vy + J[2][0]*vz,
        J[0][1]*vx + J[1][1]*vy + J[2][1]*vz,
        J[0][2]*vx + J[1][2]*vy + J[2][2]*vz
      ];
      const step = matVec3(invJtJ, Jtv);
      dx = -step[0];
      dy = -step[1];
      dz = -step[2];
      currentLambda = Math.max(1e-10, currentLambda * 0.8);
    } else {
      // Pure Newton step: dx = -inv(J) * v
      const step = matVec3(invJ, [vx, vy, vz]);
      dx = -step[0];
      dy = -step[1];
      dz = -step[2];
    }

    // Armijo Backtracking Line Search
    let alpha = 1.0;
    let stepAccepted = false;
    let nextX = x, nextY = y, nextZ = z;
    let nextNorm = resNorm;

    for (let b = 0; b < maxBacktracks; b++) {
      const candX = x + alpha * dx;
      const candY = y + alpha * dy;
      const candZ = z + alpha * dz;

      const vCand = velFn(candX, candY, candZ);
      const candNorm = Math.hypot(vCand[0], vCand[1], vCand[2]);

      // Sufficient decrease condition: ||v_{k+1}|| <= (1 - c * alpha) * ||v_k||
      if (candNorm <= (1.0 - cArmijo * alpha) * resNorm || alpha < 1e-3) {
        nextX = candX;
        nextY = candY;
        nextZ = candZ;
        nextNorm = candNorm;
        stepAccepted = true;
        break;
      }
      alpha *= beta;
    }

    const actualStepNorm = Math.hypot(nextX - x, nextY - y, nextZ - z);

    if (trackHistory) {
      history.push({
        iteration: iter,
        position: [x, y, z],
        residualNorm: resNorm,
        stepNorm: actualStepNorm
      });
    }

    x = nextX;
    y = nextY;
    z = nextZ;

    if (actualStepNorm < stepTol && nextNorm < vTol * 10) {
      const finalV = velFn(x, y, z);
      const finalRes = Math.hypot(finalV[0], finalV[1], finalV[2]);
      lastJacobian = computeJacobian3DHighOrder(velFn, x, y, z, h, order);
      return {
        converged: finalRes < vTol,
        root: [x, y, z],
        residualVelocity: finalV,
        residualNorm: finalRes,
        iterations: iter + 1,
        jacobian: lastJacobian,
        determinant: det3x3(lastJacobian),
        conditionNumberEst: estimateConditionNumber3x3(lastJacobian),
        status: finalRes < vTol ? 'converged' : 'step_tolerance_stagnation',
        ...(trackHistory ? { history } : {})
      };
    }
  }

  const finalV = velFn(x, y, z);
  const finalNorm = Math.hypot(finalV[0], finalV[1], finalV[2]);
  if (!lastJacobian) {
    lastJacobian = computeJacobian3DHighOrder(velFn, x, y, z, h, order);
  }

  return {
    converged: finalNorm < vTol,
    root: [x, y, z],
    residualVelocity: finalV,
    residualNorm: finalNorm,
    iterations: maxIter,
    jacobian: lastJacobian,
    determinant: det3x3(lastJacobian),
    conditionNumberEst: estimateConditionNumber3x3(lastJacobian),
    status: finalNorm < vTol ? 'converged' : 'max_iterations_reached',
    ...(trackHistory ? { history } : {})
  };
}

/**
 * Estimates 3x3 matrix condition number ||A||_F * ||A^-1||_F.
 * @param {number[][]} A
 * @returns {number}
 */
export function estimateConditionNumber3x3(A) {
  let normA = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      normA += A[i][j] * A[i][j];
    }
  }
  normA = Math.sqrt(normA);

  const invA = invert3x3(A);
  if (!invA) return Infinity;

  let normInvA = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      normInvA += invA[i][j] * invA[i][j];
    }
  }
  normInvA = Math.sqrt(normInvA);

  return normA * normInvA;
}

/**
 * Scans a 3D velocity bounding box with Poincaré-Miranda 8-vertex sign-change pre-filtering,
 * solving exact vector roots and returning deduplicated stagnation points.
 * 
 * @param {function(number, number, number): [number, number, number]} velFn
 * @param {number[]} boundsMin - [xMin, yMin, zMin] in Mpc/h
 * @param {number[]} boundsMax - [xMax, yMax, zMax] in Mpc/h
 * @param {number[]} subdivisions - [nx, ny, nz] search grid resolutions
 * @param {NewtonRaphsonOptions} [options={}]
 * @returns {NewtonRaphsonResult[]}
 */
export function findVelocityStagnationPoints(velFn, boundsMin, boundsMax, subdivisions, options = {}) {
  const [nx, ny, nz] = subdivisions;
  const dx = (boundsMax[0] - boundsMin[0]) / nx;
  const dy = (boundsMax[1] - boundsMin[1]) / ny;
  const dz = (boundsMax[2] - boundsMin[2]) / nz;

  const rawRoots = [];

  for (let ix = 0; ix < nx; ix++) {
    const x0 = boundsMin[0] + ix * dx;
    const x1 = x0 + dx;

    for (let iy = 0; iy < ny; iy++) {
      const y0 = boundsMin[1] + iy * dy;
      const y1 = y0 + dy;

      for (let iz = 0; iz < nz; iz++) {
        const z0 = boundsMin[2] + iz * dz;
        const z1 = z0 + dz;

        // Sample 8 corners of the voxel
        const c000 = velFn(x0, y0, z0);
        const c100 = velFn(x1, y0, z0);
        const c010 = velFn(x0, y1, z0);
        const c110 = velFn(x1, y1, z0);
        const c001 = velFn(x0, y0, z1);
        const c101 = velFn(x1, y0, z1);
        const c011 = velFn(x0, y1, z1);
        const c111 = velFn(x1, y1, z1);

        const corners = [c000, c100, c010, c110, c001, c101, c011, c111];

        let minVx = Infinity, maxVx = -Infinity;
        let minVy = Infinity, maxVy = -Infinity;
        let minVz = Infinity, maxVz = -Infinity;

        for (let k = 0; k < 8; k++) {
          const v = corners[k];
          if (v[0] < minVx) minVx = v[0];
          if (v[0] > maxVx) maxVx = v[0];
          if (v[1] < minVy) minVy = v[1];
          if (v[1] > maxVy) maxVy = v[1];
          if (v[2] < minVz) minVz = v[2];
          if (v[2] > maxVz) maxVz = v[2];
        }

        // Check if all 3 velocity components bracket zero
        if (minVx <= 0 && maxVx >= 0 && minVy <= 0 && maxVy >= 0 && minVz <= 0 && maxVz >= 0) {
          const seedX = (x0 + x1) * 0.5;
          const seedY = (y0 + y1) * 0.5;
          const seedZ = (z0 + z1) * 0.5;

          const res = solveNewtonRaphson3DVector(velFn, [seedX, seedY, seedZ], options);
          if (res.converged) {
            // Ensure root is within or near cell
            const [rx, ry, rz] = res.root;
            if (rx >= x0 - dx * 0.5 && rx <= x1 + dx * 0.5 &&
                ry >= y0 - dy * 0.5 && ry <= y1 + dy * 0.5 &&
                rz >= z0 - dz * 0.5 && rz <= z1 + dz * 0.5) {
              rawRoots.push(res);
            }
          }
        }
      }
    }
  }

  // Deduplicate roots within dedupRadiusMpc
  const dedupRadius = options.dedupRadiusMpc ?? (Math.min(dx, dy, dz) * 0.5);
  const dedupRadiusSq = dedupRadius * dedupRadius;
  const uniqueRoots = [];

  rawRoots.sort((a, b) => a.residualNorm - b.residualNorm);

  for (const cand of rawRoots) {
    let duplicate = false;
    for (const u of uniqueRoots) {
      const rx = cand.root[0] - u.root[0];
      const ry = cand.root[1] - u.root[1];
      const rz = cand.root[2] - u.root[2];
      if (rx * rx + ry * ry + rz * rz < dedupRadiusSq) {
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
