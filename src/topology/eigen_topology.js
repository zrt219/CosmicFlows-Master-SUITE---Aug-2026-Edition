/**
 * @file eigen_topology.js
 * @description Full 3x3 Strain-Rate Tensor Jacobi Diagonalizer and Hyperbolic Cosmic Web Classifier.
 * Diagonalizes symmetric strain-rate tensors to machine precision, computes divergence, curl (vorticity),
 * shear scalars, second invariants, and classifies critical points into Morse topological types:
 * - Attractor (Supercluster Megasink: λ1, λ2, λ3 < 0)
 * - Saddle-Filament (Type-2 Saddle / Cosmic Filament Hub: 2 contracting, 1 expanding)
 * - Saddle-Wall (Type-1 Saddle / Cosmic Sheet Hub: 1 contracting, 2 expanding)
 * - Repeller (Cosmic Void Center: λ1, λ2, λ3 > 0)
 * - Degenerate (Non-hyperbolic zero eigenvalues)
 */

/**
 * @typedef {Object} EigenResult
 * @property {number[]} eigenvalues - Sorted eigenvalues [λ1, λ2, λ3] (ascending: λ1 <= λ2 <= λ3)
 * @property {number[][]} eigenvectors - 3 orthonormal column/row eigenvectors [e1, e2, e3]
 * @property {number} iterations - Jacobi sweeps to convergence
 */

/**
 * @typedef {'ATTRACTOR'|'SADDLE_FILAMENT'|'SADDLE_WALL'|'REPELLER'|'DEGENERATE'} MorseClassification
 */

/**
 * @typedef {Object} TopologicalProfile
 * @property {MorseClassification} type - Hyperbolic Morse critical point classification
 * @property {number} morseIndex - Number of negative eigenvalues (0 = Repeller, 1 = Wall, 2 = Filament, 3 = Attractor)
 * @property {number[]} eigenvalues - Sorted eigenvalues [λ1, λ2, λ3]
 * @property {number[][]} eigenvectors - Orthonormal eigenvectors [e1, e2, e3]
 * @property {number} trace - Trace Tr(J) = λ1 + λ2 + λ3 = div(v)
 * @property {number} divergence - Velocity field divergence ∇·v
 * @property {number[]} curl - Vorticity vector ∇×v = [ωx, ωy, ωz]
 * @property {number} vorticityMagnitude - ||ω||
 * @property {number} strainNorm - Frobenius norm of symmetric strain-rate tensor ||S||
 * @property {number} secondInvariantQ - Q-criterion invariant Q = 0.5 * (||ω||^2 - ||S||^2)
 * @property {number[][]} strainRateTensor - S_ij = 0.5 * (J_ij + J_ji)
 * @property {number[][]} vorticityTensor - Ω_ij = 0.5 * (J_ij - J_ji)
 * @property {number[]} primarySpineVector - Primary 1D filament outflow or wall normal eigenvector
 */

/**
 * Full Jacobi Eigenvalue Algorithm for a 3x3 symmetric matrix.
 * Solves S * v = λ * v iteratively via Givens plane rotations until off-diagonal elements < tol.
 *
 * @param {number[][]} mat - 3x3 symmetric matrix
 * @param {number} [tol=1e-15] - Convergence threshold for sum of squared off-diagonals
 * @param {number} [maxSweeps=50] - Maximum Jacobi sweeps
 * @returns {EigenResult}
 */
export function jacobiDiagonalize3x3(mat, tol = 1e-15, maxSweeps = 50) {
  // Symmetrize copy of matrix
  const A = [
    [mat[0][0], 0.5 * (mat[0][1] + mat[1][0]), 0.5 * (mat[0][2] + mat[2][0])],
    [0.5 * (mat[0][1] + mat[1][0]), mat[1][1], 0.5 * (mat[1][2] + mat[2][1])],
    [0.5 * (mat[0][2] + mat[2][0]), 0.5 * (mat[1][2] + mat[2][1]), mat[2][2]]
  ];

  // Identity eigenvector matrix V
  const V = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
  ];

  let sweep = 0;
  for (sweep = 0; sweep < maxSweeps; sweep++) {
    // Compute off-diagonal norm
    const offDiag = Math.abs(A[0][1]) + Math.abs(A[0][2]) + Math.abs(A[1][2]);
    if (offDiag < tol) break;

    // Zero each off-diagonal pair (p, q) in cyclic order (0,1), (0,2), (1,2)
    const pairs = [[0, 1], [0, 2], [1, 2]];
    for (const [p, q] of pairs) {
      const app = A[p][p];
      const aqq = A[q][q];
      const apq = A[p][q];

      if (Math.abs(apq) < 1e-16) continue;

      const tau = (aqq - app) / (2.0 * apq);
      let t;
      if (tau >= 0) {
        t = 1.0 / (tau + Math.sqrt(1.0 + tau * tau));
      } else {
        t = -1.0 / (-tau + Math.sqrt(1.0 + tau * tau));
      }

      const c = 1.0 / Math.sqrt(1.0 + t * t);
      const s = t * c;
      const h = t * apq;

      A[p][p] -= h;
      A[q][q] += h;
      A[p][q] = 0.0;
      A[q][p] = 0.0;

      // Update remaining elements in row/col
      for (let r = 0; r < 3; r++) {
        if (r !== p && r !== q) {
          const arp = A[r][p];
          const arq = A[r][q];
          A[r][p] = c * arp - s * arq;
          A[p][r] = A[r][p];
          A[r][q] = s * arp + c * arq;
          A[q][r] = A[r][q];
        }
      }

      // Accumulate eigenvector transformations: V = V * R(p, q)
      for (let r = 0; r < 3; r++) {
        const vrp = V[r][p];
        const vrq = V[r][q];
        V[r][p] = c * vrp - s * vrq;
        V[r][q] = s * vrp + c * vrq;
      }
    }
  }

  // Extract eigenvalues and column eigenvectors
  let eValues = [A[0][0], A[1][1], A[2][2]];
  let eVectors = [
    [V[0][0], V[1][0], V[2][0]], // eigenvector 0
    [V[0][1], V[1][1], V[2][1]], // eigenvector 1
    [V[0][2], V[1][2], V[2][2]]  // eigenvector 2
  ];

  // Sort eigenvalues ascending: λ1 <= λ2 <= λ3
  const indices = [0, 1, 2];
  indices.sort((a, b) => eValues[a] - eValues[b]);

  const sortedValues = indices.map(i => eValues[i]);
  const sortedVectors = indices.map(i => eVectors[i]);

  return {
    eigenvalues: sortedValues,
    eigenvectors: sortedVectors,
    iterations: sweep
  };
}

/**
 * Analyzes the complete 3x3 velocity Jacobian tensor J = ∇v, performing full strain-rate
 * decomposition, Jacobi diagonalization, vorticity extraction, and Morse classification.
 *
 * @param {number[][]} J - 3x3 velocity Jacobian tensor J_ij = dv_i / dx_j
 * @param {number} [zeroTol=1e-6] - Tolerance to distinguish degenerate (zero) eigenvalues
 * @returns {TopologicalProfile}
 */
export function analyzeVelocityTensor(J, zeroTol = 1e-6) {
  // Symmetric strain-rate tensor: S_ij = 0.5 * (J_ij + J_ji)
  const S = [
    [J[0][0], 0.5 * (J[0][1] + J[1][0]), 0.5 * (J[0][2] + J[2][0])],
    [0.5 * (J[1][0] + J[0][1]), J[1][1], 0.5 * (J[1][2] + J[2][1])],
    [0.5 * (J[2][0] + J[0][2]), 0.5 * (J[2][1] + J[1][2]), J[2][2]]
  ];

  // Anti-symmetric vorticity tensor: Ω_ij = 0.5 * (J_ij - J_ji)
  const Omega = [
    [0.0, 0.5 * (J[0][1] - J[1][0]), 0.5 * (J[0][2] - J[2][0])],
    [0.5 * (J[1][0] - J[0][1]), 0.0, 0.5 * (J[1][2] - J[2][1])],
    [0.5 * (J[2][0] - J[0][2]), 0.5 * (J[2][1] - J[1][2]), 0.0]
  ];

  // Vorticity vector (curl): ω = ∇ × v = [dv_z/dy - dv_y/dz, dv_x/dz - dv_z/dx, dv_y/dx - dv_x/dy]
  const curl = [
    J[2][1] - J[1][2],
    J[0][2] - J[2][0],
    J[1][0] - J[0][1]
  ];
  const vorticityMag = Math.hypot(curl[0], curl[1], curl[2]);

  // Trace = divergence
  const trace = J[0][0] + J[1][1] + J[2][2];
  const divergence = trace;

  // Strain norm ||S|| = sqrt(sum S_ij^2)
  let strainNormSq = 0;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      strainNormSq += S[i][j] * S[i][j];
    }
  }
  const strainNorm = Math.sqrt(strainNormSq);

  // Second invariant Q = 0.5 * (||ω||^2 - ||S||^2)
  const secondInvariantQ = 0.5 * (vorticityMag * vorticityMag - strainNormSq);

  // Jacobi diagonalization on symmetric strain tensor S
  const { eigenvalues, eigenvectors } = jacobiDiagonalize3x3(S);
  const [l1, l2, l3] = eigenvalues;

  // Morse classification based on eigenvalue signs
  let negCount = 0;
  let posCount = 0;
  let zeroCount = 0;

  for (const l of eigenvalues) {
    if (Math.abs(l) < zeroTol) zeroCount++;
    else if (l < 0) negCount++;
    else posCount++;
  }

  /** @type {MorseClassification} */
  let classification;
  let morseIndex = negCount;
  let primarySpine = [...eigenvectors[2]]; // default

  if (zeroCount > 0) {
    classification = 'DEGENERATE';
  } else if (negCount === 3) {
    classification = 'ATTRACTOR'; // Node sink / supercluster core
    primarySpine = [...eigenvectors[0]]; // strongest contraction
  } else if (negCount === 2 && posCount === 1) {
    classification = 'SADDLE_FILAMENT'; // Type-2 saddle: 2 contracting (l1, l2 < 0), 1 expanding (l3 > 0)
    primarySpine = [...eigenvectors[2]]; // outflow spine direction
  } else if (negCount === 1 && posCount === 2) {
    classification = 'SADDLE_WALL'; // Type-1 saddle: 1 contracting (l1 < 0), 2 expanding (l2, l3 > 0)
    primarySpine = [...eigenvectors[0]]; // sheet normal (inflow) direction
  } else if (posCount === 3) {
    classification = 'REPELLER'; // Node source / void center
    primarySpine = [...eigenvectors[2]]; // strongest expansion
  } else {
    classification = 'DEGENERATE';
  }

  return {
    type: classification,
    morseIndex,
    eigenvalues,
    eigenvectors,
    trace,
    divergence,
    curl,
    vorticityMagnitude: vorticityMag,
    strainNorm,
    secondInvariantQ,
    strainRateTensor: S,
    vorticityTensor: Omega,
    primarySpineVector: primarySpine
  };
}
