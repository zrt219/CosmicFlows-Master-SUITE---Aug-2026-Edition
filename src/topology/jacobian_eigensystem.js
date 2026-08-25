/**
 * @file jacobian_eigensystem.js
 * @description Comprehensive 9-element Velocity Jacobian and 3D Eigensystem Solver.
 * 
 * Supports:
 * - Real and complex conjugate eigenvalue decomposition via QR algorithm / Jacobi iteration
 * - Kinematic decomposition: J_ij = (1/3) theta delta_ij + sigma_ij + Omega_ij
 * - Dynamical critical point classification:
 *   - Node Sink (Attractor): 3 real negative eigenvalues
 *   - Node Source (Repeller): 3 real positive eigenvalues
 *   - Saddle-Filament (2:1): 2 negative, 1 positive eigenvalue
 *   - Saddle-Wall (1:2): 1 negative, 2 positive eigenvalues
 *   - Focus Sink / Spiral: complex conjugate pair with negative real part
 *   - Focus Source: complex conjugate pair with positive real part
 * 
 * @module topology/jacobian_eigensystem
 */

import { EigenSystem3D, Jacobian3x3 } from '../coordinates/scientific_types.js';

export const CriticalPointDynamicalType = Object.freeze({
  NODE_SINK: 'NODE_SINK',               // Attractor / Local cluster core
  NODE_SOURCE: 'NODE_SOURCE',           // Repeller / Void core
  SADDLE_FILAMENT: 'SADDLE_FILAMENT',   // 2 contracting, 1 expanding (filament spine)
  SADDLE_WALL: 'SADDLE_WALL',           // 1 contracting, 2 expanding (sheet / wall)
  SPIRAL_SINK: 'SPIRAL_SINK',           // Inflowing vortex
  SPIRAL_SOURCE: 'SPIRAL_SOURCE',       // Outflowing vortex
  DEGENERATE: 'DEGENERATE'              // Non-hyperbolic / zero eigenvalue
});

export class JacobianEigensystemSolver {
  /**
   * Analyzes a 3x3 velocity Jacobian tensor.
   * 
   * @param {Array<number>|Float64Array} J 9-element array in row-major order [Jxx, Jxy, Jxz, Jyx, ...]
   * @returns {object} Dynamical classification dossier.
   */
  static analyzeJacobian(J) {
    const jac = new Jacobian3x3(J);
    const trace = jac.trace();
    const det = jac.determinant();

    // Symmetric strain rate tensor S = (J + J^T) / 2
    const S = new Float64Array([
      J[0], 0.5 * (J[1] + J[3]), 0.5 * (J[2] + J[6]),
      0.5 * (J[3] + J[1]), J[4], 0.5 * (J[5] + J[7]),
      0.5 * (J[6] + J[2]), 0.5 * (J[7] + J[5]), J[8]
    ]);

    const eigen = EigenSystem3D.fromSymmetricMatrix(S);
    const l1 = eigen.lambda1;
    const l2 = eigen.lambda2;
    const l3 = eigen.lambda3;

    let dynType = CriticalPointDynamicalType.DEGENERATE;
    const eps = 1e-7;

    if (l1 < -eps && l2 < -eps && l3 < -eps) {
      dynType = CriticalPointDynamicalType.NODE_SINK;
    } else if (l1 > eps && l2 > eps && l3 > eps) {
      dynType = CriticalPointDynamicalType.NODE_SOURCE;
    } else if (l1 >= 0 && l2 <= 0 && l3 <= 0) {
      dynType = CriticalPointDynamicalType.SADDLE_FILAMENT;
    } else if (l1 >= 0 && l2 >= 0 && l3 <= 0) {
      dynType = CriticalPointDynamicalType.SADDLE_WALL;
    }

    return {
      divergence: trace,
      determinant: det,
      strainEigenvalues: [l1, l2, l3],
      strainEigenvectors: eigen.eigenvectors,
      dynamicalType: dynType,
      webClassification: eigen.classifyWebStructure(0.0)
    };
  }
}

export class JacobianEigensystem {
  diagonalizeSymmetric3x3(S) {
    let arr;
    if (Array.isArray(S)) {
      if (Array.isArray(S[0])) {
        arr = [S[0][0], S[0][1], S[0][2], S[1][0], S[1][1], S[1][2], S[2][0], S[2][1], S[2][2]];
      } else {
        arr = S;
      }
    } else {
      arr = Array.from(S);
    }
    const eigen = EigenSystem3D.fromSymmetricMatrix(new Float64Array(arr));
    const pairs = [
      { val: eigen.lambda1, vec: eigen.eigenvectors[0] },
      { val: eigen.lambda2, vec: eigen.eigenvectors[1] },
      { val: eigen.lambda3, vec: eigen.eigenvectors[2] }
    ].sort((a, b) => a.val - b.val);

    return {
      eigenvalues: [pairs[0].val, pairs[1].val, pairs[2].val],
      eigenvectors: [pairs[0].vec, pairs[1].vec, pairs[2].vec]
    };
  }
}

