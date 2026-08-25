/**
 * @file poincare_section_tracer.js
 * @description Complete, rigorous scientific implementation of the Poincaré Section
 * and Orbit Recurrence Analyzer for 3D cosmological velocity fields, peculiar flow streamlines,
 * and cluster accretion basins (Laniakea, Shapley, Virgo, Coma, Great Attractor, Perseus-Pisces).
 *
 * Theoretical Foundations:
 * 1. Surface of Section (Poincaré Map):
 *    Given a continuous dynamical system dx/dt = v(x) in R^3 and a codimension-1 transversal
 *    cutting hypersurface Sigma = { x in R^3 : n · (x - x_0) = 0 }, the Poincaré first return map
 *    P: Sigma -> Sigma maps a puncture P_k to the subsequent puncture P_{k+1} with directional
 *    filter v · n > 0.
 * 2. Exact Puncture Root Interpolation:
 *    Sub-grid intersection solving via linear, cubic Hermite spline with Newton-Raphson / Brent-Dekker
 *    bisection, and dense RK output coordinate transformation.
 * 3. Poincaré Recurrence Maps & RQA:
 *    Planar coordinate projections (u_k, v_k) in orthonormal 2D frame (e_u, e_v), return time moments,
 *    symplectic return Jacobian estimation D P(u, v), and complete Recurrence Quantification Analysis:
 *    RR, DET, LAM, TT, ENTR, L_max, V_max, TREND.
 * 4. Lyapunov Exponents & Chaos Diagnostics:
 *    - Benettin two-particle shadow renormalization for Maximal Lyapunov Exponent (MLE, lambda_max),
 *      FTLE time series, and Lyapunov time horizon t_Lyap = 1 / lambda_max.
 *    - Continuous Variational Tangent Map dPhi/dt = J(x) * Phi with periodic QR Gram-Schmidt
 *      decomposition for the complete 3D spectrum (lambda_1 >= lambda_2 >= lambda_3), Kolmogorov-Sinai
 *      entropy h_KS = sum_{lambda_i > 0} lambda_i, and Kaplan-Yorke attractor dimension D_KY.
 * 5. KAM Invariant Tori & Accretion Basin Classification:
 *    - 1D closed invariant curves (2-torus cross-sections) with radial dispersion and phase continuity.
 *    - Resonant island chains (period-q island structures with angular gaps).
 *    - Fixed point classification: Elliptic centers (O-points) vs Hyperbolic saddles (X-points).
 *    - Discrete Power Spectral Density (PSD) analysis of puncture coordinates.
 *
 * Physical Units:
 * - Coordinates: Mpc/h (h^-1 Mpc)
 * - Peculiar Velocities: km/s
 * - Cosmic Time: Gyr
 * - Lyapunov Exponents: Gyr^-1 or (km/s)/(Mpc/h)
 * - Frequency: Gyr^-1
 *
 * @module streamlines/poincare_section_tracer
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { VelocityField } from '../fields/velocity_field.js';

/**
 * Directional crossing filter for surface of section punctures.
 * @readonly
 * @enum {string}
 */
export const CrossingDirection = Object.freeze({
  POSITIVE: 'POSITIVE',   // v · n > 0 (one-way forward puncture)
  NEGATIVE: 'NEGATIVE',   // v · n < 0 (one-way backward puncture)
  BOTH: 'BOTH'            // v · n != 0 (bidirectional puncture)
});

/**
 * High-precision root finding and intersection interpolation methods.
 * @readonly
 * @enum {string}
 */
export const IntersectionMethod = Object.freeze({
  LINEAR: 'LINEAR',
  HERMITE_CUBIC: 'HERMITE_CUBIC',
  BRENT_EXACT: 'BRENT_EXACT',
  NEWTON_RAPHSON: 'NEWTON_RAPHSON'
});

/**
 * Scientific classification of dynamical orbits and cosmic streamline trajectories.
 * @readonly
 * @enum {string}
 */
export const OrbitClassification = Object.freeze({
  PERIODIC_LIMIT_CYCLE: 'PERIODIC_LIMIT_CYCLE',
  QUASI_PERIODIC_INVARIANT_TORUS: 'QUASI_PERIODIC_INVARIANT_TORUS',
  ISLAND_CHAIN_RESONANCE: 'ISLAND_CHAIN_RESONANCE',
  CHAOTIC_SEA_STOCHASTIC: 'CHAOTIC_SEA_STOCHASTIC',
  OPEN_STREAMLINE_ESCAPE: 'OPEN_STREAMLINE_ESCAPE',
  DIRECT_COLLAPSE_SINK: 'DIRECT_COLLAPSE_SINK',
  ASYMPTOTIC_STALL: 'ASYMPTOTIC_STALL',
  INSUFFICIENT_DATA: 'INSUFFICIENT_DATA'
});

/**
 * Predefined cosmological and cluster-centric cutting plane orientations.
 * @readonly
 * @enum {string}
 */
export const PlaneOrientation = Object.freeze({
  SGZ_PLANE: 'SGZ_PLANE',
  SGY_PLANE: 'SGY_PLANE',
  SGX_PLANE: 'SGX_PLANE',
  CUSTOM_PLANE: 'CUSTOM_PLANE',
  PRINCIPAL_TIDAL_PLANE: 'PRINCIPAL_TIDAL_PLANE',
  CLUSTER_MERGER_AXIS: 'CLUSTER_MERGER_AXIS'
});

/**
 * Numerical ODE integrator algorithms.
 * @readonly
 * @enum {string}
 */
export const PoincareIntegratorType = Object.freeze({
  RK4: 'RK4',
  RK45_CASH_KARP: 'RK45_CASH_KARP',
  DOPRI5: 'DOPRI5'
});

/**
 * Fixed point stability classifications in 2D Poincaré section.
 * @readonly
 * @enum {string}
 */
export const FixedPointType = Object.freeze({
  ELLIPTIC_STABLE_CENTER: 'ELLIPTIC_STABLE_CENTER',
  HYPERBOLIC_SADDLE: 'HYPERBOLIC_SADDLE',
  ATTRACTING_FOCUS_SINK: 'ATTRACTING_FOCUS_SINK',
  REPELLING_FOCUS_SOURCE: 'REPELLING_FOCUS_SOURCE',
  PARABOLIC_DEGENERATE: 'PARABOLIC_DEGENERATE'
});

// Standard cosmological unit conversion factor: (km/s) / (Mpc/h) -> Gyr^-1
// 1 Mpc ≈ 3.085677581e19 km, 1 Gyr ≈ 3.15576e16 s
export const KMS_MPC_TO_GYR_INV = 1.0227121650537077e-3;

// Butcher tableau for Cash-Karp RK45
const CK_A = [
  [],
  [1.0 / 5.0],
  [3.0 / 40.0, 9.0 / 40.0],
  [3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0],
  [-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0],
  [1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0]
];
const CK_B5 = [37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0];
const CK_B4 = [2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0];

// Butcher tableau for Dormand-Prince 5(4) (DOPRI5)
const DP_A21 = 1.0 / 5.0;
const DP_A31 = 3.0 / 40.0, DP_A32 = 9.0 / 40.0;
const DP_A41 = 44.0 / 45.0, DP_A42 = -56.0 / 15.0, DP_A43 = 32.0 / 9.0;
const DP_A51 = 19372.0 / 6561.0, DP_A52 = -25360.0 / 2187.0, DP_A53 = 64448.0 / 6561.0, DP_A54 = -212.0 / 729.0;
const DP_A61 = 9017.0 / 3168.0, DP_A62 = -355.0 / 33.0, DP_A63 = 46732.0 / 5247.0, DP_A64 = 49.0 / 176.0, DP_A65 = -5103.0 / 18656.0;
const DP_B = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0];
const DP_E = [71.0 / 57600.0, 0.0, -71.0 / 16695.0, 71.0 / 1920.0, -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0];

// ============================================================================
// 1. VECTOR, MATRIX & GEOMETRIC ALGEBRA UTILITIES
// ============================================================================

/**
 * Computes Euclidean 2-norm of a 3-vector.
 * @param {Array<number>|Float64Array} v
 * @returns {number}
 */
export function vec3Norm(v) {
  return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/**
 * Returns normalized unit 3-vector.
 * @param {Array<number>|Float64Array} v
 * @returns {[number, number, number]}
 */
export function vec3Normalize(v) {
  const n = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (n < 1e-15) return [0, 0, 1];
  const inv = 1.0 / n;
  return [v[0] * inv, v[1] * inv, v[2] * inv];
}

/**
 * Computes scalar dot product of two 3-vectors.
 * @param {Array<number>|Float64Array} a
 * @param {Array<number>|Float64Array} b
 * @returns {number}
 */
export function vec3Dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Computes cross product a x b.
 * @param {Array<number>|Float64Array} a
 * @param {Array<number>|Float64Array} b
 * @returns {[number, number, number]}
 */
export function vec3Cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

/**
 * Matrix-vector product for 3x3 row-major matrix and 3-vector.
 * @param {Array<number>|Float64Array} M 3x3 row-major
 * @param {Array<number>|Float64Array} v
 * @returns {[number, number, number]}
 */
export function mat3VecMul(M, v) {
  return [
    M[0] * v[0] + M[1] * v[1] + M[2] * v[2],
    M[3] * v[0] + M[4] * v[1] + M[5] * v[2],
    M[6] * v[0] + M[7] * v[1] + M[8] * v[2]
  ];
}

/**
 * Matrix-matrix product for two 3x3 row-major matrices.
 * @param {Array<number>|Float64Array} A
 * @param {Array<number>|Float64Array} B
 * @returns {Float64Array}
 */
export function mat3Mul(A, B) {
  const out = new Float64Array(9);
  for (let r = 0; r < 3; r++) {
    const r3 = r * 3;
    for (let c = 0; c < 3; c++) {
      out[r3 + c] = A[r3 + 0] * B[0 * 3 + c] +
                   A[r3 + 1] * B[1 * 3 + c] +
                   A[r3 + 2] * B[2 * 3 + c];
    }
  }
  return out;
}

/**
 * Modified Gram-Schmidt QR decomposition of a 3x3 column-matrix.
 * Used for continuous tangent frame re-orthonormalization in Lyapunov spectrum analysis.
 *
 * @param {Float64Array|Array<number>} A Column-vector matrix [v0, v1, v2].
 * @returns {{ Q: Float64Array, R: Float64Array }} Q is orthogonal (columns q0, q1, q2), R is upper triangular.
 */
export function qrDecomposition3x3(A) {
  const q0 = [A[0], A[3], A[6]];
  const q1 = [A[1], A[4], A[7]];
  const q2 = [A[2], A[5], A[8]];

  const R = new Float64Array(9);

  // Column 0
  const r00 = vec3Norm(q0);
  R[0] = r00;
  if (r00 > 1e-15) {
    const inv00 = 1.0 / r00;
    q0[0] *= inv00; q0[1] *= inv00; q0[2] *= inv00;
  }

  // Column 1
  const r01 = vec3Dot(q0, q1);
  R[1] = r01;
  q1[0] -= r01 * q0[0];
  q1[1] -= r01 * q0[1];
  q1[2] -= r01 * q0[2];

  const r11 = vec3Norm(q1);
  R[4] = r11;
  if (r11 > 1e-15) {
    const inv11 = 1.0 / r11;
    q1[0] *= inv11; q1[1] *= inv11; q1[2] *= inv11;
  }

  // Column 2
  const r02 = vec3Dot(q0, q2);
  const r12 = vec3Dot(q1, q2);
  R[2] = r02;
  R[5] = r12;
  q2[0] -= (r02 * q0[0] + r12 * q1[0]);
  q2[1] -= (r02 * q0[1] + r12 * q1[1]);
  q2[2] -= (r02 * q0[2] + r12 * q1[2]);

  const r22 = vec3Norm(q2);
  R[8] = r22;
  if (r22 > 1e-15) {
    const inv22 = 1.0 / r22;
    q2[0] *= inv22; q2[1] *= inv22; q2[2] *= inv22;
  }

  const Q = new Float64Array(9);
  Q[0] = q0[0]; Q[1] = q1[0]; Q[2] = q2[0];
  Q[3] = q0[1]; Q[4] = q1[1]; Q[5] = q2[1];
  Q[6] = q0[2]; Q[7] = q1[2]; Q[8] = q2[2];

  return { Q, R };
}

/**
 * Jacobi eigenvalue algorithm for 3x3 real symmetric matrices.
 * @param {Array<Array<number>>|Float64Array|Array<number>} matrix
 * @returns {{ eigenvalues: [number, number, number], eigenvectors: Array<[number, number, number]> }}
 */
export function jacobiEigenvalues3x3(matrix) {
  const A = new Float64Array(9);
  if (matrix.length === 3 && Array.isArray(matrix[0])) {
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        A[r * 3 + c] = matrix[r][c];
      }
    }
  } else {
    for (let i = 0; i < 9; i++) A[i] = matrix[i];
  }

  const V = new Float64Array([
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  ]);

  const maxIter = 50;
  for (let iter = 0; iter < maxIter; iter++) {
    let p = 0, q = 1;
    let maxOff = Math.abs(A[1]);

    if (Math.abs(A[2]) > maxOff) { p = 0; q = 2; maxOff = Math.abs(A[2]); }
    if (Math.abs(A[5]) > maxOff) { p = 1; q = 2; maxOff = Math.abs(A[5]); }

    if (maxOff < 1e-14) break;

    const app = A[p * 3 + p];
    const aqq = A[q * 3 + q];
    const apq = A[p * 3 + q];

    const theta = 0.5 * (aqq - app) / apq;
    const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1.0));
    const c = 1.0 / Math.sqrt(t * t + 1.0);
    const s = t * c;
    const tau = s / (1.0 + c);

    A[p * 3 + p] = app - t * apq;
    A[q * 3 + q] = aqq + t * apq;
    A[p * 3 + q] = 0.0;
    A[q * 3 + p] = 0.0;

    for (let r = 0; r < 3; r++) {
      if (r !== p && r !== q) {
        const arp = A[r * 3 + p];
        const arq = A[r * 3 + q];
        A[r * 3 + p] = arp - s * (arq + tau * arp);
        A[p * 3 + r] = A[r * 3 + p];
        A[r * 3 + q] = arq + s * (arp - tau * arq);
        A[q * 3 + r] = A[r * 3 + q];
      }
    }

    for (let r = 0; r < 3; r++) {
      const vrp = V[r * 3 + p];
      const vrq = V[r * 3 + q];
      V[r * 3 + p] = vrp - s * (vrq + tau * vrp);
      V[r * 3 + q] = vrq + s * (vrp - tau * vrq);
    }
  }

  const eigenvalues = [A[0], A[4], A[8]];
  const eigenvectors = [
    [V[0], V[3], V[6]],
    [V[1], V[4], V[7]],
    [V[2], V[5], V[8]]
  ];

  return { eigenvalues, eigenvectors };
}

// ============================================================================
// 2. CROSS-SECTIONAL CUTTING PLANE ENGINE
// ============================================================================

/**
 * Geometric cutting plane engine defining the Poincaré section hypersurface:
 * Sigma = { x in R^3 : n · (x - x_0) = 0 }
 * and an orthonormal 2D tangent coordinate frame (e_u, e_v).
 */
export class PoincarePlane {
  /**
   * @param {Array<number>} origin Point x_0 on the cutting plane [x, y, z] in Mpc/h.
   * @param {Array<number>} normal Normal unit vector n.
   * @param {Object} [options]
   * @param {Array<number>} [options.uAxis] Optional hint vector for u-axis.
   * @param {string} [options.name] Human-readable identifier.
   * @param {string} [options.orientation] PlaneOrientation enum.
   */
  constructor(origin = [0, 0, 0], normal = [0, 0, 1], options = {}) {
    if (!Array.isArray(origin) || origin.length < 3) {
      throw new TypeError('PoincarePlane: origin must be a 3-element array.');
    }
    if (!Array.isArray(normal) || normal.length < 3) {
      throw new TypeError('PoincarePlane: normal must be a 3-element array.');
    }

    this.origin = [Number(origin[0]), Number(origin[1]), Number(origin[2])];
    this.normal = vec3Normalize([Number(normal[0]), Number(normal[1]), Number(normal[2])]);
    this.name = options.name || 'CustomPlane';
    this.orientation = options.orientation || PlaneOrientation.CUSTOM_PLANE;

    this._setupTangentBasis(options.uAxis);
  }

  /**
   * Sets up a robust orthonormal 2D basis (e_u, e_v) perpendicular to normal.
   * @private
   * @param {Array<number>} [hintU]
   */
  _setupTangentBasis(hintU) {
    let u;
    if (hintU && Array.isArray(hintU) && hintU.length >= 3) {
      const proj = vec3Dot(hintU, this.normal);
      u = [hintU[0] - proj * this.normal[0], hintU[1] - proj * this.normal[1], hintU[2] - proj * this.normal[2]];
      if (vec3Norm(u) < 1e-6) {
        u = null;
      }
    }

    if (!u) {
      const nx = Math.abs(this.normal[0]);
      const ny = Math.abs(this.normal[1]);
      const nz = Math.abs(this.normal[2]);

      let arbitrary;
      if (nz < 0.9) {
        arbitrary = [0, 0, 1];
      } else if (ny < 0.9) {
        arbitrary = [0, 1, 0];
      } else {
        arbitrary = [1, 0, 0];
      }

      u = vec3Cross(this.normal, arbitrary);
    }

    this.uAxis = vec3Normalize(u);
    this.vAxis = vec3Normalize(vec3Cross(this.normal, this.uAxis));
  }

  /**
   * Factory: Supergalactic SGZ=0 cutting plane (SGX vs SGY).
   * @param {number} [z0=0] SGZ intercept.
   * @param {Array<number>} [origin=null] Center point.
   * @returns {PoincarePlane}
   */
  static createSGZ(z0 = 0.0, origin = null) {
    const orig = origin || [0.0, 0.0, z0];
    const plane = new PoincarePlane(orig, [0.0, 0.0, 1.0], {
      uAxis: [1.0, 0.0, 0.0],
      name: 'Supergalactic SGZ Plane',
      orientation: PlaneOrientation.SGZ_PLANE
    });
    plane.uAxis = [1.0, 0.0, 0.0];
    plane.vAxis = [0.0, 1.0, 0.0];
    return plane;
  }

  /**
   * Factory: Supergalactic SGY=0 cutting plane (SGX vs SGZ).
   * @param {number} [y0=0] SGY intercept.
   * @param {Array<number>} [origin=null] Center point.
   * @returns {PoincarePlane}
   */
  static createSGY(y0 = 0.0, origin = null) {
    const orig = origin || [0.0, y0, 0.0];
    const plane = new PoincarePlane(orig, [0.0, 1.0, 0.0], {
      uAxis: [1.0, 0.0, 0.0],
      name: 'Supergalactic SGY Plane',
      orientation: PlaneOrientation.SGY_PLANE
    });
    plane.uAxis = [1.0, 0.0, 0.0];
    plane.vAxis = [0.0, 0.0, 1.0];
    return plane;
  }

  /**
   * Factory: Supergalactic SGX=0 cutting plane (SGY vs SGZ).
   * @param {number} [x0=0] SGX intercept.
   * @param {Array<number>} [origin=null] Center point.
   * @returns {PoincarePlane}
   */
  static createSGX(x0 = 0.0, origin = null) {
    const orig = origin || [x0, 0.0, 0.0];
    const plane = new PoincarePlane(orig, [1.0, 0.0, 0.0], {
      uAxis: [0.0, 1.0, 0.0],
      name: 'Supergalactic SGX Plane',
      orientation: PlaneOrientation.SGX_PLANE
    });
    plane.uAxis = [0.0, 1.0, 0.0];
    plane.vAxis = [0.0, 0.0, 1.0];
    return plane;
  }

  /**
   * Factory: Cutting plane defined by 3 non-collinear 3D points.
   * @param {Array<number>} p1
   * @param {Array<number>} p2
   * @param {Array<number>} p3
   * @param {string} [name='ThreePointPlane']
   * @returns {PoincarePlane}
   */
  static fromThreePoints(p1, p2, p3, name = 'ThreePointPlane') {
    const v12 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    const v13 = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];
    const normal = vec3Cross(v12, v13);
    if (vec3Norm(normal) < 1e-12) {
      throw new Error('PoincarePlane.fromThreePoints: Points are collinear or coincident.');
    }
    return new PoincarePlane(p1, normal, { uAxis: v12, name, orientation: PlaneOrientation.CUSTOM_PLANE });
  }

  /**
   * Factory: Plane defined by point and arbitrary normal.
   * @param {Array<number>} origin
   * @param {Array<number>} normal
   * @param {Object} [options]
   * @returns {PoincarePlane}
   */
  static fromPointAndNormal(origin, normal, options = {}) {
    return new PoincarePlane(origin, normal, options);
  }

  /**
   * Factory: Principal cluster cutting plane aligned with tidal or inertia tensor eigenvectors.
   * @param {Array<number>} clusterCenter [x, y, z] in Mpc/h.
   * @param {Array<Array<number>>|Float64Array} tensor 3x3 symmetric tensor.
   * @param {number} [principalAxisIndex=2] Principal normal axis (0=minor, 1=intermediate, 2=major).
   * @returns {PoincarePlane}
   */
  static fromClusterTidalTensor(clusterCenter, tensor, principalAxisIndex = 2) {
    const eigen = jacobiEigenvalues3x3(tensor);
    const normal = eigen.eigenvectors[principalAxisIndex];
    const uIdx = (principalAxisIndex + 1) % 3;
    const uAxis = eigen.eigenvectors[uIdx];

    return new PoincarePlane(clusterCenter, normal, {
      uAxis,
      name: `Cluster Tidal Principal Plane (Axis ${principalAxisIndex})`,
      orientation: PlaneOrientation.PRINCIPAL_TIDAL_PLANE
    });
  }

  /**
   * Signed Euclidean distance from point x to the cutting plane: d(x) = (x - x_0) · n.
   * @param {Array<number>|Float64Array} pos [x, y, z] in Mpc/h.
   * @returns {number} Distance in Mpc/h.
   */
  signedDistance(pos) {
    return (pos[0] - this.origin[0]) * this.normal[0] +
           (pos[1] - this.origin[1]) * this.normal[1] +
           (pos[2] - this.origin[2]) * this.normal[2];
  }

  /**
   * Projects 3D spatial position x onto 2D planar coordinates (u, v).
   * u = (x - x_0) · e_u, v = (x - x_0) · e_v.
   *
   * @param {Array<number>|Float64Array} pos [x, y, z] in Mpc/h.
   * @returns {[number, number]} Planar coordinates [u, v] in Mpc/h.
   */
  projectPoint(pos) {
    const dx = pos[0] - this.origin[0];
    const dy = pos[1] - this.origin[1];
    const dz = pos[2] - this.origin[2];
    const u = dx * this.uAxis[0] + dy * this.uAxis[1] + dz * this.uAxis[2];
    const v = dx * this.vAxis[0] + dy * this.vAxis[1] + dz * this.vAxis[2];
    return [u, v];
  }

  /**
   * Unprojects 2D planar coordinates (u, v) and perpendicular distance back into 3D space.
   * x = x_0 + u * e_u + v * e_v + d * n.
   *
   * @param {number} u Planar u in Mpc/h.
   * @param {number} v Planar v in Mpc/h.
   * @param {number} [distance=0] Perpendicular offset in Mpc/h.
   * @returns {[number, number, number]} 3D coordinate [x, y, z] in Mpc/h.
   */
  unprojectPoint(u, v, distance = 0.0) {
    return [
      this.origin[0] + u * this.uAxis[0] + v * this.vAxis[0] + distance * this.normal[0],
      this.origin[1] + u * this.uAxis[1] + v * this.vAxis[1] + distance * this.normal[1],
      this.origin[2] + u * this.uAxis[2] + v * this.vAxis[2] + distance * this.normal[2]
    ];
  }

  /**
   * Projects 3D velocity vector onto planar tangent frame and normal component.
   *
   * @param {Array<number>|Float64Array} vel 3D velocity vector [vx, vy, vz] in km/s.
   * @returns {{ vu: number, vv: number, vPerp: number, speed: number }}
   */
  projectVelocity(vel) {
    const vu = vel[0] * this.uAxis[0] + vel[1] * this.uAxis[1] + vel[2] * this.uAxis[2];
    const vv = vel[0] * this.vAxis[0] + vel[1] * this.vAxis[1] + vel[2] * this.vAxis[2];
    const vPerp = vel[0] * this.normal[0] + vel[1] * this.normal[1] + vel[2] * this.normal[2];
    const speed = Math.sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    return { vu, vv, vPerp, speed };
  }

  /**
   * Tests whether a trajectory segment between posA and posB intersects the plane.
   *
   * @param {Array<number>} posA Position at step k.
   * @param {Array<number>} posB Position at step k+1.
   * @param {string} [direction=CrossingDirection.POSITIVE]
   * @returns {boolean} True if crossing occurred in requested direction.
   */
  isCrossing(posA, posB, direction = CrossingDirection.POSITIVE) {
    const dA = this.signedDistance(posA);
    const dB = this.signedDistance(posB);

    if (direction === CrossingDirection.POSITIVE) {
      return (dA < 0 && dB >= 0) || (dA <= 0 && dB > 0);
    } else if (direction === CrossingDirection.NEGATIVE) {
      return (dA > 0 && dB <= 0) || (dA >= 0 && dB < 0);
    } else {
      return (dA < 0 && dB >= 0) || (dA > 0 && dB <= 0);
    }
  }

  /**
   * Serializes plane configuration.
   * @returns {Object}
   */
  toJSON() {
    return {
      name: this.name,
      orientation: this.orientation,
      origin: [...this.origin],
      normal: [...this.normal],
      uAxis: [...this.uAxis],
      vAxis: [...this.vAxis]
    };
  }
}

// ============================================================================
// 3. PUNCTURE POINT DATA STRUCTURE & EXACT INTERSECTION SOLVER
// ============================================================================

/**
 * Data structure storing an exact puncture on the Poincaré cutting plane.
 */
export class PuncturePoint {
  /**
   * @param {Object} data
   */
  constructor(data = {}) {
    this.index = data.index !== undefined ? data.index : 0;
    this.position = data.position ? [...data.position] : [0, 0, 0]; // [x, y, z] in Mpc/h
    this.velocity = data.velocity ? [...data.velocity] : [0, 0, 0]; // [vx, vy, vz] in km/s
    this.planarCoords = data.planarCoords ? [...data.planarCoords] : [0, 0]; // [u, v] in Mpc/h
    this.planarVelocity = data.planarVelocity ? [...data.planarVelocity] : [0, 0]; // [vu, vv] in km/s
    this.normalVelocity = data.normalVelocity !== undefined ? data.normalVelocity : 0.0; // km/s
    this.time = data.time !== undefined ? data.time : 0.0; // Gyr
    this.arcLength = data.arcLength !== undefined ? data.arcLength : 0.0; // Mpc/h
    this.stepIndex = data.stepIndex !== undefined ? data.stepIndex : 0;
    this.incidentSpeed = data.incidentSpeed !== undefined ? data.incidentSpeed : vec3Norm(this.velocity);
    this.crossingAngle = data.crossingAngle !== undefined ? data.crossingAngle : 0.0; // radians
  }

  /**
   * Returns formatted summary string.
   * @returns {string}
   */
  toString() {
    return `Puncture #${this.index}: (u=${this.planarCoords[0].toFixed(4)}, v=${this.planarCoords[1].toFixed(4)}) Mpc/h, ` +
           `v_perp=${this.normalVelocity.toFixed(2)} km/s, t=${this.time.toFixed(4)} Gyr`;
  }
}

/**
 * High-precision exact root finding and cutting plane intersection solver.
 */
export class ExactIntersectionSolver {
  /**
   * Linear root interpolation between two trajectory states.
   *
   * @param {Array<number>} xA Position at step k.
   * @param {Array<number>} xB Position at step k+1.
   * @param {number} tA Time at step k.
   * @param {number} tB Time at step k+1.
   * @param {Array<number>} vA Velocity at step k.
   * @param {Array<number>} vB Velocity at step k+1.
   * @param {number} sA Arc-length at step k.
   * @param {number} sB Arc-length at step k+1.
   * @param {PoincarePlane} plane
   * @param {number} [stepIndex=0]
   * @param {number} [punctureIndex=0]
   * @returns {PuncturePoint}
   */
  static solveLinear(xA, xB, tA, tB, vA, vB, sA, sB, plane, stepIndex = 0, punctureIndex = 0) {
    const dA = plane.signedDistance(xA);
    const dB = plane.signedDistance(xB);
    const denom = dB - dA;
    const alpha = Math.abs(denom) > 1e-15 ? -dA / denom : 0.5;
    const clampedAlpha = Math.max(0.0, Math.min(1.0, alpha));

    const xCross = [
      xA[0] + clampedAlpha * (xB[0] - xA[0]),
      xA[1] + clampedAlpha * (xB[1] - xA[1]),
      xA[2] + clampedAlpha * (xB[2] - xA[2])
    ];

    const vCross = [
      vA[0] + clampedAlpha * (vB[0] - vA[0]),
      vA[1] + clampedAlpha * (vB[1] - vA[1]),
      vA[2] + clampedAlpha * (vB[2] - vA[2])
    ];

    const tCross = tA + clampedAlpha * (tB - tA);
    const sCross = sA + clampedAlpha * (sB - sA);

    const [u, v] = plane.projectPoint(xCross);
    const { vu, vv, vPerp, speed } = plane.projectVelocity(vCross);
    const crossingAngle = speed > 1e-12 ? Math.asin(Math.max(-1.0, Math.min(1.0, vPerp / speed))) : 0.0;

    return new PuncturePoint({
      index: punctureIndex,
      position: xCross,
      velocity: vCross,
      planarCoords: [u, v],
      planarVelocity: [vu, vv],
      normalVelocity: vPerp,
      time: tCross,
      arcLength: sCross,
      stepIndex,
      incidentSpeed: speed,
      crossingAngle
    });
  }

  /**
   * Cubic Hermite spline exact root solver:
   * Uses positions xA, xB and velocities vA, vB to construct a C^1 continuous trajectory
   * p(theta), theta in [0, 1], and solves the scalar equation d(p(theta)) = 0 via Newton-Raphson
   * with Brent bisection fallback.
   *
   * @param {Array<number>} xA
   * @param {Array<number>} xB
   * @param {number} tA
   * @param {number} tB
   * @param {Array<number>} vA
   * @param {Array<number>} vB
   * @param {number} sA
   * @param {number} sB
   * @param {PoincarePlane} plane
   * @param {Object} [options]
   * @param {VelocityField} [options.field=null]
   * @param {number} [options.tol=1e-8]
   * @param {number} [options.maxIter=30]
   * @param {number} [options.stepIndex=0]
   * @param {number} [options.punctureIndex=0]
   * @returns {PuncturePoint}
   */
  static solveHermiteCubic(xA, xB, tA, tB, vA, vB, sA, sB, plane, options = {}) {
    const dt = tB - tA;
    if (Math.abs(dt) < 1e-15) {
      return this.solveLinear(xA, xB, tA, tB, vA, vB, sA, sB, plane, options.stepIndex, options.punctureIndex);
    }

    const tol = options.tol || 1e-8;
    const maxIter = options.maxIter || 30;

    const m0 = [vA[0] * dt, vA[1] * dt, vA[2] * dt];
    const m1 = [vB[0] * dt, vB[1] * dt, vB[2] * dt];

    const dA = plane.signedDistance(xA);
    const dB = plane.signedDistance(xB);
    const dPrimeA = plane.signedDistance(m0) - plane.signedDistance([0, 0, 0]);
    const dPrimeB = plane.signedDistance(m1) - plane.signedDistance([0, 0, 0]);

    const c0 = dA;
    const c1 = dPrimeA;
    const c2 = 3.0 * (dB - dA) - 2.0 * dPrimeA - dPrimeB;
    const c3 = 2.0 * (dA - dB) + dPrimeA + dPrimeB;

    const evalD = (th) => ((c3 * th + c2) * th + c1) * th + c0;

    // Brent's method root finding in [0, 1]
    let a = 0.0, b = 1.0;
    let fa = evalD(a), fb = evalD(b);

    let thetaRoot = 0.5;

    if (fa * fb <= 0.0) {
      let c = a, fc = fa, d = 0.0, e = 0.0;
      for (let iter = 0; iter < maxIter; iter++) {
        if (Math.abs(fb) < tol) {
          thetaRoot = b;
          break;
        }

        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
          c = a; fc = fa;
          d = b - a; e = d;
        }

        if (Math.abs(fc) < Math.abs(fb)) {
          a = b; b = c; c = a;
          fa = fb; fb = fc; fc = fa;
        }

        const tol1 = 2.0 * 1e-12 * Math.abs(b) + 0.5 * tol;
        const xm = 0.5 * (c - b);

        if (Math.abs(xm) <= tol1 || fb === 0.0) {
          thetaRoot = b;
          break;
        }

        if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
          let s = fb / fa;
          let p, q;
          if (a === c) {
            p = 2.0 * xm * s;
            q = 1.0 - s;
          } else {
            q = fa / fc;
            const r = fb / fc;
            p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
            q = (q - 1.0) * (r - 1.0) * (s - 1.0);
          }
          if (p > 0.0) q = -q;
          p = Math.abs(p);

          if (2.0 * p < Math.min(3.0 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) {
            e = d;
            d = p / q;
          } else {
            d = xm;
            e = d;
          }
        } else {
          d = xm;
          e = d;
        }

        a = b; fa = fb;
        if (Math.abs(d) > tol1) {
          b += d;
        } else {
          b += xm >= 0.0 ? tol1 : -tol1;
        }
        fb = evalD(b);
        thetaRoot = b;
      }
    } else {
      thetaRoot = Math.abs(dB - dA) > 1e-15 ? -dA / (dB - dA) : 0.5;
    }

    thetaRoot = Math.max(0.0, Math.min(1.0, thetaRoot));

    const th2 = thetaRoot * thetaRoot;
    const th3 = th2 * thetaRoot;
    const h00 = 2.0 * th3 - 3.0 * th2 + 1.0;
    const h10 = th3 - 2.0 * th2 + thetaRoot;
    const h01 = -2.0 * th3 + 3.0 * th2;
    const h11 = th3 - th2;

    const xCross = [
      h00 * xA[0] + h10 * m0[0] + h01 * xB[0] + h11 * m1[0],
      h00 * xA[1] + h10 * m0[1] + h01 * xB[1] + h11 * m1[1],
      h00 * xA[2] + h10 * m0[2] + h01 * xB[2] + h11 * m1[2]
    ];

    let vCross;
    if (options.field && typeof options.field.sampleVelocity === 'function') {
      vCross = options.field.sampleVelocity(xCross[0], xCross[1], xCross[2]);
    } else {
      const dh00 = 6.0 * th2 - 6.0 * thetaRoot;
      const dh10 = 3.0 * th2 - 4.0 * thetaRoot + 1.0;
      const dh01 = -6.0 * th2 + 6.0 * thetaRoot;
      const dh11 = 3.0 * th2 - 2.0 * thetaRoot;

      vCross = [
        (dh00 * xA[0] + dh10 * m0[0] + dh01 * xB[0] + dh11 * m1[0]) / dt,
        (dh00 * xA[1] + dh10 * m0[1] + dh01 * xB[1] + dh11 * m1[1]) / dt,
        (dh00 * xA[2] + dh10 * m0[2] + dh01 * xB[2] + dh11 * m1[2]) / dt
      ];
    }

    const tCross = tA + thetaRoot * dt;
    const sCross = sA + thetaRoot * (sB - sA);

    const [u, v] = plane.projectPoint(xCross);
    const { vu, vv, vPerp, speed } = plane.projectVelocity(vCross);
    const crossingAngle = speed > 1e-12 ? Math.asin(Math.max(-1.0, Math.min(1.0, vPerp / speed))) : 0.0;

    return new PuncturePoint({
      index: options.punctureIndex || 0,
      position: xCross,
      velocity: vCross,
      planarCoords: [u, v],
      planarVelocity: [vu, vv],
      normalVelocity: vPerp,
      time: tCross,
      arcLength: sCross,
      stepIndex: options.stepIndex || 0,
      incidentSpeed: speed,
      crossingAngle
    });
  }
}

// ============================================================================
// 4. POINCARÉ RECURRENCE MAP & RQA METRICS ENGINE
// ============================================================================

/**
 * Full Poincaré recurrence mapping, return time statistics, symplectic area tests,
 * and complete Recurrence Quantification Analysis (RQA).
 */
export class PoincareRecurrenceMap {
  /**
   * @param {Array<PuncturePoint>} punctures Array of ordered puncture points.
   * @param {PoincarePlane} plane Associated cutting plane.
   */
  constructor(punctures = [], plane = null) {
    this.punctures = punctures;
    this.plane = plane;
    this.numPunctures = punctures.length;
  }

  /**
   * Returns array of planar puncture coordinates [[u0, v0], [u1, v1], ...].
   * @returns {Array<[number, number]>}
   */
  getPlanarPoints() {
    return this.punctures.map(p => [p.planarCoords[0], p.planarCoords[1]]);
  }

  /**
   * Returns consecutive return times Delta t_k = t_{k+1} - t_k.
   * @returns {Array<number>}
   */
  getReturnTimes() {
    if (this.numPunctures < 2) return [];
    const dt = [];
    for (let i = 0; i < this.numPunctures - 1; i++) {
      dt.push(this.punctures[i + 1].time - this.punctures[i].time);
    }
    return dt;
  }

  /**
   * Returns consecutive return arc lengths Delta s_k = s_{k+1} - s_k.
   * @returns {Array<number>}
   */
  getReturnArcLengths() {
    if (this.numPunctures < 2) return [];
    const ds = [];
    for (let i = 0; i < this.numPunctures - 1; i++) {
      ds.push(this.punctures[i + 1].arcLength - this.punctures[i].arcLength);
    }
    return ds;
  }

  /**
   * Computes statistical moments of return times.
   * @returns {Object}
   */
  computeReturnTimeStats() {
    const dts = this.getReturnTimes();
    if (dts.length === 0) {
      return { mean: 0, std: 0, variance: 0, min: 0, max: 0, count: 0, skewness: 0, kurtosis: 0 };
    }

    const n = dts.length;
    let sum = 0;
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < n; i++) {
      const v = dts[i];
      sum += v;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    const mean = sum / n;

    let varSum = 0;
    let skewSum = 0;
    let kurtSum = 0;
    for (let i = 0; i < n; i++) {
      const diff = dts[i] - mean;
      const d2 = diff * diff;
      varSum += d2;
      skewSum += d2 * diff;
      kurtSum += d2 * d2;
    }
    const variance = varSum / n;
    const std = Math.sqrt(variance);
    const skewness = std > 1e-12 ? (skewSum / n) / (std * std * std) : 0.0;
    const kurtosis = std > 1e-12 ? (kurtSum / n) / (std * std * std * std) - 3.0 : 0.0;

    return { mean, std, variance, min, max, count: n, skewness, kurtosis };
  }

  /**
   * Computes centroid (u_c, v_c) of the Poincaré punctures.
   * @returns {[number, number]}
   */
  computeCentroid() {
    if (this.numPunctures === 0) return [0, 0];
    let su = 0, sv = 0;
    for (let i = 0; i < this.numPunctures; i++) {
      su += this.punctures[i].planarCoords[0];
      sv += this.punctures[i].planarCoords[1];
    }
    return [su / this.numPunctures, sv / this.numPunctures];
  }

  /**
   * Computes rotation / winding number around the orbit centroid:
   * omega = 1/(2*pi * (N-1)) * sum_{k=0}^{N-2} Delta theta_k.
   *
   * @returns {number} Winding number in turns per puncture.
   */
  computeWindingNumber() {
    if (this.numPunctures < 3) return 0.0;
    const [uc, vc] = this.computeCentroid();

    let totalAngle = 0.0;
    let prevAngle = Math.atan2(this.punctures[0].planarCoords[1] - vc, this.punctures[0].planarCoords[0] - uc);

    for (let i = 1; i < this.numPunctures; i++) {
      const curAngle = Math.atan2(this.punctures[i].planarCoords[1] - vc, this.punctures[i].planarCoords[0] - uc);
      let dTheta = curAngle - prevAngle;

      while (dTheta > Math.PI) dTheta -= 2.0 * Math.PI;
      while (dTheta < -Math.PI) dTheta += 2.0 * Math.PI;

      totalAngle += dTheta;
      prevAngle = curAngle;
    }

    return totalAngle / (2.0 * Math.PI * (this.numPunctures - 1));
  }

  /**
   * Local least-squares Jacobian matrix estimation of the Poincaré return map (u_{k+1}, v_{k+1}) = P(u_k, v_k).
   * Tests symplectic area conservation: |det(J)| ~ 1.
   *
   * @returns {{ jacobian: [[number, number], [number, number]], det: number, trace: number, isAreaPreserving: boolean }}
   */
  estimateAverageReturnJacobian() {
    if (this.numPunctures < 6) {
      return { jacobian: [[1, 0], [0, 1]], det: 1.0, trace: 2.0, isAreaPreserving: true };
    }

    const N = this.numPunctures - 1;
    let s_uu = 0, s_uv = 0, s_vv = 0;
    let s_u_next_u = 0, s_u_next_v = 0;
    let s_v_next_u = 0, s_v_next_v = 0;

    const [uc, vc] = this.computeCentroid();

    for (let i = 0; i < N; i++) {
      const u = this.punctures[i].planarCoords[0] - uc;
      const v = this.punctures[i].planarCoords[1] - vc;
      const uNext = this.punctures[i + 1].planarCoords[0] - uc;
      const vNext = this.punctures[i + 1].planarCoords[1] - vc;

      s_uu += u * u;
      s_uv += u * v;
      s_vv += v * v;

      s_u_next_u += uNext * u;
      s_u_next_v += uNext * v;
      s_v_next_u += vNext * u;
      s_v_next_v += vNext * v;
    }

    const covDet = s_uu * s_vv - s_uv * s_uv;
    if (Math.abs(covDet) < 1e-12) {
      return { jacobian: [[1, 0], [0, 1]], det: 1.0, trace: 2.0, isAreaPreserving: true };
    }

    const invDet = 1.0 / covDet;
    const j11 = (s_u_next_u * s_vv - s_u_next_v * s_uv) * invDet;
    const j12 = (s_u_next_v * s_uu - s_u_next_u * s_uv) * invDet;
    const j21 = (s_v_next_u * s_vv - s_v_next_v * s_uv) * invDet;
    const j22 = (s_v_next_v * s_uu - s_v_next_u * s_uv) * invDet;

    const det = j11 * j22 - j12 * j21;
    const trace = j11 + j22;
    const isAreaPreserving = Math.abs(Math.abs(det) - 1.0) < 0.18;

    return {
      jacobian: [[j11, j12], [j21, j22]],
      det,
      trace,
      isAreaPreserving
    };
  }

  /**
   * Recurrence Quantification Analysis (RQA):
   * Computes recurrence plot distance matrix R_{i,j} = Theta(eps - ||P_i - P_j||) and metrics:
   * - Recurrence Rate (RR)
   * - Determinism (DET)
   * - Laminarity (LAM)
   * - Trapping Time (TT)
   * - Longest diagonal line (L_max)
   * - Shannon diagonal entropy (ENTR)
   * - Recurrence Trend (TREND)
   *
   * @param {Object} [options]
   * @param {number} [options.epsilon=null] Neighborhood threshold (Mpc/h).
   * @param {number} [options.lMin=2] Minimum diagonal line length.
   * @param {number} [options.vMin=2] Minimum vertical line length.
   * @returns {Object} Complete RQA metrics dictionary.
   */
  computeRQA(options = {}) {
    const N = this.numPunctures;
    if (N < 4) {
      return {
        recurrenceRate: 0.0,
        determinism: 0.0,
        laminarity: 0.0,
        trappingTime: 0.0,
        longestDiagonalLine: 0,
        shannonEntropy: 0.0,
        divergence: 0.0,
        trend: 0.0,
        numPunctures: N,
        epsilon: 0.0
      };
    }

    const pts = this.getPlanarPoints();

    const distMat = new Float64Array(N * N);
    let maxDist = 0.0;
    for (let i = 0; i < N; i++) {
      for (let j = i; j < N; j++) {
        const du = pts[i][0] - pts[j][0];
        const dv = pts[i][1] - pts[j][1];
        const d = Math.sqrt(du * du + dv * dv);
        distMat[i * N + j] = d;
        distMat[j * N + i] = d;
        if (d > maxDist) maxDist = d;
      }
    }

    const eps = options.epsilon || Math.max(0.01, 0.15 * maxDist);
    const lMin = options.lMin || 2;
    const vMin = options.vMin || 2;

    const recMat = new Uint8Array(N * N);
    let totalRec = 0;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        if (distMat[i * N + j] <= eps) {
          recMat[i * N + j] = 1;
          totalRec++;
        }
      }
    }

    const recurrenceRate = totalRec / (N * N);

    // Diagonal lines analysis (excluding main diagonal i == j)
    const diagonalLengths = [];
    let lMax = 0;
    const diagDensities = [];

    for (let k = 1; k < N; k++) {
      let curLen = 0;
      let diagRec = 0;
      const diagTotal = N - k;

      for (let i = 0; i < diagTotal; i++) {
        const j = i + k;
        if (recMat[i * N + j] === 1) {
          curLen++;
          diagRec++;
        } else {
          if (curLen >= lMin) {
            diagonalLengths.push(curLen);
            if (curLen > lMax) lMax = curLen;
          }
          curLen = 0;
        }
      }
      if (curLen >= lMin) {
        diagonalLengths.push(curLen);
        if (curLen > lMax) lMax = curLen;
      }

      diagDensities.push({ k, density: diagRec / diagTotal });
    }

    let diagPoints = 0;
    for (let i = 0; i < diagonalLengths.length; i++) {
      diagPoints += diagonalLengths[i];
    }
    const offDiagTotal = (totalRec - N) / 2.0;
    const determinism = offDiagTotal > 0 ? Math.min(1.0, diagPoints / offDiagTotal) : 0.0;

    let shannonEntropy = 0.0;
    if (diagonalLengths.length > 0) {
      const freqMap = new Map();
      for (let i = 0; i < diagonalLengths.length; i++) {
        const len = diagonalLengths[i];
        freqMap.set(len, (freqMap.get(len) || 0) + 1);
      }
      const totalLines = diagonalLengths.length;
      for (const count of freqMap.values()) {
        const p = count / totalLines;
        shannonEntropy -= p * Math.log(p);
      }
    }

    // Vertical lines analysis
    const verticalLengths = [];
    for (let j = 0; j < N; j++) {
      let curLen = 0;
      for (let i = 0; i < N; i++) {
        if (i !== j && recMat[i * N + j] === 1) {
          curLen++;
        } else {
          if (curLen >= vMin) {
            verticalLengths.push(curLen);
          }
          curLen = 0;
        }
      }
      if (curLen >= vMin) {
        verticalLengths.push(curLen);
      }
    }

    let vertPoints = 0;
    for (let i = 0; i < verticalLengths.length; i++) {
      vertPoints += verticalLengths[i];
    }
    const laminarity = offDiagTotal > 0 ? Math.min(1.0, vertPoints / offDiagTotal) : 0.0;
    const trappingTime = verticalLengths.length > 0 ? vertPoints / verticalLengths.length : 0.0;
    const divergence = lMax > 0 ? 1.0 / lMax : 0.0;

    // Linear regression for TREND (drift of recurrence density away from diagonal)
    let trend = 0.0;
    if (diagDensities.length > 2) {
      let sumK = 0, sumD = 0, sumKD = 0, sumK2 = 0;
      const numD = diagDensities.length;
      for (let i = 0; i < numD; i++) {
        const k = diagDensities[i].k;
        const d = diagDensities[i].density;
        sumK += k;
        sumD += d;
        sumKD += k * d;
        sumK2 += k * k;
      }
      const denom = numD * sumK2 - sumK * sumK;
      if (Math.abs(denom) > 1e-12) {
        trend = (numD * sumKD - sumK * sumD) / denom;
      }
    }

    return {
      recurrenceRate,
      determinism,
      laminarity,
      trappingTime,
      longestDiagonalLine: lMax,
      shannonEntropy,
      divergence,
      trend,
      numPunctures: N,
      epsilon: eps
    };
  }
}

// ============================================================================
// 5. LYAPUNOV EXPONENT ESTIMATORS (BENETTIN SHADOW & VARIATIONAL TANGENT QR)
// ============================================================================

/**
 * High-accuracy multi-method Lyapunov exponent estimator.
 */
export class LyapunovExponentEstimator {
  /**
   * Benettin two-particle shadow trajectory tracking algorithm:
   * Evolving reference trajectory x_ref(t) and shadow x_shadow(t) with initial
   * separation d0 ~ 1e-5 Mpc/h.
   *
   * @param {VelocityField|Function} field Velocity field or function (x, y, z) => [vx, vy, vz].
   * @param {Array<number>} initialSeed [x, y, z] in Mpc/h.
   * @param {Object} [options]
   * @param {number} [options.d0=1e-5] Initial perturbation amplitude (Mpc/h).
   * @param {number} [options.totalTime=50.0] Total integration time (Gyr).
   * @param {number} [options.dt=0.05] Time step (Gyr).
   * @param {number} [options.renormInterval=1] Steps between renormalization.
   * @param {Array<number>} [options.perturbationDir] Perturbation unit vector.
   * @returns {Object} { lambdaMax, ftleHistory, lyapunovTime, isChaotic, finalSeparation }
   */
  static estimateBenettinMLE(field, initialSeed, options = {}) {
    const d0 = options.d0 || 1e-5;
    const totalTime = options.totalTime || 50.0;
    const dt = options.dt || 0.05;
    const renormInterval = options.renormInterval || 1;
    const totalSteps = Math.floor(totalTime / dt);

    const velSample = typeof field === 'function' ? field : (x, y, z) => field.sampleVelocity(x, y, z);

    let pDir = options.perturbationDir ? vec3Normalize(options.perturbationDir) : [1.0, 0.0, 0.0];

    let xRef = [initialSeed[0], initialSeed[1], initialSeed[2]];
    let xShadow = [xRef[0] + d0 * pDir[0], xRef[1] + d0 * pDir[1], xRef[2] + d0 * pDir[2]];

    let logSum = 0.0;
    let accumulatedTime = 0.0;
    const ftleHistory = [];

    for (let step = 0; step < totalSteps; step++) {
      xRef = rk4Step(xRef, dt, velSample);
      xShadow = rk4Step(xShadow, dt, velSample);

      accumulatedTime += dt;

      if ((step + 1) % renormInterval === 0) {
        const dx = [xShadow[0] - xRef[0], xShadow[1] - xRef[1], xShadow[2] - xRef[2]];
        const dCurrent = vec3Norm(dx);

        if (dCurrent > 1e-18) {
          const ratio = dCurrent / d0;
          logSum += Math.log(ratio);

          const invD = 1.0 / dCurrent;
          xShadow[0] = xRef[0] + d0 * (dx[0] * invD);
          xShadow[1] = xRef[1] + d0 * (dx[1] * invD);
          xShadow[2] = xRef[2] + d0 * (dx[2] * invD);
        }

        const currentMLE = accumulatedTime > 0 ? (logSum / accumulatedTime) : 0.0;
        ftleHistory.push({
          time: accumulatedTime,
          mle: currentMLE,
          dCurrent
        });
      }
    }

    const lambdaMax = accumulatedTime > 0 ? (logSum / accumulatedTime) : 0.0;
    const lyapunovTime = lambdaMax > 1e-5 ? (1.0 / lambdaMax) : Infinity;
    const isChaotic = lambdaMax > 0.05;

    return {
      lambdaMax,
      ftleHistory,
      lyapunovTime,
      isChaotic,
      totalTime: accumulatedTime,
      logExpansionSum: logSum
    };
  }

  /**
   * Continuous Variational Tangent Map Integration with Continuous QR Decomposition:
   * Solves dx/dt = v(x) and dPhi/dt = J(x) * Phi, where J(x) = grad v is the velocity gradient tensor.
   * Provides the complete 3D Lyapunov spectrum (lambda_1 >= lambda_2 >= lambda_3),
   * Kolmogorov-Sinai entropy, and Kaplan-Yorke attractor dimension.
   *
   * @param {VelocityField|Function} field
   * @param {Array<number>} initialSeed
   * @param {Object} [options]
   * @param {number} [options.totalTime=50.0]
   * @param {number} [options.dt=0.05]
   * @param {number} [options.qrInterval=1] Steps between Gram-Schmidt QR orthonormalizations.
   * @param {number} [options.fdH=1e-4] Finite difference step for Jacobian grad v.
   * @returns {Object} { spectrum: [lambda1, lambda2, lambda3], ksEntropy, lyapunovDimension, history }
   */
  static estimateVariationalSpectrum(field, initialSeed, options = {}) {
    const totalTime = options.totalTime || 50.0;
    const dt = options.dt || 0.05;
    const qrInterval = options.qrInterval || 1;
    const fdH = options.fdH || 1e-4;
    const totalSteps = Math.floor(totalTime / dt);

    const velSample = typeof field === 'function' ? field : (x, y, z) => field.sampleVelocity(x, y, z);

    let x = [initialSeed[0], initialSeed[1], initialSeed[2]];
    let Phi = new Float64Array([
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ]);

    const cumLogs = [0.0, 0.0, 0.0];
    let accumulatedTime = 0.0;
    const history = [];

    for (let step = 0; step < totalSteps; step++) {
      const J = computeJacobian3x3(x, velSample, fdH);

      const k1_x = velSample(x[0], x[1], x[2]);
      const k1_Phi = mat3Mul(J, Phi);

      const x2 = [x[0] + 0.5 * dt * k1_x[0], x[1] + 0.5 * dt * k1_x[1], x[2] + 0.5 * dt * k1_x[2]];
      const Phi2 = addScaledMat3(Phi, k1_Phi, 0.5 * dt);
      const J2 = computeJacobian3x3(x2, velSample, fdH);
      const k2_x = velSample(x2[0], x2[1], x2[2]);
      const k2_Phi = mat3Mul(J2, Phi2);

      const x3 = [x[0] + 0.5 * dt * k2_x[0], x[1] + 0.5 * dt * k2_x[1], x[2] + 0.5 * dt * k2_x[2]];
      const Phi3 = addScaledMat3(Phi, k2_Phi, 0.5 * dt);
      const J3 = computeJacobian3x3(x3, velSample, fdH);
      const k3_x = velSample(x3[0], x3[1], x3[2]);
      const k3_Phi = mat3Mul(J3, Phi3);

      const x4 = [x[0] + dt * k3_x[0], x[1] + dt * k3_x[1], x[2] + dt * k3_x[2]];
      const Phi4 = addScaledMat3(Phi, k3_Phi, dt);
      const J4 = computeJacobian3x3(x4, velSample, fdH);
      const k4_x = velSample(x4[0], x4[1], x4[2]);
      const k4_Phi = mat3Mul(J4, Phi4);

      x[0] += (dt / 6.0) * (k1_x[0] + 2.0 * k2_x[0] + 2.0 * k3_x[0] + k4_x[0]);
      x[1] += (dt / 6.0) * (k1_x[1] + 2.0 * k2_x[1] + 2.0 * k3_x[1] + k4_x[1]);
      x[2] += (dt / 6.0) * (k1_x[2] + 2.0 * k2_x[2] + 2.0 * k3_x[2] + k4_x[2]);

      for (let i = 0; i < 9; i++) {
        Phi[i] += (dt / 6.0) * (k1_Phi[i] + 2.0 * k2_Phi[i] + 2.0 * k3_Phi[i] + k4_Phi[i]);
      }

      accumulatedTime += dt;

      if ((step + 1) % qrInterval === 0) {
        const { Q, R } = qrDecomposition3x3(Phi);
        const r0 = Math.abs(R[0]);
        const r1 = Math.abs(R[4]);
        const r2 = Math.abs(R[8]);

        if (r0 > 1e-18) cumLogs[0] += Math.log(r0);
        if (r1 > 1e-18) cumLogs[1] += Math.log(r1);
        if (r2 > 1e-18) cumLogs[2] += Math.log(r2);

        for (let i = 0; i < 9; i++) {
          Phi[i] = Q[i];
        }

        if (accumulatedTime > 0) {
          history.push({
            time: accumulatedTime,
            lambda1: cumLogs[0] / accumulatedTime,
            lambda2: cumLogs[1] / accumulatedTime,
            lambda3: cumLogs[2] / accumulatedTime
          });
        }
      }
    }

    const tInv = accumulatedTime > 0 ? 1.0 / accumulatedTime : 1.0;
    const l1 = cumLogs[0] * tInv;
    const l2 = cumLogs[1] * tInv;
    const l3 = cumLogs[2] * tInv;

    const spectrum = [l1, l2, l3].sort((a, b) => b - a);

    let ksEntropy = 0.0;
    for (let i = 0; i < 3; i++) {
      if (spectrum[i] > 0) ksEntropy += spectrum[i];
    }

    let lyapunovDimension = 0.0;
    if (spectrum[0] <= 0) {
      lyapunovDimension = 0.0;
    } else if (spectrum[0] + spectrum[1] <= 0) {
      lyapunovDimension = 1.0 + spectrum[0] / Math.abs(spectrum[1]);
    } else if (spectrum[0] + spectrum[1] + spectrum[2] <= 0) {
      lyapunovDimension = 2.0 + (spectrum[0] + spectrum[1]) / Math.abs(spectrum[2]);
    } else {
      lyapunovDimension = 3.0;
    }

    return {
      spectrum,
      lambdaMax: spectrum[0],
      ksEntropy,
      lyapunovDimension,
      isChaotic: spectrum[0] > 0.05,
      totalTime: accumulatedTime,
      history
    };
  }
}

// ============================================================================
// 6. INVARIANT TORI, FIXED POINTS & BASIN TOPOLOGY CLASSIFIER
// ============================================================================

/**
 * Classifies the dynamical morphology of streamline orbits on the Poincaré section:
 * - KAM invariant tori (smooth 1D closed loops)
 * - Resonance island chains (q-order period-q islands with angular gaps)
 * - Elliptic O-points & Hyperbolic X-points
 * - Chaotic stochastic seas
 */
export class InvariantTorusClassifier {
  /**
   * Classifies the orbit represented by a set of Poincaré punctures.
   *
   * @param {Array<PuncturePoint>} punctures
   * @param {Object} [options]
   * @param {number} [options.lambdaMax=0.0] Precomputed MLE (Gyr^-1).
   * @param {number} [options.chaoticThreshold=0.08]
   * @param {number} [options.loopDispersalTol=0.35]
   * @returns {Object} Detailed classification report.
   */
  static classifyOrbit(punctures, options = {}) {
    const N = punctures.length;
    const lambdaMax = options.lambdaMax !== undefined ? options.lambdaMax : 0.0;
    const chaoticThresh = options.chaoticThreshold || 0.08;
    const loopDispersalTol = options.loopDispersalTol || 0.35;

    if (N < 2) {
      return {
        classification: OrbitClassification.INSUFFICIENT_DATA,
        confidence: 0.0,
        numPunctures: N,
        rotationNumber: 0.0,
        islandCount: 0,
        radialDispersion: 0.0,
        details: 'Fewer than 2 punctures detected on Poincaré section.'
      };
    }

    const pts = punctures.map(p => [p.planarCoords[0], p.planarCoords[1]]);

    let su = 0, sv = 0;
    for (let i = 0; i < N; i++) {
      su += pts[i][0];
      sv += pts[i][1];
    }
    const uc = su / N;
    const vc = sv / N;

    const radii = [];
    let rSum = 0;
    for (let i = 0; i < N; i++) {
      const r = Math.hypot(pts[i][0] - uc, pts[i][1] - vc);
      radii.push(r);
      rSum += r;
    }
    const rMean = rSum / N;

    let rVar = 0;
    for (let i = 0; i < N; i++) {
      const dr = radii[i] - rMean;
      rVar += dr * dr;
    }
    const rStd = Math.sqrt(rVar / N);
    const radialDispersion = rMean > 1e-9 ? (rStd / rMean) : 0.0;

    // 1. Strict discrete limit cycle check (small discrete point clusters)
    const discreteClusters = clusterPunctures(pts, 0.08 * (rMean + 1.0));
    if (discreteClusters.length >= 1 && discreteClusters.length <= 8 && N >= discreteClusters.length * 3) {
      let isStrictCycle = true;
      for (const cl of discreteClusters) {
        if (cl.radius > 0.15 * rMean) {
          isStrictCycle = false;
          break;
        }
      }
      if (isStrictCycle) {
        return {
          classification: OrbitClassification.PERIODIC_LIMIT_CYCLE,
          confidence: 0.95,
          numPunctures: N,
          rotationNumber: 1.0 / discreteClusters.length,
          islandCount: discreteClusters.length,
          radialDispersion,
          details: `Discrete period-${discreteClusters.length} limit cycle.`
        };
      }
    }

    // 2. High Lyapunov exponent or severe radial scattering -> Chaotic Sea
    if (lambdaMax >= chaoticThresh || radialDispersion > 0.65) {
      return {
        classification: OrbitClassification.CHAOTIC_SEA_STOCHASTIC,
        confidence: 0.90,
        numPunctures: N,
        rotationNumber: computeWindingFromPoints(pts, uc, vc),
        islandCount: 0,
        radialDispersion,
        details: `Chaotic stochastic orbital sea (lambda_max = ${lambdaMax.toFixed(3)} Gyr^-1).`
      };
    }

    // 3. Check for resonant island chains (must have true angular gaps between islands)
    const islandClusters = detectIslandChainsWithGaps(pts, uc, vc);
    if (islandClusters.length >= 2 && islandClusters.length <= 12) {
      return {
        classification: OrbitClassification.ISLAND_CHAIN_RESONANCE,
        confidence: 0.90,
        numPunctures: N,
        rotationNumber: 1.0 / islandClusters.length,
        islandCount: islandClusters.length,
        radialDispersion,
        details: `Resonance island chain with ${islandClusters.length} invariant islands.`
      };
    }

    // 4. Continuous KAM Invariant Torus (smooth 1D closed loop)
    const angularCoverage = computeAngularCoverage(pts, uc, vc);
    if (radialDispersion < loopDispersalTol && angularCoverage > 0.60 && lambdaMax < chaoticThresh) {
      return {
        classification: OrbitClassification.QUASI_PERIODIC_INVARIANT_TORUS,
        confidence: 0.92,
        numPunctures: N,
        rotationNumber: computeWindingFromPoints(pts, uc, vc),
        islandCount: 1,
        radialDispersion,
        details: 'Smooth 1D KAM invariant torus cross-section in accretion basin.'
      };
    }

    // Fallback: regular quasi-periodic or chaotic
    if (lambdaMax >= 0.05) {
      return {
        classification: OrbitClassification.CHAOTIC_SEA_STOCHASTIC,
        confidence: 0.75,
        numPunctures: N,
        rotationNumber: computeWindingFromPoints(pts, uc, vc),
        islandCount: 0,
        radialDispersion,
        details: `Weakly chaotic flow (lambda_max = ${lambdaMax.toFixed(3)} Gyr^-1).`
      };
    }

    return {
      classification: OrbitClassification.QUASI_PERIODIC_INVARIANT_TORUS,
      confidence: 0.75,
      numPunctures: N,
      rotationNumber: computeWindingFromPoints(pts, uc, vc),
      islandCount: 1,
      radialDispersion,
      details: 'Quasi-periodic regular orbit.'
    };
  }

  /**
   * Identifies candidate Elliptic (O-point) and Hyperbolic (X-point) fixed points
   * in the Poincaré section plane.
   *
   * @param {PoincareRecurrenceMap} recurrenceMap
   * @returns {{ ellipticPoints: Array<Object>, hyperbolicPoints: Array<Object> }}
   */
  static findFixedPoints(recurrenceMap) {
    const jac = recurrenceMap.estimateAverageReturnJacobian();
    const centroid = recurrenceMap.computeCentroid();
    const trace = jac.trace;
    const det = jac.det;

    const discriminant = trace * trace - 4.0 * det;

    const ellipticPoints = [];
    const hyperbolicPoints = [];

    if (discriminant < 0.0 && Math.abs(det - 1.0) < 0.3) {
      ellipticPoints.push({
        position2D: centroid,
        trace,
        det,
        stability: FixedPointType.ELLIPTIC_STABLE_CENTER,
        description: 'Elliptic O-point center of invariant nested KAM tori'
      });
    } else if (Math.abs(trace) > 2.0 && Math.abs(det - 1.0) < 0.35) {
      hyperbolicPoints.push({
        position2D: centroid,
        trace,
        det,
        stability: FixedPointType.HYPERBOLIC_SADDLE,
        description: 'Hyperbolic X-point saddle with separatrix manifold tangles'
      });
    }

    return { ellipticPoints, hyperbolicPoints };
  }
}

// ============================================================================
// 7. MAIN HIGH-LEVEL POINCARÉ SECTION TRACER ORCHESTRATOR
// ============================================================================

/**
 * High-level cosmological Poincaré Section and Orbit Recurrence Analyzer Engine.
 */
export class PoincareSectionTracer {
  /**
   * @param {VelocityField|Function} field Cosmic velocity field or vector sampler (x, y, z) => [vx, vy, vz].
   * @param {Object} [options] Configuration options.
   * @param {PoincareIntegratorType} [options.integrator=PoincareIntegratorType.DOPRI5]
   * @param {IntersectionMethod} [options.intersectionMethod=IntersectionMethod.HERMITE_CUBIC]
   * @param {CrossingDirection} [options.crossingDirection=CrossingDirection.POSITIVE]
   * @param {number} [options.maxPunctures=100] Maximum cutting plane punctures to collect.
   * @param {number} [options.maxSteps=10000] Maximum integration steps.
   * @param {number} [options.maxArcLength=5000.0] Maximum trajectory arc length (Mpc/h).
   * @param {number} [options.maxTime=100.0] Maximum cosmic time (Gyr).
   * @param {number} [options.dtInit=0.05] Initial integration time step (Gyr).
   * @param {number} [options.dtMin=1e-5] Minimum adaptive time step (Gyr).
   * @param {number} [options.dtMax=0.5] Maximum adaptive time step (Gyr).
   * @param {number} [options.relTol=1e-5] Relative tolerance for adaptive integrators.
   * @param {number} [options.absTol=1e-7] Absolute tolerance for adaptive integrators.
   * @param {boolean} [options.computeLyapunov=true] Whether to track Lyapunov exponents.
   */
  constructor(field, options = {}) {
    if (!field) {
      throw new Error('PoincareSectionTracer: A valid VelocityField or velocity sampling function must be provided.');
    }

    this.field = field;
    this.sampler = typeof field === 'function' ? field : (x, y, z) => field.sampleVelocity(x, y, z);

    this.integrator = options.integrator || PoincareIntegratorType.DOPRI5;
    this.intersectionMethod = options.intersectionMethod || IntersectionMethod.HERMITE_CUBIC;
    this.crossingDirection = options.crossingDirection || CrossingDirection.POSITIVE;

    this.maxPunctures = options.maxPunctures || 100;
    this.maxSteps = options.maxSteps || 10000;
    this.maxArcLength = options.maxArcLength || 5000.0;
    this.maxTime = options.maxTime || 100.0;

    this.dtInit = options.dtInit || 0.05;
    this.dtMin = options.dtMin || 1e-5;
    this.dtMax = options.dtMax || 0.5;
    this.relTol = options.relTol || 1e-5;
    this.absTol = options.absTol || 1e-7;

    this.computeLyapunov = options.computeLyapunov !== false;
  }

  /**
   * Traces a single streamline orbit and collects all Poincaré section punctures
   * through the specified cutting plane.
   *
   * @param {Array<number>} seedPoint [x, y, z] in Mpc/h.
   * @param {PoincarePlane} plane Poincaré cutting plane.
   * @param {Object} [traceOptions]
   * @returns {Object} Complete analysis package { punctures, recurrenceMap, classification, lyapunov, rqa, stats }
   */
  traceSection(seedPoint, plane, traceOptions = {}) {
    if (!(plane instanceof PoincarePlane)) {
      throw new TypeError('PoincareSectionTracer.traceSection: plane must be an instance of PoincarePlane.');
    }

    const maxPunctures = traceOptions.maxPunctures || this.maxPunctures;
    const maxSteps = traceOptions.maxSteps || this.maxSteps;
    const maxArcLength = traceOptions.maxArcLength || this.maxArcLength;
    const maxTime = traceOptions.maxTime || this.maxTime;
    const dir = traceOptions.crossingDirection || this.crossingDirection;
    const method = traceOptions.intersectionMethod || this.intersectionMethod;

    const punctures = [];
    const trajectory = [];

    let currentPos = [Number(seedPoint[0]), Number(seedPoint[1]), Number(seedPoint[2])];
    let currentVel = this.sampler(currentPos[0], currentPos[1], currentPos[2]);
    let currentTime = 0.0;
    let currentArcLength = 0.0;
    let dt = this.dtInit;

    trajectory.push({
      position: [...currentPos],
      velocity: [...currentVel],
      time: currentTime,
      arcLength: currentArcLength
    });

    let prevPos = [...currentPos];
    let prevVel = [...currentVel];
    let prevTime = currentTime;
    let prevArcLength = currentArcLength;

    let stepCount = 0;
    let terminationReason = 'MAX_STEPS';

    let xShadow = null;
    const d0 = 1e-5;
    let logExpSum = 0.0;
    if (this.computeLyapunov) {
      xShadow = [currentPos[0] + d0, currentPos[1], currentPos[2]];
    }

    while (stepCount < maxSteps && punctures.length < maxPunctures) {
      const stepRes = this._advanceStep(currentPos, currentVel, dt);
      const nextPos = stepRes.nextPos;
      const nextVel = stepRes.nextVel;
      const actualDt = stepRes.actualDt;
      dt = stepRes.nextDt;

      const ds = Math.hypot(nextPos[0] - currentPos[0], nextPos[1] - currentPos[1], nextPos[2] - currentPos[2]);
      const nextTime = currentTime + actualDt;
      const nextArcLength = currentArcLength + ds;

      if (plane.isCrossing(currentPos, nextPos, dir)) {
        let puncture = null;
        if (method === IntersectionMethod.HERMITE_CUBIC || method === IntersectionMethod.BRENT_EXACT) {
          puncture = ExactIntersectionSolver.solveHermiteCubic(
            currentPos, nextPos,
            currentTime, nextTime,
            currentVel, nextVel,
            currentArcLength, nextArcLength,
            plane,
            {
              field: this.field,
              stepIndex: stepCount,
              punctureIndex: punctures.length
            }
          );
        } else {
          puncture = ExactIntersectionSolver.solveLinear(
            currentPos, nextPos,
            currentTime, nextTime,
            currentVel, nextVel,
            currentArcLength, nextArcLength,
            plane,
            stepCount,
            punctures.length
          );
        }

        let acceptsCrossing = true;
        if (dir === CrossingDirection.POSITIVE && puncture.normalVelocity <= 0) {
          acceptsCrossing = false;
        } else if (dir === CrossingDirection.NEGATIVE && puncture.normalVelocity >= 0) {
          acceptsCrossing = false;
        }

        if (acceptsCrossing) {
          punctures.push(puncture);
        }
      }

      if (this.computeLyapunov && xShadow) {
        const shadowStep = this._advanceStep(xShadow, this.sampler(xShadow[0], xShadow[1], xShadow[2]), actualDt);
        xShadow = shadowStep.nextPos;
        const dCurrent = Math.hypot(xShadow[0] - nextPos[0], xShadow[1] - nextPos[1], xShadow[2] - nextPos[2]);
        if (dCurrent > 1e-18) {
          logExpSum += Math.log(dCurrent / d0);
          const invD = 1.0 / dCurrent;
          xShadow[0] = nextPos[0] + d0 * (xShadow[0] - nextPos[0]) * invD;
          xShadow[1] = nextPos[1] + d0 * (xShadow[1] - nextPos[1]) * invD;
          xShadow[2] = nextPos[2] + d0 * (xShadow[2] - nextPos[2]) * invD;
        }
      }

      prevPos = [...currentPos];
      prevVel = [...currentVel];
      prevTime = currentTime;
      prevArcLength = currentArcLength;

      currentPos = nextPos;
      currentVel = nextVel;
      currentTime = nextTime;
      currentArcLength = nextArcLength;

      stepCount++;

      if (currentArcLength >= maxArcLength) {
        terminationReason = 'MAX_ARCLENGTH';
        break;
      }
      if (currentTime >= maxTime) {
        terminationReason = 'MAX_TIME';
        break;
      }
      if (vec3Norm(currentVel) < 1e-6) {
        terminationReason = 'LOW_SPEED_STALL';
        break;
      }
    }

    if (punctures.length >= maxPunctures) {
      terminationReason = 'MAX_PUNCTURES_REACHED';
    }

    const recurrenceMap = new PoincareRecurrenceMap(punctures, plane);
    const rqa = recurrenceMap.computeRQA();

    const lambdaMax = currentTime > 0 ? (logExpSum / currentTime) : 0.0;
    const lyapunovTime = lambdaMax > 1e-5 ? (1.0 / lambdaMax) : Infinity;

    const classification = InvariantTorusClassifier.classifyOrbit(punctures, {
      lambdaMax,
      chaoticThreshold: 0.08
    });

    const fixedPoints = InvariantTorusClassifier.findFixedPoints(recurrenceMap);

    return {
      punctures,
      numPunctures: punctures.length,
      recurrenceMap,
      classification,
      fixedPoints,
      rqa,
      lyapunov: {
        lambdaMax,
        lyapunovTime,
        isChaotic: lambdaMax > 0.05,
        totalTime: currentTime
      },
      stats: {
        totalSteps: stepCount,
        totalArcLength: currentArcLength,
        totalTime: currentTime,
        terminationReason
      }
    };
  }

  /**
   * Advances a single step using the configured integrator.
   * @private
   * @param {Array<number>} pos
   * @param {Array<number>} vel
   * @param {number} dt
   * @returns {{ nextPos: Array<number>, nextVel: Array<number>, actualDt: number, nextDt: number }}
   */
  _advanceStep(pos, vel, dt) {
    if (this.integrator === PoincareIntegratorType.RK4) {
      const nextPos = rk4Step(pos, dt, this.sampler);
      const nextVel = this.sampler(nextPos[0], nextPos[1], nextPos[2]);
      return { nextPos, nextVel, actualDt: dt, nextDt: dt };
    } else if (this.integrator === PoincareIntegratorType.RK45_CASH_KARP) {
      return cashKarpStep(pos, dt, this.sampler, this.relTol, this.absTol, this.dtMin, this.dtMax);
    } else {
      return dopri5Step(pos, dt, this.sampler, this.relTol, this.absTol, this.dtMin, this.dtMax);
    }
  }

  /**
   * Batch process an ensemble of seed points.
   *
   * @param {Array<Array<number>>} seedPoints Array of [x, y, z] seeds.
   * @param {PoincarePlane} plane
   * @param {Object} [options]
   * @returns {Array<Object>}
   */
  traceBatch(seedPoints, plane, options = {}) {
    if (!Array.isArray(seedPoints)) {
      throw new TypeError('PoincareSectionTracer.traceBatch: seedPoints must be an array of seeds.');
    }
    return seedPoints.map(seed => this.traceSection(seed, plane, options));
  }
}

// ============================================================================
// 8. INTEGRATOR STEP IMPLEMENTATIONS & AUXILIARY HELPERS
// ============================================================================

/**
 * Single step of 4th-order classical Runge-Kutta.
 * @param {Array<number>} x
 * @param {number} dt
 * @param {Function} f Sampler (x, y, z) => [vx, vy, vz]
 * @returns {[number, number, number]}
 */
export function rk4Step(x, dt, f) {
  const k1 = f(x[0], x[1], x[2]);
  const x2 = [x[0] + 0.5 * dt * k1[0], x[1] + 0.5 * dt * k1[1], x[2] + 0.5 * dt * k1[2]];
  const k2 = f(x2[0], x2[1], x2[2]);
  const x3 = [x[0] + 0.5 * dt * k2[0], x[1] + 0.5 * dt * k2[1], x[2] + 0.5 * dt * k2[2]];
  const k3 = f(x3[0], x3[1], x3[2]);
  const x4 = [x[0] + dt * k3[0], x[1] + dt * k3[1], x[2] + dt * k3[2]];
  const k4 = f(x4[0], x4[1], x4[2]);

  return [
    x[0] + (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
    x[1] + (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
    x[2] + (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])
  ];
}

/**
 * Adaptive step for Cash-Karp RK45.
 */
function cashKarpStep(pos, dt, f, relTol, absTol, dtMin, dtMax) {
  let h = dt;
  const maxAttempts = 10;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    const k0 = f(pos[0], pos[1], pos[2]);
    const k1 = f(
      pos[0] + h * CK_A[1][0] * k0[0],
      pos[1] + h * CK_A[1][0] * k0[1],
      pos[2] + h * CK_A[1][0] * k0[2]
    );
    const k2 = f(
      pos[0] + h * (CK_A[2][0] * k0[0] + CK_A[2][1] * k1[0]),
      pos[1] + h * (CK_A[2][0] * k0[1] + CK_A[2][1] * k1[1]),
      pos[2] + h * (CK_A[2][0] * k0[2] + CK_A[2][1] * k1[2])
    );
    const k3 = f(
      pos[0] + h * (CK_A[3][0] * k0[0] + CK_A[3][1] * k1[0] + CK_A[3][2] * k2[0]),
      pos[1] + h * (CK_A[3][0] * k0[1] + CK_A[3][1] * k1[1] + CK_A[3][2] * k2[1]),
      pos[2] + h * (CK_A[3][0] * k0[2] + CK_A[3][1] * k1[2] + CK_A[3][2] * k2[2])
    );
    const k4 = f(
      pos[0] + h * (CK_A[4][0] * k0[0] + CK_A[4][1] * k1[0] + CK_A[4][2] * k2[0] + CK_A[4][3] * k3[0]),
      pos[1] + h * (CK_A[4][0] * k0[1] + CK_A[4][1] * k1[1] + CK_A[4][3] * k3[1]),
      pos[2] + h * (CK_A[4][0] * k0[2] + CK_A[4][1] * k1[2] + CK_A[4][2] * k2[2] + CK_A[4][3] * k3[2])
    );
    const k5 = f(
      pos[0] + h * (CK_A[5][0] * k0[0] + CK_A[5][1] * k1[0] + CK_A[5][2] * k2[0] + CK_A[5][3] * k3[0] + CK_A[5][4] * k4[0]),
      pos[1] + h * (CK_A[5][0] * k0[1] + CK_A[5][1] * k1[1] + CK_A[5][2] * k2[1] + CK_A[5][3] * k3[1] + CK_A[5][4] * k4[1]),
      pos[2] + h * (CK_A[5][0] * k0[2] + CK_A[5][1] * k1[2] + CK_A[5][2] * k2[2] + CK_A[5][3] * k3[2] + CK_A[5][4] * k4[2])
    );

    const x5 = [
      pos[0] + h * (CK_B5[0] * k0[0] + CK_B5[2] * k2[0] + CK_B5[3] * k3[0] + CK_B5[5] * k5[0]),
      pos[1] + h * (CK_B5[0] * k0[1] + CK_B5[2] * k2[1] + CK_B5[3] * k3[1] + CK_B5[5] * k5[1]),
      pos[2] + h * (CK_B5[0] * k0[2] + CK_B5[2] * k2[2] + CK_B5[3] * k3[2] + CK_B5[5] * k5[2])
    ];

    const x4 = [
      pos[0] + h * (CK_B4[0] * k0[0] + CK_B4[2] * k2[0] + CK_B4[3] * k3[0] + CK_B4[4] * k4[0] + CK_B4[5] * k5[0]),
      pos[1] + h * (CK_B4[0] * k0[1] + CK_B4[2] * k2[1] + CK_B4[3] * k3[1] + CK_B4[4] * k4[1] + CK_B4[5] * k5[1]),
      pos[2] + h * (CK_B4[0] * k0[2] + CK_B4[2] * k2[2] + CK_B4[3] * k3[2] + CK_B4[4] * k4[2] + CK_B4[5] * k5[2])
    ];

    let errNormSq = 0;
    for (let i = 0; i < 3; i++) {
      const sc = absTol + relTol * Math.max(Math.abs(pos[i]), Math.abs(x5[i]));
      const diff = (x5[i] - x4[i]) / sc;
      errNormSq += diff * diff;
    }
    const err = Math.sqrt(errNormSq / 3.0);

    const factor = Math.max(0.2, Math.min(5.0, 0.9 * Math.pow(Math.max(err, 1e-10), -0.2)));
    const nextH = Math.max(dtMin, Math.min(dtMax, h * factor));

    if (err <= 1.0 || h <= dtMin) {
      const nextVel = f(x5[0], x5[1], x5[2]);
      return { nextPos: x5, nextVel, actualDt: h, nextDt: nextH };
    }
    h = nextH;
  }

  const nextPos = rk4Step(pos, dtMin, f);
  const nextVel = f(nextPos[0], nextPos[1], nextPos[2]);
  return { nextPos, nextVel, actualDt: dtMin, nextDt: dtMin };
}

/**
 * Adaptive step for Dormand-Prince 5(4) (DOPRI5).
 */
function dopri5Step(pos, dt, f, relTol, absTol, dtMin, dtMax) {
  let h = dt;
  const maxAttempts = 10;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    const k1 = f(pos[0], pos[1], pos[2]);
    const k2 = f(
      pos[0] + h * DP_A21 * k1[0],
      pos[1] + h * DP_A21 * k1[1],
      pos[2] + h * DP_A21 * k1[2]
    );
    const k3 = f(
      pos[0] + h * (DP_A31 * k1[0] + DP_A32 * k2[0]),
      pos[1] + h * (DP_A31 * k1[1] + DP_A32 * k2[1]),
      pos[2] + h * (DP_A31 * k1[2] + DP_A32 * k2[2])
    );
    const k4 = f(
      pos[0] + h * (DP_A41 * k1[0] + DP_A42 * k2[0] + DP_A43 * k3[0]),
      pos[1] + h * (DP_A41 * k1[1] + DP_A42 * k2[1] + DP_A43 * k3[1]),
      pos[2] + h * (DP_A41 * k1[2] + DP_A42 * k2[2] + DP_A43 * k3[2])
    );
    const k5 = f(
      pos[0] + h * (DP_A51 * k1[0] + DP_A52 * k2[0] + DP_A53 * k3[0] + DP_A54 * k4[0]),
      pos[1] + h * (DP_A51 * k1[1] + DP_A52 * k2[1] + DP_A53 * k3[1] + DP_A54 * k4[1]),
      pos[2] + h * (DP_A51 * k1[2] + DP_A52 * k2[2] + DP_A53 * k3[2] + DP_A54 * k4[2])
    );
    const k6 = f(
      pos[0] + h * (DP_A61 * k1[0] + DP_A62 * k2[0] + DP_A63 * k3[0] + DP_A64 * k4[0] + DP_A65 * k5[0]),
      pos[1] + h * (DP_A61 * k1[1] + DP_A62 * k2[1] + DP_A63 * k3[1] + DP_A64 * k4[1] + DP_A65 * k5[1]),
      pos[2] + h * (DP_A61 * k1[2] + DP_A62 * k2[2] + DP_A63 * k3[2] + DP_A64 * k4[2] + DP_A65 * k5[2])
    );

    const xNext = [
      pos[0] + h * (DP_B[0] * k1[0] + DP_B[2] * k3[0] + DP_B[3] * k4[0] + DP_B[4] * k5[0] + DP_B[5] * k6[0]),
      pos[1] + h * (DP_B[0] * k1[1] + DP_B[2] * k3[1] + DP_B[3] * k4[1] + DP_B[4] * k5[1] + DP_B[5] * k6[1]),
      pos[2] + h * (DP_B[0] * k1[2] + DP_B[2] * k3[2] + DP_B[3] * k4[2] + DP_B[4] * k5[2] + DP_B[5] * k6[2])
    ];

    const k7 = f(xNext[0], xNext[1], xNext[2]);

    let errNormSq = 0;
    for (let i = 0; i < 3; i++) {
      const errEst = h * (DP_E[0] * k1[i] + DP_E[2] * k3[i] + DP_E[3] * k4[i] + DP_E[4] * k5[i] + DP_E[5] * k6[i] + DP_E[6] * k7[i]);
      const sc = absTol + relTol * Math.max(Math.abs(pos[i]), Math.abs(xNext[i]));
      const diff = errEst / sc;
      errNormSq += diff * diff;
    }
    const err = Math.sqrt(errNormSq / 3.0);

    const factor = Math.max(0.2, Math.min(5.0, 0.9 * Math.pow(Math.max(err, 1e-10), -0.2)));
    const nextH = Math.max(dtMin, Math.min(dtMax, h * factor));

    if (err <= 1.0 || h <= dtMin) {
      return { nextPos: xNext, nextVel: k7, actualDt: h, nextDt: nextH };
    }
    h = nextH;
  }

  const nextPos = rk4Step(pos, dtMin, f);
  const nextVel = f(nextPos[0], nextPos[1], nextPos[2]);
  return { nextPos, nextVel, actualDt: dtMin, nextDt: dtMin };
}

/**
 * Computes 3x3 velocity Jacobian J_{ij} = dv_i / dx_j via central finite differences.
 * @param {Array<number>} x
 * @param {Function} f
 * @param {number} [h=1e-4]
 * @returns {Float64Array} 3x3 row-major
 */
export function computeJacobian3x3(x, f, h = 1e-4) {
  const J = new Float64Array(9);
  const inv2h = 0.5 / h;

  // d/dx
  const fxP = f(x[0] + h, x[1], x[2]);
  const fxM = f(x[0] - h, x[1], x[2]);
  J[0] = (fxP[0] - fxM[0]) * inv2h;
  J[3] = (fxP[1] - fxM[1]) * inv2h;
  J[6] = (fxP[2] - fxM[2]) * inv2h;

  // d/dy
  const fyP = f(x[0], x[1] + h, x[2]);
  const fyM = f(x[0], x[1] - h, x[2]);
  J[1] = (fyP[0] - fyM[0]) * inv2h;
  J[4] = (fyP[1] - fyM[1]) * inv2h;
  J[7] = (fyP[2] - fyM[2]) * inv2h;

  // d/dz
  const fzP = f(x[0], x[1], x[2] + h);
  const fzM = f(x[0], x[1], x[2] - h);
  J[2] = (fzP[0] - fzM[0]) * inv2h;
  J[5] = (fzP[1] - fzM[1]) * inv2h;
  J[8] = (fzP[2] - fzM[2]) * inv2h;

  return J;
}

/**
 * Matrix addition A + s * B.
 * @param {Float64Array} A
 * @param {Float64Array} B
 * @param {number} s
 * @returns {Float64Array}
 */
function addScaledMat3(A, B, s) {
  const out = new Float64Array(9);
  for (let i = 0; i < 9; i++) {
    out[i] = A[i] + s * B[i];
  }
  return out;
}

/**
 * Helper: simple spatial clustering of 2D points.
 */
function clusterPunctures(points, clusterRadius) {
  const clusters = [];
  const assigned = new Uint8Array(points.length);

  for (let i = 0; i < points.length; i++) {
    if (assigned[i]) continue;
    const clPoints = [points[i]];
    assigned[i] = 1;

    for (let j = i + 1; j < points.length; j++) {
      if (assigned[j]) continue;
      const d = Math.hypot(points[i][0] - points[j][0], points[i][1] - points[j][1]);
      if (d <= clusterRadius) {
        clPoints.push(points[j]);
        assigned[j] = 1;
      }
    }

    let su = 0, sv = 0;
    for (let k = 0; k < clPoints.length; k++) {
      su += clPoints[k][0];
      sv += clPoints[k][1];
    }
    const cu = su / clPoints.length;
    const cv = sv / clPoints.length;

    let maxR = 0;
    for (let k = 0; k < clPoints.length; k++) {
      const r = Math.hypot(clPoints[k][0] - cu, clPoints[k][1] - cv);
      if (r > maxR) maxR = r;
    }

    clusters.push({ center: [cu, cv], count: clPoints.length, radius: maxR });
  }

  return clusters;
}

/**
 * Helper: Detects multi-island chains with distinct angular gaps.
 */
function detectIslandChainsWithGaps(points, uc, vc) {
  if (points.length < 20) return [];

  const numBins = 24;
  const binCounts = new Int32Array(numBins);

  for (let i = 0; i < points.length; i++) {
    let angle = Math.atan2(points[i][1] - vc, points[i][0] - uc);
    if (angle < 0) angle += 2.0 * Math.PI;
    const bin = Math.min(numBins - 1, Math.floor((angle / (2.0 * Math.PI)) * numBins));
    binCounts[bin]++;
  }

  // Count empty or near-empty gap bins
  let zeroBins = 0;
  for (let i = 0; i < numBins; i++) {
    if (binCounts[i] === 0) zeroBins++;
  }

  // An island chain MUST have significant empty gaps between islands (at least 25% of circle is empty)
  if (zeroBins < 6) {
    return [];
  }

  const peaks = [];
  const avg = points.length / numBins;

  for (let i = 0; i < numBins; i++) {
    const prev = binCounts[(i - 1 + numBins) % numBins];
    const curr = binCounts[i];
    const next = binCounts[(i + 1) % numBins];

    if (curr > prev && curr >= next && curr > 1.5 * avg) {
      peaks.push(i);
    }
  }

  return peaks;
}

/**
 * Computes angular coverage fraction in [0, 1].
 */
function computeAngularCoverage(points, uc, vc) {
  const numBins = 16;
  const bins = new Uint8Array(numBins);
  for (let i = 0; i < points.length; i++) {
    let a = Math.atan2(points[i][1] - vc, points[i][0] - uc);
    if (a < 0) a += 2.0 * Math.PI;
    const b = Math.min(numBins - 1, Math.floor((a / (2.0 * Math.PI)) * numBins));
    bins[b] = 1;
  }
  let filled = 0;
  for (let i = 0; i < numBins; i++) {
    if (bins[i] === 1) filled++;
  }
  return filled / numBins;
}

/**
 * Computes winding number directly from 2D points.
 */
function computeWindingFromPoints(points, uc, vc) {
  if (points.length < 3) return 0.0;
  let totalAngle = 0.0;
  let prevAngle = Math.atan2(points[0][1] - vc, points[0][0] - uc);

  for (let i = 1; i < points.length; i++) {
    const curAngle = Math.atan2(points[i][1] - vc, points[i][0] - uc);
    let dTheta = curAngle - prevAngle;
    while (dTheta > Math.PI) dTheta -= 2.0 * Math.PI;
    while (dTheta < -Math.PI) dTheta += 2.0 * Math.PI;
    totalAngle += dTheta;
    prevAngle = curAngle;
  }

  return totalAngle / (2.0 * Math.PI * (points.length - 1));
}
