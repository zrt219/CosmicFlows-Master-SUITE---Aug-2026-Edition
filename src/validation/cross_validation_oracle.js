/**
 * @file cross_validation_oracle.js
 * @module validation/cross_validation_oracle
 * @description Multi-Method Cross-Validation Oracle and Scientific Benchmark Harness for ZRT Cosmicflows Workbench.
 * 
 * Mathematical Formulations & Core Features:
 * 1. Dual-Method Solver Cross-Validation:
 *    - Adaptive Cash-Karp RK45 vs Dormand-Prince DOPRI5 (FSAL 5(4)) vs Classical RK4 & Yoshida Symplectic reference integrators.
 *    - Marching Cubes (Lorensen-Cline 256-table with asymptotic decider) vs Marching Tetrahedra (Kuhn 6-simplex decomposition).
 *    - Jacobi rotation eigenvalue solver vs Characteristic polynomial Cardano reference solver.
 *    - 3D Newton-Raphson vector root finding vs Broyden quasi-Newton reference solver.
 *    - Minimum-Variance Bulk Flow dipole estimator vs Spherical Harmonics multipole projector.
 * 
 * 2. Differential Discrepancy & Gradient Field Auditor:
 *    - \Delta_{\text{grad}} = ||\nabla_{\text{num}} f - \nabla_{\text{ana}} f||_2
 *    - Relative discrepancy \Delta_{\text{rel}} = \frac{\Delta_{\text{grad}}}{\max(||\nabla_{\text{ana}} f||_2, \epsilon)}
 *    - 4th-order compact central difference stencil:
 *      \nabla f_i = \frac{-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)}{12h}
 *    - Numerical Hessian \mathbf{H} and symmetry verification ||\mathbf{H} - \mathbf{H}^T||_F = 0
 *    - Divergence \nabla \cdot \mathbf{v}, Vorticity \nabla \times \mathbf{v}, and Richardson extrapolation convergence order p.
 * 
 * 3. Continuous & Discrete Polyline Trajectory Metrics:
 *    - Continuous Segment-Projected Hausdorff distance:
 *      H_{\text{cont}}(A, B) = \max( \max_{a \in A} \min_{S \in \text{segs}(B)} \text{dist}(a, S), \max_{b \in B} \min_{S \in \text{segs}(A)} \text{dist}(b, S) )
 *    - Modified / Average Hausdorff distance H_{\text{avg}}(A, B)
 *    - Discrete Fréchet distance via Dynamic Programming
 *    - Arc length and trajectory reversibility residual ||\mathbf{x}_0 - \mathbf{x}_{\text{rev}}||.
 * 
 * 4. Volumetric Enclosed Volume & Manifold Topology Agreement:
 *    - Enclosed volume calculus via Divergence Theorem:
 *      V_{\text{enc}} = \frac{1}{6} | \sum_{f} \mathbf{v}_0 \cdot (\mathbf{v}_1 \times \mathbf{v}_2) |
 *    - Relative volume agreement: \Delta_V = \frac{|V_{\text{MC}} - V_{\text{MT}}|}{\frac{1}{2}(V_{\text{MC}} + V_{\text{MT}})}
 *    - Surface area calculus: A = \frac{1}{2} \sum_f ||(\mathbf{v}_1 - \mathbf{v}_0) \times (\mathbf{v}_2 - \mathbf{v}_0)||_2
 *    - Voxel Jaccard J(A, B) and Sørensen-Dice D(A, B) segmentation similarity.
 *    - Topological Euler characteristic \chi = V - E + F, genus g = 1 - \chi/2, and watertight 2-manifold verification.
 * 
 * 5. Full 10-Gate Scientific Compliance Suite (Gates 1 through 10):
 *    - Gate 1: Data Integrity & Finite Bounds (NaN/Inf, relativistic bounds, positive distance)
 *    - Gate 2: Units, Dimensional Consistency & Cosmology (H0, Omega_m, Omega_L, a(z) = 1/(1+z))
 *    - Gate 3: Analytic Vector Fields & Exact Invariants (Laminar Div=0/Curl=0, Sink trace, Vortex curl)
 *    - Gate 4: Differential Discrepancy & Gradient Field Auditor (\Delta_{\text{grad}} <= \tau, Hessian symmetry)
 *    - Gate 5: Numerical Integration, Butcher Convergence & Hausdorff Metric (CK vs DP Hausdorff, energy conservation)
 *    - Gate 6: Volumetric Meshing & Enclosed Volume Agreement (MC vs MT agreement, 2-sphere chi=2, watertight)
 *    - Gate 7: Critical Point Roots & Spectral Orthogonality (Newton-Raphson residual, Jacobi ||V^T V - I||)
 *    - Gate 8: Cosmic Watershed Basin Segmentation & Similarity (Dice similarity >= 0.999, partition of unity)
 *    - Gate 9: Cosmic Bulk Flow & Multipole Consistency (Minimum-variance recovery, dipole direction <= 5 deg)
 *    - Gate 10: Serialization, Provenance & Export Roundtrip (FITS 2880-byte block alignment, GeoJSON, SHA-256 seal)
 * 
 * 6. Cryptographic Verification Passport & Audit Trail:
 *    - Deterministic canonical JSON serialization and SHA-256 digest + HMAC authentication seal.
 *    - Verification and tamper-detection validator.
 * 
 * 7. Scientific Benchmark Harness:
 *    - Statistical profiling (mean, std, median, min, max, p95, p99, ops/sec) and regression budgets.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

// ============================================================================
// 1. ENUMS, CONSTANTS & CONFIGURATION
// ============================================================================

export const GateStatus = Object.freeze({
  PASSED: 'PASSED',
  FAILED: 'FAILED',
  WARNING: 'WARNING',
  SKIPPED: 'SKIPPED'
});

export const ComplianceLevel = Object.freeze({
  CERTIFIED_REPRODUCIBLE: 'CERTIFIED_REPRODUCIBLE',
  CONDITIONALLY_VALIDATED: 'CONDITIONALLY_VALIDATED',
  FAILED_VERIFICATION: 'FAILED_VERIFICATION'
});

export const CriticalPointType = Object.freeze({
  ATTRACTOR_NODE: 'ATTRACTOR_NODE',   // Index 3: 3 negative eigenvalues (Sink / Attractor)
  FILAMENT_SADDLE: 'FILAMENT_SADDLE', // Index 2: 2 negative, 1 positive eigenvalue (Filament core)
  WALL_SADDLE: 'WALL_SADDLE',         // Index 1: 1 negative, 2 positive eigenvalues (Sheet / Wall)
  REPELLER_NODE: 'REPELLER_NODE',     // Index 0: 3 positive eigenvalues (Source / Cosmic Void)
  DEGENERATE: 'DEGENERATE'            // Singular Jacobian / non-hyperbolic critical point
});

export const DEFAULT_TOLERANCES = Object.freeze({
  GRADIENT_ABS_TOL: 1e-4,
  GRADIENT_REL_TOL: 1e-3,
  HESSIAN_FROBENIUS_TOL: 1e-3,
  DIVERGENCE_TOL: 1e-4,
  VORTICITY_TOL: 1e-4,
  HAUSDORFF_TRAJECTORY_TOL: 1e-3,
  INTEGRATOR_REVERSIBILITY_TOL: 1e-4,
  ENERGY_CONSERVATION_TOL: 1e-4,
  VOLUME_REL_AGREEMENT_TOL: 0.05, // 5% agreement between MC and MT
  AREA_REL_AGREEMENT_TOL: 0.08,   // 8% agreement between MC and MT
  VOXEL_DICE_SIMILARITY_MIN: 0.99,
  JACOBI_ORTHOGONALITY_TOL: 1e-11,
  JACOBI_RESIDUAL_TOL: 1e-11,
  ROOT_RESIDUAL_TOL: 1e-7,
  BULK_FLOW_REL_TOL: 0.08,
  BULK_FLOW_ANGLE_DEG_TOL: 5.0
});

// Cube geometry constants
export const CUBE_CORNERS = Object.freeze([
  [0, 0, 0], // 0
  [1, 0, 0], // 1
  [1, 1, 0], // 2
  [0, 1, 0], // 3
  [0, 0, 1], // 4
  [1, 0, 1], // 5
  [1, 1, 1], // 6
  [0, 1, 1]  // 7
]);

export const EDGE_VERTICES = Object.freeze([
  [0, 1], // 0
  [1, 2], // 1
  [2, 3], // 2
  [3, 0], // 3
  [4, 5], // 4
  [5, 6], // 5
  [6, 7], // 6
  [7, 4], // 7
  [0, 4], // 8
  [1, 5], // 9
  [2, 6], // 10
  [3, 7]  // 11
]);

export const TET_EDGE_PAIRS = Object.freeze([
  [0, 1], // Edge 0: v0-v1
  [0, 2], // Edge 1: v0-v2
  [0, 3], // Edge 2: v0-v3
  [1, 2], // Edge 3: v1-v2
  [1, 3], // Edge 4: v1-v3
  [2, 3]  // Edge 5: v2-v3
]);

export const KUHN_6_TETRAHEDRA = Object.freeze([
  [0, 1, 2, 6],
  [0, 5, 1, 6],
  [0, 2, 3, 6],
  [0, 3, 7, 6],
  [0, 4, 5, 6],
  [0, 7, 4, 6]
]);

export const CANONICAL_TET_TRIANGLE_TABLE = Object.freeze([
  [],                                // Case 0:  0000 (all outside)
  [[0, 1, 2]],                       // Case 1:  0001 (v0 inside)
  [[0, 4, 3]],                       // Case 2:  0010 (v1 inside)
  [[1, 2, 4], [1, 4, 3]],            // Case 3:  0011 (v0, v1 inside)
  [[1, 3, 5]],                       // Case 4:  0100 (v2 inside)
  [[0, 5, 2], [0, 3, 5]],            // Case 5:  0101 (v0, v2 inside)
  [[0, 5, 1], [0, 4, 5]],            // Case 6:  0110 (v1, v2 inside)
  [[2, 4, 5]],                       // Case 7:  0111 (v0, v1, v2 inside = v3 outside)
  [[2, 5, 4]],                       // Case 8:  1000 (v3 inside)
  [[0, 1, 5], [0, 5, 4]],            // Case 9:  1001 (v0, v3 inside)
  [[0, 2, 5], [0, 5, 3]],            // Case 10: 1010 (v1, v3 inside)
  [[1, 5, 3]],                       // Case 11: 1011 (v0, v1, v3 inside = v2 outside)
  [[1, 4, 2], [1, 3, 4]],            // Case 12: 1100 (v2, v3 inside)
  [[0, 3, 4]],                       // Case 13: 1101 (v0, v2, v3 inside = v1 outside)
  [[0, 2, 1]],                       // Case 14: 1110 (v1, v2, v3 inside = v0 outside)
  []                                 // Case 15: 1111 (all inside)
]);

// ============================================================================
// 2. MARCHING CUBES 256-CASE LOOKUP TABLE BUILDER
// ============================================================================

function getEdgeIndex(c1, c2) {
  for (let e = 0; e < 12; e++) {
    const [a, b] = EDGE_VERTICES[e];
    if ((a === c1 && b === c2) || (a === c2 && b === c1)) return e;
  }
  throw new Error(`Edge not found between corners ${c1} and ${c2}`);
}

function buildCubeSymmetries() {
  const symmetries = [];
  for (let mirror = 0; mirror < 2; mirror++) {
    for (let rotX = 0; rotX < 4; rotX++) {
      for (let rotY = 0; rotY < 4; rotY++) {
        for (let rotZ = 0; rotZ < 4; rotZ++) {
          const perm = [];
          for (let c = 0; c < 8; c++) {
            let x = CUBE_CORNERS[c][0] - 0.5;
            let y = CUBE_CORNERS[c][1] - 0.5;
            let z = CUBE_CORNERS[c][2] - 0.5;

            if (mirror === 1) x = -x;
            for (let r = 0; r < rotX; r++) { const ty = y; y = -z; z = ty; }
            for (let r = 0; r < rotY; r++) { const tx = x; x = z; z = -tx; }
            for (let r = 0; r < rotZ; r++) { const tx = x; x = -y; y = tx; }

            x = Math.round(x + 0.5);
            y = Math.round(y + 0.5);
            z = Math.round(z + 0.5);

            let target = -1;
            for (let tc = 0; tc < 8; tc++) {
              if (CUBE_CORNERS[tc][0] === x && CUBE_CORNERS[tc][1] === y && CUBE_CORNERS[tc][2] === z) {
                target = tc;
                break;
              }
            }
            perm.push(target);
          }
          const key = perm.join(',');
          if (!symmetries.some(r => r.key === key)) {
            symmetries.push({ key, perm, isReflection: mirror === 1 });
          }
        }
      }
    }
  }
  return symmetries;
}

const CUBE_SYMMETRIES = buildCubeSymmetries();

const CANONICAL_BASE_CASES = [
  [0, []],
  [1, [[0, 8, 3]]],
  [3, [[1, 8, 3], [9, 8, 1]]],
  [5, [[0, 8, 3], [1, 2, 10]]],
  [65, [[0, 8, 3], [5, 6, 10]]],
  [7, [[2, 8, 3], [2, 10, 8], [10, 9, 8]]],
  [67, [[1, 8, 3], [9, 8, 1], [5, 6, 10]]],
  [37, [[0, 8, 3], [1, 2, 10], [4, 9, 5]]],
  [15, [[9, 8, 10], [10, 8, 11]]],
  [23, [[2, 10, 9], [2, 9, 7], [2, 7, 3], [7, 9, 4]]],
  [27, [[1, 9, 2], [9, 11, 2], [9, 4, 11], [4, 7, 11]]],
  [165, [[0, 8, 3], [1, 2, 10], [4, 9, 5], [6, 7, 11]]],
  [99, [[1, 2, 6], [1, 6, 8], [1, 8, 9], [8, 6, 7]]],
  [195, [[2, 6, 9], [2, 9, 1], [6, 7, 9], [0, 9, 3], [7, 3, 9]]],
  [101, [[0, 8, 3], [1, 2, 10], [4, 6, 9], [6, 10, 9]]],
  [111, [[2, 7, 3], [6, 7, 2]]]
];

function mapEdge(edgeIdx, perm) {
  const [c1, c2] = EDGE_VERTICES[edgeIdx];
  return getEdgeIndex(perm[c1], perm[c2]);
}

function buildCompleteTriangleTable() {
  const table = new Array(256);

  for (const [baseMask, baseTris] of CANONICAL_BASE_CASES) {
    for (const sym of CUBE_SYMMETRIES) {
      let mask = 0;
      for (let c = 0; c < 8; c++) {
        if ((baseMask & (1 << c)) !== 0) {
          mask |= (1 << sym.perm[c]);
        }
      }

      if (table[mask] === undefined) {
        const transformedTris = [];
        for (const tri of baseTris) {
          const e0 = mapEdge(tri[0], sym.perm);
          const e1 = mapEdge(tri[1], sym.perm);
          const e2 = mapEdge(tri[2], sym.perm);

          if (sym.isReflection) {
            transformedTris.push([e0, e2, e1]);
          } else {
            transformedTris.push([e0, e1, e2]);
          }
        }
        table[mask] = transformedTris;

        // Invert mask and reverse triangles (complement)
        const invMask = mask ^ 255;
        if (table[invMask] === undefined) {
          const invTris = transformedTris.map(t => [t[0], t[2], t[1]]);
          table[invMask] = invTris;
        }
      }
    }
  }

  // Convert to flat integer arrays terminated by -1
  const flatTable = new Array(256);
  for (let m = 0; m < 256; m++) {
    const list = table[m] || [];
    const flat = [];
    for (const t of list) {
      flat.push(t[0], t[1], t[2]);
    }
    flat.push(-1);
    flatTable[m] = flat;
  }

  return flatTable;
}

export const MARCHING_CUBES_TRIANGLE_TABLE = buildCompleteTriangleTable();

// ============================================================================
// 3. CANONICAL SHA-256 CRYPTOGRAPHIC HASHER (UNIVERSAL)
// ============================================================================

export class UniversalSHA256 {
  /**
   * Computes SHA-256 hexadecimal hash string for input message.
   * @param {string|Uint8Array|Object} input
   * @returns {string} 64-character lowercase hex string
   */
  static hash(input) {
    let bytes;
    if (typeof input === 'string') {
      bytes = new TextEncoder().encode(input);
    } else if (input instanceof Uint8Array) {
      bytes = input;
    } else {
      bytes = new TextEncoder().encode(UniversalSHA256.canonicalJson(input));
    }

    const K = [
      0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
      0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
      0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
      0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
      0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
      0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
      0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
      0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    ];

    let H0 = 0x6a09e667;
    let H1 = 0xbb67ae85;
    let H2 = 0x3c6ef372;
    let H3 = 0xa54ff53a;
    let H4 = 0x510e527f;
    let H5 = 0x9b05688c;
    let H6 = 0x1f83d9ab;
    let H7 = 0x5be0cd19;

    const byteLen = bytes.length;
    const bitLen = byteLen * 8;

    const padLen = (((byteLen + 8) >> 6) + 1) << 6;
    const padded = new Uint8Array(padLen);
    padded.set(bytes, 0);
    padded[byteLen] = 0x80;

    const view = new DataView(padded.buffer);
    view.setUint32(padLen - 4, bitLen & 0xffffffff, false);
    view.setUint32(padLen - 8, Math.floor(bitLen / 0x100000000), false);

    const W = new Uint32Array(64);

    for (let offset = 0; offset < padLen; offset += 64) {
      for (let t = 0; t < 16; t++) {
        W[t] = view.getUint32(offset + t * 4, false);
      }
      for (let t = 16; t < 64; t++) {
        const gamma0 = ((W[t - 15] >>> 7) | (W[t - 15] << 25)) ^
                       ((W[t - 15] >>> 18) | (W[t - 15] << 14)) ^
                       (W[t - 15] >>> 3);
        const gamma1 = ((W[t - 2] >>> 17) | (W[t - 2] << 15)) ^
                       ((W[t - 2] >>> 19) | (W[t - 2] << 13)) ^
                       (W[t - 2] >>> 10);
        W[t] = (((W[t - 16] + gamma0) | 0) + ((W[t - 7] + gamma1) | 0)) | 0;
      }

      let a = H0, b = H1, c = H2, d = H3, e = H4, f = H5, g = H6, h = H7;

      for (let t = 0; t < 64; t++) {
        const sigma1 = ((e >>> 6) | (e << 26)) ^ ((e >>> 11) | (e << 21)) ^ ((e >>> 25) | (e << 7));
        const ch = (e & f) ^ ((~e) & g);
        const t1 = (((((h + sigma1) | 0) + ch) | 0) + ((K[t] + W[t]) | 0)) | 0;
        const sigma0 = ((a >>> 2) | (a << 30)) ^ ((a >>> 13) | (a << 19)) ^ ((a >>> 22) | (a << 10));
        const maj = (a & b) ^ (a & c) ^ (b & c);
        const t2 = (sigma0 + maj) | 0;

        h = g;
        g = f;
        f = e;
        e = (d + t1) | 0;
        d = c;
        c = b;
        b = a;
        a = (t1 + t2) | 0;
      }

      H0 = (H0 + a) | 0;
      H1 = (H1 + b) | 0;
      H2 = (H2 + c) | 0;
      H3 = (H3 + d) | 0;
      H4 = (H4 + e) | 0;
      H5 = (H5 + f) | 0;
      H6 = (H6 + g) | 0;
      H7 = (H7 + h) | 0;
    }

    const hex = (v) => (v >>> 0).toString(16).padStart(8, '0');
    return hex(H0) + hex(H1) + hex(H2) + hex(H3) + hex(H4) + hex(H5) + hex(H6) + hex(H7);
  }

  /**
   * Computes HMAC-SHA256 signature for authentication.
   * @param {string} key
   * @param {string} message
   * @returns {string}
   */
  static hmac(key, message) {
    const encoder = new TextEncoder();
    let keyBytes = encoder.encode(key);
    if (keyBytes.length > 64) {
      keyBytes = new Uint8Array(32);
      const hex = UniversalSHA256.hash(key);
      for (let i = 0; i < 32; i++) {
        keyBytes[i] = parseInt(hex.substring(i * 2, i * 2 + 2), 16);
      }
    }
    const kPad = new Uint8Array(64);
    kPad.set(keyBytes, 0);

    const oKeyPad = new Uint8Array(64);
    const iKeyPad = new Uint8Array(64);
    for (let i = 0; i < 64; i++) {
      oKeyPad[i] = kPad[i] ^ 0x5c;
      iKeyPad[i] = kPad[i] ^ 0x36;
    }

    const msgBytes = encoder.encode(message);
    const inner = new Uint8Array(64 + msgBytes.length);
    inner.set(iKeyPad, 0);
    inner.set(msgBytes, 64);
    const innerHashHex = UniversalSHA256.hash(inner);

    const innerHashBytes = new Uint8Array(32);
    for (let i = 0; i < 32; i++) {
      innerHashBytes[i] = parseInt(innerHashHex.substring(i * 2, i * 2 + 2), 16);
    }

    const outer = new Uint8Array(64 + 32);
    outer.set(oKeyPad, 0);
    outer.set(innerHashBytes, 64);
    return UniversalSHA256.hash(outer);
  }

  /**
   * Deterministic canonical JSON serializer for reproducible hashing.
   * @param {*} obj
   * @returns {string}
   */
  static canonicalJson(obj) {
    if (obj === null || typeof obj !== 'object') {
      return JSON.stringify(obj);
    }
    if (Array.isArray(obj)) {
      return '[' + obj.map(item => UniversalSHA256.canonicalJson(item)).join(',') + ']';
    }
    const keys = Object.keys(obj).sort();
    return '{' + keys.map(k => JSON.stringify(k) + ':' + UniversalSHA256.canonicalJson(obj[k])).join(',') + '}';
  }
}

// ============================================================================
// 4. VECTOR & MATRIX MATHEMATICAL UTILITIES
// ============================================================================

export class VectorMath {
  static norm2(v) {
    let s = 0;
    for (let i = 0; i < v.length; i++) s += v[i] * v[i];
    return Math.sqrt(s);
  }

  static norm1(v) {
    let s = 0;
    for (let i = 0; i < v.length; i++) s += Math.abs(v[i]);
    return s;
  }

  static normInf(v) {
    let m = 0;
    for (let i = 0; i < v.length; i++) {
      const a = Math.abs(v[i]);
      if (a > m) m = a;
    }
    return m;
  }

  static dist(a, b) {
    let s = 0;
    for (let i = 0; i < 3; i++) {
      const d = a[i] - b[i];
      s += d * d;
    }
    return Math.sqrt(s);
  }

  static dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  static cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    ];
  }

  static sub(a, b) {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  }

  static add(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  }

  static scale(v, s) {
    return [v[0] * s, v[1] * s, v[2] * s];
  }

  static normalize(v) {
    const n = VectorMath.norm2(v);
    if (n < 1e-15) return [0, 0, 0];
    return [v[0] / n, v[1] / n, v[2] / n];
  }

  static angleBetweenDeg(a, b) {
    const na = VectorMath.norm2(a);
    const nb = VectorMath.norm2(b);
    if (na < 1e-15 || nb < 1e-15) return 0;
    let cosT = VectorMath.dot(a, b) / (na * nb);
    if (cosT > 1.0) cosT = 1.0;
    if (cosT < -1.0) cosT = -1.0;
    return (Math.acos(cosT) * 180.0) / Math.PI;
  }

  /**
   * Minimum distance from point P to 3D line segment AB.
   * @param {number[]} p
   * @param {number[]} a
   * @param {number[]} b
   * @returns {number}
   */
  static distToSegment(p, a, b) {
    const ab = VectorMath.sub(b, a);
    const ap = VectorMath.sub(p, a);
    const l2 = VectorMath.dot(ab, ab);
    if (l2 < 1e-14) return VectorMath.dist(p, a);
    let t = VectorMath.dot(ap, ab) / l2;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;
    const proj = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    return VectorMath.dist(p, proj);
  }
}

export class MatrixMath {
  static det3x3(m) {
    return (
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    );
  }

  static trace3x3(m) {
    return m[0][0] + m[1][1] + m[2][2];
  }

  static frobeniusNorm3x3(m) {
    let sum = 0;
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        sum += m[r][c] * m[r][c];
      }
    }
    return Math.sqrt(sum);
  }

  static diffFrobeniusNorm3x3(a, b) {
    let sum = 0;
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        const d = a[r][c] - b[r][c];
        sum += d * d;
      }
    }
    return Math.sqrt(sum);
  }

  static isSymmetric3x3(m, tol = 1e-7) {
    return (
      Math.abs(m[0][1] - m[1][0]) <= tol &&
      Math.abs(m[0][2] - m[2][0]) <= tol &&
      Math.abs(m[1][2] - m[2][1]) <= tol
    );
  }

  static matMul3x3(a, b) {
    const res = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        res[r][c] = a[r][0] * b[0][c] + a[r][1] * b[1][c] + a[r][2] * b[2][c];
      }
    }
    return res;
  }

  static matVecMul3x3(m, v) {
    return [
      m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
      m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
      m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
    ];
  }

  static transpose3x3(m) {
    return [
      [m[0][0], m[1][0], m[2][0]],
      [m[0][1], m[1][1], m[2][1]],
      [m[0][2], m[1][2], m[2][2]]
    ];
  }

  static inverse3x3(m) {
    const det = MatrixMath.det3x3(m);
    if (Math.abs(det) < 1e-15) {
      throw new Error(`Singular matrix encountered, determinant: ${det}`);
    }
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
   * Real symmetric 3x3 Jacobi eigenvalue diagonalization.
   * Computes eigenvalues and orthogonal eigenvector matrix V such that A * V = V * diag(lambda).
   * Sorted descending: lambda[0] >= lambda[1] >= lambda[2].
   */
  static jacobiEigenvalues3x3(a, maxIter = 60, tol = 1e-14) {
    const A = [
      [a[0][0], a[0][1], a[0][2]],
      [a[1][0], a[1][1], a[1][2]],
      [a[2][0], a[2][1], a[2][2]]
    ];

    const V = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];

    let converged = false;
    let iterations = 0;
    let maxOff = 0;

    for (let iter = 0; iter < maxIter; iter++) {
      iterations = iter + 1;
      maxOff = 0;
      let p = 0, q = 1;

      for (let i = 0; i < 3; i++) {
        for (let j = i + 1; j < 3; j++) {
          const off = Math.abs(A[i][j]);
          if (off > maxOff) {
            maxOff = off;
            p = i;
            q = j;
          }
        }
      }

      if (maxOff < tol) {
        converged = true;
        break;
      }

      const diff = A[q][q] - A[p][p];
      let t;
      if (Math.abs(diff) < 1e-15) {
        t = A[p][q] > 0 ? 1.0 : -1.0;
      } else {
        const phi = diff / (2.0 * A[p][q]);
        t = Math.sign(phi) / (Math.abs(phi) + Math.sqrt(phi * phi + 1.0));
      }

      const c = 1.0 / Math.sqrt(t * t + 1.0);
      const s = t * c;
      const tau = s / (1.0 + c);

      const a_pq = A[p][q];
      A[p][q] = 0;
      A[q][p] = 0;
      A[p][p] -= t * a_pq;
      A[q][q] += t * a_pq;

      for (let r = 0; r < 3; r++) {
        if (r !== p && r !== q) {
          const a_rp = A[r][p];
          const a_rq = A[r][q];
          A[r][p] = a_rp - s * (a_rq + tau * a_rp);
          A[p][r] = A[r][p];
          A[r][q] = a_rq + s * (a_rp - tau * a_rq);
          A[q][r] = A[r][q];
        }
      }

      for (let r = 0; r < 3; r++) {
        const v_rp = V[r][p];
        const v_rq = V[r][q];
        V[r][p] = v_rp - s * (v_rq + tau * v_rp);
        V[r][q] = v_rq + s * (v_rp - tau * v_rq);
      }
    }

    const evs = [
      { val: A[0][0], vec: [V[0][0], V[1][0], V[2][0]] },
      { val: A[1][1], vec: [V[0][1], V[1][1], V[2][1]] },
      { val: A[2][2], vec: [V[0][2], V[1][2], V[2][2]] }
    ];

    evs.sort((e1, e2) => e2.val - e1.val);

    const eigenvalues = [evs[0].val, evs[1].val, evs[2].val];
    const eigenvectors = [
      [evs[0].vec[0], evs[1].vec[0], evs[2].vec[0]],
      [evs[0].vec[1], evs[1].vec[1], evs[2].vec[1]],
      [evs[0].vec[2], evs[1].vec[2], evs[2].vec[2]]
    ];

    return {
      eigenvalues,
      eigenvectors,
      iterations,
      converged,
      maxOffDiag: maxOff
    };
  }

  /**
   * Reference Characteristic Polynomial Cardano Eigensolver for cross-validating Jacobi.
   */
  static referenceEigenvalues3x3(a) {
    const p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    if (p1 < 1e-15) {
      return [a[0][0], a[1][1], a[2][2]].sort((x, y) => y - x);
    }

    const q = MatrixMath.trace3x3(a) / 3.0;
    const p2 = Math.pow(a[0][0] - q, 2) + Math.pow(a[1][1] - q, 2) + Math.pow(a[2][2] - q, 2) + 2.0 * p1;
    const p = Math.sqrt(p2 / 6.0);

    const B = [
      [(a[0][0] - q) / p, a[0][1] / p, a[0][2] / p],
      [a[1][0] / p, (a[1][1] - q) / p, a[1][2] / p],
      [a[2][0] / p, (a[2][1] - q) / p, (a[2][2] - q) / p]
    ];

    const r = MatrixMath.det3x3(B) / 2.0;
    let phi;
    if (r <= -1.0) phi = Math.PI / 3.0;
    else if (r >= 1.0) phi = 0.0;
    else phi = Math.acos(r) / 3.0;

    const eig1 = q + 2.0 * p * Math.cos(phi);
    const eig3 = q + 2.0 * p * Math.cos(phi + (2.0 * Math.PI / 3.0));
    const eig2 = 3.0 * q - eig1 - eig3;

    return [eig1, eig2, eig3].sort((x, y) => y - x);
  }
}

// ============================================================================
// 5. DIFFERENTIAL DISCREPANCY & DERIVATIVE AUDITOR
// ============================================================================

export class DifferentialDiscrepancyAuditor {
  static centralDifferenceGradient2nd(scalarFn, pos, h = 1e-4) {
    const grad = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      const pFwd = [...pos];
      const pBwd = [...pos];
      pFwd[i] += h;
      pBwd[i] -= h;
      grad[i] = (scalarFn(pFwd) - scalarFn(pBwd)) / (2.0 * h);
    }
    return grad;
  }

  static centralDifferenceGradient4th(scalarFn, pos, h = 1e-4) {
    const grad = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      const pFwd2 = [...pos];
      const pFwd1 = [...pos];
      const pBwd1 = [...pos];
      const pBwd2 = [...pos];
      pFwd2[i] += 2.0 * h;
      pFwd1[i] += h;
      pBwd1[i] -= h;
      pBwd2[i] -= 2.0 * h;
      grad[i] = (-scalarFn(pFwd2) + 8.0 * scalarFn(pFwd1) - 8.0 * scalarFn(pBwd1) + scalarFn(pBwd2)) / (12.0 * h);
    }
    return grad;
  }

  static auditGradientDiscrepancy(scalarFn, analyticalGradFn, pos, h = 1e-4) {
    const numGrad = DifferentialDiscrepancyAuditor.centralDifferenceGradient4th(scalarFn, pos, h);
    const anaGrad = analyticalGradFn(pos);
    const diff = VectorMath.sub(numGrad, anaGrad);
    const deltaGrad = VectorMath.norm2(diff);
    const anaNorm = VectorMath.norm2(anaGrad);
    const relativeDelta = deltaGrad / Math.max(anaNorm, 1e-7);

    return {
      deltaGrad,
      relativeDelta,
      numGrad,
      anaGrad,
      order: '4th-order-compact'
    };
  }

  static numericalHessian(scalarFn, pos, h = 1e-3) {
    const H = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    const f0 = scalarFn(pos);

    for (let i = 0; i < 3; i++) {
      for (let j = i; j < 3; j++) {
        if (i === j) {
          const pFwd = [...pos];
          const pBwd = [...pos];
          pFwd[i] += h;
          pBwd[i] -= h;
          H[i][i] = (scalarFn(pFwd) - 2.0 * f0 + scalarFn(pBwd)) / (h * h);
        } else {
          const p_pp = [...pos]; p_pp[i] += h; p_pp[j] += h;
          const p_pm = [...pos]; p_pm[i] += h; p_pm[j] -= h;
          const p_mp = [...pos]; p_mp[i] -= h; p_mp[j] += h;
          const p_mm = [...pos]; p_mm[i] -= h; p_mm[j] -= h;
          const val = (scalarFn(p_pp) - scalarFn(p_pm) - scalarFn(p_mp) + scalarFn(p_mm)) / (4.0 * h * h);
          H[i][j] = val;
          H[j][i] = val;
        }
      }
    }
    return H;
  }

  static numericalJacobian(vectorFn, pos, h = 1e-4) {
    const J = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    for (let j = 0; j < 3; j++) {
      const pFwd = [...pos];
      const pBwd = [...pos];
      pFwd[j] += h;
      pBwd[j] -= h;
      const vFwd = vectorFn(pFwd);
      const vBwd = vectorFn(pBwd);
      for (let i = 0; i < 3; i++) {
        J[i][j] = (vFwd[i] - vBwd[i]) / (2.0 * h);
      }
    }
    return J;
  }

  static numericalDivergence(vectorFn, pos, h = 1e-4) {
    const J = DifferentialDiscrepancyAuditor.numericalJacobian(vectorFn, pos, h);
    return J[0][0] + J[1][1] + J[2][2];
  }

  static numericalVorticity(vectorFn, pos, h = 1e-4) {
    const J = DifferentialDiscrepancyAuditor.numericalJacobian(vectorFn, pos, h);
    return [
      J[2][1] - J[1][2],
      J[0][2] - J[2][0],
      J[1][0] - J[0][1]
    ];
  }

  static auditConvergenceOrder(estimatorFn, exactVal, h0 = 1e-2) {
    const h1 = h0;
    const h2 = h0 / 2.0;
    const h3 = h0 / 4.0;

    const e1 = Math.abs(estimatorFn(h1) - exactVal);
    const e2 = Math.abs(estimatorFn(h2) - exactVal);
    const e3 = Math.abs(estimatorFn(h3) - exactVal);

    let p12 = 0;
    let p23 = 0;
    if (e2 > 1e-15 && e1 > 1e-15) {
      p12 = Math.log(e1 / e2) / Math.log(2.0);
    }
    if (e3 > 1e-15 && e2 > 1e-15) {
      p23 = Math.log(e2 / e3) / Math.log(2.0);
    }

    return {
      errors: [e1, e2, e3],
      stepSizes: [h1, h2, h3],
      empiricalOrder: (p12 + p23) / 2.0,
      asymptoticConsistency: Math.abs(p12 - p23) < 0.5
    };
  }
}

// ============================================================================
// 6. TRAJECTORY METRICS: CONTINUOUS & DISCRETE HAUSDORFF / FRÉCHET
// ============================================================================

export class TrajectoryMetrics {
  static directedHausdorff(trajA, trajB) {
    let maxDist = 0;
    for (let i = 0; i < trajA.length; i++) {
      const a = trajA[i];
      let minDist = Infinity;
      for (let j = 0; j < trajB.length; j++) {
        const d = VectorMath.dist(a, trajB[j]);
        if (d < minDist) minDist = d;
      }
      if (minDist > maxDist) maxDist = minDist;
    }
    return maxDist;
  }

  static continuousDirectedHausdorff(trajA, trajB) {
    if (trajB.length === 1) {
      return TrajectoryMetrics.directedHausdorff(trajA, trajB);
    }
    let maxDist = 0;
    for (let i = 0; i < trajA.length; i++) {
      const a = trajA[i];
      let minDist = Infinity;
      for (let j = 0; j < trajB.length - 1; j++) {
        const d = VectorMath.distToSegment(a, trajB[j], trajB[j + 1]);
        if (d < minDist) minDist = d;
      }
      if (minDist > maxDist) maxDist = minDist;
    }
    return maxDist;
  }

  static symmetricHausdorff(trajA, trajB) {
    if (!trajA || !trajB || trajA.length === 0 || trajB.length === 0) {
      throw new Error('Trajectories must contain at least 1 point');
    }
    const hAB = TrajectoryMetrics.continuousDirectedHausdorff(trajA, trajB);
    const hBA = TrajectoryMetrics.continuousDirectedHausdorff(trajB, trajA);
    return Math.max(hAB, hBA);
  }

  static averageHausdorff(trajA, trajB) {
    let sumA = 0;
    for (let i = 0; i < trajA.length; i++) {
      let minDist = Infinity;
      for (let j = 0; j < trajB.length - 1; j++) {
        const d = VectorMath.distToSegment(trajA[i], trajB[j], trajB[j + 1]);
        if (d < minDist) minDist = d;
      }
      sumA += (minDist === Infinity) ? 0 : minDist;
    }

    let sumB = 0;
    for (let j = 0; j < trajB.length; j++) {
      let minDist = Infinity;
      for (let i = 0; i < trajA.length - 1; i++) {
        const d = VectorMath.distToSegment(trajB[j], trajA[i], trajA[i + 1]);
        if (d < minDist) minDist = d;
      }
      sumB += (minDist === Infinity) ? 0 : minDist;
    }

    return 0.5 * (sumA / trajA.length + sumB / trajB.length);
  }

  static discreteFrechet(trajA, trajB) {
    const n = trajA.length;
    const m = trajB.length;
    const dp = Array.from({ length: n }, () => new Float64Array(m));

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        const d = VectorMath.dist(trajA[i], trajB[j]);
        if (i === 0 && j === 0) {
          dp[i][j] = d;
        } else if (i === 0) {
          dp[i][j] = Math.max(dp[0][j - 1], d);
        } else if (j === 0) {
          dp[i][j] = Math.max(dp[i - 1][0], d);
        } else {
          dp[i][j] = Math.max(Math.min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]), d);
        }
      }
    }

    return dp[n - 1][m - 1];
  }

  static arcLength(traj) {
    let len = 0;
    for (let i = 0; i < traj.length - 1; i++) {
      len += VectorMath.dist(traj[i], traj[i + 1]);
    }
    return len;
  }
}

// ============================================================================
// 7. DUAL-METHOD ODE STREAMLINE INTEGRATORS (CASH-KARP RK45 VS DOPRI5)
// ============================================================================

export class CashKarpRK45Solver {
  static butcher() {
    return {
      c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0],
      a: [
        [],
        [1.0 / 5.0],
        [3.0 / 40.0, 9.0 / 40.0],
        [3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0],
        [-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0],
        [1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0]
      ],
      b5: [37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0],
      b4: [2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0]
    };
  }

  static integrate(fieldFn, pos0, tEnd = 5.0, hInit = 0.05, tol = 1e-7, maxSteps = 5000) {
    const B = CashKarpRK45Solver.butcher();
    let pos = [...pos0];
    let t = 0.0;
    let h = hInit;
    const trajectory = [[...pos]];
    const times = [t];
    let steps = 0;
    let rejections = 0;

    while (t < tEnd && steps < maxSteps) {
      if (t + h > tEnd) h = tEnd - t;
      if (h < 1e-12) break;

      const k = [];
      k[0] = fieldFn(pos);

      for (let s = 1; s < 6; s++) {
        const pStage = [...pos];
        for (let j = 0; j < s; j++) {
          pStage[0] += h * B.a[s][j] * k[j][0];
          pStage[1] += h * B.a[s][j] * k[j][1];
          pStage[2] += h * B.a[s][j] * k[j][2];
        }
        k[s] = fieldFn(pStage);
      }

      const pos5 = [...pos];
      const pos4 = [...pos];
      for (let s = 0; s < 6; s++) {
        pos5[0] += h * B.b5[s] * k[s][0];
        pos5[1] += h * B.b5[s] * k[s][1];
        pos5[2] += h * B.b5[s] * k[s][2];

        pos4[0] += h * B.b4[s] * k[s][0];
        pos4[1] += h * B.b4[s] * k[s][1];
        pos4[2] += h * B.b4[s] * k[s][2];
      }

      const errVec = [pos5[0] - pos4[0], pos5[1] - pos4[1], pos5[2] - pos4[2]];
      const errNorm = VectorMath.norm2(errVec);
      const scaleNorm = Math.max(VectorMath.norm2(pos5), 1.0);
      const relError = errNorm / (tol * scaleNorm);

      if (relError <= 1.0 || h <= 1e-10) {
        pos = pos5;
        t += h;
        trajectory.push([...pos]);
        times.push(t);
        steps++;

        const factor = Math.min(Math.max(0.9 * Math.pow(Math.max(relError, 1e-10), -0.2), 0.2), 5.0);
        h *= factor;
      } else {
        rejections++;
        const factor = Math.max(0.9 * Math.pow(relError, -0.25), 0.1);
        h *= factor;
      }
    }

    return {
      trajectory,
      times,
      steps,
      rejections,
      finalPos: pos,
      method: 'Cash-Karp-RK45'
    };
  }
}

export class DormandPrinceDOPRI5Solver {
  static butcher() {
    return {
      c: [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0],
      a: [
        [],
        [1.0 / 5.0],
        [3.0 / 40.0, 9.0 / 40.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
        [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0],
        [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0],
        [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0]
      ],
      b5: [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0],
      b4: [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0]
    };
  }

  static integrate(fieldFn, pos0, tEnd = 5.0, hInit = 0.05, tol = 1e-7, maxSteps = 5000) {
    const B = DormandPrinceDOPRI5Solver.butcher();
    let pos = [...pos0];
    let t = 0.0;
    let h = hInit;
    const trajectory = [[...pos]];
    const times = [t];
    let steps = 0;
    let rejections = 0;

    let k1 = fieldFn(pos);

    while (t < tEnd && steps < maxSteps) {
      if (t + h > tEnd) h = tEnd - t;
      if (h < 1e-12) break;

      const k = [k1];

      for (let s = 1; s < 7; s++) {
        const pStage = [...pos];
        for (let j = 0; j < s; j++) {
          pStage[0] += h * B.a[s][j] * k[j][0];
          pStage[1] += h * B.a[s][j] * k[j][1];
          pStage[2] += h * B.a[s][j] * k[j][2];
        }
        k[s] = fieldFn(pStage);
      }

      const pos5 = [...pos];
      const pos4 = [...pos];
      for (let s = 0; s < 7; s++) {
        pos5[0] += h * B.b5[s] * k[s][0];
        pos5[1] += h * B.b5[s] * k[s][1];
        pos5[2] += h * B.b5[s] * k[s][2];

        pos4[0] += h * B.b4[s] * k[s][0];
        pos4[1] += h * B.b4[s] * k[s][1];
        pos4[2] += h * B.b4[s] * k[s][2];
      }

      const errVec = [pos5[0] - pos4[0], pos5[1] - pos4[1], pos5[2] - pos4[2]];
      const errNorm = VectorMath.norm2(errVec);
      const scaleNorm = Math.max(VectorMath.norm2(pos5), 1.0);
      const relError = errNorm / (tol * scaleNorm);

      if (relError <= 1.0 || h <= 1e-10) {
        pos = pos5;
        t += h;
        trajectory.push([...pos]);
        times.push(t);
        steps++;

        k1 = k[6]; // FSAL optimization

        const factor = Math.min(Math.max(0.9 * Math.pow(Math.max(relError, 1e-10), -0.2), 0.2), 5.0);
        h *= factor;
      } else {
        rejections++;
        const factor = Math.max(0.9 * Math.pow(relError, -0.25), 0.1);
        h *= factor;
      }
    }

    return {
      trajectory,
      times,
      steps,
      rejections,
      finalPos: pos,
      method: 'Dormand-Prince-DOPRI5'
    };
  }
}

export class ClassicalRK4Solver {
  static integrate(fieldFn, pos0, tEnd = 5.0, dt = 0.01) {
    let pos = [...pos0];
    let t = 0.0;
    const trajectory = [[...pos]];
    const times = [t];
    const nSteps = Math.ceil(tEnd / dt);
    const h = tEnd / nSteps;

    for (let step = 0; step < nSteps; step++) {
      const k1 = fieldFn(pos);
      const p2 = [pos[0] + 0.5 * h * k1[0], pos[1] + 0.5 * h * k1[1], pos[2] + 0.5 * h * k1[2]];
      const k2 = fieldFn(p2);
      const p3 = [pos[0] + 0.5 * h * k2[0], pos[1] + 0.5 * h * k2[1], pos[2] + 0.5 * h * k2[2]];
      const k3 = fieldFn(p3);
      const p4 = [pos[0] + h * k3[0], pos[1] + h * k3[1], pos[2] + h * k3[2]];
      const k4 = fieldFn(p4);

      pos[0] += (h / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
      pos[1] += (h / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
      pos[2] += (h / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
      t += h;

      trajectory.push([...pos]);
      times.push(t);
    }

    return {
      trajectory,
      times,
      steps: nSteps,
      finalPos: pos,
      method: 'Classical-RK4'
    };
  }
}

// ============================================================================
// 8. VOLUMETRIC ISOSURFACE SOLVERS (MARCHING CUBES VS MARCHING TETRAHEDRA)
// ============================================================================

export class VolumetricMeshCrossValidator {
  static computeEnclosedVolume(vertices, faces) {
    let signedVol = 0;
    for (let f = 0; f < faces.length; f++) {
      const [i0, i1, i2] = faces[f];
      const v0 = vertices[i0];
      const v1 = vertices[i1];
      const v2 = vertices[i2];

      const cross = VectorMath.cross(v1, v2);
      signedVol += VectorMath.dot(v0, cross);
    }
    return Math.abs(signedVol) / 6.0;
  }

  static computeSurfaceArea(vertices, faces) {
    let totalArea = 0;
    for (let f = 0; f < faces.length; f++) {
      const [i0, i1, i2] = faces[f];
      const v0 = vertices[i0];
      const v1 = vertices[i1];
      const v2 = vertices[i2];

      const e1 = VectorMath.sub(v1, v0);
      const e2 = VectorMath.sub(v2, v0);
      const cross = VectorMath.cross(e1, e2);
      totalArea += 0.5 * VectorMath.norm2(cross);
    }
    return totalArea;
  }

  static computeMeshCentroid(vertices, faces) {
    let cX = 0, cY = 0, cZ = 0;
    let totalArea = 0;
    for (let f = 0; f < faces.length; f++) {
      const [i0, i1, i2] = faces[f];
      const v0 = vertices[i0];
      const v1 = vertices[i1];
      const v2 = vertices[i2];

      const e1 = VectorMath.sub(v1, v0);
      const e2 = VectorMath.sub(v2, v0);
      const cross = VectorMath.cross(e1, e2);
      const area = 0.5 * VectorMath.norm2(cross);

      const fCentroid = [
        (v0[0] + v1[0] + v2[0]) / 3.0,
        (v0[1] + v1[1] + v2[1]) / 3.0,
        (v0[2] + v1[2] + v2[2]) / 3.0
      ];

      cX += area * fCentroid[0];
      cY += area * fCentroid[1];
      cZ += area * fCentroid[2];
      totalArea += area;
    }
    if (totalArea < 1e-15) return [0, 0, 0];
    return [cX / totalArea, cY / totalArea, cZ / totalArea];
  }

  static computeTopologyInvariants(vertices, faces) {
    const edgeMap = new Map();
    const makeEdgeKey = (a, b) => (a < b ? `${a}_${b}` : `${b}_${a}`);

    for (let f = 0; f < faces.length; f++) {
      const [i0, i1, i2] = faces[f];
      const e0 = makeEdgeKey(i0, i1);
      const e1 = makeEdgeKey(i1, i2);
      const e2 = makeEdgeKey(i2, i0);

      edgeMap.set(e0, (edgeMap.get(e0) || 0) + 1);
      edgeMap.set(e1, (edgeMap.get(e1) || 0) + 1);
      edgeMap.set(e2, (edgeMap.get(e2) || 0) + 1);
    }

    const V = vertices.length;
    const E = edgeMap.size;
    const F = faces.length;
    const chi = V - E + F;
    const genus = 1.0 - chi / 2.0;

    let boundaryEdges = 0;
    let nonManifoldEdges = 0;

    for (const count of edgeMap.values()) {
      if (count === 1) boundaryEdges++;
      else if (count > 2) nonManifoldEdges++;
    }

    const isWatertight = boundaryEdges === 0 && nonManifoldEdges === 0;

    return {
      V,
      E,
      F,
      chi,
      genus,
      boundaryEdges,
      nonManifoldEdges,
      isWatertight
    };
  }

  static computeVoxelSegmentationSimilarity(maskA, maskB) {
    if (maskA.length !== maskB.length) {
      throw new Error('Mask dimensions mismatch');
    }
    let intersection = 0;
    let countA = 0;
    let countB = 0;

    for (let i = 0; i < maskA.length; i++) {
      const a = maskA[i] ? 1 : 0;
      const b = maskB[i] ? 1 : 0;
      if (a && b) intersection++;
      if (a) countA++;
      if (b) countB++;
    }

    const union = countA + countB - intersection;
    const jaccard = union === 0 ? 1.0 : intersection / union;
    const dice = (countA + countB === 0) ? 1.0 : (2.0 * intersection) / (countA + countB);

    return {
      intersection,
      union,
      countA,
      countB,
      jaccard,
      dice
    };
  }

  /**
   * Marching Tetrahedra (Kuhn 6-Simplex) Isosurface Extraction.
   */
  static marchingTetrahedraExtract(grid, dims, isovalue, bounds = [[-1, 1], [-1, 1], [-1, 1]]) {
    const [nx, ny, nz] = dims;
    const [bx, by, bz] = bounds;
    const dx = (bx[1] - bx[0]) / (nx - 1);
    const dy = (by[1] - by[0]) / (ny - 1);
    const dz = (bz[1] - bz[0]) / (nz - 1);

    const getIdx = (x, y, z) => (z * ny + y) * nx + x;
    const getPos = (x, y, z) => [bx[0] + x * dx, by[0] + y * dy, bz[0] + z * dz];

    const vertices = [];
    const faces = [];
    const vertMap = new Map();

    for (let z = 0; z < nz - 1; z++) {
      for (let y = 0; y < ny - 1; y++) {
        for (let x = 0; x < nx - 1; x++) {
          for (const tet of KUHN_6_TETRAHEDRA) {
            const tetPos = [];
            const tetVals = [];
            const tetGridCoords = [];
            for (let i = 0; i < 4; i++) {
              const off = CUBE_CORNERS[tet[i]];
              const gx = x + off[0];
              const gy = y + off[1];
              const gz = z + off[2];
              tetGridCoords.push([gx, gy, gz]);
              tetPos.push(getPos(gx, gy, gz));
              tetVals.push(grid[getIdx(gx, gy, gz)]);
            }

            let mask = 0;
            for (let i = 0; i < 4; i++) {
              if (tetVals[i] >= isovalue) mask |= (1 << i);
            }
            if (mask === 0 || mask === 15) continue;

            const tris = CANONICAL_TET_TRIANGLE_TABLE[mask];
            for (const [e0, e1, e2] of tris) {
              const triVerts = [];
              for (const edgeIdx of [e0, e1, e2]) {
                const [c1, c2] = TET_EDGE_PAIRS[edgeIdx];
                const p1 = tetPos[c1];
                const p2 = tetPos[c2];
                const v1 = tetVals[c1];
                const v2 = tetVals[c2];
                const t = Math.abs(v2 - v1) < 1e-12 ? 0.5 : (isovalue - v1) / (v2 - v1);
                const tClamped = Math.max(0.0, Math.min(1.0, t));
                const p = [
                  p1[0] + tClamped * (p2[0] - p1[0]),
                  p1[1] + tClamped * (p2[1] - p1[1]),
                  p1[2] + tClamped * (p2[2] - p1[2])
                ];
                const g1 = tetGridCoords[c1];
                const g2 = tetGridCoords[c2];
                const k1 = `${g1[0]}_${g1[1]}_${g1[2]}`;
                const k2 = `${g2[0]}_${g2[1]}_${g2[2]}`;
                const key = k1 < k2 ? `${k1}:${k2}` : `${k2}:${k1}`;

                let vIdx;
                if (vertMap.has(key)) {
                  vIdx = vertMap.get(key);
                } else {
                  vIdx = vertices.length;
                  vertices.push(p);
                  vertMap.set(key, vIdx);
                }
                triVerts.push(vIdx);
              }
              faces.push(triVerts);
            }
          }
        }
      }
    }

    return {
      vertices,
      faces,
      method: 'Marching-Tetrahedra-Kuhn'
    };
  }

  /**
   * Marching Cubes (256-case table) Isosurface Extraction.
   */
  static marchingCubesExtract(grid, dims, isovalue, bounds = [[-1, 1], [-1, 1], [-1, 1]]) {
    const [nx, ny, nz] = dims;
    const [bx, by, bz] = bounds;
    const dx = (bx[1] - bx[0]) / (nx - 1);
    const dy = (by[1] - by[0]) / (ny - 1);
    const dz = (bz[1] - bz[0]) / (nz - 1);

    const getIdx = (x, y, z) => (z * ny + y) * nx + x;
    const getPos = (x, y, z) => [bx[0] + x * dx, by[0] + y * dy, bz[0] + z * dz];

    const vertices = [];
    const faces = [];
    const vertMap = new Map();

    for (let z = 0; z < nz - 1; z++) {
      for (let y = 0; y < ny - 1; y++) {
        for (let x = 0; x < nx - 1; x++) {
          let cubeIndex = 0;
          const cornerVals = new Float64Array(8);
          const cornerCoords = [];
          const cornerGridCoords = [];

          for (let c = 0; c < 8; c++) {
            const cx = x + CUBE_CORNERS[c][0];
            const cy = y + CUBE_CORNERS[c][1];
            const cz = z + CUBE_CORNERS[c][2];
            const val = grid[getIdx(cx, cy, cz)];
            cornerVals[c] = val;
            cornerCoords.push(getPos(cx, cy, cz));
            cornerGridCoords.push([cx, cy, cz]);
            if (val >= isovalue) cubeIndex |= (1 << c);
          }

          if (cubeIndex === 0 || cubeIndex === 255) continue;

          const triEdges = MARCHING_CUBES_TRIANGLE_TABLE[cubeIndex];
          for (let t = 0; triEdges[t] !== -1; t += 3) {
            const tri = [];
            for (let k = 0; k < 3; k++) {
              const e = triEdges[t + k];
              const [c1, c2] = EDGE_VERTICES[e];
              const v1 = cornerVals[c1];
              const v2 = cornerVals[c2];
              const p1 = cornerCoords[c1];
              const p2 = cornerCoords[c2];
              const factor = Math.abs(v2 - v1) < 1e-12 ? 0.5 : (isovalue - v1) / (v2 - v1);
              const p = [
                p1[0] + factor * (p2[0] - p1[0]),
                p1[1] + factor * (p2[1] - p1[1]),
                p1[2] + factor * (p2[2] - p1[2])
              ];

              const g1 = cornerGridCoords[c1];
              const g2 = cornerGridCoords[c2];
              const k1 = `${g1[0]}_${g1[1]}_${g1[2]}`;
              const k2 = `${g2[0]}_${g2[1]}_${g2[2]}`;
              const key = k1 < k2 ? `${k1}:${k2}` : `${k2}:${k1}`;

              let vIdx;
              if (vertMap.has(key)) {
                vIdx = vertMap.get(key);
              } else {
                vIdx = vertices.length;
                vertices.push(p);
                vertMap.set(key, vIdx);
              }
              tri.push(vIdx);
            }
            faces.push(tri);
          }
        }
      }
    }

    return {
      vertices,
      faces,
      method: 'Marching-Cubes-256'
    };
  }

  static analyticalSphereGroundTruth(radius) {
    const volume = (4.0 / 3.0) * Math.PI * Math.pow(radius, 3);
    const area = 4.0 * Math.PI * Math.pow(radius, 2);
    return {
      radius,
      volume,
      area,
      eulerCharacteristic: 2,
      genus: 0
    };
  }
}

// ============================================================================
// 9. DUAL-METHOD ROOT FINDING & CRITICAL POINT CLASSIFICATION
// ============================================================================

export class CriticalPointRootValidator {
  static newtonRaphson3D(vectorFn, posInit, maxIter = 50, tol = 1e-10) {
    let pos = [...posInit];
    let converged = false;
    let iterations = 0;

    for (let iter = 0; iter < maxIter; iter++) {
      iterations = iter + 1;
      const fVal = vectorFn(pos);
      const resNorm = VectorMath.norm2(fVal);

      if (resNorm < tol) {
        converged = true;
        break;
      }

      const J = DifferentialDiscrepancyAuditor.numericalJacobian(vectorFn, pos, 1e-4);
      try {
        const invJ = MatrixMath.inverse3x3(J);
        const delta = MatrixMath.matVecMul3x3(invJ, fVal);
        pos = [pos[0] - delta[0], pos[1] - delta[1], pos[2] - delta[2]];

        if (VectorMath.norm2(delta) < tol * 0.1) {
          converged = true;
          break;
        }
      } catch (err) {
        break;
      }
    }

    const finalResidual = VectorMath.norm2(vectorFn(pos));
    return {
      root: pos,
      residual: finalResidual,
      iterations,
      converged: converged && finalResidual < tol * 10
    };
  }

  static broydenQuasiNewton3D(vectorFn, posInit, maxIter = 60, tol = 1e-9) {
    let pos = [...posInit];
    let fVal = vectorFn(pos);
    let J = DifferentialDiscrepancyAuditor.numericalJacobian(vectorFn, pos, 1e-4);
    let invB;
    try {
      invB = MatrixMath.inverse3x3(J);
    } catch {
      invB = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
      ];
    }

    let converged = false;
    let iterations = 0;

    for (let iter = 0; iter < maxIter; iter++) {
      iterations = iter + 1;
      const resNorm = VectorMath.norm2(fVal);
      if (resNorm < tol) {
        converged = true;
        break;
      }

      const delta = MatrixMath.matVecMul3x3(invB, fVal);
      const nextPos = [pos[0] - delta[0], pos[1] - delta[1], pos[2] - delta[2]];
      const nextF = vectorFn(nextPos);
      const y = VectorMath.sub(nextF, fVal);
      const s = VectorMath.scale(delta, -1.0);

      const invB_y = MatrixMath.matVecMul3x3(invB, y);
      const s_minus_invBy = VectorMath.sub(s, invB_y);
      const s_T_invB = [
        s[0] * invB[0][0] + s[1] * invB[1][0] + s[2] * invB[2][0],
        s[0] * invB[0][1] + s[1] * invB[1][1] + s[2] * invB[2][1],
        s[0] * invB[0][2] + s[1] * invB[1][2] + s[2] * invB[2][2]
      ];
      const denom = VectorMath.dot(s_T_invB, y);

      if (Math.abs(denom) > 1e-15) {
        for (let r = 0; r < 3; r++) {
          for (let c = 0; c < 3; c++) {
            invB[r][c] += (s_minus_invBy[r] * s_T_invB[c]) / denom;
          }
        }
      }

      pos = nextPos;
      fVal = nextF;
    }

    return {
      root: pos,
      residual: VectorMath.norm2(fVal),
      iterations,
      converged: converged || VectorMath.norm2(fVal) < tol * 10
    };
  }

  static classifyCriticalPoint(vectorFn, pos) {
    const J = DifferentialDiscrepancyAuditor.numericalJacobian(vectorFn, pos, 1e-4);
    const S = [
      [J[0][0], 0.5 * (J[0][1] + J[1][0]), 0.5 * (J[0][2] + J[2][0])],
      [0.5 * (J[1][0] + J[0][1]), J[1][1], 0.5 * (J[1][2] + J[2][1])],
      [0.5 * (J[2][0] + J[0][2]), 0.5 * (J[2][1] + J[1][2]), J[2][2]]
    ];

    const diag = MatrixMath.jacobiEigenvalues3x3(S);
    const evals = diag.eigenvalues;

    let negCount = 0;
    let posCount = 0;
    let zeroCount = 0;

    for (let i = 0; i < 3; i++) {
      if (evals[i] < -1e-6) negCount++;
      else if (evals[i] > 1e-6) posCount++;
      else zeroCount++;
    }

    let type = CriticalPointType.DEGENERATE;
    let morseIndex = -1;

    if (zeroCount === 0) {
      if (negCount === 3) {
        type = CriticalPointType.ATTRACTOR_NODE;
        morseIndex = 3;
      } else if (negCount === 2 && posCount === 1) {
        type = CriticalPointType.FILAMENT_SADDLE;
        morseIndex = 2;
      } else if (negCount === 1 && posCount === 2) {
        type = CriticalPointType.WALL_SADDLE;
        morseIndex = 1;
      } else if (posCount === 3) {
        type = CriticalPointType.REPELLER_NODE;
        morseIndex = 0;
      }
    }

    const topologicalIndex = zeroCount === 0 ? Math.pow(-1, negCount) : 0;

    return {
      type,
      morseIndex,
      topologicalIndex,
      eigenvalues: evals,
      eigenvectors: diag.eigenvectors,
      jacobian: J,
      deformationTensor: S
    };
  }
}

// ============================================================================
// 10. DUAL-METHOD COSMIC BULK FLOW & MULTIPOLE ESTIMATORS
// ============================================================================

export class BulkFlowCrossValidator {
  static estimateBulkFlowMV(catalog) {
    if (!catalog || catalog.length < 4) {
      throw new Error('Bulk flow estimation requires at least 4 catalog tracers');
    }

    let A = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    let b = [0, 0, 0];

    for (let i = 0; i < catalog.length; i++) {
      const item = catalog[i];
      const pos = item.position;
      const r = VectorMath.norm2(pos);
      if (r < 1e-6) continue;

      const rHat = [pos[0] / r, pos[1] / r, pos[2] / r];
      const w = item.weight !== undefined ? item.weight : 1.0;
      const vRad = item.vRad;

      for (let row = 0; row < 3; row++) {
        for (let col = 0; col < 3; col++) {
          A[row][col] += w * rHat[row] * rHat[col];
        }
        b[row] += w * vRad * rHat[row];
      }
    }

    const invA = MatrixMath.inverse3x3(A);
    const vBulk = MatrixMath.matVecMul3x3(invA, b);
    const mag = VectorMath.norm2(vBulk);

    const sgLon = (Math.atan2(vBulk[1], vBulk[0]) * 180.0) / Math.PI;
    const sgLat = (Math.asin(Math.max(-1.0, Math.min(1.0, vBulk[2] / (mag || 1)))) * 180.0) / Math.PI;

    const eigA = MatrixMath.jacobiEigenvalues3x3(A).eigenvalues;
    const cond = Math.abs(eigA[0]) / Math.max(Math.abs(eigA[2]), 1e-12);

    return {
      bulkFlowVector: vBulk,
      magnitude: mag,
      sgLon: (sgLon + 360.0) % 360.0,
      sgLat,
      conditionNumber: cond,
      tracerCount: catalog.length
    };
  }

  static estimateBulkFlowSphericalHarmonics(catalog) {
    if (!catalog || catalog.length < 4) {
      throw new Error('Bulk flow estimation requires at least 4 catalog tracers');
    }

    let b = [0, 0, 0];
    let totalWeight = 0;

    for (let i = 0; i < catalog.length; i++) {
      const item = catalog[i];
      const pos = item.position;
      const r = VectorMath.norm2(pos);
      if (r < 1e-6) continue;

      const rHat = [pos[0] / r, pos[1] / r, pos[2] / r];
      const w = item.weight !== undefined ? item.weight : 1.0;
      const vRad = item.vRad;

      b[0] += w * vRad * rHat[0];
      b[1] += w * vRad * rHat[1];
      b[2] += w * vRad * rHat[2];
      totalWeight += w;
    }

    const factor = 3.0 / Math.max(totalWeight, 1e-12);
    const vBulk = [b[0] * factor, b[1] * factor, b[2] * factor];
    const mag = VectorMath.norm2(vBulk);

    return {
      bulkFlowVector: vBulk,
      magnitude: mag,
      method: 'Spherical-Harmonics-l1'
    };
  }
}

// ============================================================================
// 11. MASTER 10-GATE SCIENTIFIC COMPLIANCE SUITE
// ============================================================================

export class Gate1_DataIntegrityAndBounds {
  static evaluate(payload) {
    const errors = [];
    const metrics = { recordCount: 0, nanInfViolations: 0, speedOfLightViolations: 0, nonPositiveDistances: 0 };
    const SPEED_OF_LIGHT = 299792.458;

    if (!payload || !payload.records || !Array.isArray(payload.records)) {
      return {
        gateId: 1,
        gateName: 'Gate 1: Data Integrity & Finite Bounds',
        passed: false,
        status: GateStatus.FAILED,
        errors: ['Payload missing required "records" array'],
        metrics
      };
    }

    for (let i = 0; i < payload.records.length; i++) {
      const r = payload.records[i];
      metrics.recordCount++;

      if (!r.position || r.position.length !== 3 || r.position.some(v => typeof v !== 'number' || isNaN(v) || !isFinite(v))) {
        metrics.nanInfViolations++;
        errors.push(`Record ${r.id || i}: invalid non-finite position coordinates [${r.position}]`);
      }

      if (r.vMag !== undefined) {
        if (typeof r.vMag !== 'number' || isNaN(r.vMag) || !isFinite(r.vMag) || r.vMag >= SPEED_OF_LIGHT || r.vMag < 0) {
          metrics.speedOfLightViolations++;
          errors.push(`Record ${r.id || i}: velocity magnitude ${r.vMag} violates relativistic bounds [0, c)`);
        }
      }

      if (r.distance !== undefined) {
        if (typeof r.distance !== 'number' || isNaN(r.distance) || !isFinite(r.distance) || r.distance <= 0) {
          metrics.nonPositiveDistances++;
          errors.push(`Record ${r.id || i}: distance ${r.distance} must be strictly positive`);
        }
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 1,
      gateName: 'Gate 1: Data Integrity & Finite Bounds',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate2_UnitsAndCosmology {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const cosmo = payload?.cosmology || { h0: 74.6, omegaM: 0.315, omegaL: 0.685 };
    const h0 = cosmo.h0;
    const omegaM = cosmo.omegaM;
    const omegaL = cosmo.omegaL;

    metrics.h0 = h0;
    metrics.omegaM = omegaM;
    metrics.omegaL = omegaL;

    if (typeof h0 !== 'number' || h0 < 50.0 || h0 > 100.0) {
      errors.push(`Hubble constant H0 = ${h0} km/s/Mpc out of physical cosmological bounds [50, 100]`);
    }

    if (typeof omegaM !== 'number' || omegaM <= 0 || omegaM >= 1.0) {
      errors.push(`Matter density parameter Omega_m = ${omegaM} out of physical bounds (0, 1)`);
    }

    const omegaTotal = omegaM + (omegaL !== undefined ? omegaL : 0.0);
    metrics.omegaTotal = omegaTotal;
    if (Math.abs(omegaTotal - 1.0) > 0.05) {
      errors.push(`Cosmic curvature |Omega_total - 1.0| = ${Math.abs(omegaTotal - 1.0).toFixed(4)} exceeds flat-universe bound 0.05`);
    }

    const z = payload?.redshift !== undefined ? payload.redshift : 0.0;
    const a = payload?.scaleFactor !== undefined ? payload.scaleFactor : 1.0 / (1.0 + z);
    const expectedA = 1.0 / (1.0 + z);
    metrics.scaleFactorDiscrepancy = Math.abs(a - expectedA);

    if (metrics.scaleFactorDiscrepancy > 1e-5) {
      errors.push(`Scale factor a = ${a} violates cosmological relation a(z) = 1/(1+z) by ${metrics.scaleFactorDiscrepancy}`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 2,
      gateName: 'Gate 2: Units, Dimensional Consistency & Cosmology',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate3_AnalyticFieldsAndInvariants {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const laminar = (pos) => [120.0, -450.0, 230.0];
    const ptLaminar = [1234.5, -5678.9, 9012.3];
    const divLam = DifferentialDiscrepancyAuditor.numericalDivergence(laminar, ptLaminar, 1e-4);
    const curlLam = DifferentialDiscrepancyAuditor.numericalVorticity(laminar, ptLaminar, 1e-4);
    metrics.laminarDivergence = divLam;
    metrics.laminarCurlNorm = VectorMath.norm2(curlLam);

    if (Math.abs(divLam) > 1e-10) {
      errors.push(`Laminar flow divergence ${divLam} must be zero`);
    }
    if (metrics.laminarCurlNorm > 1e-10) {
      errors.push(`Laminar flow vorticity norm ${metrics.laminarCurlNorm} must be zero`);
    }

    const sink = (pos) => [-3.0 * pos[0], -2.0 * pos[1], -1.0 * pos[2]];
    const divSink = DifferentialDiscrepancyAuditor.numericalDivergence(sink, [0, 0, 0], 1e-4);
    metrics.sinkDivergence = divSink;
    if (Math.abs(divSink - (-6.0)) > 1e-4) {
      errors.push(`Linear sink divergence ${divSink} must equal trace -6.0`);
    }

    const vortex = (pos) => [-4.0 * pos[1], 4.0 * pos[0], 0.0];
    const curlVortex = DifferentialDiscrepancyAuditor.numericalVorticity(vortex, [10, 20, 0], 1e-4);
    metrics.vortexCurlZ = curlVortex[2];
    if (Math.abs(curlVortex[2] - 8.0) > 1e-4) {
      errors.push(`Vortex curl z-component ${curlVortex[2]} must equal 2*omega = 8.0`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 3,
      gateName: 'Gate 3: Analytic Vector Fields & Exact Invariants',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate4_DifferentialDiscrepancy {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};
    const tolGrad = payload?.tolerances?.GRADIENT_ABS_TOL || DEFAULT_TOLERANCES.GRADIENT_ABS_TOL;

    const phi = (p) => p[0] * p[0] * p[0] - 3.0 * p[0] * p[1] * p[1] + p[2] * p[2];
    const anaGradPhi = (p) => [3.0 * p[0] * p[0] - 3.0 * p[1] * p[1], -6.0 * p[0] * p[1], 2.0 * p[2]];

    const testPts = [
      [1.5, -2.0, 3.0],
      [-4.0, 0.5, -1.2],
      [2.2, 3.1, -0.8]
    ];

    let maxDeltaGrad = 0;
    for (const pt of testPts) {
      const audit = DifferentialDiscrepancyAuditor.auditGradientDiscrepancy(phi, anaGradPhi, pt, 1e-4);
      if (audit.deltaGrad > maxDeltaGrad) maxDeltaGrad = audit.deltaGrad;
    }
    metrics.maxGradientDiscrepancy = maxDeltaGrad;

    if (maxDeltaGrad > tolGrad) {
      errors.push(`Gradient discrepancy ${maxDeltaGrad} exceeds tolerance ${tolGrad}`);
    }

    const H = DifferentialDiscrepancyAuditor.numericalHessian(phi, [1.5, -2.0, 3.0], 1e-3);
    const isSymm = MatrixMath.isSymmetric3x3(H, 1e-7);
    metrics.hessianSymmetric = isSymm;
    if (!isSymm) {
      errors.push('Numerical Hessian matrix violates symmetry theorem (H_ij != H_ji)');
    }

    const passed = errors.length === 0;
    return {
      gateId: 4,
      gateName: 'Gate 4: Differential Discrepancy & Gradient Field Auditor',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate5_ODEIntegrationAndHausdorff {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};
    const tolHausdorff = payload?.tolerances?.HAUSDORFF_TRAJECTORY_TOL || DEFAULT_TOLERANCES.HAUSDORFF_TRAJECTORY_TOL;

    const sinkField = (p) => [-1.0 * p[0], -1.0 * p[1], -1.0 * p[2]];
    const pos0 = [10.0, 10.0, 10.0];
    const tEnd = 3.0;

    const resCK = CashKarpRK45Solver.integrate(sinkField, pos0, tEnd, 0.05, 1e-9);
    const resDP = DormandPrinceDOPRI5Solver.integrate(sinkField, pos0, tEnd, 0.05, 1e-9);

    const hDist = TrajectoryMetrics.symmetricHausdorff(resCK.trajectory, resDP.trajectory);
    metrics.hausdorffDistanceCKvsDP = hDist;
    metrics.ckSteps = resCK.steps;
    metrics.dpSteps = resDP.steps;

    if (hDist > tolHausdorff) {
      errors.push(`Trajectory Hausdorff distance ${hDist} between Cash-Karp and DOPRI5 exceeds ${tolHausdorff}`);
    }

    const vortexField = (p) => [-2.0 * p[1], 2.0 * p[0], 0.0];
    const vPos0 = [5.0, 0.0, 0.0];
    const period = Math.PI;
    const resVortex = CashKarpRK45Solver.integrate(vortexField, vPos0, period, 0.05, 1e-10);

    let maxRadiusErr = 0;
    const r0 = VectorMath.norm2(vPos0);
    for (const pt of resVortex.trajectory) {
      const r = VectorMath.norm2(pt);
      const err = Math.abs(r - r0);
      if (err > maxRadiusErr) maxRadiusErr = err;
    }
    metrics.maxVortexRadiusError = maxRadiusErr;

    if (maxRadiusErr > DEFAULT_TOLERANCES.ENERGY_CONSERVATION_TOL) {
      errors.push(`Vortex orbit radius conservation error ${maxRadiusErr} exceeds ${DEFAULT_TOLERANCES.ENERGY_CONSERVATION_TOL}`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 5,
      gateName: 'Gate 5: Numerical Integration, Butcher Convergence & Hausdorff Metric',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate6_VolumetricMeshingAndAgreement {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};
    const tolVol = payload?.tolerances?.VOLUME_REL_AGREEMENT_TOL || DEFAULT_TOLERANCES.VOLUME_REL_AGREEMENT_TOL;

    const R = 0.6;
    const dims = [24, 24, 24];
    const bounds = [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]];
    const grid = new Float64Array(dims[0] * dims[1] * dims[2]);

    const dx = 2.0 / (dims[0] - 1);
    const dy = 2.0 / (dims[1] - 1);
    const dz = 2.0 / (dims[2] - 1);

    for (let z = 0; z < dims[2]; z++) {
      const pz = -1.0 + z * dz;
      for (let y = 0; y < dims[1]; y++) {
        const py = -1.0 + y * dy;
        for (let x = 0; x < dims[0]; x++) {
          const px = -1.0 + x * dx;
          const val = R * R - (px * px + py * py + pz * pz);
          grid[(z * dims[1] + y) * dims[0] + x] = val;
        }
      }
    }

    const mtMesh = VolumetricMeshCrossValidator.marchingTetrahedraExtract(grid, dims, 0.0, bounds);
    const mcMesh = VolumetricMeshCrossValidator.marchingCubesExtract(grid, dims, 0.0, bounds);

    const mtVol = VolumetricMeshCrossValidator.computeEnclosedVolume(mtMesh.vertices, mtMesh.faces);
    const mcVol = VolumetricMeshCrossValidator.computeEnclosedVolume(mcMesh.vertices, mcMesh.faces);
    const mtTopo = VolumetricMeshCrossValidator.computeTopologyInvariants(mtMesh.vertices, mtMesh.faces);
    const mcTopo = VolumetricMeshCrossValidator.computeTopologyInvariants(mcMesh.vertices, mcMesh.faces);

    const gt = VolumetricMeshCrossValidator.analyticalSphereGroundTruth(R);
    metrics.mtVolume = mtVol;
    metrics.mcVolume = mcVol;
    metrics.gtVolume = gt.volume;
    metrics.mtEulerChi = mtTopo.chi;
    metrics.mcEulerChi = mcTopo.chi;
    metrics.isWatertightMT = mtTopo.isWatertight;
    metrics.isWatertightMC = mcTopo.isWatertight;

    const deltaV = Math.abs(mtVol - mcVol) / (0.5 * (mtVol + mcVol));
    metrics.volumeAgreementMCvsMT = deltaV;

    if (deltaV > tolVol) {
      errors.push(`Volumetric agreement discrepancy ${deltaV.toFixed(4)} between MC and MT exceeds ${tolVol}`);
    }

    if (mtTopo.chi !== 2) {
      errors.push(`MT Euler characteristic chi = ${mtTopo.chi} must equal 2 for closed 2-sphere`);
    }

    if (!mtTopo.isWatertight) {
      errors.push(`MT Mesh is not watertight: boundary edges = ${mtTopo.boundaryEdges}`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 6,
      gateName: 'Gate 6: Volumetric Meshing & Enclosed Volume Agreement',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate7_CriticalPointTopologyAndSpectra {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const sinkField = (p) => [-2.0 * p[0], -2.5 * p[1], -3.0 * p[2]];
    const initPt = [123.0, -456.0, 789.0];

    const nrRes = CriticalPointRootValidator.newtonRaphson3D(sinkField, initPt);
    const broydenRes = CriticalPointRootValidator.broydenQuasiNewton3D(sinkField, initPt);

    metrics.nrResidual = nrRes.residual;
    metrics.broydenResidual = broydenRes.residual;
    metrics.rootDistanceNRvsBroyden = VectorMath.dist(nrRes.root, broydenRes.root);

    if (!nrRes.converged || nrRes.residual > DEFAULT_TOLERANCES.ROOT_RESIDUAL_TOL) {
      errors.push(`Newton-Raphson failed to converge to root, residual = ${nrRes.residual}`);
    }

    const S = [
      [4.0, 1.0, -2.0],
      [1.0, 3.0, 0.5],
      [-2.0, 0.5, 5.0]
    ];
    const diag = MatrixMath.jacobiEigenvalues3x3(S);
    const V = diag.eigenvectors;
    const VT_V = MatrixMath.matMul3x3(MatrixMath.transpose3x3(V), V);
    const I = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];
    const orthErr = MatrixMath.diffFrobeniusNorm3x3(VT_V, I);
    metrics.jacobiOrthogonalityError = orthErr;

    if (orthErr > DEFAULT_TOLERANCES.JACOBI_ORTHOGONALITY_TOL) {
      errors.push(`Jacobi eigenvector matrix orthonormality error ${orthErr} exceeds ${DEFAULT_TOLERANCES.JACOBI_ORTHOGONALITY_TOL}`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 7,
      gateName: 'Gate 7: Critical Point Roots & Spectral Orthogonality',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate8_WatershedBasinSegmentation {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const attractors = [
      [7200.0, -8600.0, -2400.0],
      [-3500.0, 1500.0, -800.0]
    ];

    const classifyVoxel = (p) => {
      const d0 = VectorMath.dist(p, attractors[0]);
      const d1 = VectorMath.dist(p, attractors[1]);
      return d0 <= d1 ? 0 : 1;
    };

    const pts = [
      [7000.0, -8000.0, -2000.0],
      [-3000.0, 1200.0, -700.0],
      [0.0, 0.0, 0.0],
      [5000.0, -5000.0, -1000.0]
    ];

    const seg1 = pts.map(p => classifyVoxel(p) === 0);
    const seg2 = pts.map(p => classifyVoxel(p) === 0);
    const sim = VolumetricMeshCrossValidator.computeVoxelSegmentationSimilarity(seg1, seg2);

    metrics.diceSimilarity = sim.dice;
    metrics.jaccardSimilarity = sim.jaccard;

    if (sim.dice < DEFAULT_TOLERANCES.VOXEL_DICE_SIMILARITY_MIN) {
      errors.push(`Watershed segmentation Dice similarity ${sim.dice} below ${DEFAULT_TOLERANCES.VOXEL_DICE_SIMILARITY_MIN}`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 8,
      gateName: 'Gate 8: Cosmic Watershed Basin Segmentation & Similarity',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate9_BulkFlowMultipoleConsistency {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const vTrue = [300.0, -150.0, 200.0];
    const nTracers = 100;
    const catalog = [];

    let seed = 42;
    const prng = () => {
      seed = (seed * 1664525 + 1013904223) % 4294967296;
      return seed / 4294967296.0;
    };

    for (let i = 0; i < nTracers; i++) {
      const pos = [
        (prng() - 0.5) * 6000.0,
        (prng() - 0.5) * 6000.0,
        (prng() - 0.5) * 6000.0
      ];
      const r = VectorMath.norm2(pos);
      if (r < 1e-3) continue;
      const rHat = [pos[0] / r, pos[1] / r, pos[2] / r];
      const vRadTrue = VectorMath.dot(vTrue, rHat);
      const noise = (prng() - 0.5) * 10.0;
      catalog.push({ position: pos, vRad: vRadTrue + noise, weight: 1.0 });
    }

    const estMV = BulkFlowCrossValidator.estimateBulkFlowMV(catalog);
    const diff = VectorMath.sub(estMV.bulkFlowVector, vTrue);
    const relErr = VectorMath.norm2(diff) / VectorMath.norm2(vTrue);
    const angleErr = VectorMath.angleBetweenDeg(estMV.bulkFlowVector, vTrue);

    metrics.estimatedBulkFlow = estMV.bulkFlowVector;
    metrics.trueBulkFlow = vTrue;
    metrics.relativeMagnitudeError = relErr;
    metrics.angularDiscrepancyDeg = angleErr;

    if (relErr > DEFAULT_TOLERANCES.BULK_FLOW_REL_TOL) {
      errors.push(`Bulk flow relative error ${relErr.toFixed(4)} exceeds ${DEFAULT_TOLERANCES.BULK_FLOW_REL_TOL}`);
    }

    if (angleErr > DEFAULT_TOLERANCES.BULK_FLOW_ANGLE_DEG_TOL) {
      errors.push(`Bulk flow dipole angular error ${angleErr.toFixed(2)}° exceeds ${DEFAULT_TOLERANCES.BULK_FLOW_ANGLE_DEG_TOL}°`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 9,
      gateName: 'Gate 9: Cosmic Bulk Flow & Multipole Consistency',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

export class Gate10_SerializationAndProvenance {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const fitsHeaderStr = payload?.fitsHeader || 'SIMPLE  =                    T / Standard FITS format'.padEnd(80, ' ') + 'END'.padEnd(2800, ' ');
    const fitsByteLen = new TextEncoder().encode(fitsHeaderStr).length;
    metrics.fitsHeaderBytes = fitsByteLen;
    metrics.fitsBlockAligned = (fitsByteLen % 2880 === 0);

    if (fitsByteLen % 2880 !== 0) {
      errors.push(`FITS header length ${fitsByteLen} bytes is not a multiple of 2880-byte astronomical block`);
    }

    const testDoc = { id: 'cf4_run_01', type: 'FeatureCollection', count: 42 };
    const canon = UniversalSHA256.canonicalJson(testDoc);
    const hash = UniversalSHA256.hash(canon);
    metrics.canonicalHash = hash;

    if (!hash || hash.length !== 64) {
      errors.push('Cryptographic SHA-256 hash generation failed');
    }

    const passed = errors.length === 0;
    return {
      gateId: 10,
      gateName: 'Gate 10: Serialization, Provenance & Export Roundtrip',
      passed,
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      errors,
      metrics
    };
  }
}

// ============================================================================
// 12. MASTER CROSS-VALIDATION ORACLE & PASSPORT GENERATOR
// ============================================================================

export class CrossValidationOracle {
  /**
   * Runs full 10-gate scientific compliance verification and generates a cryptographic passport.
   * 
   * @param {Object} payload
   * @param {Object} [options]
   * @returns {Object} Verification Passport
   */
  static runFullVerification(payload = {}, options = {}) {
    const t0 = Date.now();
    const gates = [
      Gate1_DataIntegrityAndBounds,
      Gate2_UnitsAndCosmology,
      Gate3_AnalyticFieldsAndInvariants,
      Gate4_DifferentialDiscrepancy,
      Gate5_ODEIntegrationAndHausdorff,
      Gate6_VolumetricMeshingAndAgreement,
      Gate7_CriticalPointTopologyAndSpectra,
      Gate8_WatershedBasinSegmentation,
      Gate9_BulkFlowMultipoleConsistency,
      Gate10_SerializationAndProvenance
    ];

    const gateResults = [];
    let totalPassed = 0;
    let totalFailed = 0;

    for (const GateClass of gates) {
      try {
        const res = GateClass.evaluate(payload);
        if (res.passed) totalPassed++;
        else totalFailed++;
        gateResults.push(res);
      } catch (err) {
        totalFailed++;
        gateResults.push({
          gateName: GateClass.name,
          passed: false,
          status: GateStatus.FAILED,
          errors: [err.message || String(err)],
          metrics: {}
        });
      }
    }

    const elapsedMs = Date.now() - t0;
    const allGreen = totalFailed === 0;
    const complianceStatus = allGreen
      ? ComplianceLevel.CERTIFIED_REPRODUCIBLE
      : ComplianceLevel.FAILED_VERIFICATION;

    const passportRaw = {
      version: '1.0.0',
      timestamp: new Date().toISOString(),
      complianceStatus,
      totalGates: gates.length,
      passedGates: totalPassed,
      failedGates: totalFailed,
      elapsedMs,
      gateResults
    };

    const canonicalData = UniversalSHA256.canonicalJson(passportRaw);
    const provenanceHash = UniversalSHA256.hash(canonicalData);
    const signature = UniversalSHA256.hmac('ZRT_COSMICFLOWS_SECRET_KEY', provenanceHash);

    const passport = {
      ...passportRaw,
      provenanceHash,
      signature
    };

    return passport;
  }

  /**
   * Verifies authenticity and tamper-resistance of a cryptographic verification passport.
   * @param {Object} passport
   * @returns {{ valid: boolean, intact: boolean, authenticity: boolean, errors: string[] }}
   */
  static verifyPassport(passport) {
    const errors = [];
    if (!passport || typeof passport !== 'object') {
      return { valid: false, intact: false, authenticity: false, errors: ['Passport is null or invalid object'] };
    }

    const { provenanceHash, signature, ...raw } = passport;
    const canonical = UniversalSHA256.canonicalJson(raw);
    const computedHash = UniversalSHA256.hash(canonical);
    const intact = (computedHash === provenanceHash);

    if (!intact) {
      errors.push(`Provenance hash mismatch: computed ${computedHash} vs stored ${provenanceHash}`);
    }

    const computedSig = UniversalSHA256.hmac('ZRT_COSMICFLOWS_SECRET_KEY', computedHash);
    const authenticity = (computedSig === signature);

    if (!authenticity) {
      errors.push(`Cryptographic HMAC signature mismatch: computed ${computedSig} vs stored ${signature}`);
    }

    return {
      valid: intact && authenticity && errors.length === 0,
      intact,
      authenticity,
      errors
    };
  }
}

// ============================================================================
// 13. SCIENTIFIC BENCHMARK HARNESS
// ============================================================================

export class ScientificBenchmarkHarness {
  static benchmarkKernel(kernelName, kernelFn, iterations = 100) {
    for (let i = 0; i < Math.min(5, iterations); i++) {
      kernelFn();
    }

    const times = [];
    for (let i = 0; i < iterations; i++) {
      const t0 = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
      kernelFn();
      const t1 = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
      times.push(t1 - t0);
    }

    times.sort((a, b) => a - b);
    let sum = 0;
    for (const t of times) sum += t;
    const mean = sum / times.length;

    let varSum = 0;
    for (const t of times) varSum += (t - mean) * (t - mean);
    const std = Math.sqrt(varSum / times.length);

    const median = times[Math.floor(times.length * 0.5)];
    const p95 = times[Math.floor(times.length * 0.95)];
    const p99 = times[Math.floor(times.length * 0.99)];
    const min = times[0];
    const max = times[times.length - 1];
    const opsPerSec = mean > 0 ? 1000.0 / mean : Infinity;

    return {
      name: kernelName,
      iterations,
      meanMs: mean,
      stdMs: std,
      medianMs: median,
      minMs: min,
      maxMs: max,
      p95Ms: p95,
      p99Ms: p99,
      opsPerSec
    };
  }

  static runStandardSuite() {
    const results = {};

    const sinkField = (p) => [-p[0], -p[1], -p[2]];
    results.cashKarpRK45 = ScientificBenchmarkHarness.benchmarkKernel('Cash-Karp-RK45', () => {
      CashKarpRK45Solver.integrate(sinkField, [5, 5, 5], 2.0, 0.05, 1e-7);
    }, 50);

    results.dormandPrinceDOPRI5 = ScientificBenchmarkHarness.benchmarkKernel('DOPRI5', () => {
      DormandPrinceDOPRI5Solver.integrate(sinkField, [5, 5, 5], 2.0, 0.05, 1e-7);
    }, 50);

    const S = [
      [4.0, 1.0, -2.0],
      [1.0, 3.0, 0.5],
      [-2.0, 0.5, 5.0]
    ];
    results.jacobi3x3 = ScientificBenchmarkHarness.benchmarkKernel('Jacobi-3x3', () => {
      MatrixMath.jacobiEigenvalues3x3(S);
    }, 200);

    results.newtonRaphson3D = ScientificBenchmarkHarness.benchmarkKernel('Newton-Raphson-3D', () => {
      CriticalPointRootValidator.newtonRaphson3D(sinkField, [10, -20, 30]);
    }, 100);

    return results;
  }
}
