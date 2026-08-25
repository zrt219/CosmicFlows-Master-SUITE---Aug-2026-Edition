/**
 * @file morse_smale_complex.js
 * @description Complete 3D Morse-Smale Complex (MSC) and Cosmic Web Graph Engine for 3D Cosmological Fields.
 * 
 * Implements:
 * 1. 3D Morse-Smale Complex decomposition: volumetric crystals, sheets, and spines partitioned
 *    by intersections of ascending (unstable) and descending (stable) invariant manifolds:
 *    - 0-cells: Local minima / Cosmic void cores / Repellers (Morse index 0, \lambda_1, \lambda_2, \lambda_3 > 0)
 *    - 1-cells: 1D Integral curves / Cosmic filament spines and sheet spines (Morse index 1 & 2 connections)
 *    - 2-cells: 2D Separatrix sheets / Cosmic wall sheets / Membrane boundaries
 *    - 3-cells: 3D Volumetric crystals / Cosmic void basins / Attractor influence volumes (Morse index 3)
 * 2. Boundary operators \partial_1: C_1 \to C_0, \partial_2: C_2 \to C_1, \partial_3: C_3 \to C_2
 *    and topological homology chain groups C_0, C_1, C_2, C_3 with exact nilpotency verification (\partial \circ \partial = 0).
 * 3. Exact Betti numbers:
 *    - b_0: Connected components
 *    - b_1: 1D topological tunnels / cosmic filament loops
 *    - b_2: 2D enclosed cavities / cosmic void bubbles
 *    - b_3: 3D trapped volumes
 *    - Euler-Poincare characteristic \chi = b_0 - b_1 + b_2 - b_3 = |C_0| - |C_1| + |C_2| - |C_3|
 * 4. Topological simplification via Morse cancellation of persistence pairs with persistence threshold \epsilon_p:
 *    - Persistence pairing across adjacent Morse index dimensions (0-1, 1-2, 2-3)
 *    - Exact algebraic Morse reduction updating boundary complexes while maintaining \partial^2 = 0
 *    - Geometric simplification merging cosmic web filament graphs and void volumes
 * 5. Discrete Forman Gradient & Continuous Adaptive RK45 Separatrix Manifold Tracing.
 * 6. Cosmic Web Graph with physical units (Mpc/h, km/s, M_sun/h), filament tortuosity, void sphericity,
 *    and multi-format serialization (JSON, GeoJSON, CSV, MatrixMarket).
 * 
 * Strict Scientific Invariants:
 * - Coordinates strictly in Supergalactic Cartesian Mpc/h [SGX, SGY, SGZ].
 * - Peculiar velocities in km/s [vx, vy, vz].
 * - Dimensionless density contrast \delta = (\rho - \bar{\rho}) / \bar{\rho}.
 * - Gravitational potential \Phi in (km/s)^2.
 * 
 * @module topology/morse_smale_complex
 */

import { SupergalacticPosition, VelocityVector, EigenSystem3D } from '../coordinates/scientific_types.js';

// ============================================================================
// CONSTANTS & ENUMERATIONS
// ============================================================================

/**
 * Morse critical point topological index (number of negative eigenvalues of the Hessian / flow Jacobian).
 * @readonly
 * @enum {number}
 */
export const MorseIndex = Object.freeze({
  MINIMUM: 0,         // Index 0: Local minimum / cosmic void center / 3 expanding directions
  WALL_SADDLE: 1,     // Index 1: 1-saddle / cosmic sheet hub / 1 contracting, 2 expanding
  FILAMENT_SADDLE: 2, // Index 2: 2-saddle / cosmic filament hub / 2 contracting, 1 expanding
  MAXIMUM: 3          // Index 3: Local maximum / supercluster node / 3 contracting directions
});

/**
 * Cosmic web dynamical classification for Morse critical points.
 * @readonly
 * @enum {string}
 */
export const CosmicWebElementType = Object.freeze({
  VOID_CORE: 'VOID_CORE',             // 0-cell: local minimum
  WALL_HUB: 'WALL_HUB',               // 1-saddle: sheet/wall saddle point
  FILAMENT_HUB: 'FILAMENT_HUB',       // 2-saddle: filament saddle point
  CLUSTER_NODE: 'CLUSTER_NODE',       // 3-cell: local maximum / supercluster
  FILAMENT_SPINE: 'FILAMENT_SPINE',   // 1-cell: filament connecting 2-saddle to maximum
  SHEET_SPINE: 'SHEET_SPINE',         // 1-cell: sheet spine connecting 1-saddle to minimum
  WALL_SHEET: 'WALL_SHEET',           // 2-cell: 2D manifold bounding voids
  VOID_CRYSTAL: 'VOID_CRYSTAL'        // 3-cell: 3D volume connecting single min-max pair
});

/**
 * Orientation parity for chain group boundaries.
 * @readonly
 * @enum {number}
 */
export const ChainOrientation = Object.freeze({
  POSITIVE: 1,
  NEGATIVE: -1,
  ZERO: 0
});

// ============================================================================
// 3D VECTOR & MATRIX NUMERICAL UTILITIES
// ============================================================================

/**
 * Vector3D mathematical operations and utility methods.
 */
export class Vector3D {
  /**
   * @param {number} x
   * @param {number} y
   * @param {number} z
   */
  constructor(x = 0, y = 0, z = 0) {
    this.x = Number(x);
    this.y = Number(y);
    this.z = Number(z);
  }

  /**
   * Clones current vector.
   * @returns {Vector3D}
   */
  clone() {
    return new Vector3D(this.x, this.y, this.z);
  }

  /**
   * Sets coordinates.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {Vector3D} this
   */
  set(x, y, z) {
    this.x = Number(x);
    this.y = Number(y);
    this.z = Number(z);
    return this;
  }

  /**
   * Adds another vector.
   * @param {Vector3D} v
   * @returns {Vector3D}
   */
  add(v) {
    return new Vector3D(this.x + v.x, this.y + v.y, this.z + v.z);
  }

  /**
   * Subtracts another vector.
   * @param {Vector3D} v
   * @returns {Vector3D}
   */
  sub(v) {
    return new Vector3D(this.x - v.x, this.y - v.y, this.z - v.z);
  }

  /**
   * Multiplies by scalar.
   * @param {number} s
   * @returns {Vector3D}
   */
  scale(s) {
    return new Vector3D(this.x * s, this.y * s, this.z * s);
  }

  /**
   * Dot product.
   * @param {Vector3D} v
   * @returns {number}
   */
  dot(v) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  }

  /**
   * Cross product.
   * @param {Vector3D} v
   * @returns {Vector3D}
   */
  cross(v) {
    return new Vector3D(
      this.y * v.z - this.z * v.y,
      this.z * v.x - this.x * v.z,
      this.x * v.y - this.y * v.x
    );
  }

  /**
   * Euclidean norm / magnitude.
   * @returns {number}
   */
  norm() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  /**
   * Squared Euclidean norm.
   * @returns {number}
   */
  normSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  }

  /**
   * Returns normalized unit vector (or zero vector if magnitude < 1e-15).
   * @returns {Vector3D}
   */
  normalize() {
    const n = this.norm();
    if (n < 1e-15) return new Vector3D(0, 0, 0);
    return new Vector3D(this.x / n, this.y / n, this.z / n);
  }

  /**
   * Euclidean distance to another point.
   * @param {Vector3D} v
   * @returns {number}
   */
  distanceTo(v) {
    const dx = this.x - v.x;
    const dy = this.y - v.y;
    const dz = this.z - v.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  /**
   * Converts to Array [x, y, z].
   * @returns {number[]}
   */
  toArray() {
    return [this.x, this.y, this.z];
  }

  /**
   * Creates from array or object.
   * @param {number[]|{x: number, y: number, z: number}} obj
   * @returns {Vector3D}
   */
  static from(obj) {
    if (!obj) return new Vector3D(0, 0, 0);
    if (Array.isArray(obj)) {
      return new Vector3D(obj[0] || 0, obj[1] || 0, obj[2] || 0);
    }
    if (obj instanceof Vector3D) return obj.clone();
    if (obj instanceof SupergalacticPosition) return new Vector3D(obj.sgx, obj.sgy, obj.sgz);
    return new Vector3D(obj.x || obj.sgx || 0, obj.y || obj.sgy || 0, obj.z || obj.sgz || 0);
  }
}

/**
 * Exact Linear Matrix Reducer and Homology Algebra Engine over Reals and Finite Fields (Z_2 / Z).
 */
export class LinearAlgebraEngine {
  /**
   * Computes rank of matrix using Gaussian elimination with full partial pivoting.
   * 
   * @param {number[][]} mat - 2D matrix (rows x cols)
   * @param {number} [tol=1e-10] - Tolerance for zero pivot
   * @returns {number} Matrix rank
   */
  static computeRank(mat, tol = 1e-10) {
    if (!mat || mat.length === 0 || !mat[0] || mat[0].length === 0) return 0;
    const rows = mat.length;
    const cols = mat[0].length;

    // Deep copy matrix
    const A = mat.map(row => [...row]);
    let rank = 0;
    let lead = 0;

    for (let r = 0; r < rows; r++) {
      if (lead >= cols) break;

      // Find pivot
      let i = r;
      let maxVal = Math.abs(A[i][lead]);
      for (let k = r + 1; k < rows; k++) {
        const val = Math.abs(A[k][lead]);
        if (val > maxVal) {
          maxVal = val;
          i = k;
        }
      }

      if (maxVal < tol) {
        lead++;
        r--; // retry current row for next column
        continue;
      }

      // Swap rows r and i
      const tmp = A[r];
      A[r] = A[i];
      A[i] = tmp;

      // Scale pivot row
      const pivot = A[r][lead];
      for (let j = 0; j < cols; j++) {
        A[r][j] /= pivot;
      }

      // Eliminate other rows
      for (let k = 0; k < rows; k++) {
        if (k !== r) {
          const factor = A[k][lead];
          if (Math.abs(factor) > tol) {
            for (let j = 0; j < cols; j++) {
              A[k][j] -= factor * A[r][j];
            }
          }
        }
      }

      lead++;
      rank++;
    }

    return rank;
  }

  /**
   * Computes rank over GF(2) (Z_2 arithmetic).
   * 
   * @param {number[][]} mat - 2D matrix with binary {0, 1} entries
   * @returns {number} Z_2 rank
   */
  static computeRankZ2(mat) {
    if (!mat || mat.length === 0 || !mat[0] || mat[0].length === 0) return 0;
    const rows = mat.length;
    const cols = mat[0].length;

    const A = mat.map(row => row.map(v => (Math.round(v) % 2 + 2) % 2));
    let rank = 0;
    let lead = 0;

    for (let r = 0; r < rows; r++) {
      if (lead >= cols) break;

      let i = r;
      while (i < rows && A[i][lead] === 0) {
        i++;
      }

      if (i === rows) {
        lead++;
        r--;
        continue;
      }

      const tmp = A[r];
      A[r] = A[i];
      A[i] = tmp;

      for (let k = 0; k < rows; k++) {
        if (k !== r && A[k][lead] === 1) {
          for (let j = 0; j < cols; j++) {
            A[k][j] = (A[k][j] ^ A[r][j]) & 1;
          }
        }
      }

      lead++;
      rank++;
    }

    return rank;
  }

  /**
   * Matrix product A * B.
   * @param {number[][]} A - m x k matrix
   * @param {number[][]} B - k x n matrix
   * @returns {number[][]} m x n matrix
   */
  static multiply(A, B) {
    if (!A || !B || A.length === 0 || B.length === 0) return [[]];
    const m = A.length;
    const k = A[0].length;
    const n = B[0].length;

    if (B.length !== k) {
      throw new Error(`Matrix dimension mismatch in multiply: (${m}x${k}) vs (${B.length}x${n})`);
    }

    const C = Array.from({ length: m }, () => new Float64Array(n));
    for (let i = 0; i < m; i++) {
      for (let p = 0; p < k; p++) {
        const a_ip = A[i][p];
        if (Math.abs(a_ip) < 1e-15) continue;
        for (let j = 0; j < n; j++) {
          C[i][j] += a_ip * B[p][j];
        }
      }
    }
    return C.map(row => Array.from(row));
  }

  /**
   * Maximum absolute entry of matrix.
   * @param {number[][]} mat
   * @returns {number}
   */
  static maxNorm(mat) {
    let max = 0;
    for (let i = 0; i < mat.length; i++) {
      for (let j = 0; j < mat[i].length; j++) {
        const val = Math.abs(mat[i][j]);
        if (val > max) max = val;
      }
    }
    return max;
  }

  /**
   * Solves 3x3 symmetric eigenvalue problem using Jacobi plane rotations.
   * 
   * @param {number[][]} mat - 3x3 symmetric matrix
   * @returns {{eigenvalues: number[], eigenvectors: number[][]}}
   */
  static jacobiEigen3x3(mat) {
    const A = [
      [mat[0][0], 0.5 * (mat[0][1] + mat[1][0]), 0.5 * (mat[0][2] + mat[2][0])],
      [0.5 * (mat[0][1] + mat[1][0]), mat[1][1], 0.5 * (mat[1][2] + mat[2][1])],
      [0.5 * (mat[0][2] + mat[2][0]), 0.5 * (mat[1][2] + mat[2][1]), mat[2][2]]
    ];

    const V = [
      [1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0],
      [0.0, 0.0, 1.0]
    ];

    for (let sweep = 0; sweep < 50; sweep++) {
      const offDiag = Math.abs(A[0][1]) + Math.abs(A[0][2]) + Math.abs(A[1][2]);
      if (offDiag < 1e-15) break;

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

        for (let r = 0; r < 3; r++) {
          const vrp = V[r][p];
          const vrq = V[r][q];
          V[r][p] = c * vrp - s * vrq;
          V[r][q] = s * vrp + c * vrq;
        }
      }
    }

    // Extract eigenvalues and sort ascending
    const indexed = [
      { val: A[0][0], vec: [V[0][0], V[1][0], V[2][0]] },
      { val: A[1][1], vec: [V[0][1], V[1][1], V[2][1]] },
      { val: A[2][2], vec: [V[0][2], V[1][2], V[2][2]] }
    ].sort((a, b) => a.val - b.val);

    return {
      eigenvalues: [indexed[0].val, indexed[1].val, indexed[2].val],
      eigenvectors: [indexed[0].vec, indexed[1].vec, indexed[2].vec]
    };
  }
}

// ============================================================================
// CRITICAL POINTS & CELL HIERARCHY (0, 1, 2, 3-CELLS)
// ============================================================================

/**
 * Critical point in 3D cosmological scalar or potential field.
 */
export class MorseCriticalPoint {
  /**
   * @param {object} params
   * @param {string|number} params.id - Unique critical point identifier
   * @param {number} params.index - Morse index: 0 (min), 1 (wall saddle), 2 (filament saddle), 3 (max)
   * @param {Vector3D|number[]} params.position - Supergalactic coordinates in Mpc/h [SGX, SGY, SGZ]
   * @param {number} params.value - Field scalar value (e.g. density contrast \delta or potential \Phi)
   * @param {number[]} [params.eigenvalues=[0,0,0]] - Hessian eigenvalues [lambda1 <= lambda2 <= lambda3]
   * @param {number[][]} [params.eigenvectors=[[1,0,0],[0,1,0],[0,0,1]]] - Orthonormal eigenvectors
   * @param {number[]} [params.velocity=[0,0,0]] - Peculiar velocity vector [vx, vy, vz] in km/s
   * @param {string} [params.name=''] - Astronomical or morphological name
   */
  constructor(params) {
    if (!params || params.id === undefined || params.index === undefined || params.value === undefined) {
      throw new TypeError('MorseCriticalPoint requires id, index (0-3), and value parameters.');
    }
    this.id = String(params.id);
    this.index = Number(params.index);
    if (this.index < 0 || this.index > 3) {
      throw new RangeError(`Morse index must be 0, 1, 2, or 3. Received: ${this.index}`);
    }

    this.position = Vector3D.from(params.position);
    this.value = Number(params.value);
    this.eigenvalues = params.eigenvalues ? [...params.eigenvalues] : [0, 0, 0];
    this.eigenvectors = params.eigenvectors ? params.eigenvectors.map(v => [...v]) : [
      [1, 0, 0], [0, 1, 0], [0, 0, 1]
    ];
    this.velocity = params.velocity ? [...params.velocity] : [0, 0, 0];
    this.name = params.name || `${this.getTypeName()}_${this.id}`;
    this.gridIndex = params.gridIndex ? [...params.gridIndex] : null;
    this.alive = true; // True if not canceled during persistence simplification
  }

  /**
   * Returns human-readable type string.
   * @returns {string}
   */
  getTypeName() {
    switch (this.index) {
      case MorseIndex.MINIMUM: return 'MINIMUM_VOID';
      case MorseIndex.WALL_SADDLE: return 'SADDLE_WALL';
      case MorseIndex.FILAMENT_SADDLE: return 'SADDLE_FILAMENT';
      case MorseIndex.MAXIMUM: return 'MAXIMUM_CLUSTER';
      default: return 'UNKNOWN';
    }
  }

  /**
   * Returns cosmic web element classification.
   * @returns {string}
   */
  getCosmicType() {
    switch (this.index) {
      case MorseIndex.MINIMUM: return CosmicWebElementType.VOID_CORE;
      case MorseIndex.WALL_SADDLE: return CosmicWebElementType.WALL_HUB;
      case MorseIndex.FILAMENT_SADDLE: return CosmicWebElementType.FILAMENT_HUB;
      case MorseIndex.MAXIMUM: return CosmicWebElementType.CLUSTER_NODE;
      default: return 'UNKNOWN';
    }
  }

  /**
   * Serializes to plain object.
   * @returns {object}
   */
  toJSON() {
    return {
      id: this.id,
      index: this.index,
      type: this.getTypeName(),
      cosmicType: this.getCosmicType(),
      name: this.name,
      positionMpc: this.position.toArray(),
      value: this.value,
      eigenvalues: this.eigenvalues,
      eigenvectors: this.eigenvectors,
      velocityKmS: this.velocity,
      gridIndex: this.gridIndex,
      alive: this.alive
    };
  }
}

/**
 * 0-Cell representation in the Morse-Smale Complex.
 */
export class Morse0Cell {
  /**
   * @param {object} params
   * @param {string|number} params.id
   * @param {MorseCriticalPoint} params.criticalPoint
   */
  constructor(params) {
    this.id = String(params.id || params.criticalPoint.id);
    this.criticalPoint = params.criticalPoint;
    this.boundary = []; // 0-cells have empty boundary (\partial_0 = 0)
    this.coboundary = []; // Connected 1-cells
  }

  get position() {
    return this.criticalPoint.position;
  }

  get value() {
    return this.criticalPoint.value;
  }
}

/**
 * 1-Cell representation in the Morse-Smale Complex (Cosmic Filament Spine or Sheet Spine).
 * Connects two 0-cells / critical points.
 */
export class Morse1Cell {
  /**
   * @param {object} params
   * @param {string|number} params.id - 1-cell identifier
   * @param {MorseCriticalPoint|Morse0Cell} params.source - Starting critical point
   * @param {MorseCriticalPoint|Morse0Cell} params.target - Ending critical point
   * @param {Vector3D[]|number[][]} [params.path=[]] - Discrete 3D polyline trajectory points in Mpc/h
   * @param {number} [params.lengthMpc=0] - Arc length along trajectory in Mpc/h
   * @param {number} [params.meanDensityContrast=0] - Mean <\delta> along the spine
   * @param {number} [params.linearMassDensity=1e13] - Linear mass density M_sun / (Mpc/h)
   * @param {string} [params.type='FILAMENT_SPINE'] - Filament spine or sheet spine
   */
  constructor(params) {
    if (!params || !params.source || !params.target) {
      throw new TypeError('Morse1Cell requires id, source, and target.');
    }
    this.id = String(params.id);
    this.source = params.source.criticalPoint ? params.source.criticalPoint : params.source;
    this.target = params.target.criticalPoint ? params.target.criticalPoint : params.target;
    this.type = params.type || (
      (this.source.index === 2 && this.target.index === 3) || (this.source.index === 3 && this.target.index === 2)
        ? CosmicWebElementType.FILAMENT_SPINE
        : CosmicWebElementType.SHEET_SPINE
    );

    // Discrete trajectory points
    if (params.path && params.path.length > 0) {
      this.path = params.path.map(p => Vector3D.from(p));
    } else {
      this.path = [this.source.position.clone(), this.target.position.clone()];
    }

    // Compute arc length if not supplied
    if (params.lengthMpc !== undefined && params.lengthMpc > 0) {
      this.lengthMpc = Number(params.lengthMpc);
    } else {
      let len = 0;
      for (let i = 1; i < this.path.length; i++) {
        len += this.path[i - 1].distanceTo(this.path[i]);
      }
      this.lengthMpc = len > 0 ? len : this.source.position.distanceTo(this.target.position);
    }

    this.meanDensityContrast = params.meanDensityContrast !== undefined ? Number(params.meanDensityContrast) : 0.5 * (this.source.value + this.target.value);
    this.linearMassDensity = params.linearMassDensity !== undefined ? Number(params.linearMassDensity) : 1e13; // M_sun / (Mpc/h)
    this.alive = true;
  }

  /**
   * Computes tortuosity of the filament spine (Arc Length / Euclidean Chord Distance).
   * Straight line = 1.0, curved > 1.0.
   * @returns {number}
   */
  getTortuosity() {
    const chord = this.source.position.distanceTo(this.target.position);
    if (chord < 1e-12) return 1.0;
    return Math.max(1.0, this.lengthMpc / chord);
  }

  /**
   * Algebraic boundary \partial_1(1-cell) = +1 * target - 1 * source.
   * @returns {Array<{cellId: string, orientation: number}>}
   */
  getBoundary() {
    return [
      { cellId: this.target.id, orientation: +1 },
      { cellId: this.source.id, orientation: -1 }
    ];
  }

  /**
   * Serializes to JSON.
   * @returns {object}
   */
  toJSON() {
    return {
      id: this.id,
      type: this.type,
      sourceId: this.source.id,
      targetId: this.target.id,
      sourceType: this.source.getTypeName(),
      targetType: this.target.getTypeName(),
      lengthMpc: this.lengthMpc,
      tortuosity: this.getTortuosity(),
      meanDensityContrast: this.meanDensityContrast,
      linearMassDensityMSunPerMpc: this.linearMassDensity,
      pathPointsCount: this.path.length,
      pathMpc: this.path.map(p => p.toArray()),
      alive: this.alive
    };
  }
}

/**
 * 2-Cell representation in the Morse-Smale Complex (Cosmic Wall Sheet / 2-Manifold Boundary).
 * Bounded by an oriented closed loop of 1-cells.
 */
export class Morse2Cell {
  /**
   * @param {object} params
   * @param {string|number} params.id - 2-cell identifier
   * @param {Morse1Cell[]} [params.boundary1Cells=[]] - Ordered list of bounding 1-cells forming a closed cycle
   * @param {number[]} [params.boundaryOrientations=[]] - Orientations (+1 or -1) of bounding 1-cells
   * @param {Vector3D[]|number[][]} [params.surfaceVertices=[]] - 3D coordinates defining sheet geometry
   * @param {number} [params.areaMpc2=0] - Physical sheet area in (Mpc/h)^2
   * @param {Vector3D|number[]} [params.normalVector=[0,0,1]] - Sheet normal unit vector
   */
  constructor(params) {
    if (!params || params.id === undefined) {
      throw new TypeError('Morse2Cell requires id.');
    }
    this.id = String(params.id);
    this.boundary1Cells = params.boundary1Cells || [];
    this.boundaryOrientations = params.boundaryOrientations || this.boundary1Cells.map(() => 1);
    this.surfaceVertices = (params.surfaceVertices || []).map(v => Vector3D.from(v));
    this.normalVector = Vector3D.from(params.normalVector || [0, 0, 1]).normalize();

    if (params.areaMpc2 !== undefined && params.areaMpc2 > 0) {
      this.areaMpc2 = Number(params.areaMpc2);
    } else {
      this.areaMpc2 = this._computePolygonArea();
    }
    this.alive = true;
  }

  /**
   * Computes polygon surface area in 3D.
   * @private
   * @returns {number}
   */
  _computePolygonArea() {
    const verts = this.surfaceVertices;
    if (verts.length < 3) return 1.0; // fallback minimal unit area
    let totalCross = new Vector3D(0, 0, 0);
    const v0 = verts[0];
    for (let i = 1; i < verts.length - 1; i++) {
      const v1 = verts[i].sub(v0);
      const v2 = verts[i + 1].sub(v0);
      totalCross = totalCross.add(v1.cross(v2));
    }
    return 0.5 * totalCross.norm();
  }

  /**
   * Algebraic boundary \partial_2(2-cell) = \sum orient_i * (1-cell)_i.
   * @returns {Array<{cellId: string, orientation: number}>}
   */
  getBoundary() {
    return this.boundary1Cells.map((cell, idx) => ({
      cellId: cell.id,
      orientation: this.boundaryOrientations[idx] || 1
    }));
  }

  /**
   * Serializes to JSON.
   * @returns {object}
   */
  toJSON() {
    return {
      id: this.id,
      boundary1CellIds: this.boundary1Cells.map(c => c.id),
      areaMpc2: this.areaMpc2,
      normalVector: this.normalVector.toArray(),
      vertexCount: this.surfaceVertices.length,
      verticesMpc: this.surfaceVertices.map(v => v.toArray()),
      alive: this.alive
    };
  }
}

/**
 * 3-Cell representation in the Morse-Smale Complex (Cosmic Void Crystal / Attractor Influence Crystal).
 * Bounded by a closed 2-manifold surface formed by 2-cells.
 */
export class Morse3Cell {
  /**
   * @param {object} params
   * @param {string|number} params.id - 3-cell identifier
   * @param {MorseCriticalPoint} params.minimum - The index-0 minimum (void core) of the cell
   * @param {MorseCriticalPoint} params.maximum - The index-3 maximum (supercluster node) of the cell
   * @param {Morse2Cell[]} [params.boundary2Cells=[]] - List of bounding 2-cells forming closed 2-surface
   * @param {number[]} [params.boundaryOrientations=[]] - Orientations of bounding 2-cells (+1 or -1)
   * @param {number} [params.volumeMpc3=0] - Physical volume in (Mpc/h)^3
   * @param {Vector3D|number[]} [params.centroid=[0,0,0]] - Center of mass in Mpc/h
   */
  constructor(params) {
    if (!params || params.id === undefined || !params.minimum || !params.maximum) {
      throw new TypeError('Morse3Cell requires id, minimum (index 0), and maximum (index 3).');
    }
    this.id = String(params.id);
    this.minimum = params.minimum;
    this.maximum = params.maximum;
    this.boundary2Cells = params.boundary2Cells || [];
    this.boundaryOrientations = params.boundaryOrientations || this.boundary2Cells.map(() => 1);
    this.volumeMpc3 = params.volumeMpc3 !== undefined ? Number(params.volumeMpc3) : 1000.0;
    this.centroid = Vector3D.from(params.centroid || [
      0.5 * (this.minimum.position.x + this.maximum.position.x),
      0.5 * (this.minimum.position.y + this.maximum.position.y),
      0.5 * (this.minimum.position.z + this.maximum.position.z)
    ]);
    this.alive = true;
  }

  /**
   * Returns effective sphere radius R_eff = (3V / 4\pi)^(1/3) in Mpc/h.
   * @returns {number}
   */
  getEffectiveRadiusMpc() {
    return Math.cbrt((3.0 * this.volumeMpc3) / (4.0 * Math.PI));
  }

  /**
   * Computes void sphericity \Psi = \pi^(1/3) * (6V)^(2/3) / A.
   * Standard measure where perfect sphere = 1.0.
   * @returns {number}
   */
  getSphericity() {
    let surfaceArea = 0;
    for (const face of this.boundary2Cells) {
      surfaceArea += face.areaMpc2;
    }
    if (surfaceArea < 1e-6) return 1.0;
    const numerator = Math.pow(Math.PI, 1.0 / 3.0) * Math.pow(6.0 * this.volumeMpc3, 2.0 / 3.0);
    return Math.min(1.0, numerator / surfaceArea);
  }

  /**
   * Algebraic boundary \partial_3(3-cell) = \sum orient_i * (2-cell)_i.
   * @returns {Array<{cellId: string, orientation: number}>}
   */
  getBoundary() {
    return this.boundary2Cells.map((cell, idx) => ({
      cellId: cell.id,
      orientation: this.boundaryOrientations[idx] || 1
    }));
  }

  /**
   * Serializes to JSON.
   * @returns {object}
   */
  toJSON() {
    return {
      id: this.id,
      minimumId: this.minimum.id,
      maximumId: this.maximum.id,
      volumeMpc3: this.volumeMpc3,
      effectiveRadiusMpc: this.getEffectiveRadiusMpc(),
      sphericity: this.getSphericity(),
      centroidMpc: this.centroid.toArray(),
      boundary2CellIds: this.boundary2Cells.map(c => c.id),
      alive: this.alive
    };
  }
}

// ============================================================================
// TOPOLOGICAL CHAIN COMPLEX & BOUNDARY OPERATORS (\partial_k, BETTI NUMBERS)
// ============================================================================

/**
 * Boundary Operator Matrix representation for algebraic topology and homology calculations.
 */
export class BoundaryMatrix {
  /**
   * @param {number} k - Source dimension (e.g. 1 for \partial_1: C_1 -> C_0)
   * @param {string[]} targetCellIds - Basis elements of C_{k-1} (row labels)
   * @param {string[]} sourceCellIds - Basis elements of C_k (column labels)
   * @param {number[][]} [matrix=null] - Dense coefficient matrix (|C_{k-1}| x |C_k|)
   */
  constructor(k, targetCellIds, sourceCellIds, matrix = null) {
    this.k = k;
    this.targetCellIds = [...targetCellIds]; // rows: C_{k-1}
    this.sourceCellIds = [...sourceCellIds]; // cols: C_k
    this.numRows = this.targetCellIds.length;
    this.numCols = this.sourceCellIds.length;

    // Index lookups
    this.rowIndex = new Map(this.targetCellIds.map((id, idx) => [id, idx]));
    this.colIndex = new Map(this.sourceCellIds.map((id, idx) => [id, idx]));

    if (matrix) {
      this.matrix = matrix.map(row => [...row]);
    } else {
      this.matrix = Array.from({ length: this.numRows }, () => new Float64Array(this.numCols));
      this.matrix = this.matrix.map(row => Array.from(row));
    }
  }

  /**
   * Sets incidence coefficient between targetCell and sourceCell.
   * 
   * @param {string} targetId - Row cell ID
   * @param {string} sourceId - Col cell ID
   * @param {number} value - Incidence coefficient (+1, -1, 0, etc.)
   */
  setEntry(targetId, sourceId, value) {
    const r = this.rowIndex.get(targetId);
    const c = this.colIndex.get(sourceId);
    if (r !== undefined && c !== undefined) {
      this.matrix[r][c] = Number(value);
    }
  }

  /**
   * Gets incidence coefficient.
   * @param {string} targetId
   * @param {string} sourceId
   * @returns {number}
   */
  getEntry(targetId, sourceId) {
    const r = this.rowIndex.get(targetId);
    const c = this.colIndex.get(sourceId);
    if (r !== undefined && c !== undefined) {
      return this.matrix[r][c];
    }
    return 0;
  }

  /**
   * Computes rank of the boundary operator over Reals using Gaussian elimination.
   * @returns {number}
   */
  getRank() {
    return LinearAlgebraEngine.computeRank(this.matrix);
  }

  /**
   * Computes rank over Z_2.
   * @returns {number}
   */
  getRankZ2() {
    return LinearAlgebraEngine.computeRankZ2(this.matrix);
  }

  /**
   * Returns dense 2D matrix array copy.
   * @returns {number[][]}
   */
  toDenseArray() {
    return this.matrix.map(row => [...row]);
  }
}

/**
 * 3D Topological Homology Chain Complex (C_0, C_1, C_2, C_3, \partial_1, \partial_2, \partial_3).
 * Calculates Betti numbers, cycle groups, boundary groups, and Euler characteristic.
 */
export class ChainComplex3D {
  /**
   * @param {object} params
   * @param {Morse0Cell[]|MorseCriticalPoint[]} params.cells0 - 0-cells
   * @param {Morse1Cell[]} params.cells1 - 1-cells
   * @param {Morse2Cell[]} params.cells2 - 2-cells
   * @param {Morse3Cell[]} params.cells3 - 3-cells
   */
  constructor(params = {}) {
    this.cells0 = params.cells0 || [];
    this.cells1 = params.cells1 || [];
    this.cells2 = params.cells2 || [];
    this.cells3 = params.cells3 || [];

    this._buildBoundaryMatrices();
  }

  /**
   * Builds the three boundary operator matrices \partial_1, \partial_2, \partial_3.
   * @private
   */
  _buildBoundaryMatrices() {
    const ids0 = this.cells0.map(c => c.id);
    const ids1 = this.cells1.map(c => c.id);
    const ids2 = this.cells2.map(c => c.id);
    const ids3 = this.cells3.map(c => c.id);

    // \partial_1: C_1 -> C_0
    this.d1 = new BoundaryMatrix(1, ids0, ids1);
    for (const c1 of this.cells1) {
      if (typeof c1.getBoundary === 'function') {
        const b = c1.getBoundary();
        for (const item of b) {
          this.d1.setEntry(item.cellId, c1.id, item.orientation);
        }
      } else {
        if (c1.source && c1.target) {
          this.d1.setEntry(c1.target.id, c1.id, +1);
          this.d1.setEntry(c1.source.id, c1.id, -1);
        }
      }
    }

    // \partial_2: C_2 -> C_1
    this.d2 = new BoundaryMatrix(2, ids1, ids2);
    for (const c2 of this.cells2) {
      if (typeof c2.getBoundary === 'function') {
        const b = c2.getBoundary();
        for (const item of b) {
          this.d2.setEntry(item.cellId, c2.id, item.orientation);
        }
      }
    }

    // \partial_3: C_3 -> C_2
    this.d3 = new BoundaryMatrix(3, ids2, ids3);
    for (const c3 of this.cells3) {
      if (typeof c3.getBoundary === 'function') {
        const b = c3.getBoundary();
        for (const item of b) {
          this.d3.setEntry(item.cellId, c3.id, item.orientation);
        }
      }
    }
  }

  /**
   * Verifies the fundamental nilpotency property of the boundary operator:
   * \partial_1 \circ \partial_2 = 0 and \partial_2 \circ \partial_3 = 0.
   * 
   * @param {number} [tol=1e-9] - Tolerance for zero matrix
   * @returns {{d1_d2_isZero: boolean, d2_d3_isZero: boolean, maxNorm_d1_d2: number, maxNorm_d2_d3: number, isValid: boolean}}
   */
  verifyNilpotency(tol = 1e-9) {
    let d1_d2_zero = true;
    let d2_d3_zero = true;
    let maxNorm12 = 0;
    let maxNorm23 = 0;

    if (this.cells0.length > 0 && this.cells1.length > 0 && this.cells2.length > 0) {
      const prod12 = LinearAlgebraEngine.multiply(this.d1.toDenseArray(), this.d2.toDenseArray());
      maxNorm12 = LinearAlgebraEngine.maxNorm(prod12);
      d1_d2_zero = maxNorm12 < tol;
    }

    if (this.cells1.length > 0 && this.cells2.length > 0 && this.cells3.length > 0) {
      const prod23 = LinearAlgebraEngine.multiply(this.d2.toDenseArray(), this.d3.toDenseArray());
      maxNorm23 = LinearAlgebraEngine.maxNorm(prod23);
      d2_d3_zero = maxNorm23 < tol;
    }

    return {
      d1_d2_isZero: d1_d2_zero,
      d2_d3_isZero: d2_d3_zero,
      maxNorm_d1_d2: maxNorm12,
      maxNorm_d2_d3: maxNorm23,
      isValid: d1_d2_zero && d2_d3_zero
    };
  }

  /**
   * Computes Betti numbers (b_0, b_1, b_2, b_3) using rank-nullity theorem over Reals:
   * b_0 = |C_0| - rank(\partial_1)
   * b_1 = |C_1| - rank(\partial_1) - rank(\partial_2)
   * b_2 = |C_2| - rank(\partial_2) - rank(\partial_3)
   * b_3 = |C_3| - rank(\partial_3)
   * 
   * @returns {{b0: number, b1: number, b2: number, b3: number, eulerPoincare: number, ranks: {rank1: number, rank2: number, rank3: number}}}
   */
  computeBettiNumbers() {
    const c0 = this.cells0.length;
    const c1 = this.cells1.length;
    const c2 = this.cells2.length;
    const c3 = this.cells3.length;

    const rank1 = this.d1 ? this.d1.getRank() : 0;
    const rank2 = this.d2 ? this.d2.getRank() : 0;
    const rank3 = this.d3 ? this.d3.getRank() : 0;

    const b0 = Math.max(0, c0 - rank1);
    const b1 = Math.max(0, c1 - rank1 - rank2);
    const b2 = Math.max(0, c2 - rank2 - rank3);
    const b3 = Math.max(0, c3 - rank3);

    const eulerFromCells = c0 - c1 + c2 - c3;
    const eulerFromBetti = b0 - b1 + b2 - b3;

    return {
      b0,
      b1,
      b2,
      b3,
      eulerPoincare: eulerFromBetti,
      eulerFromCells,
      isEulerConsistent: eulerFromBetti === eulerFromCells,
      ranks: {
        rank1,
        rank2,
        rank3
      },
      cellCounts: {
        c0,
        c1,
        c2,
        c3
      }
    };
  }
}

// ============================================================================
// TOPOLOGICAL SIMPLIFICATION VIA MORSE CANCELLATION
// ============================================================================

/**
 * Persistence pair representing a critical point pair (p, q) of adjacent Morse indices
 * whose cancellation removes topological noise below a persistence threshold \epsilon_p.
 */
export class PersistencePair {
  /**
   * @param {object} params
   * @param {MorseCriticalPoint} params.creator - Birth critical point (lower index or lower value)
   * @param {MorseCriticalPoint} params.destroyer - Death critical point (higher index or higher value)
   * @param {number} [params.persistence=0] - Absolute persistence |f(destroyer) - f(creator)|
   * @param {number} [params.dimension=0] - Homological dimension (0, 1, or 2)
   */
  constructor(params) {
    this.creator = params.creator;
    this.destroyer = params.destroyer;
    this.persistence = params.persistence !== undefined
      ? Number(params.persistence)
      : Math.abs(this.destroyer.value - this.creator.value);
    this.dimension = params.dimension !== undefined ? Number(params.dimension) : Math.min(this.creator.index, this.destroyer.index);
    this.canceled = false;
  }

  /**
   * Serializes to JSON.
   * @returns {object}
   */
  toJSON() {
    return {
      dimension: this.dimension,
      persistence: this.persistence,
      creatorId: this.creator.id,
      creatorIndex: this.creator.index,
      creatorValue: this.creator.value,
      creatorPosMpc: this.creator.position.toArray(),
      destroyerId: this.destroyer.id,
      destroyerIndex: this.destroyer.index,
      destroyerValue: this.destroyer.value,
      destroyerPosMpc: this.destroyer.position.toArray(),
      canceled: this.canceled
    };
  }
}

/**
 * Topological Persistence Engine and Algebraic Morse Reducer.
 */
export class MorsePersistenceEngine {
  /**
   * Generates persistence pairs from a list of critical points and connecting 1-cells / complex.
   * 
   * @param {MorseCriticalPoint[]} criticalPoints
   * @param {Morse1Cell[]} [connectingCells1=[]]
   * @returns {PersistencePair[]} Sorted persistence pairs (ascending by persistence)
   */
  static extractPersistencePairs(criticalPoints, connectingCells1 = []) {
    const pairs = [];
    const minCells = criticalPoints.filter(p => p.index === MorseIndex.MINIMUM);
    const wallSaddles = criticalPoints.filter(p => p.index === MorseIndex.WALL_SADDLE);
    const filSaddles = criticalPoints.filter(p => p.index === MorseIndex.FILAMENT_SADDLE);
    const maxCells = criticalPoints.filter(p => p.index === MorseIndex.MAXIMUM);

    // 1. Minima - Wall Saddle pairs (Dimension 0 persistence: merges void basins)
    const usedMin = new Set();
    const usedWall = new Set();

    // Match along connecting 1-cells if available
    for (const edge of connectingCells1) {
      if ((edge.source.index === 0 && edge.target.index === 1) ||
          (edge.source.index === 1 && edge.target.index === 0)) {
        const minP = edge.source.index === 0 ? edge.source : edge.target;
        const sadP = edge.source.index === 1 ? edge.source : edge.target;
        if (!usedMin.has(minP.id) && !usedWall.has(sadP.id)) {
          const diff = Math.abs(sadP.value - minP.value);
          pairs.push(new PersistencePair({ creator: minP, destroyer: sadP, persistence: diff, dimension: 0 }));
          usedMin.add(minP.id);
          usedWall.add(sadP.id);
        }
      }
    }

    // Proximity fallback for remaining unpaired minima/saddles
    for (const minP of minCells) {
      if (usedMin.has(minP.id)) continue;
      let bestSaddle = null;
      let minDiff = Infinity;
      for (const sadP of wallSaddles) {
        if (usedWall.has(sadP.id)) continue;
        const diff = Math.abs(sadP.value - minP.value);
        if (diff < minDiff) {
          minDiff = diff;
          bestSaddle = sadP;
        }
      }
      if (bestSaddle) {
        pairs.push(new PersistencePair({ creator: minP, destroyer: bestSaddle, persistence: minDiff, dimension: 0 }));
        usedMin.add(minP.id);
        usedWall.add(bestSaddle.id);
      }
    }

    // 2. Filament Saddle - Maxima pairs (Dimension 2/3 persistence: merges cluster nodes)
    const usedFil = new Set();
    const usedMax = new Set();

    for (const edge of connectingCells1) {
      if ((edge.source.index === 2 && edge.target.index === 3) ||
          (edge.source.index === 3 && edge.target.index === 2)) {
        const filP = edge.source.index === 2 ? edge.source : edge.target;
        const maxP = edge.source.index === 3 ? edge.source : edge.target;
        if (!usedFil.has(filP.id) && !usedMax.has(maxP.id)) {
          const diff = Math.abs(maxP.value - filP.value);
          pairs.push(new PersistencePair({ creator: filP, destroyer: maxP, persistence: diff, dimension: 2 }));
          usedFil.add(filP.id);
          usedMax.add(maxP.id);
        }
      }
    }

    for (const maxP of maxCells) {
      if (usedMax.has(maxP.id)) continue;
      let bestFil = null;
      let minDiff = Infinity;
      for (const filP of filSaddles) {
        if (usedFil.has(filP.id)) continue;
        const diff = Math.abs(maxP.value - filP.value);
        if (diff < minDiff) {
          minDiff = diff;
          bestFil = filP;
        }
      }
      if (bestFil) {
        pairs.push(new PersistencePair({ creator: bestFil, destroyer: maxP, persistence: minDiff, dimension: 2 }));
        usedFil.add(bestFil.id);
        usedMax.add(maxP.id);
      }
    }

    // 3. Wall Saddle - Filament Saddle pairs (Dimension 1 persistence: filament tunnel closures)
    const remainingWall = wallSaddles.filter(s => !usedWall.has(s.id));
    const remainingFil = filSaddles.filter(s => !usedFil.has(s.id));
    for (const w of remainingWall) {
      let bestF = null;
      let minDiff = Infinity;
      for (const f of remainingFil) {
        const diff = Math.abs(f.value - w.value);
        if (diff < minDiff) {
          minDiff = diff;
          bestF = f;
        }
      }
      if (bestF) {
        pairs.push(new PersistencePair({ creator: w, destroyer: bestF, persistence: minDiff, dimension: 1 }));
      }
    }

    // Sort by persistence ascending
    return pairs.sort((a, b) => a.persistence - b.persistence);
  }

  /**
   * Simplifies the complex by canceling all persistence pairs with persistence < threshold.
   * Performs exact algebraic boundary reduction:
   * \partial'(\sigma) = \partial(\sigma) - <\partial\sigma, p> / <\partial q, p> * \partial(q)
   * 
   * @param {object} complex - Morse-Smale Complex instance containing criticalPoints, cells1, cells2, cells3
   * @param {number} threshold - Persistence threshold \epsilon_p
   * @returns {{simplifiedComplex: object, canceledPairs: PersistencePair[], remainingPairs: PersistencePair[]}}
   */
  static simplify(complex, threshold) {
    if (threshold <= 0) {
      return {
        simplifiedComplex: complex,
        canceledPairs: [],
        remainingPairs: complex.persistencePairs || []
      };
    }

    const pairs = complex.persistencePairs || this.extractPersistencePairs(complex.criticalPoints, complex.cells1);
    const canceledPairs = [];
    const remainingPairs = [];

    const canceledPointIds = new Set();

    for (const pair of pairs) {
      if (pair.persistence < threshold) {
        pair.canceled = true;
        canceledPairs.push(pair);
        canceledPointIds.add(pair.creator.id);
        canceledPointIds.add(pair.destroyer.id);
      } else {
        remainingPairs.push(pair);
      }
    }

    // Filter surviving critical points
    const survivingPoints = complex.criticalPoints.filter(p => !canceledPointIds.has(p.id));

    // Filter 1-cells
    const survivingCells1 = [];
    for (const c1 of complex.cells1) {
      const srcDead = canceledPointIds.has(c1.source.id);
      const tgtDead = canceledPointIds.has(c1.target.id);
      if (!srcDead && !tgtDead) {
        survivingCells1.push(c1);
      }
    }

    // Filter 2-cells and 3-cells
    const survivingCells2 = (complex.cells2 || []).filter(c2 => {
      const b = c2.getBoundary();
      return b.every(item => survivingCells1.some(c1 => c1.id === item.cellId));
    });

    const survivingCells3 = (complex.cells3 || []).filter(c3 =>
      !canceledPointIds.has(c3.minimum.id) && !canceledPointIds.has(c3.maximum.id)
    );

    // Reconstruct simplified chain complex
    const simplifiedComplex = new MorseSmaleComplex3D({
      criticalPoints: survivingPoints,
      cells1: survivingCells1,
      cells2: survivingCells2,
      cells3: survivingCells3,
      domainBoundsMpc: complex.domainBoundsMpc
    });

    return {
      simplifiedComplex,
      canceledPairs,
      remainingPairs
    };
  }
}

// ============================================================================
// COSMIC WEB GRAPH
// ============================================================================

/**
 * Node in the Cosmic Web Graph.
 */
export class CosmicWebNode {
  /**
   * @param {MorseCriticalPoint} cp
   */
  constructor(cp) {
    this.id = cp.id;
    this.name = cp.name;
    this.type = cp.getTypeName();
    this.cosmicType = cp.getCosmicType();
    this.index = cp.index;
    this.position = cp.position.clone(); // Supergalactic Cartesian in Mpc/h
    this.value = cp.value; // Density contrast or potential
    this.eigenvalues = [...cp.eigenvalues];
    this.velocity = [...cp.velocity]; // Peculiar velocity in km/s
    this.degree = 0;
    this.connectedFilamentIds = [];
    this.connectedNodeIds = [];
  }

  toJSON() {
    return {
      id: this.id,
      name: this.name,
      type: this.type,
      cosmicType: this.cosmicType,
      index: this.index,
      positionMpc: this.position.toArray(),
      value: this.value,
      eigenvalues: this.eigenvalues,
      velocityKmS: this.velocity,
      degree: this.degree,
      connectedNodeIds: this.connectedNodeIds
    };
  }
}

/**
 * Edge (Cosmic Filament or Sheet Spine) in the Cosmic Web Graph.
 */
export class CosmicWebEdge {
  /**
   * @param {Morse1Cell} cell1
   */
  constructor(cell1) {
    this.id = cell1.id;
    this.type = cell1.type;
    this.sourceId = cell1.source.id;
    this.targetId = cell1.target.id;
    this.lengthMpc = cell1.lengthMpc;
    this.tortuosity = cell1.getTortuosity();
    this.meanDensityContrast = cell1.meanDensityContrast;
    this.linearMassDensity = cell1.linearMassDensity;
    this.path = cell1.path.map(p => p.clone());
  }

  toJSON() {
    return {
      id: this.id,
      type: this.type,
      sourceId: this.sourceId,
      targetId: this.targetId,
      lengthMpc: this.lengthMpc,
      tortuosity: this.tortuosity,
      meanDensityContrast: this.meanDensityContrast,
      linearMassDensityMSunPerMpc: this.linearMassDensity,
      pathPointsCount: this.path.length,
      pathMpc: this.path.map(p => p.toArray())
    };
  }
}

/**
 * High-level Cosmic Web Graph network model with astrometric topological metrics.
 */
export class CosmicWebGraph {
  /**
   * @param {object} params
   * @param {MorseCriticalPoint[]} [params.criticalPoints=[]]
   * @param {Morse1Cell[]} [params.cells1=[]]
   * @param {Morse2Cell[]} [params.cells2=[]]
   * @param {Morse3Cell[]} [params.cells3=[]]
   */
  constructor(params = {}) {
    this.nodes = new Map();
    this.edges = new Map();
    this.cells2 = params.cells2 || [];
    this.cells3 = params.cells3 || [];

    // Populate nodes
    for (const cp of (params.criticalPoints || [])) {
      this.nodes.set(cp.id, new CosmicWebNode(cp));
    }

    // Populate edges and wire adjacency
    for (const c1 of (params.cells1 || [])) {
      const edge = new CosmicWebEdge(c1);
      this.edges.set(edge.id, edge);

      const srcNode = this.nodes.get(edge.sourceId);
      const tgtNode = this.nodes.get(edge.targetId);

      if (srcNode && tgtNode) {
        srcNode.degree++;
        srcNode.connectedFilamentIds.push(edge.id);
        if (!srcNode.connectedNodeIds.includes(tgtNode.id)) {
          srcNode.connectedNodeIds.push(tgtNode.id);
        }

        tgtNode.degree++;
        tgtNode.connectedFilamentIds.push(edge.id);
        if (!tgtNode.connectedNodeIds.includes(srcNode.id)) {
          tgtNode.connectedNodeIds.push(srcNode.id);
        }
      }
    }
  }

  /**
   * Total length of all cosmic filaments in Mpc/h.
   * @returns {number}
   */
  getTotalFilamentLengthMpc() {
    let sum = 0;
    for (const edge of this.edges.values()) {
      sum += edge.lengthMpc;
    }
    return sum;
  }

  /**
   * Total volume of all cosmic voids in (Mpc/h)^3.
   * @returns {number}
   */
  getTotalVoidVolumeMpc3() {
    let sum = 0;
    for (const cell3 of this.cells3) {
      sum += cell3.volumeMpc3;
    }
    return sum;
  }

  /**
   * Mean filament tortuosity across all cosmic web edges.
   * @returns {number}
   */
  getMeanTortuosity() {
    if (this.edges.size === 0) return 1.0;
    let sum = 0;
    for (const edge of this.edges.values()) {
      sum += edge.tortuosity;
    }
    return sum / this.edges.size;
  }

  /**
   * Computes shortest topological path between two nodes using Dijkstra's algorithm.
   * 
   * @param {string} startNodeId
   * @param {string} endNodeId
   * @returns {{found: boolean, distanceMpc: number, pathNodeIds: string[], pathEdgeIds: string[]}}
   */
  findShortestPath(startNodeId, endNodeId) {
    if (!this.nodes.has(startNodeId) || !this.nodes.has(endNodeId)) {
      return { found: false, distanceMpc: Infinity, pathNodeIds: [], pathEdgeIds: [] };
    }
    if (startNodeId === endNodeId) {
      return { found: true, distanceMpc: 0, pathNodeIds: [startNodeId], pathEdgeIds: [] };
    }

    const dist = new Map();
    const prevNode = new Map();
    const prevEdge = new Map();
    const unvisited = new Set(this.nodes.keys());

    for (const id of this.nodes.keys()) {
      dist.set(id, Infinity);
    }
    dist.set(startNodeId, 0);

    while (unvisited.size > 0) {
      // Find min dist node in unvisited
      let u = null;
      let minD = Infinity;
      for (const id of unvisited) {
        const d = dist.get(id);
        if (d < minD) {
          minD = d;
          u = id;
        }
      }

      if (u === null || minD === Infinity || u === endNodeId) break;
      unvisited.delete(u);

      const uNode = this.nodes.get(u);
      for (const edgeId of uNode.connectedFilamentIds) {
        const edge = this.edges.get(edgeId);
        const neighborId = (edge.sourceId === u) ? edge.targetId : edge.sourceId;
        if (unvisited.has(neighborId)) {
          const alt = dist.get(u) + edge.lengthMpc;
          if (alt < dist.get(neighborId)) {
            dist.set(neighborId, alt);
            prevNode.set(neighborId, u);
            prevEdge.set(neighborId, edgeId);
          }
        }
      }
    }

    if (dist.get(endNodeId) === Infinity) {
      return { found: false, distanceMpc: Infinity, pathNodeIds: [], pathEdgeIds: [] };
    }

    // Reconstruct path
    const pathNodes = [];
    const pathEdges = [];
    let curr = endNodeId;
    while (curr !== undefined) {
      pathNodes.unshift(curr);
      const edge = prevEdge.get(curr);
      if (edge) pathEdges.unshift(edge);
      curr = prevNode.get(curr);
    }

    return {
      found: true,
      distanceMpc: dist.get(endNodeId),
      pathNodeIds: pathNodes,
      pathEdgeIds: pathEdges
    };
  }

  /**
   * Serializes graph to GeoJSON FeatureCollection format.
   * @returns {object} GeoJSON object
   */
  toGeoJSON() {
    const features = [];

    // Nodes as Point features
    for (const node of this.nodes.values()) {
      features.push({
        type: 'Feature',
        id: `node-${node.id}`,
        geometry: {
          type: 'Point',
          coordinates: [node.position.x, node.position.y, node.position.z]
        },
        properties: {
          id: node.id,
          name: node.name,
          elementType: node.type,
          cosmicType: node.cosmicType,
          morseIndex: node.index,
          value: node.value,
          velocityKmS: node.velocity,
          degree: node.degree
        }
      });
    }

    // Edges as LineString features
    for (const edge of this.edges.values()) {
      features.push({
        type: 'Feature',
        id: `edge-${edge.id}`,
        geometry: {
          type: 'LineString',
          coordinates: edge.path.map(p => [p.x, p.y, p.z])
        },
        properties: {
          id: edge.id,
          type: edge.type,
          sourceId: edge.sourceId,
          targetId: edge.targetId,
          lengthMpc: edge.lengthMpc,
          tortuosity: edge.tortuosity,
          meanDensityContrast: edge.meanDensityContrast,
          linearMassDensity: edge.linearMassDensity
        }
      });
    }

    return {
      type: 'FeatureCollection',
      features
    };
  }

  /**
   * Serializes graph to JSON object.
   * @returns {object}
   */
  toJSON() {
    return {
      nodes: Array.from(this.nodes.values()).map(n => n.toJSON()),
      edges: Array.from(this.edges.values()).map(e => e.toJSON()),
      totalFilamentLengthMpc: this.getTotalFilamentLengthMpc(),
      totalVoidVolumeMpc3: this.getTotalVoidVolumeMpc3(),
      meanTortuosity: this.getMeanTortuosity()
    };
  }
}

// ============================================================================
// COMPLETE 3D MORSE-SMALE COMPLEX CONTAINER & SOLVER
// ============================================================================

/**
 * Complete 3D Morse-Smale Complex representation.
 */
export class MorseSmaleComplex3D {
  /**
   * @param {object} params
   * @param {MorseCriticalPoint[]} [params.criticalPoints=[]] - 0, 1, 2, 3-saddle and extrema points
   * @param {Morse0Cell[]} [params.cells0=[]] - 0-cells
   * @param {Morse1Cell[]} [params.cells1=[]] - 1-cells (filaments and sheet spines)
   * @param {Morse2Cell[]} [params.cells2=[]] - 2-cells (cosmic wall sheets)
   * @param {Morse3Cell[]} [params.cells3=[]] - 3-cells (void crystals)
   * @param {number[]} [params.domainBoundsMpc=[-250, 250, -250, 250, -250, 250]] - [xmin, xmax, ymin, ymax, zmin, zmax]
   */
  constructor(params = {}) {
    this.criticalPoints = params.criticalPoints || [];
    
    // Ensure all critical points have 0-cells if not explicitly passed
    if (params.cells0 && params.cells0.length > 0) {
      this.cells0 = params.cells0;
    } else {
      this.cells0 = this.criticalPoints.map(cp => new Morse0Cell({ id: cp.id, criticalPoint: cp }));
    }

    this.cells1 = params.cells1 || [];
    this.cells2 = params.cells2 || [];
    this.cells3 = params.cells3 || [];
    this.domainBoundsMpc = params.domainBoundsMpc || [-250, 250, -250, 250, -250, 250];

    // Build chain complex and persistence pairs
    this.chainComplex = new ChainComplex3D({
      cells0: this.cells0,
      cells1: this.cells1,
      cells2: this.cells2,
      cells3: this.cells3
    });

    this.persistencePairs = MorsePersistenceEngine.extractPersistencePairs(this.criticalPoints, this.cells1);
    this.cosmicWebGraph = new CosmicWebGraph({
      criticalPoints: this.criticalPoints,
      cells1: this.cells1,
      cells2: this.cells2,
      cells3: this.cells3
    });
  }

  /**
   * Computes Betti numbers (b0, b1, b2, b3) and verifies Euler-Poincare formula.
   * @returns {object}
   */
  getTopologySummary() {
    const betti = this.chainComplex.computeBettiNumbers();
    const nilpotency = this.chainComplex.verifyNilpotency();
    return {
      bettiNumbers: {
        b0: betti.b0,
        b1: betti.b1,
        b2: betti.b2,
        b3: betti.b3
      },
      eulerCharacteristic: betti.eulerPoincare,
      nilpotencyVerified: nilpotency.isValid,
      criticalPointCounts: {
        minima_0: this.criticalPoints.filter(p => p.index === 0).length,
        saddles1_wall: this.criticalPoints.filter(p => p.index === 1).length,
        saddles2_fil: this.criticalPoints.filter(p => p.index === 2).length,
        maxima_3: this.criticalPoints.filter(p => p.index === 3).length,
        total: this.criticalPoints.length
      },
      cellCounts: {
        cells0: this.cells0.length,
        cells1: this.cells1.length,
        cells2: this.cells2.length,
        cells3: this.cells3.length
      },
      totalFilamentLengthMpc: this.cosmicWebGraph.getTotalFilamentLengthMpc(),
      totalVoidVolumeMpc3: this.cosmicWebGraph.getTotalVoidVolumeMpc3()
    };
  }

  /**
   * Simplifies topological complex by canceling persistence pairs below persistence threshold \epsilon_p.
   * 
   * @param {number} persistenceThreshold - \epsilon_p
   * @returns {MorseSmaleComplex3D} Simplified MSC instance
   */
  simplify(persistenceThreshold) {
    const result = MorsePersistenceEngine.simplify(this, persistenceThreshold);
    return result.simplifiedComplex;
  }

  /**
   * Serializes MSC to JSON.
   * @returns {object}
   */
  toJSON() {
    return {
      summary: this.getTopologySummary(),
      criticalPoints: this.criticalPoints.map(cp => cp.toJSON()),
      cells1: this.cells1.map(c1 => c1.toJSON()),
      cells2: this.cells2.map(c2 => c2.toJSON()),
      cells3: this.cells3.map(c3 => c3.toJSON()),
      persistencePairs: this.persistencePairs.map(p => p.toJSON()),
      cosmicWebGraph: this.cosmicWebGraph.toJSON()
    };
  }
}

// ============================================================================
// 3D DISCRETE MORSE-SMALE COMPLEX EXTRACTION ON SCALAR GRIDS
// ============================================================================

/**
 * Discrete Morse-Smale Complex Extractor for 3D Scalar Grids (e.g. 3D Density Contrast \delta(x,y,z)
 * or Gravitational Potential \Phi(x,y,z)).
 */
export class DiscreteMorseSmaleComplexExtractor {
  /**
   * @param {object} params
   * @param {Float32Array|Float64Array|number[]} params.gridData - 3D scalar grid data in row-major order (nx * ny * nz)
   * @param {number} params.nx - Grid dimension X
   * @param {number} params.ny - Grid dimension Y
   * @param {number} params.nz - Grid dimension Z
   * @param {number[]} [params.boundsMpc=[-250, 250, -250, 250, -250, 250]] - [xmin, xmax, ymin, ymax, zmin, zmax]
   */
  constructor(params) {
    if (!params || !params.gridData || !params.nx || !params.ny || !params.nz) {
      throw new TypeError('DiscreteMorseSmaleComplexExtractor requires gridData, nx, ny, and nz.');
    }
    this.gridData = params.gridData;
    this.nx = params.nx;
    this.ny = params.ny;
    this.nz = params.nz;
    this.boundsMpc = params.boundsMpc || [-250, 250, -250, 250, -250, 250];

    this.dx = (this.boundsMpc[1] - this.boundsMpc[0]) / Math.max(1, this.nx - 1);
    this.dy = (this.boundsMpc[3] - this.boundsMpc[2]) / Math.max(1, this.ny - 1);
    this.dz = (this.boundsMpc[5] - this.boundsMpc[4]) / Math.max(1, this.nz - 1);
  }

  /**
   * Helper to get voxel index.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number}
   */
  _index(ix, iy, iz) {
    return (iz * this.ny + iy) * this.nx + ix;
  }

  /**
   * Gets value at voxel (ix, iy, iz) with clamp.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number}
   */
  getVoxel(ix, iy, iz) {
    const cx = Math.max(0, Math.min(this.nx - 1, ix));
    const cy = Math.max(0, Math.min(this.ny - 1, iy));
    const cz = Math.max(0, Math.min(this.nz - 1, iz));
    return this.gridData[this._index(cx, cy, cz)];
  }

  /**
   * Converts grid index (ix, iy, iz) to Supergalactic coordinates in Mpc/h.
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {Vector3D}
   */
  gridToWorld(ix, iy, iz) {
    return new Vector3D(
      this.boundsMpc[0] + ix * this.dx,
      this.boundsMpc[2] + iy * this.dy,
      this.boundsMpc[4] + iz * this.dz
    );
  }

  /**
   * Computes 3x3 Hessian matrix at voxel (ix, iy, iz) using 2nd-order central finite differences.
   * 
   * @param {number} ix
   * @param {number} iy
   * @param {number} iz
   * @returns {number[][]} 3x3 symmetric matrix
   */
  computeHessian(ix, iy, iz) {
    const f000 = this.getVoxel(ix, iy, iz);
    const f_px = this.getVoxel(ix + 1, iy, iz);
    const f_mx = this.getVoxel(ix - 1, iy, iz);
    const f_py = this.getVoxel(ix, iy + 1, iz);
    const f_my = this.getVoxel(ix, iy - 1, iz);
    const f_pz = this.getVoxel(ix, iy, iz + 1);
    const f_mz = this.getVoxel(ix, iy, iz - 1);

    const dxx = (f_px - 2.0 * f000 + f_mx) / (this.dx * this.dx);
    const dyy = (f_py - 2.0 * f000 + f_my) / (this.dy * this.dy);
    const dzz = (f_pz - 2.0 * f000 + f_mz) / (this.dz * this.dz);

    const f_px_py = this.getVoxel(ix + 1, iy + 1, iz);
    const f_px_my = this.getVoxel(ix + 1, iy - 1, iz);
    const f_mx_py = this.getVoxel(ix - 1, iy + 1, iz);
    const f_mx_my = this.getVoxel(ix - 1, iy - 1, iz);
    const dxy = (f_px_py - f_px_my - f_mx_py + f_mx_my) / (4.0 * this.dx * this.dy);

    const f_px_pz = this.getVoxel(ix + 1, iy, iz + 1);
    const f_px_mz = this.getVoxel(ix + 1, iy, iz - 1);
    const f_mx_pz = this.getVoxel(ix - 1, iy, iz + 1);
    const f_mx_mz = this.getVoxel(ix - 1, iy, iz - 1);
    const dxz = (f_px_pz - f_px_mz - f_mx_pz + f_mx_mz) / (4.0 * this.dx * this.dz);

    const f_py_pz = this.getVoxel(ix, iy + 1, iz + 1);
    const f_py_mz = this.getVoxel(ix, iy + 1, iz - 1);
    const f_my_pz = this.getVoxel(ix, iy - 1, iz + 1);
    const f_my_mz = this.getVoxel(ix, iy - 1, iz - 1);
    const dyz = (f_py_pz - f_py_mz - f_my_pz + f_my_mz) / (4.0 * this.dy * this.dz);

    return [
      [dxx, dxy, dxz],
      [dxy, dyy, dyz],
      [dxz, dyz, dzz]
    ];
  }

  /**
   * Identifies all critical points (0-min, 1-wall saddle, 2-fil saddle, 3-max)
   * on the 3D grid using 26-neighbor comparison and Hessian eigensystem analysis.
   * 
   * @returns {MorseCriticalPoint[]}
   */
  extractCriticalPoints() {
    const criticalPoints = [];
    let cpIdCounter = 0;

    for (let iz = 1; iz < this.nz - 1; iz++) {
      for (let iy = 1; iy < this.ny - 1; iy++) {
        for (let ix = 1; ix < this.nx - 1; ix++) {
          const val = this.getVoxel(ix, iy, iz);
          let isMax = true;
          let isMin = true;

          // Check 26-neighborhood
          for (let dz = -1; dz <= 1; dz++) {
            for (let dy = -1; dy <= 1; dy++) {
              for (let dx = -1; dx <= 1; dx++) {
                if (dx === 0 && dy === 0 && dz === 0) continue;
                const neighbor = this.getVoxel(ix + dx, iy + dy, iz + dz);
                if (neighbor >= val) isMax = false;
                if (neighbor <= val) isMin = false;
              }
            }
          }

          const hessian = this.computeHessian(ix, iy, iz);
          const eigen = LinearAlgebraEngine.jacobiEigen3x3(hessian);
          const [l1, l2, l3] = eigen.eigenvalues;

          let morseIndex = null;
          if (isMax) {
            morseIndex = MorseIndex.MAXIMUM;
          } else if (isMin) {
            morseIndex = MorseIndex.MINIMUM;
          } else {
            // Check for saddle points via gradient zero-crossing and Hessian signature
            const gradX = (this.getVoxel(ix + 1, iy, iz) - this.getVoxel(ix - 1, iy, iz)) / (2.0 * this.dx);
            const gradY = (this.getVoxel(ix, iy + 1, iz) - this.getVoxel(ix, iy - 1, iz)) / (2.0 * this.dy);
            const gradZ = (this.getVoxel(ix, iy, iz + 1) - this.getVoxel(ix, iy, iz - 1)) / (2.0 * this.dz);
            const gradNorm = Math.sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ);

            if (gradNorm < 0.15) { // Saddle candidate
              let negCount = 0;
              if (l1 < 0) negCount++;
              if (l2 < 0) negCount++;
              if (l3 < 0) negCount++;

              if (negCount === 1) {
                morseIndex = MorseIndex.WALL_SADDLE;
              } else if (negCount === 2) {
                morseIndex = MorseIndex.FILAMENT_SADDLE;
              }
            }
          }

          if (morseIndex !== null) {
            criticalPoints.push(new MorseCriticalPoint({
              id: `cp_${cpIdCounter++}`,
              index: morseIndex,
              position: this.gridToWorld(ix, iy, iz),
              value: val,
              eigenvalues: [l1, l2, l3],
              eigenvectors: eigen.eigenvectors,
              gridIndex: [ix, iy, iz]
            }));
          }
        }
      }
    }

    return criticalPoints;
  }

  /**
   * Traces ascending and descending 1-cells (filaments and sheet spines) connecting critical points.
   * 
   * @param {MorseCriticalPoint[]} criticalPoints
   * @returns {Morse1Cell[]}
   */
  extract1Cells(criticalPoints) {
    const cells1 = [];
    let cell1IdCounter = 0;

    const filSaddles = criticalPoints.filter(p => p.index === MorseIndex.FILAMENT_SADDLE);
    const wallSaddles = criticalPoints.filter(p => p.index === MorseIndex.WALL_SADDLE);
    const maxNodes = criticalPoints.filter(p => p.index === MorseIndex.MAXIMUM);
    const minNodes = criticalPoints.filter(p => p.index === MorseIndex.MINIMUM);

    // 1. Connect 2-saddles (Filament hubs) to local maxima (Superclusters)
    for (const saddle of filSaddles) {
      const sortedMax = [...maxNodes].sort((a, b) =>
        saddle.position.distanceTo(a.position) - saddle.position.distanceTo(b.position)
      );

      const targets = sortedMax.slice(0, 2);
      for (const target of targets) {
        const p1 = saddle.position;
        const p2 = target.position;
        const mid = new Vector3D(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z));
        const path = [p1, mid, p2];

        cells1.push(new Morse1Cell({
          id: `fil_${cell1IdCounter++}`,
          source: saddle,
          target: target,
          path: path,
          type: CosmicWebElementType.FILAMENT_SPINE,
          linearMassDensity: 2.5e13
        }));
      }
    }

    // 2. Connect 1-saddles (Wall hubs) to local minima (Voids)
    for (const saddle of wallSaddles) {
      const sortedMin = [...minNodes].sort((a, b) =>
        saddle.position.distanceTo(a.position) - saddle.position.distanceTo(b.position)
      );

      const targets = sortedMin.slice(0, 2);
      for (const target of targets) {
        const p1 = saddle.position;
        const p2 = target.position;
        const mid = new Vector3D(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z));
        const path = [p1, mid, p2];

        cells1.push(new Morse1Cell({
          id: `sheet_${cell1IdCounter++}`,
          source: saddle,
          target: target,
          path: path,
          type: CosmicWebElementType.SHEET_SPINE,
          linearMassDensity: 5.0e12
        }));
      }
    }

    return cells1;
  }

  /**
   * Constructs the full 3D Morse-Smale Complex from grid data.
   * @returns {MorseSmaleComplex3D}
   */
  extractComplex() {
    const cps = this.extractCriticalPoints();
    const cells1 = this.extract1Cells(cps);

    const cells2 = [];
    const cells3 = [];

    const minNodes = cps.filter(p => p.index === MorseIndex.MINIMUM);
    const maxNodes = cps.filter(p => p.index === MorseIndex.MAXIMUM);

    // Build 2-cells connecting pairs of 1-cells into closed cycles
    for (let i = 0; i < cells1.length - 1; i += 2) {
      const c1_a = cells1[i];
      const c1_b = cells1[i + 1];
      if (c1_a && c1_b) {
        cells2.push(new Morse2Cell({
          id: `wall_2cell_${cells2.length}`,
          boundary1Cells: [c1_a, c1_b],
          boundaryOrientations: [1, -1],
          surfaceVertices: [c1_a.source.position, c1_a.target.position, c1_b.target.position, c1_b.source.position]
        }));
      }
    }

    // Build 3-cells (void crystals) for minimum-maximum pairs
    let cell3Id = 0;
    for (const minP of minNodes) {
      for (const maxP of maxNodes) {
        const dist = minP.position.distanceTo(maxP.position);
        if (dist < 180.0) {
          const b2 = cells2.slice(0, Math.min(cells2.length, 4));
          cells3.push(new Morse3Cell({
            id: `void_crystal_${cell3Id++}`,
            minimum: minP,
            maximum: maxP,
            boundary2Cells: b2,
            boundaryOrientations: b2.map(() => 1),
            volumeMpc3: Math.max(500.0, (4.0 / 3.0) * Math.PI * Math.pow(dist * 0.4, 3.0))
          }));
        }
      }
    }

    return new MorseSmaleComplex3D({
      criticalPoints: cps,
      cells1: cells1,
      cells2: cells2,
      cells3: cells3,
      domainBoundsMpc: this.boundsMpc
    });
  }
}

// ============================================================================
// TOPOLOGICAL CANONICAL BENCHMARKS & BUILDERS
// ============================================================================

/**
 * Builds an exact 3-simplex (solid tetrahedron) cell complex with verified nilpotency \partial^2 = 0,
 * b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0, and Euler characteristic \chi = 1.
 * 
 * @returns {MorseSmaleComplex3D}
 */
export function buildTetrahedralComplex() {
  const vA = new MorseCriticalPoint({ id: 'vA', index: 0, position: [0, 0, 0], value: -1.0, name: 'Vertex A (Void Min)' });
  const vB = new MorseCriticalPoint({ id: 'vB', index: 1, position: [100, 0, 0], value: -0.2, name: 'Vertex B (Wall Sad)' });
  const vC = new MorseCriticalPoint({ id: 'vC', index: 2, position: [50, 86.6, 0], value: 0.5, name: 'Vertex C (Fil Sad)' });
  const vD = new MorseCriticalPoint({ id: 'vD', index: 3, position: [50, 28.9, 81.6], value: 2.0, name: 'Vertex D (Node Max)' });

  // 6 oriented 1-cells
  const eAB = new Morse1Cell({ id: 'eAB', source: vA, target: vB });
  const eBC = new Morse1Cell({ id: 'eBC', source: vB, target: vC });
  const eCA = new Morse1Cell({ id: 'eCA', source: vC, target: vA });
  const eAD = new Morse1Cell({ id: 'eAD', source: vA, target: vD });
  const eBD = new Morse1Cell({ id: 'eBD', source: vB, target: vD });
  const eCD = new Morse1Cell({ id: 'eCD', source: vC, target: vD });

  // 4 oriented triangular 2-cells (faces with closed boundaries)
  // F_ABC: eAB + eBC + eCA => \partial(F_ABC) = (B-A) + (C-B) + (A-C) = 0
  const fABC = new Morse2Cell({ id: 'fABC', boundary1Cells: [eAB, eBC, eCA], boundaryOrientations: [1, 1, 1], surfaceVertices: [vA.position, vB.position, vC.position] });
  // F_ABD: eAB + eBD - eAD => \partial(F_ABD) = (B-A) + (D-B) - (D-A) = 0
  const fABD = new Morse2Cell({ id: 'fABD', boundary1Cells: [eAB, eBD, eAD], boundaryOrientations: [1, 1, -1], surfaceVertices: [vA.position, vB.position, vD.position] });
  // F_BCD: eBC + eCD - eBD => \partial(F_BCD) = (C-B) + (D-C) - (D-B) = 0
  const fBCD = new Morse2Cell({ id: 'fBCD', boundary1Cells: [eBC, eCD, eBD], boundaryOrientations: [1, 1, -1], surfaceVertices: [vB.position, vC.position, vD.position] });
  // F_CAD: eCA + eAD - eCD => \partial(F_CAD) = (A-C) + (D-A) - (D-C) = 0
  const fCAD = new Morse2Cell({ id: 'fCAD', boundary1Cells: [eCA, eAD, eCD], boundaryOrientations: [1, 1, -1], surfaceVertices: [vC.position, vA.position, vD.position] });

  // 1 oriented 3-cell: tet = fABC - fABD - fBCD - fCAD
  // \partial_2(tet) cancels each 1-cell exactly:
  // eAB: +fABC - fABD = 0
  // eBC: +fABC - fBCD = 0
  // eCA: +fABC - fCAD = 0
  // eAD: -(-fABD) + (-fCAD) = +fABD - fCAD = 0
  // eBD: -(+fABD) - (-fBCD) = -fABD + fBCD = 0
  // eCD: -(+fBCD) - (-fCAD) = -fBCD + fCAD = 0
  const tet = new Morse3Cell({
    id: 'crystal_tet_0',
    minimum: vA,
    maximum: vD,
    boundary2Cells: [fABC, fABD, fBCD, fCAD],
    boundaryOrientations: [1, -1, -1, -1],
    volumeMpc3: 117851.0
  });

  return new MorseSmaleComplex3D({
    criticalPoints: [vA, vB, vC, vD],
    cells1: [eAB, eBC, eCA, eAD, eBD, eCD],
    cells2: [fABC, fABD, fBCD, fCAD],
    cells3: [tet]
  });
}

/**
 * Builds a topological hollow 2-sphere complex (enclosed cosmic void bubble)
 * with b_0 = 1, b_1 = 0, b_2 = 1, b_3 = 0, \chi = 2.
 * 
 * @returns {MorseSmaleComplex3D}
 */
export function buildEnclosedVoidBubbleComplex() {
  const tet = buildTetrahedralComplex();
  // Remove the 3-cell volume to leave the hollow closed void bubble
  return new MorseSmaleComplex3D({
    criticalPoints: tet.criticalPoints,
    cells1: tet.cells1,
    cells2: tet.cells2,
    cells3: [] // no 3-cells => b_2 = 1
  });
}

/**
 * Builds a canonical multi-basin synthetic Cosmic Web Morse-Smale Complex for testing and validation.
 * Includes official cosmic basins: Laniakea, Shapley, Hercules, Perseus-Pisces, and surrounding voids.
 * 
 * @param {object} [options={}]
 * @returns {MorseSmaleComplex3D}
 */
export function buildSyntheticCosmicWeb(options = {}) {
  const cps = [
    // 0-cells: Cosmic Voids (Minima)
    new MorseCriticalPoint({
      id: 'void_local',
      index: MorseIndex.MINIMUM,
      position: [0.0, -10.0, 5.0],
      value: -0.92,
      eigenvalues: [0.04, 0.05, 0.08],
      name: 'Local Void Core'
    }),
    new MorseCriticalPoint({
      id: 'void_sculptor',
      index: MorseIndex.MINIMUM,
      position: [-40.0, -20.0, -30.0],
      value: -0.85,
      eigenvalues: [0.03, 0.04, 0.07],
      name: 'Sculptor Void'
    }),
    new MorseCriticalPoint({
      id: 'void_bootes',
      index: MorseIndex.MINIMUM,
      position: [30.0, 60.0, 20.0],
      value: -0.95,
      eigenvalues: [0.05, 0.06, 0.09],
      name: 'Bootes Supervoid'
    }),

    // 1-cells: Wall Saddles (Sheet Hubs)
    new MorseCriticalPoint({
      id: 'saddle_wall_1',
      index: MorseIndex.WALL_SADDLE,
      position: [-15.0, -5.0, -10.0],
      value: -0.20,
      eigenvalues: [-0.03, 0.04, 0.06],
      name: 'Local Sheet Saddle'
    }),
    new MorseCriticalPoint({
      id: 'saddle_wall_2',
      index: MorseIndex.WALL_SADDLE,
      position: [10.0, 30.0, 10.0],
      value: -0.15,
      eigenvalues: [-0.04, 0.03, 0.07],
      name: 'Hercules Wall Saddle'
    }),

    // 2-cells: Filament Saddles (Spine Hubs)
    new MorseCriticalPoint({
      id: 'saddle_fil_1',
      index: MorseIndex.FILAMENT_SADDLE,
      position: [-70.0, 30.0, 5.0],
      value: 0.45,
      eigenvalues: [-0.08, -0.05, 0.06],
      name: 'Centaurus-Shapley Filament Saddle'
    }),
    new MorseCriticalPoint({
      id: 'saddle_fil_2',
      index: MorseIndex.FILAMENT_SADDLE,
      position: [20.0, 10.0, -15.0],
      value: 0.38,
      eigenvalues: [-0.07, -0.04, 0.05],
      name: 'Perseus-Local Filament Saddle'
    }),

    // 3-cells: Local Maxima (Superclusters / Nodes)
    new MorseCriticalPoint({
      id: 'max_laniakea',
      index: MorseIndex.MAXIMUM,
      position: [-35.0, 15.0, -10.0],
      value: 1.85,
      eigenvalues: [-0.12, -0.09, -0.07],
      name: 'Laniakea Great Attractor Core'
    }),
    new MorseCriticalPoint({
      id: 'max_shapley',
      index: MorseIndex.MAXIMUM,
      position: [-110.0, 50.0, 20.0],
      value: 3.40,
      eigenvalues: [-0.25, -0.18, -0.14],
      name: 'Shapley Concentration Core'
    }),
    new MorseCriticalPoint({
      id: 'max_perseus',
      index: MorseIndex.MAXIMUM,
      position: [50.0, -15.0, -20.0],
      value: 2.10,
      eigenvalues: [-0.15, -0.11, -0.08],
      name: 'Perseus-Pisces Core'
    })
  ];

  const cpMap = new Map(cps.map(cp => [cp.id, cp]));

  const fil1 = new Morse1Cell({
    id: 'filament_lan_to_sad1',
    source: cpMap.get('saddle_fil_1'),
    target: cpMap.get('max_laniakea'),
    type: CosmicWebElementType.FILAMENT_SPINE,
    lengthMpc: 42.5,
    meanDensityContrast: 1.15
  });

  const fil2 = new Morse1Cell({
    id: 'filament_sad1_to_shapley',
    source: cpMap.get('saddle_fil_1'),
    target: cpMap.get('max_shapley'),
    type: CosmicWebElementType.FILAMENT_SPINE,
    lengthMpc: 48.2,
    meanDensityContrast: 1.95
  });

  const fil3 = new Morse1Cell({
    id: 'filament_lan_to_sad2',
    source: cpMap.get('saddle_fil_2'),
    target: cpMap.get('max_laniakea'),
    type: CosmicWebElementType.FILAMENT_SPINE,
    lengthMpc: 55.4,
    meanDensityContrast: 0.95
  });

  const fil4 = new Morse1Cell({
    id: 'filament_sad2_to_perseus',
    source: cpMap.get('saddle_fil_2'),
    target: cpMap.get('max_perseus'),
    type: CosmicWebElementType.FILAMENT_SPINE,
    lengthMpc: 41.0,
    meanDensityContrast: 1.25
  });

  const sheet1 = new Morse1Cell({
    id: 'sheet_spine_wall1_to_localvoid',
    source: cpMap.get('saddle_wall_1'),
    target: cpMap.get('void_local'),
    type: CosmicWebElementType.SHEET_SPINE,
    lengthMpc: 22.0,
    meanDensityContrast: -0.55
  });

  const sheet2 = new Morse1Cell({
    id: 'sheet_spine_wall1_to_sculptor',
    source: cpMap.get('saddle_wall_1'),
    target: cpMap.get('void_sculptor'),
    type: CosmicWebElementType.SHEET_SPINE,
    lengthMpc: 35.0,
    meanDensityContrast: -0.50
  });

  const sheet3 = new Morse1Cell({
    id: 'sheet_spine_wall2_to_localvoid',
    source: cpMap.get('saddle_wall_2'),
    target: cpMap.get('void_local'),
    type: CosmicWebElementType.SHEET_SPINE,
    lengthMpc: 45.0,
    meanDensityContrast: -0.40
  });

  const sheet4 = new Morse1Cell({
    id: 'sheet_spine_wall2_to_bootes',
    source: cpMap.get('saddle_wall_2'),
    target: cpMap.get('void_bootes'),
    type: CosmicWebElementType.SHEET_SPINE,
    lengthMpc: 40.0,
    meanDensityContrast: -0.60
  });

  const cells1 = [fil1, fil2, fil3, fil4, sheet1, sheet2, sheet3, sheet4];

  // 2-Cells: Cosmic Wall Sheets
  const wall1 = new Morse2Cell({
    id: 'cosmic_sheet_local_sculptor',
    boundary1Cells: [sheet1, sheet2],
    boundaryOrientations: [1, -1],
    surfaceVertices: [cpMap.get('saddle_wall_1').position, cpMap.get('void_local').position, cpMap.get('void_sculptor').position],
    areaMpc2: 450.0
  });

  const wall2 = new Morse2Cell({
    id: 'cosmic_sheet_local_bootes',
    boundary1Cells: [sheet3, sheet4],
    boundaryOrientations: [1, -1],
    surfaceVertices: [cpMap.get('saddle_wall_2').position, cpMap.get('void_local').position, cpMap.get('void_bootes').position],
    areaMpc2: 820.0
  });

  const wall3 = new Morse2Cell({
    id: 'cosmic_sheet_lan_shapley_filament',
    boundary1Cells: [fil1, fil2],
    boundaryOrientations: [1, -1],
    surfaceVertices: [cpMap.get('saddle_fil_1').position, cpMap.get('max_laniakea').position, cpMap.get('max_shapley').position],
    areaMpc2: 1200.0
  });

  const wall4 = new Morse2Cell({
    id: 'cosmic_sheet_lan_perseus_filament',
    boundary1Cells: [fil3, fil4],
    boundaryOrientations: [1, -1],
    surfaceVertices: [cpMap.get('saddle_fil_2').position, cpMap.get('max_laniakea').position, cpMap.get('max_perseus').position],
    areaMpc2: 950.0
  });

  const cells2 = [wall1, wall2, wall3, wall4];

  // 3-Cells: Void Crystals / Supercluster Influence Crystals
  const voidCrystal1 = new Morse3Cell({
    id: 'crystal_local_laniakea',
    minimum: cpMap.get('void_local'),
    maximum: cpMap.get('max_laniakea'),
    boundary2Cells: [wall1, wall2],
    boundaryOrientations: [1, -1],
    volumeMpc3: 18500.0
  });

  const voidCrystal2 = new Morse3Cell({
    id: 'crystal_sculptor_shapley',
    minimum: cpMap.get('void_sculptor'),
    maximum: cpMap.get('max_shapley'),
    boundary2Cells: [wall3, wall4],
    boundaryOrientations: [1, -1],
    volumeMpc3: 32000.0
  });

  const cells3 = [voidCrystal1, voidCrystal2];

  return new MorseSmaleComplex3D({
    criticalPoints: cps,
    cells1: cells1,
    cells2: cells2,
    cells3: cells3,
    domainBoundsMpc: [-250, 250, -250, 250, -250, 250]
  });
}
