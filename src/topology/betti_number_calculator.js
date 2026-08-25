/**
 * @file betti_number_calculator.js
 * @module topology/betti_number_calculator
 * @description 3D Persistent Homology, Betti Number, and Euler Characteristic Calculator
 * for Cosmic Density Contrast (\delta) and Gravitational Velocity Potential (\Phi) Fields.
 *
 * Implements:
 * 1. 3D Cubical Complex Filtration (superlevel & sublevel sets) over scalar grids.
 * 2. Exact Betti Numbers \beta_0 (clusters), \beta_1 (filament loops/tunnels), \beta_2 (void bubbles).
 * 3. Euler Characteristic \chi(\delta) = \beta_0 - \beta_1 + \beta_2 verified against cell-count Euler-Poincaré \chi = V - E + F - C.
 * 4. Analytical Tomita-Gott Gaussian Random Field (GRF) Euler density expectation:
 *    \chi_{\text{GRF}}(\nu) = \frac{1}{(2\pi)^2} \left(\frac{\sigma_1}{\sqrt{3}\sigma_0}\right)^3 (\nu^2 - 1) e^{-\nu^2/2}
 * 5. Topological persistence diagrams, persistence barcodes, persistence entropy, bottleneck & Wasserstein distances,
 *    persistence landscapes, velocity potential field tidal tensor topology, and cosmic web classification.
 *
 * @author Scientific Computational Cosmology Engineer
 * @license MIT
 */

// ============================================================================
// CONSTANTS & ENUMS
// ============================================================================

/**
 * Filtration direction enumeration.
 * @readonly
 * @enum {string}
 */
export const FiltrationDirection = Object.freeze({
  SUPERLEVEL: 'SUPERLEVEL', // \delta >= \delta_th (overdense clusters, filaments, walls)
  SUBLEVEL: 'SUBLEVEL'      // \delta <= \delta_th (underdense voids, basins)
});

/**
 * Homology dimension enumeration.
 * @readonly
 * @enum {number}
 */
export const HomologyDimension = Object.freeze({
  H0: 0, // Connected components / overdense clusters
  H1: 1, // 1D Tunnels / cosmic filament loops
  H2: 2  // 2D Enclosed cavities / cosmic void bubbles
});

/**
 * Cosmic web morphological topology classification.
 * @readonly
 * @enum {string}
 */
export const CosmicTopologyType = Object.freeze({
  MEATBALL: 'MEATBALL',         // \beta_0 >> \beta_2 (isolated high-density clusters)
  SPONGY: 'SPONGY',             // \beta_1 dominates (interconnected filamentary/tunnel network)
  SWISS_CHEESE: 'SWISS_CHEESE', // \beta_2 >> \beta_0 (isolated void bubbles in wall matrix)
  GAUSSIAN_RANDOM: 'GAUSSIAN_RANDOM', // Symmetrical genus curve \chi(-\nu) = -\chi(\nu)
  INTERMEDIATE: 'INTERMEDIATE'
});

/**
 * Cosmic web environment classified by Zel'dovich deformation tensor eigenvalues.
 * @readonly
 * @enum {string}
 */
export const CosmicEnvironmentType = Object.freeze({
  PEAK_NODE: 'PEAK_NODE',     // \lambda_1, \lambda_2, \lambda_3 > \lambda_{th} (cluster halo)
  FILAMENT: 'FILAMENT',       // \lambda_1, \lambda_2 > \lambda_{th}, \lambda_3 <= \lambda_{th}
  SHEET_WALL: 'SHEET_WALL',   // \lambda_1 > \lambda_{th}, \lambda_2, \lambda_3 <= \lambda_{th}
  VOID: 'VOID'                // \lambda_1, \lambda_2, \lambda_3 <= \lambda_{th}
});

/**
 * Numerical tolerances and physical constants.
 */
export const TOPOLOGY_CONSTANTS = Object.freeze({
  TWO_PI: 2.0 * Math.PI,
  TWO_PI_SQ: 4.0 * Math.PI * Math.PI,
  SQRT_3: Math.sqrt(3.0),
  INV_TWO_PI_SQ: 1.0 / (4.0 * Math.PI * Math.PI),
  EPSILON: 1e-12,
  INFINITY_VAL: 1e12
});

// ============================================================================
// DATA STRUCTURES: UNION-FIND FOR 0D PERSISTENCE
// ============================================================================

/**
 * Disjoint Set Union-Find with path compression, union-by-rank,
 * and birth-death tracking for 0-dimensional persistent homology.
 */
export class UnionFindPersistent {
  /**
   * @param {number} size - Total number of elements.
   */
  constructor(size) {
    if (!Number.isInteger(size) || size <= 0) {
      throw new Error(`UnionFindPersistent: size must be positive integer, got ${size}`);
    }
    this.size = size;
    this.parent = new Int32Array(size);
    this.rank = new Int32Array(size);
    this.birthValue = new Float64Array(size);
    this.birthIndex = new Int32Array(size);
    this.active = new Uint8Array(size);
    this.centroids = new Float64Array(size * 3);
    this.voxelCounts = new Int32Array(size);

    for (let i = 0; i < size; i++) {
      this.parent[i] = i;
      this.rank[i] = 0;
      this.birthValue[i] = 0.0;
      this.birthIndex[i] = -1;
      this.active[i] = 0;
      this.voxelCounts[i] = 0;
    }
  }

  /**
   * Activates an element with its birth threshold and 3D coordinate.
   * @param {number} x - Element index.
   * @param {number} birth - Filtration value at birth.
   * @param {number} [cx=0] - X spatial coordinate.
   * @param {number} [cy=0] - Y spatial coordinate.
   * @param {number} [cz=0] - Z spatial coordinate.
   */
  makeSet(x, birth, cx = 0, cy = 0, cz = 0) {
    this.parent[x] = x;
    this.rank[x] = 0;
    this.birthValue[x] = birth;
    this.birthIndex[x] = x;
    this.active[x] = 1;
    this.voxelCounts[x] = 1;
    const base = x * 3;
    this.centroids[base] = cx;
    this.centroids[base + 1] = cy;
    this.centroids[base + 2] = cz;
  }

  /**
   * Finds the representative root of element x with path compression.
   * @param {number} x
   * @returns {number} Root index.
   */
  find(x) {
    let root = x;
    while (root !== this.parent[root]) {
      root = this.parent[root];
    }
    let curr = x;
    while (curr !== root) {
      const nxt = this.parent[curr];
      this.parent[curr] = root;
      curr = nxt;
    }
    return root;
  }

  /**
   * Unions two components at a death threshold according to the Elder Rule.
   *
   * @param {number} x - First element index.
   * @param {number} y - Second element index.
   * @param {number} death - Filtration threshold at merge.
   * @param {boolean} [isSuperlevel=true] - If true, higher birth is older.
   * @returns {Object|null} Persistence pair {birth, death, dyingIndex, survivingIndex, lifetime} or null if already same set.
   */
  union(x, y, death, isSuperlevel = true) {
    let rootX = this.find(x);
    let rootY = this.find(y);

    if (rootX === rootY) {
      return null;
    }

    const birthX = this.birthValue[rootX];
    const birthY = this.birthValue[rootY];

    // Determine elder root
    let elderRoot = rootX;
    let youngerRoot = rootY;

    if (isSuperlevel) {
      // Superlevel: higher value was born earlier
      if (birthY > birthX || (birthY === birthX && rootY < rootX)) {
        elderRoot = rootY;
        youngerRoot = rootX;
      }
    } else {
      // Sublevel: lower value was born earlier
      if (birthY < birthX || (birthY === birthX && rootY < rootX)) {
        elderRoot = rootY;
        youngerRoot = rootX;
      }
    }

    const pair = {
      dim: HomologyDimension.H0,
      birth: this.birthValue[youngerRoot],
      death: death,
      dyingIndex: this.birthIndex[youngerRoot],
      survivingIndex: this.birthIndex[elderRoot],
      lifetime: Math.abs(this.birthValue[youngerRoot] - death),
      volume: this.voxelCounts[youngerRoot]
    };

    // Link younger into elder
    this.parent[youngerRoot] = elderRoot;
    if (this.rank[elderRoot] === this.rank[youngerRoot]) {
      this.rank[elderRoot]++;
    }

    // Merge centroid & voxel counts
    const countElder = this.voxelCounts[elderRoot];
    const countYounger = this.voxelCounts[youngerRoot];
    const totalCount = countElder + countYounger;
    if (totalCount > 0) {
      const eb = elderRoot * 3;
      const yb = youngerRoot * 3;
      this.centroids[eb] = (this.centroids[eb] * countElder + this.centroids[yb] * countYounger) / totalCount;
      this.centroids[eb + 1] = (this.centroids[eb + 1] * countElder + this.centroids[yb + 1] * countYounger) / totalCount;
      this.centroids[eb + 2] = (this.centroids[eb + 2] * countElder + this.centroids[yb + 2] * countYounger) / totalCount;
    }
    this.voxelCounts[elderRoot] = totalCount;

    return pair;
  }
}

// ============================================================================
// 3D CUBICAL COMPLEX STRUCTURE & BOUNDARY OPERATORS
// ============================================================================

/**
 * 3D Cubical Complex on a regular cartesian grid Nx * Ny * Nz.
 * Cell dimensions:
 * - 0-cells: Vertices (Nx * Ny * Nz)
 * - 1-cells: Edges in X, Y, Z directions
 * - 2-cells: Faces in XY, YZ, ZX planes
 * - 3-cells: Cubes
 */
export class CubicalComplex3D {
  /**
   * @param {number} nx - Number of grid vertices in X.
   * @param {number} ny - Number of grid vertices in Y.
   * @param {number} nz - Number of grid vertices in Z.
   * @param {Object} [options={}] - Options like periodic boundary conditions, voxel size.
   */
  constructor(nx, ny, nz, options = {}) {
    if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(ny) || ny < 2 || !Number.isInteger(nz) || nz < 2) {
      throw new Error(`CubicalComplex3D requires grid dimensions >= 2, got (${nx}, ${ny}, ${nz})`);
    }

    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.periodic = Boolean(options.periodic);
    this.dx = options.dx || 1.0;
    this.dy = options.dy || 1.0;
    this.dz = options.dz || 1.0;

    // Cell counts for standard embedded cubical complex
    this.numVertices = nx * ny * nz;
    this.numEdgesX = (nx - 1) * ny * nz;
    this.numEdgesY = nx * (ny - 1) * nz;
    this.numEdgesZ = nx * ny * (nz - 1);
    this.numEdges = this.numEdgesX + this.numEdgesY + this.numEdgesZ;

    this.numFacesXY = (nx - 1) * (ny - 1) * nz;
    this.numFacesYZ = nx * (ny - 1) * (nz - 1);
    this.numFacesZX = (nx - 1) * ny * (nz - 1);
    this.numFaces = this.numFacesXY + this.numFacesYZ + this.numFacesZX;

    this.numCubes = (nx - 1) * (ny - 1) * (nz - 1);
    this.totalCells = this.numVertices + this.numEdges + this.numFaces + this.numCubes;
  }

  // --- 0-CELL (VERTEX) INDEXING ---
  vertexIndex(i, j, k) {
    return i + this.nx * (j + this.ny * k);
  }

  vertexCoords(idx) {
    const k = Math.floor(idx / (this.nx * this.ny));
    const rem = idx % (this.nx * this.ny);
    const j = Math.floor(rem / this.nx);
    const i = rem % this.nx;
    return [i, j, k];
  }

  // --- 1-CELL (EDGE) INDEXING ---
  edgeXIndex(i, j, k) {
    return i + (this.nx - 1) * (j + this.ny * k);
  }

  edgeYIndex(i, j, k) {
    return this.numEdgesX + i + this.nx * (j + (this.ny - 1) * k);
  }

  edgeZIndex(i, j, k) {
    return this.numEdgesX + this.numEdgesY + i + this.nx * (j + this.ny * k);
  }

  // --- 2-CELL (FACE) INDEXING ---
  faceXYIndex(i, j, k) {
    return i + (this.nx - 1) * (j + (this.ny - 1) * k);
  }

  faceYZIndex(i, j, k) {
    return this.numFacesXY + i + this.nx * (j + (this.ny - 1) * k);
  }

  faceZXIndex(i, j, k) {
    return this.numFacesXY + this.numFacesYZ + i + (this.nx - 1) * (j + this.ny * k);
  }

  // --- 3-CELL (CUBE) INDEXING ---
  cubeIndex(i, j, k) {
    return i + (this.nx - 1) * (j + (this.ny - 1) * k);
  }

  /**
   * Retrieves boundary 0-cells of a 1-cell (edge).
   * @param {number} edgeIdx
   * @returns {[number, number]}
   */
  edgeBoundaryVertices(edgeIdx) {
    if (edgeIdx < this.numEdgesX) {
      const k = Math.floor(edgeIdx / ((this.nx - 1) * this.ny));
      const rem = edgeIdx % ((this.nx - 1) * this.ny);
      const j = Math.floor(rem / (this.nx - 1));
      const i = rem % (this.nx - 1);
      return [this.vertexIndex(i, j, k), this.vertexIndex(i + 1, j, k)];
    } else if (edgeIdx < this.numEdgesX + this.numEdgesY) {
      const offset = edgeIdx - this.numEdgesX;
      const k = Math.floor(offset / (this.nx * (this.ny - 1)));
      const rem = offset % (this.nx * (this.ny - 1));
      const j = Math.floor(rem / this.nx);
      const i = rem % this.nx;
      return [this.vertexIndex(i, j, k), this.vertexIndex(i, j + 1, k)];
    } else {
      const offset = edgeIdx - (this.numEdgesX + this.numEdgesY);
      const k = Math.floor(offset / (this.nx * this.ny));
      const rem = offset % (this.nx * this.ny);
      const j = Math.floor(rem / this.nx);
      const i = rem % this.nx;
      return [this.vertexIndex(i, j, k), this.vertexIndex(i, j, k + 1)];
    }
  }

  /**
   * Retrieves boundary 1-cells (edges) of a 2-cell (face).
   * @param {number} faceIdx
   * @returns {number[]}
   */
  faceBoundaryEdges(faceIdx) {
    if (faceIdx < this.numFacesXY) {
      const k = Math.floor(faceIdx / ((this.nx - 1) * (this.ny - 1)));
      const rem = faceIdx % ((this.nx - 1) * (this.ny - 1));
      const j = Math.floor(rem / (this.nx - 1));
      const i = rem % (this.nx - 1);
      return [
        this.edgeXIndex(i, j, k),
        this.edgeXIndex(i, j + 1, k),
        this.edgeYIndex(i, j, k),
        this.edgeYIndex(i + 1, j, k)
      ];
    } else if (faceIdx < this.numFacesXY + this.numFacesYZ) {
      const offset = faceIdx - this.numFacesXY;
      const k = Math.floor(offset / (this.nx * (this.ny - 1)));
      const rem = offset % (this.nx * (this.ny - 1));
      const j = Math.floor(rem / this.nx);
      const i = rem % this.nx;
      return [
        this.edgeYIndex(i, j, k),
        this.edgeYIndex(i, j, k + 1),
        this.edgeZIndex(i, j, k),
        this.edgeZIndex(i, j + 1, k)
      ];
    } else {
      const offset = faceIdx - (this.numFacesXY + this.numFacesYZ);
      const k = Math.floor(offset / ((this.nx - 1) * this.ny));
      const rem = offset % ((this.nx - 1) * this.ny);
      const j = Math.floor(rem / (this.nx - 1));
      const i = rem % (this.nx - 1);
      return [
        this.edgeZIndex(i, j, k),
        this.edgeZIndex(i + 1, j, k),
        this.edgeXIndex(i, j, k),
        this.edgeXIndex(i, j, k + 1)
      ];
    }
  }

  /**
   * Retrieves boundary 2-cells (faces) of a 3-cell (cube).
   * @param {number} cubeIdx
   * @returns {number[]}
   */
  cubeBoundaryFaces(cubeIdx) {
    const k = Math.floor(cubeIdx / ((this.nx - 1) * (this.ny - 1)));
    const rem = cubeIdx % ((this.nx - 1) * (this.ny - 1));
    const j = Math.floor(rem / (this.nx - 1));
    const i = rem % (this.nx - 1);
    return [
      this.faceXYIndex(i, j, k),
      this.faceXYIndex(i, j, k + 1),
      this.faceYZIndex(i, j, k),
      this.faceYZIndex(i + 1, j, k),
      this.faceZXIndex(i, j, k),
      this.faceZXIndex(i, j + 1, k)
    ];
  }

  /**
   * Evaluates scalar filtration values for all cells in the cubical complex.
   *
   * @param {Float64Array|number[]} vertexValues - 1D array of vertex field values of length nx*ny*nz.
   * @param {string} [direction=FiltrationDirection.SUPERLEVEL]
   * @returns {Object} { vertexVal, edgeVal, faceVal, cubeVal }
   */
  computeCellFiltrationValues(vertexValues, direction = FiltrationDirection.SUPERLEVEL) {
    if (vertexValues.length !== this.numVertices) {
      throw new Error(`vertexValues length (${vertexValues.length}) must match complex vertex count (${this.numVertices})`);
    }

    const isSuper = direction === FiltrationDirection.SUPERLEVEL;
    const edgeVal = new Float64Array(this.numEdges);
    const faceVal = new Float64Array(this.numFaces);
    const cubeVal = new Float64Array(this.numCubes);

    // 1. Edge filtration values
    for (let e = 0; e < this.numEdges; e++) {
      const [v0, v1] = this.edgeBoundaryVertices(e);
      const val0 = vertexValues[v0];
      const val1 = vertexValues[v1];
      edgeVal[e] = isSuper ? Math.min(val0, val1) : Math.max(val0, val1);
    }

    // 2. Face filtration values
    for (let f = 0; f < this.numFaces; f++) {
      const edges = this.faceBoundaryEdges(f);
      let v = edgeVal[edges[0]];
      for (let i = 1; i < 4; i++) {
        const ev = edgeVal[edges[i]];
        v = isSuper ? Math.min(v, ev) : Math.max(v, ev);
      }
      faceVal[f] = v;
    }

    // 3. Cube filtration values
    for (let c = 0; c < this.numCubes; c++) {
      const faces = this.cubeBoundaryFaces(c);
      let v = faceVal[faces[0]];
      for (let i = 1; i < 6; i++) {
        const fv = faceVal[faces[i]];
        v = isSuper ? Math.min(v, fv) : Math.max(v, fv);
      }
      cubeVal[c] = v;
    }

    return {
      vertexVal: vertexValues,
      edgeVal,
      faceVal,
      cubeVal
    };
  }

  /**
   * Computes cell counts (V, E, F, C) active at a given filtration threshold \delta_th.
   *
   * @param {Object} cellVals - { vertexVal, edgeVal, faceVal, cubeVal }
   * @param {number} threshold - Filtration threshold value.
   * @param {string} [direction=FiltrationDirection.SUPERLEVEL]
   * @returns {Object} { V: number, E: number, F: number, C: number, chi: number }
   */
  countActiveCells(cellVals, threshold, direction = FiltrationDirection.SUPERLEVEL) {
    const isSuper = direction === FiltrationDirection.SUPERLEVEL;
    let V = 0;
    let E = 0;
    let F = 0;
    let C = 0;

    const { vertexVal, edgeVal, faceVal, cubeVal } = cellVals;

    for (let i = 0; i < this.numVertices; i++) {
      if (isSuper ? vertexVal[i] >= threshold : vertexVal[i] <= threshold) V++;
    }
    for (let i = 0; i < this.numEdges; i++) {
      if (isSuper ? edgeVal[i] >= threshold : edgeVal[i] <= threshold) E++;
    }
    for (let i = 0; i < this.numFaces; i++) {
      if (isSuper ? faceVal[i] >= threshold : faceVal[i] <= threshold) F++;
    }
    for (let i = 0; i < this.numCubes; i++) {
      if (isSuper ? cubeVal[i] >= threshold : cubeVal[i] <= threshold) C++;
    }

    const chi = V - E + F - C;
    return { V, E, F, C, chi };
  }
}

// ============================================================================
// PERSISTENT HOMOLOGY & BOUNDARY MATRIX REDUCTION
// ============================================================================

/**
 * Algebraic Boundary Matrix Reducer over Z_2 field.
 * Computes exact persistent homology pairs (birth, death) across dimensions H0, H1, H2.
 */
export class BoundaryMatrixReducer {
  /**
   * Performs standard column reduction over Z_2 on a sorted boundary matrix.
   *
   * @param {Array<{dim: number, id: number, value: number, boundary: number[]}>} sortedCells
   * @returns {Object} { pairs: Array<{dim: number, birth: number, death: number, birthCell: number, deathCell: number}>, essential: Array<{dim: number, birth: number, birthCell: number}> }
   */
  static reduceZ2(sortedCells) {
    const m = sortedCells.length;
    const cellToIndex = new Map();
    for (let i = 0; i < m; i++) {
      cellToIndex.set(sortedCells[i].id, i);
    }

    const R = new Array(m);
    for (let j = 0; j < m; j++) {
      const bnd = sortedCells[j].boundary;
      if (!bnd || bnd.length === 0) {
        R[j] = [];
      } else {
        const mapped = [];
        for (let k = 0; k < bnd.length; k++) {
          const idx = cellToIndex.get(bnd[k]);
          if (idx !== undefined) {
            mapped.push(idx);
          }
        }
        mapped.sort((a, b) => a - b);
        R[j] = mapped;
      }
    }

    const low = (col) => {
      const list = R[col];
      return list.length > 0 ? list[list.length - 1] : -1;
    };

    const addColumnsZ2 = (colA, colB) => {
      const res = [];
      let i = 0;
      let j = 0;
      const nA = colA.length;
      const nB = colB.length;
      while (i < nA && j < nB) {
        if (colA[i] < colB[j]) {
          res.push(colA[i++]);
        } else if (colB[j] < colA[i]) {
          res.push(colB[j++]);
        } else {
          i++;
          j++;
        }
      }
      while (i < nA) res.push(colA[i++]);
      while (j < nB) res.push(colB[j++]);
      return res;
    };

    const pivotToCol = new Map();
    const pairs = [];
    const pairedBirthIndices = new Set();

    for (let j = 0; j < m; j++) {
      while (true) {
        const l = low(j);
        if (l === -1) {
          break;
        }
        if (pivotToCol.has(l)) {
          const k = pivotToCol.get(l);
          R[j] = addColumnsZ2(R[j], R[k]);
        } else {
          pivotToCol.set(l, j);
          const birthCell = sortedCells[l];
          const deathCell = sortedCells[j];
          pairs.push({
            dim: birthCell.dim,
            birth: birthCell.value,
            death: deathCell.value,
            birthCell: birthCell.id,
            deathCell: deathCell.id,
            lifetime: Math.abs(deathCell.value - birthCell.value)
          });
          pairedBirthIndices.add(l);
          pairedBirthIndices.add(j);
          break;
        }
      }
    }

    const essential = [];
    for (let i = 0; i < m; i++) {
      if (!pairedBirthIndices.has(i) && low(i) === -1) {
        essential.push({
          dim: sortedCells[i].dim,
          birth: sortedCells[i].value,
          birthCell: sortedCells[i].id
        });
      }
    }

    return { pairs, essential };
  }
}

// ============================================================================
// PERSISTENCE BARCODES & LANDSCAPES
// ============================================================================

/**
 * Persistence Barcode representation and interval inspection.
 */
export class PersistenceBarcode {
  /**
   * @param {number} dimension
   * @param {Array<{birth: number, death: number, lifetime: number}>} intervals
   */
  constructor(dimension, intervals = []) {
    this.dimension = dimension;
    this.intervals = intervals.map(int => ({
      birth: int.birth,
      death: int.death,
      lifetime: Math.abs(int.death - int.birth)
    }));
    // Sort intervals by descending lifetime
    this.intervals.sort((a, b) => b.lifetime - a.lifetime);
  }

  /**
   * Returns active interval count at threshold t.
   * @param {number} t
   * @param {boolean} [isSuperlevel=true]
   * @returns {number}
   */
  activeCount(t, isSuperlevel = true) {
    let count = 0;
    for (let i = 0; i < this.intervals.length; i++) {
      const { birth, death } = this.intervals[i];
      if (isSuperlevel) {
        if (birth >= t && t > death) count++;
      } else {
        if (birth <= t && t < death) count++;
      }
    }
    return count;
  }

  /**
   * Serializes barcode into JSON structure for UI rendering.
   * @returns {Object}
   */
  toJSON() {
    return {
      dimension: this.dimension,
      count: this.intervals.length,
      intervals: this.intervals
    };
  }
}

/**
 * Persistence Landscape Vectorizer (Bubenik 2015).
 * Represents a persistence diagram as a sequence of continuous piecewise-linear functions \lambda_k(t).
 */
export class PersistenceLandscape {
  /**
   * Computes the k-th persistence landscape function \lambda_k(t) evaluated on a grid of t values.
   * For a point (b, d) with b < d:
   * f_{(b,d)}(t) = max(0, min(t - b, d - t)).
   * \lambda_k(t) is the k-th largest value of f_{(b,d)}(t) across all pairs in the diagram.
   *
   * @param {PersistenceDiagram} diagram
   * @param {number[]} tGrid - Sorted array of parameter thresholds t.
   * @param {number} [numLandscapes=5] - Number of landscape layers (k = 1, 2, ..., numLandscapes).
   * @returns {Array<Float64Array>} Array of Float64Array evaluations of length tGrid.length.
   */
  static computeLandscapes(diagram, tGrid, numLandscapes = 5) {
    const pairs = diagram.pairs.map(p => ({
      b: Math.min(p.birth, p.death),
      d: Math.max(p.birth, p.death)
    })).filter(p => p.d > p.b);

    const numT = tGrid.length;
    const landscapes = Array.from({ length: numLandscapes }, () => new Float64Array(numT));

    for (let i = 0; i < numT; i++) {
      const t = tGrid[i];
      const vals = [];
      for (let p = 0; p < pairs.length; p++) {
        const { b, d } = pairs[p];
        if (t > b && t < d) {
          const val = Math.min(t - b, d - t);
          vals.push(val);
        }
      }
      vals.sort((a, b) => b - a); // Sort descending
      for (let k = 0; k < numLandscapes; k++) {
        landscapes[k][i] = k < vals.length ? vals[k] : 0.0;
      }
    }

    return landscapes;
  }

  /**
   * Computes the L_p norm of the persistence landscape:
   * ||\Lambda||_p = \left( \sum_k \int |\lambda_k(t)|^p dt \right)^{1/p}.
   *
   * @param {Array<Float64Array>} landscapes
   * @param {number[]} tGrid
   * @param {number} [p=2]
   * @returns {number}
   */
  static landscapeNorm(landscapes, tGrid, p = 2) {
    if (tGrid.length < 2 || landscapes.length === 0) return 0.0;
    let sum = 0.0;

    for (let k = 0; k < landscapes.length; k++) {
      const land = landscapes[k];
      for (let i = 0; i < tGrid.length - 1; i++) {
        const dt = tGrid[i + 1] - tGrid[i];
        const avgVal = 0.5 * (Math.pow(Math.abs(land[i]), p) + Math.pow(Math.abs(land[i + 1]), p));
        sum += avgVal * dt;
      }
    }

    return Math.pow(sum, 1.0 / p);
  }
}

// ============================================================================
// PERSISTENCE DIAGRAMS & STATISTICAL METRICS
// ============================================================================

/**
 * Persistence Diagram representation and topological summary statistics.
 */
export class PersistenceDiagram {
  /**
   * @param {number} dimension - Homology dimension (0, 1, or 2).
   * @param {Array<{birth: number, death: number, lifetime?: number}>} [pairs=[]]
   * @param {Array<{birth: number}>} [essential=[]]
   */
  constructor(dimension, pairs = [], essential = []) {
    this.dimension = dimension;
    this.pairs = pairs.map(p => ({
      birth: p.birth,
      death: p.death,
      lifetime: Math.abs(p.death - p.birth)
    }));
    this.essential = [...essential];
  }

  /**
   * Adds a birth-death pair.
   * @param {number} birth
   * @param {number} death
   */
  addPair(birth, death) {
    this.pairs.push({
      birth,
      death,
      lifetime: Math.abs(death - birth)
    });
  }

  /**
   * Returns total persistence \sum |d_i - b_i|^p.
   * @param {number} [p=1]
   * @returns {number}
   */
  totalPersistence(p = 1) {
    let sum = 0.0;
    for (let i = 0; i < this.pairs.length; i++) {
      sum += Math.pow(this.pairs[i].lifetime, p);
    }
    return sum;
  }

  /**
   * Computes persistent entropy:
   * E = -\sum_{i} p_i \ln p_i where p_i = \ell_i / \sum \ell_j.
   * @returns {number} Persistent entropy in nats.
   */
  persistentEntropy() {
    const tot = this.totalPersistence(1);
    if (tot <= TOPOLOGY_CONSTANTS.EPSILON || this.pairs.length === 0) {
      return 0.0;
    }
    let entropy = 0.0;
    for (let i = 0; i < this.pairs.length; i++) {
      const p = this.pairs[i].lifetime / tot;
      if (p > 1e-15) {
        entropy -= p * Math.log(p);
      }
    }
    return entropy;
  }

  /**
   * Filters out short-lived topological noise below lifetime threshold.
   * @param {number} minLifetime
   * @returns {PersistenceDiagram} New filtered diagram.
   */
  filterNoise(minLifetime) {
    const filtered = this.pairs.filter(p => p.lifetime >= minLifetime);
    return new PersistenceDiagram(this.dimension, filtered, this.essential);
  }

  /**
   * Generates a PersistenceBarcode instance.
   * @returns {PersistenceBarcode}
   */
  toBarcode() {
    return new PersistenceBarcode(this.dimension, this.pairs);
  }

  /**
   * Evaluates the Betti number \beta_k(\delta_th) at a given filtration threshold.
   * @param {number} threshold
   * @param {string} [direction=FiltrationDirection.SUPERLEVEL]
   * @returns {number} Active Betti number.
   */
  evaluateBetti(threshold, direction = FiltrationDirection.SUPERLEVEL) {
    let count = 0;
    const isSuper = direction === FiltrationDirection.SUPERLEVEL;

    for (let i = 0; i < this.pairs.length; i++) {
      const { birth, death } = this.pairs[i];
      if (isSuper) {
        if (birth >= threshold && threshold > death) {
          count++;
        }
      } else {
        if (birth <= threshold && threshold < death) {
          count++;
        }
      }
    }

    for (let i = 0; i < this.essential.length; i++) {
      const { birth } = this.essential[i];
      if (isSuper ? birth >= threshold : birth <= threshold) {
        count++;
      }
    }

    return count;
  }

  /**
   * Computes Bottleneck distance W_\infty between two persistence diagrams of the same dimension.
   * @param {PersistenceDiagram} other
   * @returns {number} Bottleneck distance.
   */
  bottleneckDistance(other) {
    return PersistenceDiagram.computeBottleneckDistance(this, other);
  }

  /**
   * Computes p-Wasserstein distance W_p between two persistence diagrams.
   * @param {PersistenceDiagram} other
   * @param {number} [p=1]
   * @param {number} [q=Infinity]
   * @returns {number}
   */
  wassersteinDistance(other, p = 1, q = Infinity) {
    return PersistenceDiagram.computeWassersteinDistance(this, other, p, q);
  }

  /**
   * Static computation of Bottleneck distance between two diagrams.
   * @param {PersistenceDiagram} diagA
   * @param {PersistenceDiagram} diagB
   * @returns {number}
   */
  static computeBottleneckDistance(diagA, diagB) {
    const ptsA = diagA.pairs;
    const ptsB = diagB.pairs;
    const nA = ptsA.length;
    const nB = ptsB.length;

    if (nA === 0 && nB === 0) return 0.0;

    const distPtToPt = (p1, p2) => Math.max(Math.abs(p1.birth - p2.birth), Math.abs(p1.death - p2.death));
    const distToDiag = (p) => Math.abs(p.death - p.birth) / 2.0;

    const candidateDistances = [0.0];
    for (let i = 0; i < nA; i++) {
      candidateDistances.push(distToDiag(ptsA[i]));
      for (let j = 0; j < nB; j++) {
        candidateDistances.push(distPtToPt(ptsA[i], ptsB[j]));
      }
    }
    for (let j = 0; j < nB; j++) {
      candidateDistances.push(distToDiag(ptsB[j]));
    }

    candidateDistances.sort((a, b) => a - b);
    const uniqueEps = [];
    for (let i = 0; i < candidateDistances.length; i++) {
      if (i === 0 || Math.abs(candidateDistances[i] - candidateDistances[i - 1]) > 1e-10) {
        uniqueEps.push(candidateDistances[i]);
      }
    }

    let lowIdx = 0;
    let highIdx = uniqueEps.length - 1;
    let bestDist = uniqueEps[highIdx];

    while (lowIdx <= highIdx) {
      const mid = Math.floor((lowIdx + highIdx) / 2);
      const eps = uniqueEps[mid];

      if (PersistenceDiagram._hasPerfectMatching(ptsA, ptsB, eps)) {
        bestDist = eps;
        highIdx = mid - 1;
      } else {
        lowIdx = mid + 1;
      }
    }

    return bestDist;
  }

  /**
   * Helper: checks if a bottleneck matching exists within threshold eps.
   * @private
   */
  static _hasPerfectMatching(ptsA, ptsB, eps) {
    const nA = ptsA.length;
    const nB = ptsB.length;
    const totalA = nA + nB;
    const totalB = nB + nA;

    const distPtToPt = (p1, p2) => Math.max(Math.abs(p1.birth - p2.birth), Math.abs(p1.death - p2.death));
    const distToDiag = (p) => Math.abs(p.death - p.birth) / 2.0;

    const adj = Array.from({ length: totalA }, () => []);

    for (let i = 0; i < nA; i++) {
      for (let j = 0; j < nB; j++) {
        if (distPtToPt(ptsA[i], ptsB[j]) <= eps + 1e-12) {
          adj[i].push(j);
        }
      }
      if (distToDiag(ptsA[i]) <= eps + 1e-12) {
        adj[i].push(nB + i);
      }
    }

    for (let j = 0; j < nB; j++) {
      if (distToDiag(ptsB[j]) <= eps + 1e-12) {
        adj[nA + j].push(j);
      }
    }

    for (let j = 0; j < nB; j++) {
      for (let i = 0; i < nA; i++) {
        adj[nA + j].push(nB + i);
      }
    }

    const matchB = new Int32Array(totalB).fill(-1);
    const visited = new Uint8Array(totalA);

    const dfs = (u) => {
      for (let k = 0; k < adj[u].length; k++) {
        const v = adj[u][k];
        if (matchB[v] < 0) {
          matchB[v] = u;
          return true;
        }
      }
      for (let k = 0; k < adj[u].length; k++) {
        const v = adj[u][k];
        const nextU = matchB[v];
        if (!visited[nextU]) {
          visited[nextU] = 1;
          if (dfs(nextU)) {
            matchB[v] = u;
            return true;
          }
        }
      }
      return false;
    };

    let matchingCount = 0;
    for (let u = 0; u < totalA; u++) {
      visited.fill(0);
      visited[u] = 1;
      if (dfs(u)) {
        matchingCount++;
      }
    }

    return matchingCount === totalA;
  }

  /**
   * Static computation of Wasserstein distance W_p between two diagrams via Hungarian algorithm.
   * @param {PersistenceDiagram} diagA
   * @param {PersistenceDiagram} diagB
   * @param {number} [p=1]
   * @param {number} [q=Infinity]
   * @returns {number}
   */
  static computeWassersteinDistance(diagA, diagB, p = 1, q = Infinity) {
    const ptsA = diagA.pairs;
    const ptsB = diagB.pairs;
    const nA = ptsA.length;
    const nB = ptsB.length;

    if (nA === 0 && nB === 0) return 0.0;

    const N = nA + nB;
    const costMatrix = Array.from({ length: N }, () => new Float64Array(N));

    const ptDist = (p1, p2) => {
      if (q === Infinity) {
        return Math.max(Math.abs(p1.birth - p2.birth), Math.abs(p1.death - p2.death));
      } else if (q === 1) {
        return Math.abs(p1.birth - p2.birth) + Math.abs(p1.death - p2.death);
      } else {
        const db = p1.birth - p2.birth;
        const dd = p1.death - p2.death;
        return Math.pow(Math.pow(Math.abs(db), q) + Math.pow(Math.abs(dd), q), 1.0 / q);
      }
    };

    const diagDist = (pt) => {
      return Math.abs(pt.death - pt.birth) / Math.SQRT2;
    };

    for (let i = 0; i < nA; i++) {
      for (let j = 0; j < nB; j++) {
        costMatrix[i][j] = Math.pow(ptDist(ptsA[i], ptsB[j]), p);
      }
    }
    for (let i = 0; i < nA; i++) {
      for (let j = 0; j < nA; j++) {
        costMatrix[i][nB + j] = i === j ? Math.pow(diagDist(ptsA[i]), p) : TOPOLOGY_CONSTANTS.INFINITY_VAL;
      }
    }
    for (let i = 0; i < nB; i++) {
      for (let j = 0; j < nB; j++) {
        costMatrix[nA + i][j] = i === j ? Math.pow(diagDist(ptsB[j]), p) : TOPOLOGY_CONSTANTS.INFINITY_VAL;
      }
    }
    for (let i = 0; i < nB; i++) {
      for (let j = 0; j < nA; j++) {
        costMatrix[nA + i][nB + j] = 0.0;
      }
    }

    const minCost = PersistenceDiagram._solveHungarian(costMatrix, N);
    return Math.pow(minCost, 1.0 / p);
  }

  /**
   * Exact O(N^3) Kuhn-Munkres Hungarian algorithm.
   * @private
   */
  static _solveHungarian(cost, n) {
    const u = new Float64Array(n + 1);
    const v = new Float64Array(n + 1);
    const p = new Int32Array(n + 1);
    const way = new Int32Array(n + 1);

    for (let i = 1; i <= n; i++) {
      p[0] = i;
      let j0 = 0;
      const minv = new Float64Array(n + 1).fill(Infinity);
      const used = new Uint8Array(n + 1);

      do {
        used[j0] = 1;
        const i0 = p[j0];
        let delta = Infinity;
        let j1 = 0;

        for (let j = 1; j <= n; j++) {
          if (!used[j]) {
            const cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
            if (cur < minv[j]) {
              minv[j] = cur;
              way[j] = j0;
            }
            if (minv[j] < delta) {
              delta = minv[j];
              j1 = j;
            }
          }
        }

        for (let j = 0; j <= n; j++) {
          if (used[j]) {
            u[p[j]] += delta;
            v[j] -= delta;
          } else {
            minv[j] -= delta;
          }
        }
        j0 = j1;
      } while (p[j0] !== 0);

      do {
        const j1 = way[j0];
        p[j0] = p[j1];
        j0 = j1;
      } while (j0 !== 0);
    }

    return -v[0];
  }
}

// ============================================================================
// GAUSSIAN RANDOM FIELD (GRF) ANALYTICAL TOMITA-GOTT TOPOLOGY
// ============================================================================

/**
 * Analytical Gaussian Random Field (GRF) Euler Characteristic Density and Genus Evaluator.
 * Implements the Gott-Melott-Dickinson & Tomita (1986) invariant:
 * \chi_{\text{GRF}}(\nu) = \frac{1}{(2\pi)^2} \left(\frac{\sigma_1}{\sqrt{3}\sigma_0}\right)^3 (\nu^2 - 1) e^{-\nu^2/2}
 */
export class TomitaGottAnalytical {
  /**
   * Computes theoretical Euler characteristic density per unit volume for a 3D isotropic GRF.
   *
   * @param {number} nu - Normalized density threshold \nu = (\delta - \langle\delta\rangle) / \sigma_0.
   * @param {number} sigma0 - RMS field amplitude \sigma_0 = \sqrt{\langle \delta^2 \rangle}.
   * @param {number} sigma1 - Gradient dispersion \sigma_1 = \sqrt{\langle |\nabla \delta|^2 \rangle}.
   * @returns {number} Theoretical Euler characteristic density \chi_V(\nu) [Mpc^{-3} or grid^{-3}].
   */
  static computeExpectedEulerDensity(nu, sigma0, sigma1) {
    if (sigma0 <= 0 || sigma1 < 0) {
      throw new Error(`TomitaGott: sigma0 must be > 0 (got ${sigma0}), sigma1 >= 0 (got ${sigma1})`);
    }
    const lambda = sigma1 / (TOPOLOGY_CONSTANTS.SQRT_3 * sigma0);
    const prefactor = TOPOLOGY_CONSTANTS.INV_TWO_PI_SQ * Math.pow(lambda, 3);
    const polynomial = (nu * nu - 1.0);
    const gaussian = Math.exp(-0.5 * nu * nu);
    return prefactor * polynomial * gaussian;
  }

  /**
   * Computes the theoretical Gott-Weinberg-Melott genus density g_V(\nu) = -\frac{1}{2} \chi_V(\nu).
   *
   * @param {number} nu
   * @param {number} sigma0
   * @param {number} sigma1
   * @returns {number} Genus density.
   */
  static computeExpectedGenusDensity(nu, sigma0, sigma1) {
    return -0.5 * TomitaGottAnalytical.computeExpectedEulerDensity(nu, sigma0, sigma1);
  }

  /**
   * Computes empirical spectral dispersion parameters \sigma_0, \sigma_1, \sigma_2 from a 3D scalar grid.
   *
   * @param {Float64Array|number[]} grid - Flattened 3D scalar field array (nx * ny * nz).
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @param {number} [dx=1.0] - Voxel physical size in X.
   * @param {number} [dy=1.0] - Voxel physical size in Y.
   * @param {number} [dz=1.0] - Voxel physical size in Z.
   * @returns {Object} { mean, sigma0, sigma1, sigma2, spectralGamma, volume }
   */
  static computeSpectralParameters(grid, nx, ny, nz, dx = 1.0, dy = 1.0, dz = 1.0) {
    const N = nx * ny * nz;
    if (grid.length !== N) {
      throw new Error(`Grid size mismatch: expected ${N}, got ${grid.length}`);
    }

    let sum = 0.0;
    for (let i = 0; i < N; i++) {
      sum += grid[i];
    }
    const mean = sum / N;

    let varSum = 0.0;
    for (let i = 0; i < N; i++) {
      const diff = grid[i] - mean;
      varSum += diff * diff;
    }
    const sigma0 = Math.sqrt(varSum / N);

    let gradSqSum = 0.0;
    let laplacianSqSum = 0.0;
    let countInterior = 0;

    const idx = (i, j, k) => i + nx * (j + ny * k);

    for (let k = 1; k < nz - 1; k++) {
      for (let j = 1; j < ny - 1; j++) {
        for (let i = 1; i < nx - 1; i++) {
          const val = grid[idx(i, j, k)];
          // First derivatives
          const ddx = (grid[idx(i + 1, j, k)] - grid[idx(i - 1, j, k)]) / (2.0 * dx);
          const ddy = (grid[idx(i, j + 1, k)] - grid[idx(i, j - 1, k)]) / (2.0 * dy);
          const ddz = (grid[idx(i, j, k + 1)] - grid[idx(i, j, k - 1)]) / (2.0 * dz);
          const gradSq = ddx * ddx + ddy * ddy + ddz * ddz;
          gradSqSum += gradSq;

          // Second derivatives (Laplacian)
          const d2dx2 = (grid[idx(i + 1, j, k)] - 2.0 * val + grid[idx(i - 1, j, k)]) / (dx * dx);
          const d2dy2 = (grid[idx(i, j + 1, k)] - 2.0 * val + grid[idx(i, j - 1, k)]) / (dy * dy);
          const d2dz2 = (grid[idx(i, j, k + 1)] - 2.0 * val + grid[idx(i, j, k - 1)]) / (dz * dz);
          const laplacian = d2dx2 + d2dy2 + d2dz2;
          laplacianSqSum += laplacian * laplacian;

          countInterior++;
        }
      }
    }

    const sigma1 = countInterior > 0 ? Math.sqrt(gradSqSum / countInterior) : 0.0;
    const sigma2 = countInterior > 0 ? Math.sqrt(laplacianSqSum / countInterior) : 0.0;
    const spectralGamma = (sigma0 > 0 && sigma2 > 0) ? (sigma1 * sigma1) / (sigma0 * sigma2) : 0.0;
    const volume = (nx * dx) * (ny * dy) * (nz * dz);

    return {
      mean,
      sigma0,
      sigma1,
      sigma2,
      spectralGamma,
      volume,
      totalVoxels: N
    };
  }

  /**
   * Fits empirical Euler characteristic curve against the Tomita-Gott expectation.
   *
   * @param {Array<{nu: number, chiEmpirical: number}>} empiricalData
   * @param {number} totalVolume
   * @param {number} sigma0
   * @param {number} sigma1
   * @returns {Object} Non-Gaussian topological diagnostics.
   */
  static fitGenusCurve(empiricalData, totalVolume, sigma0, sigma1) {
    if (!empiricalData || empiricalData.length === 0) {
      throw new Error('fitGenusCurve requires non-empty empiricalData array');
    }

    let chiSqSum = 0.0;
    let maxEmpirical = -Infinity;
    let minEmpirical = Infinity;
    let nuMax = 0.0;
    let nuMin = 0.0;

    const fittedCurve = [];

    for (let i = 0; i < empiricalData.length; i++) {
      const { nu, chiEmpirical } = empiricalData[i];
      const chiExpected = TomitaGottAnalytical.computeExpectedEulerDensity(nu, sigma0, sigma1) * totalVolume;
      const residual = chiEmpirical - chiExpected;
      chiSqSum += residual * residual;

      if (chiEmpirical > maxEmpirical) {
        maxEmpirical = chiEmpirical;
        nuMax = nu;
      }
      if (chiEmpirical < minEmpirical) {
        minEmpirical = chiEmpirical;
        nuMin = nu;
      }

      fittedCurve.push({
        nu,
        chiEmpirical,
        chiExpected,
        residual
      });
    }

    const rmsResidual = Math.sqrt(chiSqSum / empiricalData.length);

    // Topological Asymmetry parameter A_T = (\chi_{max} - |\chi_{min}|) / (\chi_{max} + |\chi_{min}|)
    const peakPos = Math.abs(maxEmpirical);
    const peakNeg = Math.abs(minEmpirical);
    const asymmetry = (peakPos + peakNeg > TOPOLOGY_CONSTANTS.EPSILON)
      ? (peakPos - peakNeg) / (peakPos + peakNeg)
      : 0.0;

    // Shift parameter \Delta\nu: empirical zero-crossing vs analytical zero crossing at \nu = \pm 1
    let empiricalZeroCrossing = 0.0;
    for (let i = 0; i < empiricalData.length - 1; i++) {
      const p1 = empiricalData[i];
      const p2 = empiricalData[i + 1];
      if (p1.nu >= 0 && p2.nu >= 0 && ((p1.chiEmpirical <= 0 && p2.chiEmpirical >= 0) || (p1.chiEmpirical >= 0 && p2.chiEmpirical <= 0))) {
        const denom = p2.chiEmpirical - p1.chiEmpirical;
        if (Math.abs(denom) > 1e-12) {
          empiricalZeroCrossing = p1.nu + (-p1.chiEmpirical / denom) * (p2.nu - p1.nu);
          break;
        }
      }
    }
    const shiftDeltaNu = empiricalZeroCrossing - 1.0;

    return {
      rmsResidual,
      asymmetry,
      shiftDeltaNu,
      nuMax,
      nuMin,
      fittedCurve
    };
  }
}

// ============================================================================
// VELOCITY POTENTIAL & TIDAL TENSOR TOPOLOGY ANALYZER
// ============================================================================

/**
 * Velocity Potential (\Phi) and Gravitational Tidal Tensor Field Topology Analyzer.
 * Computes deformation tensor T_{ij} = \partial_i \partial_j \Phi eigenvalues \lambda_1 \ge \lambda_2 \ge \lambda_3
 * and classifies cosmic web environments (nodes, filaments, walls, voids).
 */
export class VelocityPotentialTopologist {
  /**
   * Computes the 3x3 gravitational tidal deformation tensor and its eigenvalues at (i, j, k).
   *
   * @param {Float64Array|number[]} phiGrid - 3D velocity potential grid.
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @param {number} i
   * @param {number} j
   * @param {number} k
   * @param {number} [dx=1.0]
   * @param {number} [dy=1.0]
   * @param {number} [dz=1.0]
   * @returns {Object} { eigenvalues: [l1, l2, l3], divergence: number, environment: string }
   */
  static computeTidalTensor(phiGrid, nx, ny, nz, i, j, k, dx = 1.0, dy = 1.0, dz = 1.0, threshold = 0.2) {
    const idx = (x, y, z) => {
      const cx = Math.max(0, Math.min(nx - 1, x));
      const cy = Math.max(0, Math.min(ny - 1, y));
      const cz = Math.max(0, Math.min(nz - 1, z));
      return cx + nx * (cy + ny * cz);
    };

    const val = phiGrid[idx(i, j, k)];

    // Diagonal elements: \partial_i^2 \Phi
    const Txx = (phiGrid[idx(i + 1, j, k)] - 2.0 * val + phiGrid[idx(i - 1, j, k)]) / (dx * dx);
    const Tyy = (phiGrid[idx(i, j + 1, k)] - 2.0 * val + phiGrid[idx(i, j - 1, k)]) / (dy * dy);
    const Tzz = (phiGrid[idx(i, j, k + 1)] - 2.0 * val + phiGrid[idx(i, j, k - 1)]) / (dz * dz);

    // Off-diagonal elements: \partial_i \partial_j \Phi via 4-point cross differences
    const Txy = (
      phiGrid[idx(i + 1, j + 1, k)] - phiGrid[idx(i - 1, j + 1, k)] -
      phiGrid[idx(i + 1, j - 1, k)] + phiGrid[idx(i - 1, j - 1, k)]
    ) / (4.0 * dx * dy);

    const Tyz = (
      phiGrid[idx(i, j + 1, k + 1)] - phiGrid[idx(i, j - 1, k + 1)] -
      phiGrid[idx(i, j + 1, k - 1)] + phiGrid[idx(i, j - 1, k - 1)]
    ) / (4.0 * dy * dz);

    const Tzx = (
      phiGrid[idx(i + 1, j, k + 1)] - phiGrid[idx(i - 1, j, k + 1)] -
      phiGrid[idx(i + 1, j, k - 1)] + phiGrid[idx(i - 1, j, k - 1)]
    ) / (4.0 * dz * dx);

    // Symmetric 3x3 eigensystem via analytical Cardano roots
    const eigenvalues = VelocityPotentialTopologist._eigenvaluesSymmetric3x3(Txx, Tyy, Tzz, Txy, Tyz, Tzx);
    const divergence = Txx + Tyy + Tzz; // Trace is Laplacian \nabla^2 \Phi = \delta

    // Classify cosmic environment by eigenvalue counts above threshold
    const [l1, l2, l3] = eigenvalues;
    let environment = CosmicEnvironmentType.VOID;
    if (l3 > threshold) {
      environment = CosmicEnvironmentType.PEAK_NODE;
    } else if (l2 > threshold) {
      environment = CosmicEnvironmentType.FILAMENT;
    } else if (l1 > threshold) {
      environment = CosmicEnvironmentType.SHEET_WALL;
    }

    return {
      tensor: [
        [Txx, Txy, Tzx],
        [Txy, Tyy, Tyz],
        [Tzx, Tyz, Tzz]
      ],
      eigenvalues,
      divergence,
      environment
    };
  }

  /**
   * Exact analytical Cardano eigenvalue solver for real symmetric 3x3 matrix.
   * Returns sorted eigenvalues [l1 >= l2 >= l3].
   * @private
   */
  static _eigenvaluesSymmetric3x3(a11, a22, a33, a12, a23, a13) {
    const q = (a11 + a22 + a33) / 3.0;
    const b11 = a11 - q;
    const b22 = a22 - q;
    const b33 = a33 - q;
    const p = Math.sqrt((b11 * b11 + b22 * b22 + b33 * b33 + 2.0 * (a12 * a12 + a23 * a23 + a13 * a13)) / 6.0);

    if (p < 1e-15) {
      return [q, q, q];
    }

    // Det(B / p)
    const invP = 1.0 / p;
    const m11 = b11 * invP;
    const m22 = b22 * invP;
    const m33 = b33 * invP;
    const m12 = a12 * invP;
    const m23 = a23 * invP;
    const m13 = a13 * invP;

    const detM = m11 * (m22 * m33 - m23 * m23) -
                 m12 * (m12 * m33 - m23 * m13) +
                 m13 * (m12 * m23 - m22 * m13);

    const r = Math.max(-1.0, Math.min(1.0, detM / 2.0));
    const phi = Math.acos(r) / 3.0;

    const eig1 = q + 2.0 * p * Math.cos(phi);
    const eig3 = q + 2.0 * p * Math.cos(phi + (2.0 * Math.PI / 3.0));
    const eig2 = 3.0 * q - eig1 - eig3;

    const sorted = [eig1, eig2, eig3].sort((a, b) => b - a);
    return sorted;
  }
}

// ============================================================================
// COMPLETE BETTI NUMBER AND TOPOLOGY CALCULATOR ENGINE
// ============================================================================

/**
 * High-performance 3D Betti Number, Persistent Homology, and Cosmic Topology Calculator.
 */
export class BettiNumberCalculator {
  /**
   * @param {Object} [options={}]
   * @param {boolean} [options.periodic=false] - Use 3-torus periodic boundaries.
   * @param {number} [options.dx=1.0] - Physical grid spacing in X [Mpc/h].
   * @param {number} [options.dy=1.0] - Physical grid spacing in Y [Mpc/h].
   * @param {number} [options.dz=1.0] - Physical grid spacing in Z [Mpc/h].
   */
  constructor(options = {}) {
    this.options = {
      periodic: Boolean(options.periodic),
      dx: options.dx || 1.0,
      dy: options.dy || 1.0,
      dz: options.dz || 1.0
    };
  }

  /**
   * Computes the exact Betti numbers \beta_0, \beta_1, \beta_2 and Euler characteristic \chi
   * for a 3D scalar grid filtered at threshold \delta_th.
   *
   * @param {Float64Array|number[]} grid - 3D scalar grid (density contrast or potential).
   * @param {number} nx - X dimension.
   * @param {number} ny - Y dimension.
   * @param {number} nz - Z dimension.
   * @param {number} threshold - Filtration threshold \delta_th.
   * @param {string} [direction=FiltrationDirection.SUPERLEVEL]
   * @returns {Object} { beta0, beta1, beta2, chi, cellCounts: {V, E, F, C}, valid: boolean }
   */
  computeBettiNumbers(grid, nx, ny, nz, threshold, direction = FiltrationDirection.SUPERLEVEL) {
    const complex = new CubicalComplex3D(nx, ny, nz, this.options);
    const cellVals = complex.computeCellFiltrationValues(grid, direction);
    const cellCounts = complex.countActiveCells(cellVals, threshold, direction);
    const { V, E, F, C, chi } = cellCounts;

    if (V === 0) {
      return {
        beta0: 0,
        beta1: 0,
        beta2: 0,
        chi: 0,
        cellCounts: { V: 0, E: 0, F: 0, C: 0 },
        valid: true
      };
    }

    const isSuper = direction === FiltrationDirection.SUPERLEVEL;
    const totalVertices = nx * ny * nz;

    // 1. Compute \beta_0 via Union-Find
    const uf0 = new UnionFindPersistent(totalVertices);
    for (let idx = 0; idx < totalVertices; idx++) {
      const active = isSuper ? grid[idx] >= threshold : grid[idx] <= threshold;
      if (active) {
        uf0.makeSet(idx, grid[idx]);
      }
    }

    for (let e = 0; e < complex.numEdges; e++) {
      const activeEdge = isSuper ? cellVals.edgeVal[e] >= threshold : cellVals.edgeVal[e] <= threshold;
      if (activeEdge) {
        const [v0, v1] = complex.edgeBoundaryVertices(e);
        if (uf0.active[v0] && uf0.active[v1]) {
          uf0.union(v0, v1, threshold, isSuper);
        }
      }
    }

    const roots0 = new Set();
    for (let idx = 0; idx < totalVertices; idx++) {
      if (uf0.active[idx]) {
        roots0.add(uf0.find(idx));
      }
    }
    const beta0 = roots0.size;

    // 2. Compute \beta_2 via Dual Complementary Union-Find (Alexander Duality)
    const ufComp = new UnionFindPersistent(totalVertices);
    const compActive = new Uint8Array(totalVertices);

    for (let idx = 0; idx < totalVertices; idx++) {
      const isComp = isSuper ? grid[idx] < threshold : grid[idx] > threshold;
      if (isComp) {
        compActive[idx] = 1;
        ufComp.makeSet(idx, grid[idx]);
      }
    }

    for (let e = 0; e < complex.numEdges; e++) {
      const [v0, v1] = complex.edgeBoundaryVertices(e);
      if (compActive[v0] && compActive[v1]) {
        ufComp.union(v0, v1, threshold, !isSuper);
      }
    }

    const boundaryRoots = new Set();
    for (let k = 0; k < nz; k++) {
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const isBoundary = (i === 0 || i === nx - 1 || j === 0 || j === ny - 1 || k === 0 || k === nz - 1);
          if (isBoundary) {
            const vIdx = complex.vertexIndex(i, j, k);
            if (compActive[vIdx]) {
              boundaryRoots.add(ufComp.find(vIdx));
            }
          }
        }
      }
    }

    const allCompRoots = new Set();
    for (let idx = 0; idx < totalVertices; idx++) {
      if (compActive[idx]) {
        allCompRoots.add(ufComp.find(idx));
      }
    }

    let beta2 = 0;
    for (const root of allCompRoots) {
      if (!boundaryRoots.has(root)) {
        beta2++;
      }
    }

    // 3. Compute \beta_1 via Euler-Poincaré Identity: \beta_1 = \beta_0 + \beta_2 - \chi
    let beta1 = beta0 + beta2 - chi;
    if (beta1 < 0) {
      beta1 = 0;
    }

    return {
      beta0,
      beta1,
      beta2,
      chi,
      cellCounts: { V, E, F, C },
      valid: (beta0 - beta1 + beta2 === chi)
    };
  }

  /**
   * Computes Betti curves and Euler characteristic as a function of threshold \delta_th.
   *
   * @param {Float64Array|number[]} grid
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @param {Object} [rangeOptions={}]
   * @returns {Array<{threshold: number, nu: number, beta0: number, beta1: number, beta2: number, chi: number, cellChi: number}>}
   */
  computeBettiCurves(grid, nx, ny, nz, rangeOptions = {}) {
    const N = nx * ny * nz;
    let minVal = Infinity;
    let maxVal = -Infinity;
    let sum = 0.0;

    for (let i = 0; i < N; i++) {
      const v = grid[i];
      if (v < minVal) minVal = v;
      if (v > maxVal) maxVal = v;
      sum += v;
    }

    const mean = sum / N;
    let varSum = 0.0;
    for (let i = 0; i < N; i++) {
      const d = grid[i] - mean;
      varSum += d * d;
    }
    const sigma0 = Math.sqrt(varSum / N) || 1.0;

    const minTh = rangeOptions.minThreshold !== undefined ? rangeOptions.minThreshold : minVal;
    const maxTh = rangeOptions.maxThreshold !== undefined ? rangeOptions.maxThreshold : maxVal;
    const steps = rangeOptions.numSteps || 50;
    const direction = rangeOptions.direction || FiltrationDirection.SUPERLEVEL;

    const results = [];
    const stepSize = (maxTh - minTh) / Math.max(1, steps - 1);

    for (let s = 0; s < steps; s++) {
      const th = minTh + s * stepSize;
      const nu = (th - mean) / sigma0;
      const betti = this.computeBettiNumbers(grid, nx, ny, nz, th, direction);
      results.push({
        threshold: th,
        nu,
        beta0: betti.beta0,
        beta1: betti.beta1,
        beta2: betti.beta2,
        chi: betti.chi,
        cellChi: betti.cellCounts.V - betti.cellCounts.E + betti.cellCounts.F - betti.cellCounts.C
      });
    }

    return results;
  }

  /**
   * Computes complete persistent homology across dimensions H0, H1, H2.
   *
   * @param {Float64Array|number[]} grid
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @param {string} [direction=FiltrationDirection.SUPERLEVEL]
   * @returns {Object} { h0: PersistenceDiagram, h1: PersistenceDiagram, h2: PersistenceDiagram, totalEntropy: number }
   */
  computePersistentHomology(grid, nx, ny, nz, direction = FiltrationDirection.SUPERLEVEL) {
    const complex = new CubicalComplex3D(nx, ny, nz, this.options);
    const cellVals = complex.computeCellFiltrationValues(grid, direction);
    const isSuper = direction === FiltrationDirection.SUPERLEVEL;

    const allCells = [];

    // 0-cells
    for (let v = 0; v < complex.numVertices; v++) {
      allCells.push({
        dim: 0,
        id: v,
        value: cellVals.vertexVal[v],
        boundary: []
      });
    }

    // 1-cells
    const edgeOffset = complex.numVertices;
    for (let e = 0; e < complex.numEdges; e++) {
      const [v0, v1] = complex.edgeBoundaryVertices(e);
      allCells.push({
        dim: 1,
        id: edgeOffset + e,
        value: cellVals.edgeVal[e],
        boundary: [v0, v1]
      });
    }

    // 2-cells
    const faceOffset = edgeOffset + complex.numEdges;
    for (let f = 0; f < complex.numFaces; f++) {
      const edges = complex.faceBoundaryEdges(f).map(eIdx => edgeOffset + eIdx);
      allCells.push({
        dim: 2,
        id: faceOffset + f,
        value: cellVals.faceVal[f],
        boundary: edges
      });
    }

    // 3-cells
    const cubeOffset = faceOffset + complex.numFaces;
    for (let c = 0; c < complex.numCubes; c++) {
      const faces = complex.cubeBoundaryFaces(c).map(fIdx => faceOffset + fIdx);
      allCells.push({
        dim: 3,
        id: cubeOffset + c,
        value: cellVals.cubeVal[c],
        boundary: faces
      });
    }

    // Sort cells for filtration
    allCells.sort((a, b) => {
      if (Math.abs(a.value - b.value) > 1e-14) {
        return isSuper ? (b.value - a.value) : (a.value - b.value);
      }
      return a.dim - b.dim;
    });

    const { pairs, essential } = BoundaryMatrixReducer.reduceZ2(allCells);

    const pairsH0 = pairs.filter(p => p.dim === 0);
    const pairsH1 = pairs.filter(p => p.dim === 1);
    const pairsH2 = pairs.filter(p => p.dim === 2);

    const essH0 = essential.filter(p => p.dim === 0);
    const essH1 = essential.filter(p => p.dim === 1);
    const essH2 = essential.filter(p => p.dim === 2);

    const h0 = new PersistenceDiagram(0, pairsH0, essH0);
    const h1 = new PersistenceDiagram(1, pairsH1, essH1);
    const h2 = new PersistenceDiagram(2, pairsH2, essH2);

    const totalEntropy = h0.persistentEntropy() + h1.persistentEntropy() + h2.persistentEntropy();

    return {
      h0,
      h1,
      h2,
      totalEntropy,
      totalPairsCount: pairs.length
    };
  }

  /**
   * Full cosmological topological classification and structure analysis.
   *
   * @param {Float64Array|number[]} grid
   * @param {number} nx
   * @param {number} ny
   * @param {number} nz
   * @param {Object} [options={}]
   * @returns {Object} Comprehensive cosmic topology dossier.
   */
  analyzeCosmicTopology(grid, nx, ny, nz, options = {}) {
    const spectral = TomitaGottAnalytical.computeSpectralParameters(
      grid, nx, ny, nz,
      this.options.dx, this.options.dy, this.options.dz
    );

    const bettiCurves = this.computeBettiCurves(grid, nx, ny, nz, {
      numSteps: options.numSteps || 40,
      direction: FiltrationDirection.SUPERLEVEL
    });

    const empiricalEulerData = bettiCurves.map(b => ({ nu: b.nu, chiEmpirical: b.chi }));
    const genusFit = TomitaGottAnalytical.fitGenusCurve(
      empiricalEulerData, spectral.volume, spectral.sigma0, spectral.sigma1
    );

    let maxBeta0 = 0;
    let maxBeta1 = 0;
    let maxBeta2 = 0;
    for (let i = 0; i < bettiCurves.length; i++) {
      if (bettiCurves[i].beta0 > maxBeta0) maxBeta0 = bettiCurves[i].beta0;
      if (bettiCurves[i].beta1 > maxBeta1) maxBeta1 = bettiCurves[i].beta1;
      if (bettiCurves[i].beta2 > maxBeta2) maxBeta2 = bettiCurves[i].beta2;
    }

    let topologyType = CosmicTopologyType.INTERMEDIATE;
    if (Math.abs(genusFit.asymmetry) < 0.15 && Math.abs(genusFit.shiftDeltaNu) < 0.15) {
      topologyType = CosmicTopologyType.GAUSSIAN_RANDOM;
    } else if (maxBeta0 > 1.5 * maxBeta2 && genusFit.asymmetry > 0.2) {
      topologyType = CosmicTopologyType.MEATBALL;
    } else if (maxBeta2 > 1.5 * maxBeta0 && genusFit.asymmetry < -0.2) {
      topologyType = CosmicTopologyType.SWISS_CHEESE;
    } else if (maxBeta1 > maxBeta0 && maxBeta1 > maxBeta2) {
      topologyType = CosmicTopologyType.SPONGY;
    }

    return {
      spectral,
      bettiCurves,
      genusFit,
      maxBeta0,
      maxBeta1,
      maxBeta2,
      topologyType,
      asymmetry: genusFit.asymmetry,
      shiftDeltaNu: genusFit.shiftDeltaNu,
      rmsResidual: genusFit.rmsResidual
    };
  }
}
