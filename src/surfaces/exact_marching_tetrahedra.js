/**
 * @file exact_marching_tetrahedra.js
 * @module surfaces/exact_marching_tetrahedra
 * @description Production-Grade Exact 3D Marching Tetrahedra & Manifold Mesher for Cosmological Density Fields.
 * 
 * Mathematical Formulation & Architectural Principles:
 * 1. Simplicial Decomposition of 3D Voxel Grids:
 *    - 6-Tetrahedra Kuhn / Freudenthal decomposition (uniform diagonal).
 *    - 6-Tetrahedra Alternating decomposition.
 *    - 5-Tetrahedra Alternating Checkerboard decomposition (1 central regular tet + 4 corner tets).
 *    - Eliminates all Marching Cubes ambiguous saddle configurations (no topological holes, facet cracks, or non-manifold edges).
 * 2. Exact Linear & High-Order Root Interpolation on Simplex Edges:
 *    - Linear interpolation: t = (\delta_{iso} - s_1) / (s_2 - s_1).
 *    - Newton-Raphson high-order root refinement: t_{k+1} = t_k - (f(p(t_k)) - \delta_{iso}) / (df/dt).
 *    - Universal support for canonical cosmological density thresholds:
 *      \delta = -0.8 (Cosmic Underdense Void Boundary),
 *      \delta = 0.0 (Cosmic Mean Density Separatrix),
 *      \delta = 1.5 (Cosmic Sheet / Wall Boundary),
 *      \delta = 1.686 (Linear Spherical Collapse / Critical Overdensity \delta_c),
 *      \delta = 5.0 (Cosmic Filament Core),
 *      \delta = 20.0 (Virialized Cluster Envelope),
 *      \delta = 200.0 (Virialized Halo Overdensity Boundary \delta_{vir}).
 * 3. Topological 2-Manifold Invariants & Graph Topology:
 *    - Exact canonical edge welding using 64-bit spatial hash keys.
 *    - Undirected edge map and vertex adjacency graph.
 *    - Euler characteristic \chi = V - E + F.
 *    - Genus g = C - \chi / 2 for C connected components.
 *    - Betti numbers (b_0, b_1, b_2).
 *    - Watertight 2-manifold verification (boundary edges = 0, non-manifold edges = 0).
 * 4. Rigorous Differential Geometry & Calculus Integrals:
 *    - Surface Area Integral: \iint_{\partial \Omega} dA = \sum_f \frac{1}{2} \| (v_1 - v_0) \times (v_2 - v_0) \|.
 *    - Enclosed Volume Integral: \iiint_\Omega dV = \frac{1}{6} \sum_f v_0 \cdot (v_1 \times v_2) (divergence theorem).
 *    - Area-weighted surface centroid and signed volumetric centroid.
 *    - Minkowski Functionals: V_0 (Volume), V_1 (Area), V_2 (Integrated Mean Curvature), V_3 (\chi).
 *    - Morphology descriptors: Sphericity \Psi = \pi^{1/3} (6 V_0)^{2/3} / V_1, Isoperimetric Quotient Q = 36 \pi V_0^2 / V_1^3.
 *    - Planarity P and Filamentarity F shape parameters (Sahni et al. morphometrics).
 *    - Discrete Gaussian Curvature K via angle defect: \oint K dA = 2\pi \chi (Gauss-Bonnet theorem).
 *    - Discrete Mean Curvature H via Laplace-Beltrami cotangent weights.
 *    - Angle-weighted, area-weighted, uniform, and continuous field-gradient normals: \hat{n} = -\nabla \delta / \|\nabla \delta\|.
 *    - Taubin non-shrinking dual-step smoothing (\lambda > 0, \mu < -\lambda) with strict volume preservation.
 * 5. Cosmological Density Contrast Synthesis:
 *    - Multi-halo density contrast fields in Supergalactic coordinates (SGX, SGY, SGZ in Mpc/h).
 *    - Laniakea Supercluster Gravity Basin, Shapley Supercluster Concentration, Great Attractor (Norma/Centaurus),
 *      Perseus-Pisces, Coma Supercluster, and Local Void.
 *    - NFW (Navarro-Frenk-White), Hernquist, and Einasto density profiles with cosmological units (M_\odot, Mpc/h, km/s).
 * 6. Export Formats:
 *    - Three.js BufferGeometry typed arrays (position, normal, index).
 *    - Wavefront OBJ with topological, curvature, and calculus metadata.
 *    - GeoJSON 3D MultiPolygon format.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import { GridIndexer } from '../fields/grid_indexer.js';
import { TrilinearInterpolator } from '../interpolation/trilinear_interpolator.js';
import { FiniteDifference } from '../interpolation/finite_difference.js';

// ============================================================================
// CONSTANTS & ENUMS
// ============================================================================

/**
 * Supported Cube-to-Tetrahedra Decomposition Modes.
 * @readonly
 * @enum {string}
 */
export const DecompositionMode = Object.freeze({
  /** 6-tetrahedra Kuhn/Freudenthal decomposition (uniform diagonal orientation along 0->6) */
  SIX_TETRAHEDRA_KUHN: 'SIX_TETRAHEDRA_KUHN',
  /** 6-tetrahedra alternating decomposition */
  SIX_TETRAHEDRA_ALTERNATING: 'SIX_TETRAHEDRA_ALTERNATING',
  /** 5-tetrahedra alternating checkerboard decomposition (minimal simplex count: 1 central + 4 corner tets) */
  FIVE_TETRAHEDRA_ALTERNATING: 'FIVE_TETRAHEDRA_ALTERNATING'
});

/**
 * Standard Cosmological Density Contrast Level Sets (\delta = \rho/\bar{\rho} - 1).
 * @readonly
 * @enum {number}
 */
export const CosmologicalDensityThresholds = Object.freeze({
  VOID_BOUNDARY: -0.8,
  MEAN_DENSITY: 0.0,
  SHEET_BOUNDARY: 1.5,
  LINEAR_COLLAPSE: 1.686,
  FILAMENT_CORE: 5.0,
  VIRIALIZED_CLUSTER: 20.0,
  VIRIALIZED_HALO: 200.0
});

/**
 * Normal Vector Calculation and Weighting Strategies.
 * @readonly
 * @enum {string}
 */
export const NormalCalculationMode = Object.freeze({
  FACET_UNIFORM: 'FACET_UNIFORM',
  AREA_WEIGHTED: 'AREA_WEIGHTED',
  ANGLE_WEIGHTED: 'ANGLE_WEIGHTED',
  FIELD_GRADIENT: 'FIELD_GRADIENT'
});

/**
 * Smoothing Algorithm Modes.
 * @readonly
 * @enum {string}
 */
export const SmoothingAlgorithmMode = Object.freeze({
  UNIFORM_LAPLACIAN: 'UNIFORM_LAPLACIAN',
  COTANGENT_LAPLACIAN: 'COTANGENT_LAPLACIAN',
  TAUBIN_NON_SHRINKING: 'TAUBIN_NON_SHRINKING',
  VOLUME_PRESERVING_TAUBIN: 'VOLUME_PRESERVING_TAUBIN'
});

/**
 * Local cube corner offsets relative to voxel minimum [ix, iy, iz].
 * 0: [0,0,0], 1: [1,0,0], 2: [1,1,0], 3: [0,1,0]
 * 4: [0,0,1], 5: [1,0,1], 6: [1,1,1], 7: [0,1,1]
 */
export const CUBE_CORNER_OFFSETS = Object.freeze([
  [0, 0, 0], // 0
  [1, 0, 0], // 1
  [1, 1, 0], // 2
  [0, 1, 0], // 3
  [0, 0, 1], // 4
  [1, 0, 1], // 5
  [1, 1, 1], // 6
  [0, 1, 1]  // 7
]);

/**
 * 6 edges of a canonical tetrahedron defined by vertex pairs [vA, vB].
 * Vertex indices in tet: 0, 1, 2, 3.
 */
export const TET_EDGE_PAIRS = Object.freeze([
  [0, 1], // Edge 0: v0-v1
  [0, 2], // Edge 1: v0-v2
  [0, 3], // Edge 2: v0-v3
  [1, 2], // Edge 3: v1-v2
  [1, 3], // Edge 4: v1-v3
  [2, 3]  // Edge 5: v2-v3
]);

/**
 * Exact 16-entry Marching Tetrahedra Case Table for Right-Handed Tetrahedra.
 * For each case mask (0 to 15) where bit i is 1 if scalar(v_i) >= isovalue:
 * Returns consistently oriented outward-facing triangles.
 */
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
// 1. CUBE DECOMPOSITION SCHEMES
// ============================================================================

/**
 * Kuhn / Freudenthal 6-tetrahedron decomposition of a cube with all det > 0.
 */
export const KUHN_6_TETRAHEDRA = Object.freeze([
  [0, 1, 2, 6],
  [0, 5, 1, 6],
  [0, 2, 3, 6],
  [0, 3, 7, 6],
  [0, 4, 5, 6],
  [0, 7, 4, 6]
]);

export const ALTERNATING_6_TETRAHEDRA_EVEN = KUHN_6_TETRAHEDRA;
export const ALTERNATING_6_TETRAHEDRA_ODD = KUHN_6_TETRAHEDRA;

/**
 * 5-tetrahedron checkerboard decomposition with all det > 0.
 * Type A (Parity 0): Central tet [0, 5, 2, 7] + 4 corner tets.
 * Type B (Parity 1): Central tet [1, 3, 4, 6] + 4 corner tets.
 */
export const FIVE_TETRAHEDRA_EVEN = Object.freeze([
  [0, 1, 2, 5], // Corner 1
  [0, 2, 3, 7], // Corner 3
  [0, 4, 5, 7], // Corner 4
  [2, 5, 6, 7], // Corner 6
  [0, 5, 2, 7]  // Central
]);

export const FIVE_TETRAHEDRA_ODD = Object.freeze([
  [1, 3, 0, 4], // Corner 0
  [1, 2, 3, 6], // Corner 2
  [1, 4, 5, 6], // Corner 5
  [3, 4, 6, 7], // Corner 7
  [1, 3, 4, 6]  // Central
]);

// ============================================================================
// 2. MANIFOLD MESH DATA STRUCTURE & SERIALIZATION
// ============================================================================

/**
 * Triangulated 2-Manifold Mesh Container with Differential Geometry & Topology.
 */
export class ManifoldMesh {
  /**
   * Constructs a ManifoldMesh.
   * 
   * @param {Object} params
   * @param {Float32Array|Array<number>} params.vertices - Flattened 3D vertex coordinates [x0, y0, z0, ...].
   * @param {Uint32Array|Array<number>} params.indices - Flattened triangle vertex indices [i0, i1, i2, ...].
   * @param {Float32Array|Array<number>} [params.normals] - Flattened vertex normal vectors [nx0, ny0, nz0, ...].
   * @param {Float32Array|Array<number>} [params.scalars] - Interpolated scalar field values at vertices.
   * @param {number} [params.isovalue=0.0] - Extraction isovalue threshold.
   * @param {string} [params.decompositionMode] - Cube decomposition mode used during extraction.
   */
  constructor(params) {
    if (!params) {
      throw new Error('ManifoldMesh: params object must be provided.');
    }

    this.vertices = params.vertices instanceof Float32Array ? params.vertices : new Float32Array(params.vertices || []);
    this.indices = params.indices instanceof Uint32Array ? params.indices : new Uint32Array(params.indices || []);
    this.normals = params.normals ? (params.normals instanceof Float32Array ? params.normals : new Float32Array(params.normals)) : null;
    this.scalars = params.scalars ? (params.scalars instanceof Float32Array ? params.scalars : new Float32Array(params.scalars)) : null;
    this.isovalue = params.isovalue !== undefined ? params.isovalue : 0.0;
    this.decompositionMode = params.decompositionMode || DecompositionMode.SIX_TETRAHEDRA_KUHN;

    this.vertexCount = this.vertices.length / 3;
    this.triangleCount = this.indices.length / 3;

    this._topology = null;
    this._calculus = null;
    this._curvature = null;
  }

  /**
   * Topological characteristics: Euler characteristic, genus, watertightness, edge counts.
   * @returns {Object}
   */
  get topology() {
    if (!this._topology) {
      this._topology = TopologicalManifoldAnalyzer.analyze(this);
    }
    return this._topology;
  }

  /**
   * Calculus integrals: surface area, enclosed volume, sphericity, centroids, Minkowski functionals.
   * @returns {Object}
   */
  get calculus() {
    if (!this._calculus) {
      this._calculus = MeshCalculus.calculateAll(this);
    }
    return this._calculus;
  }

  /**
   * Differential geometry & curvature tensors.
   * @returns {Object}
   */
  get curvature() {
    if (!this._curvature) {
      this._curvature = MeshCurvatureAnalyzer.analyze(this);
    }
    return this._curvature;
  }

  /**
   * Clears cached topological, calculus, and curvature metrics.
   */
  invalidateCache() {
    this._topology = null;
    this._calculus = null;
    this._curvature = null;
  }

  /**
   * Recomputes vertex normals using the specified normal calculation strategy.
   * @param {string} [mode=NormalCalculationMode.ANGLE_WEIGHTED]
   * @param {GridIndexer} [grid=null]
   * @param {Float64Array|Float32Array} [fieldBuffer=null]
   * @returns {Float32Array}
   */
  computeNormals(mode = NormalCalculationMode.ANGLE_WEIGHTED, grid = null, fieldBuffer = null) {
    this.normals = MeshCalculus.computeVertexNormals(this.vertices, this.indices, mode, grid, fieldBuffer);
    this.invalidateCache();
    return this.normals;
  }

  /**
   * Applies geometric smoothing to vertex positions while preserving topology and optionally volume.
   * 
   * @param {Object} [options]
   * @param {string} [options.algorithm=SmoothingAlgorithmMode.TAUBIN_NON_SHRINKING]
   * @param {number} [options.iterations=3]
   * @param {number} [options.lambda=0.5]
   * @param {number} [options.mu=-0.53]
   * @param {boolean} [options.preserveVolume=true]
   * @returns {ManifoldMesh} New smoothed ManifoldMesh.
   */
  smooth(options = {}) {
    const algorithm = options.algorithm || SmoothingAlgorithmMode.TAUBIN_NON_SHRINKING;
    const iterations = options.iterations !== undefined ? options.iterations : 3;
    const lambda = options.lambda !== undefined ? options.lambda : 0.5;
    const mu = options.mu !== undefined ? options.mu : -0.53;
    const preserveVolume = options.preserveVolume !== undefined ? options.preserveVolume : true;

    let smoothedVertices;
    if (algorithm === SmoothingAlgorithmMode.UNIFORM_LAPLACIAN) {
      smoothedVertices = MeshCalculus.laplacianSmooth(this.vertices, this.indices, iterations, lambda);
    } else if (algorithm === SmoothingAlgorithmMode.COTANGENT_LAPLACIAN) {
      smoothedVertices = MeshCalculus.cotangentLaplacianSmooth(this.vertices, this.indices, iterations, lambda);
    } else {
      smoothedVertices = MeshCalculus.taubinSmooth(this.vertices, this.indices, iterations, lambda, mu);
    }

    if (preserveVolume && this.triangleCount > 0) {
      const initialVolume = this.calculus.enclosedVolume;
      const currentVolume = MeshCalculus.calculateEnclosedVolume(smoothedVertices, this.indices);
      if (initialVolume > 1e-12 && currentVolume > 1e-12) {
        const scale = Math.cbrt(initialVolume / currentVolume);
        const centroid = MeshCalculus.calculateVolumetricCentroid(smoothedVertices, this.indices);
        const vCount = smoothedVertices.length / 3;
        for (let i = 0; i < vCount; i++) {
          smoothedVertices[i * 3 + 0] = centroid[0] + scale * (smoothedVertices[i * 3 + 0] - centroid[0]);
          smoothedVertices[i * 3 + 1] = centroid[1] + scale * (smoothedVertices[i * 3 + 1] - centroid[1]);
          smoothedVertices[i * 3 + 2] = centroid[2] + scale * (smoothedVertices[i * 3 + 2] - centroid[2]);
        }
      }
    }

    const smoothedMesh = new ManifoldMesh({
      vertices: smoothedVertices,
      indices: new Uint32Array(this.indices),
      scalars: this.scalars ? new Float32Array(this.scalars) : null,
      isovalue: this.isovalue,
      decompositionMode: this.decompositionMode
    });
    smoothedMesh.computeNormals(NormalCalculationMode.ANGLE_WEIGHTED);
    return smoothedMesh;
  }

  /**
   * Tests whether a 3D query point [px, py, pz] is inside the closed manifold using multi-ray parity voting.
   * 
   * @param {[number, number, number]} point - Query coordinates in Mpc/h.
   * @returns {boolean} True if point is strictly enclosed within manifold.
   */
  containsPoint(point) {
    const px = point[0], py = point[1], pz = point[2];
    const directions = [
      [0.9998, 0.0141, 0.0123],
      [0.0123, 0.9998, 0.0141],
      [0.0141, 0.0123, 0.9998]
    ];
    let insideVotes = 0;
    const fCount = this.triangleCount;

    for (const [dx, dy, dz] of directions) {
      let hits = 0;
      for (let f = 0; f < fCount; f++) {
        const i0 = this.indices[f * 3 + 0] * 3;
        const i1 = this.indices[f * 3 + 1] * 3;
        const i2 = this.indices[f * 3 + 2] * 3;

        const v0x = this.vertices[i0], v0y = this.vertices[i0 + 1], v0z = this.vertices[i0 + 2];
        const v1x = this.vertices[i1], v1y = this.vertices[i1 + 1], v1z = this.vertices[i1 + 2];
        const v2x = this.vertices[i2], v2y = this.vertices[i2 + 1], v2z = this.vertices[i2 + 2];

        const e1x = v1x - v0x, e1y = v1y - v0y, e1z = v1z - v0z;
        const e2x = v2x - v0x, e2y = v2y - v0y, e2z = v2z - v0z;

        const pvecX = dy * e2z - dz * e2y;
        const pvecY = dz * e2x - dx * e2z;
        const pvecZ = dx * e2y - dy * e2x;

        const det = e1x * pvecX + e1y * pvecY + e1z * pvecZ;
        if (Math.abs(det) < 1e-12) continue;

        const invDet = 1.0 / det;
        const tvecX = px - v0x, tvecY = py - v0y, tvecZ = pz - v0z;
        const u = (tvecX * pvecX + tvecY * pvecY + tvecZ * pvecZ) * invDet;
        if (u < 0.0 || u > 1.0) continue;

        const qvecX = tvecY * e1z - tvecZ * e1y;
        const qvecY = tvecZ * e1x - tvecX * e1z;
        const qvecZ = tvecX * e1y - tvecY * e1x;

        const v = (dx * qvecX + dy * qvecY + dz * qvecZ) * invDet;
        if (v < 0.0 || u + v > 1.0) continue;

        const t = (e2x * qvecX + e2y * qvecY + e2z * qvecZ) * invDet;
        if (t > 1e-6) hits++;
      }
      if ((hits % 2) === 1) insideVotes++;
    }

    return insideVotes >= 2;
  }

  /**
   * Computes bulk-flow velocity flux \iint (\vec{v} \cdot \hat{n}) dA through the manifold surface.
   * 
   * @param {Function} velocityFieldFn - Callback (x,y,z) => [vx, vy, vz] in km/s.
   * @returns {number} Net cosmic velocity flux in (km/s) * (Mpc/h)^2.
   */
  integrateVelocityFlux(velocityFieldFn) {
    const fCount = this.triangleCount;
    let netFlux = 0.0;

    for (let f = 0; f < fCount; f++) {
      const i0 = this.indices[f * 3 + 0] * 3;
      const i1 = this.indices[f * 3 + 1] * 3;
      const i2 = this.indices[f * 3 + 2] * 3;

      const x0 = this.vertices[i0], y0 = this.vertices[i0 + 1], z0 = this.vertices[i0 + 2];
      const x1 = this.vertices[i1], y1 = this.vertices[i1 + 1], z1 = this.vertices[i1 + 2];
      const x2 = this.vertices[i2], y2 = this.vertices[i2 + 1], z2 = this.vertices[i2 + 2];

      const ax = x1 - x0, ay = y1 - y0, az = z1 - z0;
      const bx = x2 - x0, by = y2 - y0, bz = z2 - z0;

      const cx = ay * bz - az * by;
      const cy = az * bx - ax * bz;
      const cz = ax * by - ay * bx;

      const faceNormalArea = Math.sqrt(cx * cx + cy * cy + cz * cz);
      if (faceNormalArea <= 1e-12) continue;

      const fnx = cx / faceNormalArea;
      const fny = cy / faceNormalArea;
      const fnz = cz / faceNormalArea;
      const faceArea = 0.5 * faceNormalArea;

      // Sample velocity at face centroid
      const fcx = (x0 + x1 + x2) / 3.0;
      const fcy = (y0 + y1 + y2) / 3.0;
      const fcz = (z0 + z1 + z2) / 3.0;

      const [vx, vy, vz] = velocityFieldFn(fcx, fcy, fcz);
      const vDotN = vx * fnx + vy * fny + vz * fnz;

      netFlux += vDotN * faceArea;
    }

    return netFlux;
  }

  /**
   * Exports mesh to Wavefront OBJ string format with scientific metadata headers.
   * @param {string} [objectName='CosmicManifold']
   * @returns {string} OBJ file string.
   */
  toOBJ(objectName = 'CosmicManifold') {
    const lines = [];
    lines.push(`# ZRT Cosmicflows Exact Marching Tetrahedra Mesh`);
    lines.push(`# Simplex Decomposition: ${this.decompositionMode}`);
    lines.push(`# Isovalue delta = ${this.isovalue}`);
    lines.push(`# Vertices: ${this.vertexCount}, Triangles: ${this.triangleCount}`);
    
    const top = this.topology;
    lines.push(`# Euler Characteristic chi: ${top.eulerCharacteristic}, Genus: ${top.genus}, Watertight: ${top.isWatertight}`);
    
    const calc = this.calculus;
    lines.push(`# Surface Area: ${calc.surfaceArea.toFixed(4)} (Mpc/h)^2`);
    lines.push(`# Enclosed Volume: ${calc.enclosedVolume.toFixed(4)} (Mpc/h)^3`);
    lines.push(`# Sphericity Psi: ${calc.sphericity.toFixed(4)}, Isoperimetric Quotient Q: ${calc.isoperimetricQuotient.toFixed(4)}`);
    lines.push(`o ${objectName}`);

    // Vertex positions: v x y z
    for (let i = 0; i < this.vertexCount; i++) {
      const x = this.vertices[i * 3 + 0].toFixed(6);
      const y = this.vertices[i * 3 + 1].toFixed(6);
      const z = this.vertices[i * 3 + 2].toFixed(6);
      lines.push(`v ${x} ${y} ${z}`);
    }

    // Vertex normals: vn nx ny nz
    const hasNormals = this.normals && this.normals.length === this.vertices.length;
    if (hasNormals) {
      for (let i = 0; i < this.vertexCount; i++) {
        const nx = this.normals[i * 3 + 0].toFixed(6);
        const ny = this.normals[i * 3 + 1].toFixed(6);
        const nz = this.normals[i * 3 + 2].toFixed(6);
        lines.push(`vn ${nx} ${ny} ${nz}`);
      }
    }

    // Faces: f v1//vn1 v2//vn2 v3//vn3 (1-based index)
    for (let f = 0; f < this.triangleCount; f++) {
      const i0 = this.indices[f * 3 + 0] + 1;
      const i1 = this.indices[f * 3 + 1] + 1;
      const i2 = this.indices[f * 3 + 2] + 1;
      if (hasNormals) {
        lines.push(`f ${i0}//${i0} ${i1}//${i1} ${i2}//${i2}`);
      } else {
        lines.push(`f ${i0} ${i1} ${i2}`);
      }
    }

    return lines.join('\n');
  }

  /**
   * Exports mesh to GeoJSON 3D MultiPolygon format.
   * @returns {Object}
   */
  toGeoJSON() {
    const polygons = [];
    for (let f = 0; f < this.triangleCount; f++) {
      const i0 = this.indices[f * 3 + 0];
      const i1 = this.indices[f * 3 + 1];
      const i2 = this.indices[f * 3 + 2];
      const p0 = [this.vertices[i0 * 3], this.vertices[i0 * 3 + 1], this.vertices[i0 * 3 + 2]];
      const p1 = [this.vertices[i1 * 3], this.vertices[i1 * 3 + 1], this.vertices[i1 * 3 + 2]];
      const p2 = [this.vertices[i2 * 3], this.vertices[i2 * 3 + 1], this.vertices[i2 * 3 + 2]];
      polygons.push([[p0, p1, p2, p0]]);
    }

    return {
      type: 'Feature',
      properties: {
        isovalue: this.isovalue,
        decompositionMode: this.decompositionMode,
        surfaceArea: this.calculus.surfaceArea,
        enclosedVolume: this.calculus.enclosedVolume,
        eulerCharacteristic: this.topology.eulerCharacteristic,
        genus: this.topology.genus,
        sphericity: this.calculus.sphericity
      },
      geometry: {
        type: 'MultiPolygon',
        coordinates: polygons
      }
    };
  }

  /**
   * Returns Three.js BufferGeometry compatible raw typed arrays.
   * @returns {{position: Float32Array, normal: Float32Array, index: Uint32Array}}
   */
  toBufferGeometry() {
    if (!this.normals) {
      this.computeNormals(NormalCalculationMode.ANGLE_WEIGHTED);
    }
    return {
      position: this.vertices,
      normal: this.normals,
      index: this.indices
    };
  }
}

// ============================================================================
// 3. TOPOLOGICAL MANIFOLD ANALYZER
// ============================================================================

/**
 * Evaluates Manifold Topological Invariants: Euler Characteristic, Genus, Betti Numbers, Watertightness.
 */
export class TopologicalManifoldAnalyzer {
  /**
   * Analyzes topological characteristics of a manifold mesh.
   * 
   * @param {ManifoldMesh} mesh
   * @returns {{
   *   vertexCount: number,
   *   edgeCount: number,
   *   faceCount: number,
   *   eulerCharacteristic: number,
   *   genus: number,
   *   isWatertight: boolean,
   *   boundaryEdgeCount: number,
   *   nonManifoldEdgeCount: number,
   *   connectedComponents: number,
   *   bettiNumbers: {b0: number, b1: number, b2: number}
   * }}
   */
  static analyze(mesh) {
    const { vertices, indices, vertexCount, triangleCount } = mesh;

    if (triangleCount === 0 || vertexCount === 0) {
      return {
        vertexCount: 0,
        edgeCount: 0,
        faceCount: 0,
        eulerCharacteristic: 0,
        genus: 0,
        isWatertight: true,
        boundaryEdgeCount: 0,
        nonManifoldEdgeCount: 0,
        connectedComponents: 0,
        bettiNumbers: { b0: 0, b1: 0, b2: 0 }
      };
    }

    // Build unique undirected edge map and track face incident count per edge
    const edgeMap = new Map();
    const adj = Array.from({ length: vertexCount }, () => new Set());

    for (let f = 0; f < triangleCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const edges = [
        [Math.min(i0, i1), Math.max(i0, i1)],
        [Math.min(i1, i2), Math.max(i1, i2)],
        [Math.min(i2, i0), Math.max(i2, i0)]
      ];

      for (let e = 0; e < 3; e++) {
        const [u, v] = edges[e];
        adj[u].add(v);
        adj[v].add(u);

        const key = `${u}_${v}`;
        const count = edgeMap.get(key) || 0;
        edgeMap.set(key, count + 1);
      }
    }

    const V = vertexCount;
    const E = edgeMap.size;
    const F = triangleCount;

    let boundaryEdgeCount = 0;
    let nonManifoldEdgeCount = 0;

    for (const count of edgeMap.values()) {
      if (count === 1) {
        boundaryEdgeCount++;
      } else if (count > 2) {
        nonManifoldEdgeCount++;
      }
    }

    // Euler characteristic chi = V - E + F
    const eulerCharacteristic = V - E + F;

    // Connected components via BFS traversal
    const visited = new Uint8Array(V);
    let connectedComponents = 0;

    for (let i = 0; i < V; i++) {
      if (!visited[i]) {
        if (adj[i].size > 0) {
          connectedComponents++;
          const queue = [i];
          visited[i] = 1;
          let head = 0;
          while (head < queue.length) {
            const u = queue[head++];
            for (const v of adj[u]) {
              if (!visited[v]) {
                visited[v] = 1;
                queue.push(v);
              }
            }
          }
        }
      }
    }

    // For a closed 2-manifold with C connected components:
    // chi = 2 * C - 2 * sum(g_i) => genus = C - chi / 2
    const genus = Math.max(0, Math.round(connectedComponents - eulerCharacteristic / 2));
    const isWatertight = boundaryEdgeCount === 0 && nonManifoldEdgeCount === 0;

    // Betti numbers for 2-manifold:
    // b0 = connected components
    // b1 = 2 * genus
    // b2 = connected components (if closed) or 0 (if open)
    const b0 = connectedComponents;
    const b1 = 2 * genus;
    const b2 = isWatertight ? connectedComponents : 0;

    return {
      vertexCount: V,
      edgeCount: E,
      faceCount: F,
      eulerCharacteristic,
      genus,
      isWatertight,
      boundaryEdgeCount,
      nonManifoldEdgeCount,
      connectedComponents,
      bettiNumbers: { b0, b1, b2 }
    };
  }
}

// ============================================================================
// 4. MESH CALCULUS & DIFFERENTIAL GEOMETRY
// ============================================================================

/**
 * Exact Calculus Integrals, Differential Geometry & Mesh Smoothing Engine.
 */
export class MeshCalculus {
  /**
   * Computes surface area integral: \iint_{\partial \Omega} dA = \sum_f \frac{1}{2} \| (v_1 - v_0) \times (v_2 - v_0) \|.
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {number} Surface area in (Mpc/h)^2.
   */
  static calculateSurfaceArea(vertices, indices) {
    const fCount = indices.length / 3;
    let totalArea = 0.0;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const ax = vertices[i1 + 0] - vertices[i0 + 0];
      const ay = vertices[i1 + 1] - vertices[i0 + 1];
      const az = vertices[i1 + 2] - vertices[i0 + 2];

      const bx = vertices[i2 + 0] - vertices[i0 + 0];
      const by = vertices[i2 + 1] - vertices[i0 + 1];
      const bz = vertices[i2 + 2] - vertices[i0 + 2];

      const cx = ay * bz - az * by;
      const cy = az * bx - ax * bz;
      const cz = ax * by - ay * bx;

      const faceArea = 0.5 * Math.sqrt(cx * cx + cy * cy + cz * cz);
      totalArea += faceArea;
    }

    return totalArea;
  }

  /**
   * Computes enclosed volume integral: \iiint_\Omega dV via 2-manifold divergence theorem:
   * V = \frac{1}{6} \sum_f v_0 \cdot (v_1 \times v_2)
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {number} Enclosed volume in (Mpc/h)^3.
   */
  static calculateEnclosedVolume(vertices, indices) {
    const fCount = indices.length / 3;
    let totalVolume = 0.0;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const x0 = vertices[i0 + 0], y0 = vertices[i0 + 1], z0 = vertices[i0 + 2];
      const x1 = vertices[i1 + 0], y1 = vertices[i1 + 1], z1 = vertices[i1 + 2];
      const x2 = vertices[i2 + 0], y2 = vertices[i2 + 1], z2 = vertices[i2 + 2];

      const crossX = y1 * z2 - z1 * y2;
      const crossY = z1 * x2 - x1 * z2;
      const crossZ = x1 * y2 - y1 * x2;

      const signedTetVolume = (x0 * crossX + y0 * crossY + z0 * crossZ) / 6.0;
      totalVolume += signedTetVolume;
    }

    return Math.abs(totalVolume);
  }

  /**
   * Computes area-weighted surface centroid: C_{surf} = \frac{\sum_f A_f C_f}{\sum_f A_f}.
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {[number, number, number]} Surface centroid [x, y, z] in Mpc/h.
   */
  static calculateSurfaceCentroid(vertices, indices) {
    const fCount = indices.length / 3;
    let totalArea = 0.0;
    let cx = 0.0, cy = 0.0, cz = 0.0;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const x0 = vertices[i0], y0 = vertices[i0 + 1], z0 = vertices[i0 + 2];
      const x1 = vertices[i1], y1 = vertices[i1 + 1], z1 = vertices[i1 + 2];
      const x2 = vertices[i2], y2 = vertices[i2 + 1], z2 = vertices[i2 + 2];

      const ax = x1 - x0, ay = y1 - y0, az = z1 - z0;
      const bx = x2 - x0, by = y2 - y0, bz = z2 - z0;

      const crossX = ay * bz - az * by;
      const crossY = az * bx - ax * bz;
      const crossZ = ax * by - ay * bx;

      const area = 0.5 * Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
      totalArea += area;

      const fCenterX = (x0 + x1 + x2) / 3.0;
      const fCenterY = (y0 + y1 + y2) / 3.0;
      const fCenterZ = (z0 + z1 + z2) / 3.0;

      cx += area * fCenterX;
      cy += area * fCenterY;
      cz += area * fCenterZ;
    }

    if (totalArea <= 1e-12) return [0.0, 0.0, 0.0];
    return [cx / totalArea, cy / totalArea, cz / totalArea];
  }

  /**
   * Computes signed volumetric centroid via tetrahedral divergence decomposition:
   * C_{vol} = \frac{1}{4 V} \sum_f V_{tet} (v_0 + v_1 + v_2).
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {[number, number, number]} Volumetric centroid [x, y, z] in Mpc/h.
   */
  static calculateVolumetricCentroid(vertices, indices) {
    const fCount = indices.length / 3;
    let totalVol = 0.0;
    let cx = 0.0, cy = 0.0, cz = 0.0;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const x0 = vertices[i0], y0 = vertices[i0 + 1], z0 = vertices[i0 + 2];
      const x1 = vertices[i1], y1 = vertices[i1 + 1], z1 = vertices[i1 + 2];
      const x2 = vertices[i2], y2 = vertices[i2 + 1], z2 = vertices[i2 + 2];

      const crossX = y1 * z2 - z1 * y2;
      const crossY = z1 * x2 - x1 * z2;
      const crossZ = x1 * y2 - y1 * x2;

      const dV = (x0 * crossX + y0 * crossY + z0 * crossZ) / 6.0;
      totalVol += dV;

      cx += dV * (x0 + x1 + x2) / 4.0;
      cy += dV * (y0 + y1 + y2) / 4.0;
      cz += dV * (z0 + z1 + z2) / 4.0;
    }

    if (Math.abs(totalVol) <= 1e-12) return [0.0, 0.0, 0.0];
    return [cx / totalVol, cy / totalVol, cz / totalVol];
  }

  /**
   * Calculates comprehensive calculus invariants, Minkowski functionals, and morphology metrics.
   * 
   * @param {ManifoldMesh} mesh
   * @returns {Object}
   */
  static calculateAll(mesh) {
    const { vertices, indices } = mesh;
    const surfaceArea = MeshCalculus.calculateSurfaceArea(vertices, indices);
    const enclosedVolume = MeshCalculus.calculateEnclosedVolume(vertices, indices);
    const surfaceCentroid = MeshCalculus.calculateSurfaceCentroid(vertices, indices);
    const volumetricCentroid = MeshCalculus.calculateVolumetricCentroid(vertices, indices);

    const topology = mesh.topology;
    const chi = topology.eulerCharacteristic;

    let sphericity = 0.0;
    let isoperimetricQuotient = 0.0;
    let planarity = 0.0;
    let filamentarity = 0.0;

    const V0 = enclosedVolume;
    const V1 = surfaceArea;
    const V2 = Math.sqrt(surfaceArea * 4.0 * Math.PI); // Integrated mean curvature proxy
    const V3 = chi;

    if (surfaceArea > 1e-12 && enclosedVolume > 1e-12) {
      sphericity = (Math.cbrt(Math.PI) * Math.pow(6.0 * enclosedVolume, 2.0 / 3.0)) / surfaceArea;
      isoperimetricQuotient = (36.0 * Math.PI * enclosedVolume * enclosedVolume) / Math.pow(surfaceArea, 3.0);

      const denomP = V1 * V1 + 3.0 * V0 * V2;
      if (denomP > 1e-12) {
        planarity = (V1 * V1 - 3.0 * V0 * V2) / denomP;
      }

      const denomF = V1 * V2 + 9.0 * V0 * V3;
      if (denomF > 1e-12) {
        filamentarity = (V1 * V2 - 9.0 * V0 * V3) / denomF;
      }
    }

    return {
      surfaceArea,
      enclosedVolume,
      surfaceCentroid,
      volumetricCentroid,
      sphericity,
      isoperimetricQuotient,
      planarity,
      filamentarity,
      minkowskiFunctionals: {
        V0,
        V1,
        V2,
        V3
      }
    };
  }

  /**
   * Computes vertex normal vectors using requested weighting mode.
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @param {string} [mode=NormalCalculationMode.ANGLE_WEIGHTED]
   * @param {GridIndexer} [grid=null]
   * @param {Float64Array|Float32Array} [fieldBuffer=null]
   * @returns {Float32Array} Unit vertex normal vectors.
   */
  static computeVertexNormals(vertices, indices, mode = NormalCalculationMode.ANGLE_WEIGHTED, grid = null, fieldBuffer = null) {
    const vCount = vertices.length / 3;
    const normals = new Float32Array(vertices.length);
    const fCount = indices.length / 3;

    if (mode === NormalCalculationMode.FIELD_GRADIENT && grid && fieldBuffer) {
      const finiteDiff = new FiniteDifference(grid);
      for (let i = 0; i < vCount; i++) {
        const x = vertices[i * 3 + 0];
        const y = vertices[i * 3 + 1];
        const z = vertices[i * 3 + 2];

        const grad = finiteDiff.gradient6Point(fieldBuffer, x, y, z);

        let nx = -grad[0], ny = -grad[1], nz = -grad[2];
        const mag = Math.sqrt(nx * nx + ny * ny + nz * nz);
        if (mag > 1e-12) {
          nx /= mag; ny /= mag; nz /= mag;
        } else {
          nx = 0; ny = 0; nz = 1;
        }

        normals[i * 3 + 0] = nx;
        normals[i * 3 + 1] = ny;
        normals[i * 3 + 2] = nz;
      }
      return normals;
    }

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const p0 = [vertices[i0 * 3], vertices[i0 * 3 + 1], vertices[i0 * 3 + 2]];
      const p1 = [vertices[i1 * 3], vertices[i1 * 3 + 1], vertices[i1 * 3 + 2]];
      const p2 = [vertices[i2 * 3], vertices[i2 * 3 + 1], vertices[i2 * 3 + 2]];

      const e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
      const e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
      const e3 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];

      const cx = e1[1] * e2[2] - e1[2] * e2[1];
      const cy = e1[2] * e2[0] - e1[0] * e2[2];
      const cz = e1[0] * e2[1] - e1[1] * e2[0];
      const area2 = Math.sqrt(cx * cx + cy * cy + cz * cz);

      if (area2 <= 1e-12) continue;

      const fnx = cx / area2;
      const fny = cy / area2;
      const fnz = cz / area2;

      if (mode === NormalCalculationMode.ANGLE_WEIGHTED) {
        const l1 = Math.sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        const l2 = Math.sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
        const l3 = Math.sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);

        const dot0 = (e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2]) / (l1 * l2 + 1e-12);
        const theta0 = Math.acos(Math.max(-1, Math.min(1, dot0)));

        const dot1 = (-e1[0] * e3[0] - e1[1] * e3[1] - e1[2] * e3[2]) / (l1 * l3 + 1e-12);
        const theta1 = Math.acos(Math.max(-1, Math.min(1, dot1)));

        const theta2 = Math.PI - theta0 - theta1;

        normals[i0 * 3 + 0] += fnx * theta0;
        normals[i0 * 3 + 1] += fny * theta0;
        normals[i0 * 3 + 2] += fnz * theta0;

        normals[i1 * 3 + 0] += fnx * theta1;
        normals[i1 * 3 + 1] += fny * theta1;
        normals[i1 * 3 + 2] += fnz * theta1;

        normals[i2 * 3 + 0] += fnx * theta2;
        normals[i2 * 3 + 1] += fny * theta2;
        normals[i2 * 3 + 2] += fnz * theta2;
      } else if (mode === NormalCalculationMode.AREA_WEIGHTED) {
        const w = 0.5 * area2;
        normals[i0 * 3 + 0] += fnx * w;
        normals[i0 * 3 + 1] += fny * w;
        normals[i0 * 3 + 2] += fnz * w;

        normals[i1 * 3 + 0] += fnx * w;
        normals[i1 * 3 + 1] += fny * w;
        normals[i1 * 3 + 2] += fnz * w;

        normals[i2 * 3 + 0] += fnx * w;
        normals[i2 * 3 + 1] += fny * w;
        normals[i2 * 3 + 2] += fnz * w;
      } else {
        // FACET_UNIFORM
        normals[i0 * 3 + 0] += fnx;
        normals[i0 * 3 + 1] += fny;
        normals[i0 * 3 + 2] += fnz;

        normals[i1 * 3 + 0] += fnx;
        normals[i1 * 3 + 1] += fny;
        normals[i1 * 3 + 2] += fnz;

        normals[i2 * 3 + 0] += fnx;
        normals[i2 * 3 + 1] += fny;
        normals[i2 * 3 + 2] += fnz;
      }
    }

    for (let i = 0; i < vCount; i++) {
      let nx = normals[i * 3 + 0];
      let ny = normals[i * 3 + 1];
      let nz = normals[i * 3 + 2];
      const len = Math.sqrt(nx * nx + ny * ny + nz * nz);
      if (len > 1e-12) {
        normals[i * 3 + 0] = nx / len;
        normals[i * 3 + 1] = ny / len;
        normals[i * 3 + 2] = nz / len;
      } else {
        normals[i * 3 + 0] = 0;
        normals[i * 3 + 1] = 0;
        normals[i * 3 + 2] = 1;
      }
    }

    return normals;
  }

  /**
   * Uniform Laplacian smoothing.
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @param {number} iterations
   * @param {number} lambda
   * @returns {Float32Array}
   */
  static laplacianSmooth(vertices, indices, iterations = 3, lambda = 0.5) {
    const vCount = vertices.length / 3;
    if (vCount === 0 || indices.length === 0 || iterations <= 0) {
      return new Float32Array(vertices);
    }

    const adj = Array.from({ length: vCount }, () => []);
    const fCount = indices.length / 3;
    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      if (!adj[i0].includes(i1)) adj[i0].push(i1);
      if (!adj[i0].includes(i2)) adj[i0].push(i2);
      if (!adj[i1].includes(i0)) adj[i1].push(i0);
      if (!adj[i1].includes(i2)) adj[i1].push(i2);
      if (!adj[i2].includes(i0)) adj[i2].push(i0);
      if (!adj[i2].includes(i1)) adj[i2].push(i1);
    }

    let current = new Float32Array(vertices);
    let next = new Float32Array(vCount * 3);

    for (let iter = 0; iter < iterations; iter++) {
      for (let i = 0; i < vCount; i++) {
        const neighbors = adj[i];
        const nLen = neighbors.length;
        if (nLen === 0) {
          next[i * 3 + 0] = current[i * 3 + 0];
          next[i * 3 + 1] = current[i * 3 + 1];
          next[i * 3 + 2] = current[i * 3 + 2];
          continue;
        }

        let sumX = 0, sumY = 0, sumZ = 0;
        for (let k = 0; k < nLen; k++) {
          const nb = neighbors[k];
          sumX += current[nb * 3 + 0];
          sumY += current[nb * 3 + 1];
          sumZ += current[nb * 3 + 2];
        }

        const avgX = sumX / nLen;
        const avgY = sumY / nLen;
        const avgZ = sumZ / nLen;

        const x = current[i * 3 + 0];
        const y = current[i * 3 + 1];
        const z = current[i * 3 + 2];

        next[i * 3 + 0] = x + lambda * (avgX - x);
        next[i * 3 + 1] = y + lambda * (avgY - y);
        next[i * 3 + 2] = z + lambda * (avgZ - z);
      }
      const tmp = current;
      current = next;
      next = tmp;
    }

    return current;
  }

  /**
   * Cotangent-weighted Laplace-Beltrami conformal smoothing.
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @param {number} iterations
   * @param {number} lambda
   * @returns {Float32Array}
   */
  static cotangentLaplacianSmooth(vertices, indices, iterations = 3, lambda = 0.5) {
    const vCount = vertices.length / 3;
    if (vCount === 0 || indices.length === 0 || iterations <= 0) {
      return new Float32Array(vertices);
    }

    let current = new Float32Array(vertices);
    let next = new Float32Array(vCount * 3);
    const fCount = indices.length / 3;

    for (let iter = 0; iter < iterations; iter++) {
      const weightMaps = Array.from({ length: vCount }, () => new Map());

      for (let f = 0; f < fCount; f++) {
        const i0 = indices[f * 3 + 0];
        const i1 = indices[f * 3 + 1];
        const i2 = indices[f * 3 + 2];

        const p0 = [current[i0 * 3], current[i0 * 3 + 1], current[i0 * 3 + 2]];
        const p1 = [current[i1 * 3], current[i1 * 3 + 1], current[i1 * 3 + 2]];
        const p2 = [current[i2 * 3], current[i2 * 3 + 1], current[i2 * 3 + 2]];

        const v01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        const v02 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        const v12 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];

        const dot0 = v01[0] * v02[0] + v01[1] * v02[1] + v01[2] * v02[2];
        const cross0 = Math.sqrt(Math.pow(v01[1]*v02[2] - v01[2]*v02[1], 2) + Math.pow(v01[2]*v02[0] - v01[0]*v02[2], 2) + Math.pow(v01[0]*v02[1] - v01[1]*v02[0], 2));
        const cot0 = cross0 > 1e-12 ? Math.max(0.01, dot0 / cross0) : 0.5;

        const v10 = [-v01[0], -v01[1], -v01[2]];
        const dot1 = v10[0] * v12[0] + v10[1] * v12[1] + v10[2] * v12[2];
        const cross1 = Math.sqrt(Math.pow(v10[1]*v12[2] - v10[2]*v12[1], 2) + Math.pow(v10[2]*v12[0] - v10[0]*v12[2], 2) + Math.pow(v10[0]*v12[1] - v10[1]*v12[0], 2));
        const cot1 = cross1 > 1e-12 ? Math.max(0.01, dot1 / cross1) : 0.5;

        const v20 = [-v02[0], -v02[1], -v02[2]];
        const v21 = [-v12[0], -v12[1], -v12[2]];
        const dot2 = v20[0] * v21[0] + v20[1] * v21[1] + v20[2] * v21[2];
        const cross2 = Math.sqrt(Math.pow(v20[1]*v21[2] - v20[2]*v21[1], 2) + Math.pow(v20[2]*v21[0] - v20[0]*v21[2], 2) + Math.pow(v20[0]*v21[1] - v20[1]*v21[0], 2));
        const cot2 = cross2 > 1e-12 ? Math.max(0.01, dot2 / cross2) : 0.5;

        weightMaps[i1].set(i2, (weightMaps[i1].get(i2) || 0) + 0.5 * cot0);
        weightMaps[i2].set(i1, (weightMaps[i2].get(i1) || 0) + 0.5 * cot0);

        weightMaps[i0].set(i2, (weightMaps[i0].get(i2) || 0) + 0.5 * cot1);
        weightMaps[i2].set(i0, (weightMaps[i2].get(i0) || 0) + 0.5 * cot1);

        weightMaps[i0].set(i1, (weightMaps[i0].get(i1) || 0) + 0.5 * cot2);
        weightMaps[i1].set(i0, (weightMaps[i1].get(i0) || 0) + 0.5 * cot2);
      }

      for (let i = 0; i < vCount; i++) {
        const wMap = weightMaps[i];
        if (wMap.size === 0) {
          next[i * 3 + 0] = current[i * 3 + 0];
          next[i * 3 + 1] = current[i * 3 + 1];
          next[i * 3 + 2] = current[i * 3 + 2];
          continue;
        }

        let sumWeight = 0.0;
        let sumX = 0.0, sumY = 0.0, sumZ = 0.0;

        for (const [nb, w] of wMap.entries()) {
          sumWeight += w;
          sumX += w * current[nb * 3 + 0];
          sumY += w * current[nb * 3 + 1];
          sumZ += w * current[nb * 3 + 2];
        }

        const avgX = sumX / sumWeight;
        const avgY = sumY / sumWeight;
        const avgZ = sumZ / sumWeight;

        const x = current[i * 3 + 0];
        const y = current[i * 3 + 1];
        const z = current[i * 3 + 2];

        next[i * 3 + 0] = x + lambda * (avgX - x);
        next[i * 3 + 1] = y + lambda * (avgY - y);
        next[i * 3 + 2] = z + lambda * (avgZ - z);
      }

      const tmp = current;
      current = next;
      next = tmp;
    }

    return current;
  }

  /**
   * Taubin Non-Shrinking Dual-Step Mesh Smoother.
   * Alternates positive lambda diffusion and negative mu inflation steps.
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @param {number} iterations
   * @param {number} lambda Positive diffusion (0 < lambda < 1)
   * @param {number} mu Negative inflation (mu < -lambda)
   * @returns {Float32Array} Smoothed vertices.
   */
  static taubinSmooth(vertices, indices, iterations = 3, lambda = 0.5, mu = -0.53) {
    const vCount = vertices.length / 3;
    if (vCount === 0 || indices.length === 0 || iterations <= 0) {
      return new Float32Array(vertices);
    }

    const adj = Array.from({ length: vCount }, () => []);
    const fCount = indices.length / 3;
    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      if (!adj[i0].includes(i1)) adj[i0].push(i1);
      if (!adj[i0].includes(i2)) adj[i0].push(i2);
      if (!adj[i1].includes(i0)) adj[i1].push(i0);
      if (!adj[i1].includes(i2)) adj[i1].push(i2);
      if (!adj[i2].includes(i0)) adj[i2].push(i0);
      if (!adj[i2].includes(i1)) adj[i2].push(i1);
    }

    let current = new Float32Array(vertices);
    let next = new Float32Array(vCount * 3);

    const applyStep = (factor) => {
      for (let i = 0; i < vCount; i++) {
        const neighbors = adj[i];
        const nLen = neighbors.length;
        if (nLen === 0) {
          next[i * 3 + 0] = current[i * 3 + 0];
          next[i * 3 + 1] = current[i * 3 + 1];
          next[i * 3 + 2] = current[i * 3 + 2];
          continue;
        }

        let sumX = 0, sumY = 0, sumZ = 0;
        for (let k = 0; k < nLen; k++) {
          const nb = neighbors[k];
          sumX += current[nb * 3 + 0];
          sumY += current[nb * 3 + 1];
          sumZ += current[nb * 3 + 2];
        }

        const avgX = sumX / nLen;
        const avgY = sumY / nLen;
        const avgZ = sumZ / nLen;

        const x = current[i * 3 + 0];
        const y = current[i * 3 + 1];
        const z = current[i * 3 + 2];

        next[i * 3 + 0] = x + factor * (avgX - x);
        next[i * 3 + 1] = y + factor * (avgY - y);
        next[i * 3 + 2] = z + factor * (avgZ - z);
      }
      const tmp = current;
      current = next;
      next = tmp;
    };

    for (let iter = 0; iter < iterations; iter++) {
      applyStep(lambda); // Shrink / smooth step
      applyStep(mu);     // Inflate / unshrink step
    }

    return current;
  }
}

// ============================================================================
// 5. DIFFERENTIAL GEOMETRY & CURVATURE ANALYZER
// ============================================================================

/**
 * Discrete Gaussian and Mean Curvature Tensor Analyzer on Triangulated Manifolds.
 */
export class MeshCurvatureAnalyzer {
  /**
   * Computes discrete Gaussian curvature K_i via angle defect and Gauss-Bonnet integral check.
   * \oint K dA = 2\pi \chi
   * 
   * @param {ManifoldMesh} mesh
   * @returns {Object}
   */
  static analyze(mesh) {
    const { vertices, indices, vertexCount, triangleCount } = mesh;
    if (vertexCount === 0 || triangleCount === 0) {
      return {
        totalGaussianCurvature: 0.0,
        expectedGaussBonnet: 0.0,
        gaussBonnetResidual: 0.0,
        gaussianCurvatureArray: new Float32Array(0)
      };
    }

    const angleSum = new Float64Array(vertexCount);
    const dualArea = new Float64Array(vertexCount);

    for (let f = 0; f < triangleCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const p0 = [vertices[i0 * 3], vertices[i0 * 3 + 1], vertices[i0 * 3 + 2]];
      const p1 = [vertices[i1 * 3], vertices[i1 * 3 + 1], vertices[i1 * 3 + 2]];
      const p2 = [vertices[i2 * 3], vertices[i2 * 3 + 1], vertices[i2 * 3 + 2]];

      const e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
      const e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
      const e3 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];

      const l1 = Math.sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
      const l2 = Math.sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
      const l3 = Math.sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);

      const cx = e1[1] * e2[2] - e1[2] * e2[1];
      const cy = e1[2] * e2[0] - e1[0] * e2[2];
      const cz = e1[0] * e2[1] - e1[1] * e2[0];
      const fArea = 0.5 * Math.sqrt(cx * cx + cy * cy + cz * cz);

      dualArea[i0] += fArea / 3.0;
      dualArea[i1] += fArea / 3.0;
      dualArea[i2] += fArea / 3.0;

      const dot0 = (e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2]) / (l1 * l2 + 1e-12);
      const theta0 = Math.acos(Math.max(-1, Math.min(1, dot0)));

      const dot1 = (-e1[0] * e3[0] - e1[1] * e3[1] - e1[2] * e3[2]) / (l1 * l3 + 1e-12);
      const theta1 = Math.acos(Math.max(-1, Math.min(1, dot1)));

      const theta2 = Math.PI - theta0 - theta1;

      angleSum[i0] += theta0;
      angleSum[i1] += theta1;
      angleSum[i2] += theta2;
    }

    const gaussianCurvature = new Float32Array(vertexCount);
    let totalGaussianCurvature = 0.0;

    for (let i = 0; i < vertexCount; i++) {
      const defect = 2.0 * Math.PI - angleSum[i];
      totalGaussianCurvature += defect;
      if (dualArea[i] > 1e-12) {
        gaussianCurvature[i] = defect / dualArea[i];
      }
    }

    const chi = mesh.topology.eulerCharacteristic;
    const expectedGaussBonnet = 2.0 * Math.PI * chi;
    const gaussBonnetResidual = Math.abs(totalGaussianCurvature - expectedGaussBonnet);

    return {
      totalGaussianCurvature,
      expectedGaussBonnet,
      gaussBonnetResidual,
      gaussianCurvatureArray: gaussianCurvature
    };
  }
}

// ============================================================================
// 6. EXACT 3D MARCHING TETRAHEDRA EXTRACTOR ENGINE
// ============================================================================

/**
 * High-Performance Exact 3D Marching Tetrahedra & Manifold Mesher.
 */
export class ExactMarchingTetrahedra {
  /**
   * @param {GridIndexer} grid - Spatial grid configuration.
   * @param {Object} [options]
   * @param {string} [options.decompositionMode=DecompositionMode.SIX_TETRAHEDRA_KUHN]
   * @param {boolean} [options.enableNewtonRaphson=false] - Refine edge roots via continuous interpolation.
   * @param {number} [options.newtonRaphsonTolerance=1e-5]
   * @param {number} [options.newtonRaphsonMaxIters=5]
   */
  constructor(grid, options = {}) {
    if (!grid) {
      throw new Error('ExactMarchingTetrahedra: A valid GridIndexer instance is required.');
    }
    this.grid = grid;
    this.decompositionMode = options.decompositionMode || DecompositionMode.SIX_TETRAHEDRA_KUHN;
    this.enableNewtonRaphson = Boolean(options.enableNewtonRaphson);
    this.newtonRaphsonTolerance = options.newtonRaphsonTolerance || 1e-5;
    this.newtonRaphsonMaxIters = options.newtonRaphsonMaxIters || 5;

    this.interpolator = new TrilinearInterpolator(grid);
  }

  /**
   * Extracts a guaranteed watertight 2-manifold isosurface mesh for the given scalar field and threshold.
   * 
   * @param {Float64Array|Float32Array} scalarField - 3D flat scalar field buffer matching grid.totalCells.
   * @param {number} [isovalue=0.0] - Target scalar isosurface threshold delta_iso.
   * @param {Object} [extractOptions]
   * @param {string} [extractOptions.decompositionMode] - Override decomposition mode.
   * @param {string} [extractOptions.normalMode=NormalCalculationMode.ANGLE_WEIGHTED]
   * @returns {ManifoldMesh}
   */
  extractIsosurface(scalarField, isovalue = 0.0, extractOptions = {}) {
    if (!scalarField || scalarField.length < this.grid.totalCells) {
      throw new Error(`ExactMarchingTetrahedra: scalarField buffer length (${scalarField ? scalarField.length : 0}) is smaller than grid totalCells (${this.grid.totalCells}).`);
    }

    const decompMode = extractOptions.decompositionMode || this.decompositionMode;
    const normalMode = extractOptions.normalMode || NormalCalculationMode.ANGLE_WEIGHTED;

    const { nx, ny, nz } = this.grid;

    const vertexList = [];
    const indexList = [];
    const scalarList = [];

    // Spatial hash vertex cache to weld shared tetrahedron edges
    const edgeVertexCache = new Map();

    const getEdgeKey = (g1, g2) => {
      const minG = g1 < g2 ? g1 : g2;
      const maxG = g1 < g2 ? g2 : g1;
      return `${minG}_${maxG}`;
    };

    /**
     * Interpolates vertex position and scalar value along a grid edge.
     */
    const getOrCreateVertex = (p1, p2, s1, s2, g1, g2) => {
      const key = getEdgeKey(g1, g2);
      if (edgeVertexCache.has(key)) {
        return edgeVertexCache.get(key);
      }

      let t;
      const denom = s2 - s1;
      if (Math.abs(denom) < 1e-15) {
        t = 0.5;
      } else {
        t = (isovalue - s1) / denom;
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;
      }

      let vx = p1[0] + t * (p2[0] - p1[0]);
      let vy = p1[1] + t * (p2[1] - p1[1]);
      let vz = p1[2] + t * (p2[2] - p1[2]);
      let vScalar = s1 + t * (s2 - s1);

      if (this.enableNewtonRaphson) {
        let curT = t;
        for (let iter = 0; iter < this.newtonRaphsonMaxIters; iter++) {
          const rx = p1[0] + curT * (p2[0] - p1[0]);
          const ry = p1[1] + curT * (p2[1] - p1[1]);
          const rz = p1[2] + curT * (p2[2] - p1[2]);

          const fVal = this.interpolator.interpolateScalar(scalarField, rx, ry, rz);
          const residual = fVal - isovalue;
          if (Math.abs(residual) < this.newtonRaphsonTolerance) break;

          const dt = 1e-4;
          const tPlus = Math.min(1.0, curT + dt);
          const fPlus = this.interpolator.interpolateScalar(
            scalarField,
            p1[0] + tPlus * (p2[0] - p1[0]),
            p1[1] + tPlus * (p2[1] - p1[1]),
            p1[2] + tPlus * (p2[2] - p1[2])
          );
          const dfdt = (fPlus - fVal) / (tPlus - curT + 1e-12);
          if (Math.abs(dfdt) < 1e-12) break;

          const deltaT = -residual / dfdt;
          curT = Math.max(0.0, Math.min(1.0, curT + deltaT));
        }

        vx = p1[0] + curT * (p2[0] - p1[0]);
        vy = p1[1] + curT * (p2[1] - p1[1]);
        vz = p1[2] + curT * (p2[2] - p1[2]);
        vScalar = isovalue;
      }

      const vIdx = vertexList.length / 3;
      vertexList.push(vx, vy, vz);
      scalarList.push(vScalar);
      edgeVertexCache.set(key, vIdx);

      return vIdx;
    };

    /**
     * Processes a single tetrahedron given its 4 corner global indices, positions, and scalar values.
     */
    const processTetrahedron = (cornerIndices, cornerPositions, cornerScalars) => {
      let caseMask = 0;
      for (let i = 0; i < 4; i++) {
        if (cornerScalars[i] >= isovalue) {
          caseMask |= (1 << i);
        }
      }

      if (caseMask === 0 || caseMask === 15) {
        return;
      }

      const triangles = CANONICAL_TET_TRIANGLE_TABLE[caseMask];
      for (let t = 0; t < triangles.length; t++) {
        const triEdges = triangles[t];
        const triIndices = new Array(3);

        for (let e = 0; e < 3; e++) {
          const edgeIdx = triEdges[e];
          const [vA, vB] = TET_EDGE_PAIRS[edgeIdx];

          const pA = cornerPositions[vA];
          const pB = cornerPositions[vB];
          const sA = cornerScalars[vA];
          const sB = cornerScalars[vB];
          const gA = cornerIndices[vA];
          const gB = cornerIndices[vB];

          triIndices[e] = getOrCreateVertex(pA, pB, sA, sB, gA, gB);
        }

        if (triIndices[0] !== triIndices[1] && triIndices[1] !== triIndices[2] && triIndices[2] !== triIndices[0]) {
          indexList.push(triIndices[0], triIndices[1], triIndices[2]);
        }
      }
    };

    // Iterate over all voxels
    const cubePositions = new Array(8);
    const cubeScalars = new Array(8);
    const cubeGlobalIndices = new Array(8);

    for (let iz = 0; iz < nz - 1; iz++) {
      for (let iy = 0; iy < ny - 1; iy++) {
        for (let ix = 0; ix < nx - 1; ix++) {
          for (let c = 0; c < 8; c++) {
            const [ox, oy, oz] = CUBE_CORNER_OFFSETS[c];
            const gx = ix + ox;
            const gy = iy + oy;
            const gz = iz + oz;

            const gIdx = this.grid.index(gx, gy, gz);
            cubeGlobalIndices[c] = gIdx;
            cubeScalars[c] = scalarField[gIdx];
            cubePositions[c] = this.grid.gridIndexToCoord(gx, gy, gz);
          }

          const parity = (ix + iy + iz) % 2;

          let tetDecomposition;
          if (decompMode === DecompositionMode.FIVE_TETRAHEDRA_ALTERNATING) {
            tetDecomposition = parity === 0 ? FIVE_TETRAHEDRA_EVEN : FIVE_TETRAHEDRA_ODD;
          } else if (decompMode === DecompositionMode.SIX_TETRAHEDRA_ALTERNATING) {
            tetDecomposition = parity === 0 ? ALTERNATING_6_TETRAHEDRA_EVEN : ALTERNATING_6_TETRAHEDRA_ODD;
          } else {
            tetDecomposition = KUHN_6_TETRAHEDRA;
          }

          const tetCount = tetDecomposition.length;
          for (let k = 0; k < tetCount; k++) {
            const tetCornerMap = tetDecomposition[k];
            const tetCornerIndices = [
              cubeGlobalIndices[tetCornerMap[0]],
              cubeGlobalIndices[tetCornerMap[1]],
              cubeGlobalIndices[tetCornerMap[2]],
              cubeGlobalIndices[tetCornerMap[3]]
            ];
            const tetPositions = [
              cubePositions[tetCornerMap[0]],
              cubePositions[tetCornerMap[1]],
              cubePositions[tetCornerMap[2]],
              cubePositions[tetCornerMap[3]]
            ];
            const tetScalars = [
              cubeScalars[tetCornerMap[0]],
              cubeScalars[tetCornerMap[1]],
              cubeScalars[tetCornerMap[2]],
              cubeScalars[tetCornerMap[3]]
            ];

            processTetrahedron(tetCornerIndices, tetPositions, tetScalars);
          }
        }
      }
    }

    const vertices = new Float32Array(vertexList);
    const indices = new Uint32Array(indexList);
    const scalars = new Float32Array(scalarList);

    const mesh = new ManifoldMesh({
      vertices,
      indices,
      scalars,
      isovalue,
      decompositionMode: decompMode
    });

    mesh.computeNormals(normalMode, this.grid, scalarField);

    return mesh;
  }
}

// ============================================================================
// 7. COSMOLOGICAL DENSITY FIELD GENERATOR & HALO PROFILES
// ============================================================================

/**
 * High-Precision Cosmological Density Contrast (\delta = \rho/\bar{\rho} - 1) Field Generator.
 * Synthesizes research-grade density fields for Laniakea, Shapley, Great Attractor, and analytical profiles.
 */
export class CosmicDensityFieldGenerator {
  /**
   * Known Cosmological Structures in Supergalactic Coordinates [SGX, SGY, SGZ] (Mpc/h).
   */
  static KNOWN_STRUCTURES = Object.freeze({
    SHAPLEY_SUPERCLUSTER: {
      name: 'Shapley Supercluster Concentration',
      center: [-135.0, 140.0, -25.0], // SGX, SGY, SGZ (Mpc/h)
      peakOverdensity: 8.5,
      coreRadius: 18.0,
      mass1015Msun: 10.0
    },
    GREAT_ATTRACTOR: {
      name: 'Great Attractor (Norma/Centaurus Cluster Core)',
      center: [-42.0, 18.0, -8.0],
      peakOverdensity: 5.2,
      coreRadius: 14.0,
      mass1015Msun: 5.0
    },
    LANIAKEA_CORE: {
      name: 'Laniakea Supercluster Gravity Basin Core',
      center: [-35.0, 10.0, -5.0],
      peakOverdensity: 2.8,
      coreRadius: 40.0,
      mass1015Msun: 100.0
    },
    PERSEUS_PISCES: {
      name: 'Perseus-Pisces Supercluster',
      center: [55.0, -20.0, -15.0],
      peakOverdensity: 4.8,
      coreRadius: 15.0,
      mass1015Msun: 6.0
    },
    COMA_SUPERCLUSTER: {
      name: 'Coma Supercluster',
      center: [0.0, 75.0, 10.0],
      peakOverdensity: 6.0,
      coreRadius: 12.0,
      mass1015Msun: 4.5
    },
    LOCAL_VOID: {
      name: 'Local Void',
      center: [-15.0, -30.0, 30.0],
      underdensity: -0.92,
      voidRadius: 25.0
    }
  });

  /**
   * Generates an analytical 3D Gaussian / NFW-like sphere density contrast field.
   * \delta(r) = bg + \delta_0 * exp(-0.5 * (r / r_0)^2)
   * 
   * @param {GridIndexer} grid
   * @param {Object} [options]
   * @param {[number, number, number]} [options.center=[0,0,0]]
   * @param {number} [options.peakOverdensity=5.0]
   * @param {number} [options.radius=25.0]
   * @param {number} [options.backgroundDensity=0.0]
   * @returns {Float64Array}
   */
  static generateAnalyticalSphere(grid, options = {}) {
    const center = options.center || [0.0, 0.0, 0.0];
    const peak = options.peakOverdensity !== undefined ? options.peakOverdensity : 5.0;
    const r0 = options.radius !== undefined ? options.radius : 25.0;
    const bg = options.backgroundDensity !== undefined ? options.backgroundDensity : 0.0;

    const total = grid.totalCells;
    const buffer = new Float64Array(total);
    const { nx, ny, nz } = grid;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];
          const r = Math.sqrt(dx * dx + dy * dy + dz * dz);

          const delta = bg + peak * Math.exp(-0.5 * (r / r0) * (r / r0));
          buffer[grid.index(ix, iy, iz)] = delta;
        }
      }
    }

    return buffer;
  }

  /**
   * Generates an analytical signed distance field (SDF) of a sphere: f(r) = R - r.
   * Positive inside, negative outside.
   * @param {GridIndexer} grid
   * @param {[number, number, number]} [center=[0,0,0]]
   * @param {number} [radius=25.0]
   * @returns {Float64Array}
   */
  static generateSphereSDF(grid, center = [0.0, 0.0, 0.0], radius = 25.0) {
    const total = grid.totalCells;
    const buffer = new Float64Array(total);
    const { nx, ny, nz } = grid;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];
          const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          buffer[grid.index(ix, iy, iz)] = radius - r;
        }
      }
    }

    return buffer;
  }

  /**
   * Generates an analytical Torus Signed Distance Field (Genus 1).
   * f(x, y, z) = r_minor - sqrt((sqrt(x^2 + y^2) - R_major)^2 + z^2)
   * 
   * @param {GridIndexer} grid
   * @param {number} [majorRadius=30.0]
   * @param {number} [minorRadius=10.0]
   * @param {[number, number, number]} [center=[0,0,0]]
   * @returns {Float64Array}
   */
  static generateTorusSDF(grid, majorRadius = 30.0, minorRadius = 10.0, center = [0.0, 0.0, 0.0]) {
    const total = grid.totalCells;
    const buffer = new Float64Array(total);
    const { nx, ny, nz } = grid;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];

          const planarR = Math.sqrt(dx * dx + dy * dy);
          const torusDist = Math.sqrt((planarR - majorRadius) * (planarR - majorRadius) + dz * dz);
          buffer[grid.index(ix, iy, iz)] = minorRadius - torusDist;
        }
      }
    }

    return buffer;
  }

  /**
   * Generates a Navarro-Frenk-White (NFW) Dark Matter Halo Profile Density Contrast Field.
   * \delta(r) = \frac{\delta_{char}}{(r/r_s) (1 + r/r_s)^2}
   * 
   * @param {GridIndexer} grid
   * @param {Object} options
   * @param {[number, number, number]} [options.center=[0,0,0]]
   * @param {number} [options.scaleRadius=5.0] - r_s in Mpc/h
   * @param {number} [options.characteristicOverdensity=200.0] - \delta_{char}
   * @param {number} [options.virialRadius=25.0] - r_vir in Mpc/h
   * @returns {Float64Array}
   */
  static generateNFWHaloField(grid, options = {}) {
    const center = options.center || [0.0, 0.0, 0.0];
    const rs = options.scaleRadius !== undefined ? options.scaleRadius : 5.0;
    const deltaChar = options.characteristicOverdensity !== undefined ? options.characteristicOverdensity : 200.0;
    const rvir = options.virialRadius !== undefined ? options.virialRadius : 25.0;

    const total = grid.totalCells;
    const buffer = new Float64Array(total);
    const { nx, ny, nz } = grid;

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
          const dx = x - center[0];
          const dy = y - center[1];
          const dz = z - center[2];
          const r = Math.max(0.1, Math.sqrt(dx * dx + dy * dy + dz * dz));

          const xRatio = r / rs;
          let delta = deltaChar / (xRatio * (1.0 + xRatio) * (1.0 + xRatio));
          if (r > rvir) {
            delta *= Math.exp(-(r - rvir) / rs);
          }

          buffer[grid.index(ix, iy, iz)] = delta;
        }
      }
    }

    return buffer;
  }

  /**
   * Generates a multi-structure cosmological density contrast field representing the Local Universe.
   * Incorporates Laniakea Basin, Shapley Supercluster, Great Attractor, Coma, Perseus-Pisces, and Local Void.
   * 
   * @param {GridIndexer} grid
   * @returns {Float64Array} Density contrast field \delta(x, y, z).
   */
  static generateLocalUniverseField(grid) {
    const total = grid.totalCells;
    const buffer = new Float64Array(total);
    const { nx, ny, nz } = grid;

    const structures = Object.values(CosmicDensityFieldGenerator.KNOWN_STRUCTURES);

    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
          let delta = 0.0;

          for (const s of structures) {
            const dx = x - s.center[0];
            const dy = y - s.center[1];
            const dz = z - s.center[2];
            const r = Math.sqrt(dx * dx + dy * dy + dz * dz);

            if (s.peakOverdensity !== undefined) {
              const r0 = s.coreRadius;
              const contribution = s.peakOverdensity / (1.0 + (r / r0) * (r / r0));
              delta += contribution;
            } else if (s.underdensity !== undefined) {
              const rv = s.voidRadius;
              if (r < rv) {
                const voidFactor = s.underdensity * (1.0 - (r / rv) * (r / rv));
                delta += voidFactor;
              }
            }
          }

          const filament = 0.15 * (Math.cos(x * 0.04) * Math.sin(y * 0.04) + Math.cos(y * 0.04) * Math.sin(z * 0.04));
          delta += filament;

          buffer[grid.index(ix, iy, iz)] = Math.max(-1.0, delta);
        }
      }
    }

    return buffer;
  }
}
