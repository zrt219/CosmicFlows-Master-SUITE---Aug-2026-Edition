/**
 * @file mesh_smoother.js
 * @module surfaces/mesh_smoother
 * @description Advanced Mesh Smoothing and Surface Fairing Algorithms for 3D Cosmological Manifolds.
 * 
 * Implements:
 * 1. Uniform Umbrella Laplacian Smoothing (iterative local neighbor centroid relaxation).
 * 2. Cotangent-Weighted Laplace-Beltrami Operator (conformal geometric curvature-weighted smoothing).
 * 3. Taubin Non-Shrinking Dual-Step Filter (alternating positive lambda and negative mu steps to prevent volume shrinkage).
 * 4. Feature-Preserving Bilateral Mesh Denoising Filter (spatial and normal variance kernels).
 * 5. Strict Volume-Preservation Invariant Enforcement (rescaling by (V_init / V_curr)^(1/3)).
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Smoothing Algorithm Types.
 * @enum {string}
 */
export const SmoothingAlgorithm = Object.freeze({
  UNIFORM_LAPLACIAN: 'UNIFORM_LAPLACIAN',
  COTANGENT_LAPLACIAN: 'COTANGENT_LAPLACIAN',
  TAUBIN_NON_SHRINKING: 'TAUBIN_NON_SHRINKING',
  BILATERAL_FEATURE_PRESERVING: 'BILATERAL_FEATURE_PRESERVING'
});

/**
 * 3D Mesh Smoothing and Geometric Fairing Engine.
 */
export class MeshSmoother {
  /**
   * Smooth a triangulated 3D mesh.
   * 
   * @param {Float32Array} vertices - Interleaved 3D coordinates [x0, y0, z0, x1, ...].
   * @param {Uint32Array} indices - Triangle index triplets [i0, i1, i2, ...].
   * @param {Object} [options]
   * @param {string} [options.algorithm=SmoothingAlgorithm.TAUBIN_NON_SHRINKING]
   * @param {number} [options.iterations=5]
   * @param {number} [options.lambda=0.5] - Positive diffusion parameter (0 < lambda < 1).
   * @param {number} [options.mu=-0.53] - Negative inflation parameter (mu < -lambda).
   * @param {boolean} [options.preserveVolume=true] - Enforce volume preservation.
   * @returns {{vertices: Float32Array, normals: Float32Array}} Smoothed mesh vertices and recomputed normals.
   */
  static smooth(vertices, indices, options = {}) {
    const algorithm = options.algorithm || SmoothingAlgorithm.TAUBIN_NON_SHRINKING;
    const iterations = options.iterations !== undefined ? options.iterations : 5;
    const lambda = options.lambda !== undefined ? options.lambda : 0.5;
    const mu = options.mu !== undefined ? options.mu : -0.53;
    const preserveVolume = options.preserveVolume !== undefined ? options.preserveVolume : true;

    const vCount = vertices.length / 3;
    if (vCount === 0 || indices.length === 0 || iterations <= 0) {
      return {
        vertices: new Float32Array(vertices),
        normals: MeshSmoother.computeVertexNormals(vertices, indices)
      };
    }

    // Build vertex adjacency graph
    const adjacency = MeshSmoother._buildAdjacency(vertices, indices);

    // Initial enclosed volume
    const initialVolume = preserveVolume ? MeshSmoother._computeEnclosedVolume(vertices, indices) : 0.0;

    let currentVertices = new Float32Array(vertices);

    for (let iter = 0; iter < iterations; iter++) {
      switch (algorithm) {
        case SmoothingAlgorithm.UNIFORM_LAPLACIAN:
          currentVertices = MeshSmoother._laplacianStep(currentVertices, adjacency, lambda);
          break;

        case SmoothingAlgorithm.TAUBIN_NON_SHRINKING:
          // Step 1: Positive shrink step
          currentVertices = MeshSmoother._laplacianStep(currentVertices, adjacency, lambda);
          // Step 2: Negative inflate step
          currentVertices = MeshSmoother._laplacianStep(currentVertices, adjacency, mu);
          break;

        case SmoothingAlgorithm.COTANGENT_LAPLACIAN:
          currentVertices = MeshSmoother._cotangentStep(currentVertices, indices, adjacency, lambda);
          break;

        case SmoothingAlgorithm.BILATERAL_FEATURE_PRESERVING:
          currentVertices = MeshSmoother._bilateralStep(currentVertices, indices, adjacency, options);
          break;

        default:
          currentVertices = MeshSmoother._laplacianStep(currentVertices, adjacency, lambda);
          break;
      }
    }

    // Restore exact volume if required
    if (preserveVolume && initialVolume > 1e-8) {
      const currentVolume = MeshSmoother._computeEnclosedVolume(currentVertices, indices);
      if (currentVolume > 1e-8) {
        const scaleFactor = Math.cbrt(initialVolume / currentVolume);
        const centroid = MeshSmoother._computeCentroid(currentVertices);

        for (let i = 0; i < vCount; i++) {
          const idx = i * 3;
          currentVertices[idx + 0] = centroid[0] + (currentVertices[idx + 0] - centroid[0]) * scaleFactor;
          currentVertices[idx + 1] = centroid[1] + (currentVertices[idx + 1] - centroid[1]) * scaleFactor;
          currentVertices[idx + 2] = centroid[2] + (currentVertices[idx + 2] - centroid[2]) * scaleFactor;
        }
      }
    }

    const recomputedNormals = MeshSmoother.computeVertexNormals(currentVertices, indices);

    return {
      vertices: currentVertices,
      normals: recomputedNormals
    };
  }

  /**
   * Build vertex 1-ring neighbor adjacency list.
   * @private
   */
  static _buildAdjacency(vertices, indices) {
    const vCount = vertices.length / 3;
    const adj = Array.from({ length: vCount }, () => new Set());
    const fCount = indices.length / 3;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      adj[i0].add(i1);
      adj[i0].add(i2);
      adj[i1].add(i0);
      adj[i1].add(i2);
      adj[i2].add(i0);
      adj[i2].add(i1);
    }

    return adj.map(set => Array.from(set));
  }

  /**
   * Single Laplacian smoothing step: v_i' = v_i + factor * (sum(v_j)/k - v_i).
   * @private
   */
  static _laplacianStep(vertices, adjacency, factor) {
    const vCount = vertices.length / 3;
    const nextVertices = new Float32Array(vertices.length);

    for (let i = 0; i < vCount; i++) {
      const neighbors = adjacency[i];
      const k = neighbors.length;

      if (k === 0) {
        nextVertices[i * 3 + 0] = vertices[i * 3 + 0];
        nextVertices[i * 3 + 1] = vertices[i * 3 + 1];
        nextVertices[i * 3 + 2] = vertices[i * 3 + 2];
        continue;
      }

      let sumX = 0.0;
      let sumY = 0.0;
      let sumZ = 0.0;

      for (let n = 0; n < k; n++) {
        const nIdx = neighbors[n] * 3;
        sumX += vertices[nIdx + 0];
        sumY += vertices[nIdx + 1];
        sumZ += vertices[nIdx + 2];
      }

      const avgX = sumX / k;
      const avgY = sumY / k;
      const avgZ = sumZ / k;

      const vx = vertices[i * 3 + 0];
      const vy = vertices[i * 3 + 1];
      const vz = vertices[i * 3 + 2];

      nextVertices[i * 3 + 0] = vx + factor * (avgX - vx);
      nextVertices[i * 3 + 1] = vy + factor * (avgY - vy);
      nextVertices[i * 3 + 2] = vz + factor * (avgZ - vz);
    }

    return nextVertices;
  }

  /**
   * Cotangent-weighted Laplace-Beltrami operator step.
   * @private
   */
  static _cotangentStep(vertices, indices, adjacency, lambda) {
    const vCount = vertices.length / 3;
    const nextVertices = new Float32Array(vertices.length);

    // Build map of directed edges to opposite angles
    const fCount = indices.length / 3;
    const edgeWeights = new Map(); // key `${a}_${b}` -> cotangent weight sum

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const p0 = [vertices[i0 * 3], vertices[i0 * 3 + 1], vertices[i0 * 3 + 2]];
      const p1 = [vertices[i1 * 3], vertices[i1 * 3 + 1], vertices[i1 * 3 + 2]];
      const p2 = [vertices[i2 * 3], vertices[i2 * 3 + 1], vertices[i2 * 3 + 2]];

      const calcCot = (a, b, c) => {
        // Cotangent of angle at vertex c opposite edge (a, b)
        const v1 = [a[0] - c[0], a[1] - c[1], a[2] - c[2]];
        const v2 = [b[0] - c[0], b[1] - c[1], b[2] - c[2]];
        const dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        const cross = Math.hypot(
          v1[1] * v2[2] - v1[2] * v2[1],
          v1[2] * v2[0] - v1[0] * v2[2],
          v1[0] * v2[1] - v1[1] * v2[0]
        );
        return cross > 1e-8 ? dot / cross : 0.0;
      };

      const cot0 = calcCot(p1, p2, p0); // opp to edge 1-2
      const cot1 = calcCot(p2, p0, p1); // opp to edge 2-0
      const cot2 = calcCot(p0, p1, p2); // opp to edge 0-1

      const addWeight = (u, v, w) => {
        const key = `${Math.min(u, v)}_${Math.max(u, v)}`;
        edgeWeights.set(key, (edgeWeights.get(key) || 0.0) + Math.max(0.0, w));
      };

      addWeight(i1, i2, cot0);
      addWeight(i2, i0, cot1);
      addWeight(i0, i1, cot2);
    }

    for (let i = 0; i < vCount; i++) {
      const neighbors = adjacency[i];
      let totalW = 0.0;
      let dispX = 0.0, dispY = 0.0, dispZ = 0.0;

      const vx = vertices[i * 3 + 0];
      const vy = vertices[i * 3 + 1];
      const vz = vertices[i * 3 + 2];

      for (const n of neighbors) {
        const key = `${Math.min(i, n)}_${Math.max(i, n)}`;
        const w = edgeWeights.get(key) || 1.0;
        totalW += w;

        dispX += w * (vertices[n * 3 + 0] - vx);
        dispY += w * (vertices[n * 3 + 1] - vy);
        dispZ += w * (vertices[n * 3 + 2] - vz);
      }

      if (totalW > 1e-8) {
        nextVertices[i * 3 + 0] = vx + (lambda * dispX) / totalW;
        nextVertices[i * 3 + 1] = vy + (lambda * dispY) / totalW;
        nextVertices[i * 3 + 2] = vz + (lambda * dispZ) / totalW;
      } else {
        nextVertices[i * 3 + 0] = vx;
        nextVertices[i * 3 + 1] = vy;
        nextVertices[i * 3 + 2] = vz;
      }
    }

    return nextVertices;
  }

  /**
   * Bilateral feature-preserving mesh smoothing step.
   * @private
   */
  static _bilateralStep(vertices, indices, adjacency, options) {
    const vCount = vertices.length / 3;
    const normals = MeshSmoother.computeVertexNormals(vertices, indices);
    const nextVertices = new Float32Array(vertices.length);

    const sigmaS = options.spatialSigma || 5.0;
    const sigmaS2 = 2.0 * sigmaS * sigmaS;
    const sigmaN = options.normalSigma || 0.3;
    const sigmaN2 = 2.0 * sigmaN * sigmaN;

    for (let i = 0; i < vCount; i++) {
      const vx = vertices[i * 3 + 0];
      const vy = vertices[i * 3 + 1];
      const vz = vertices[i * 3 + 2];

      const nx = normals[i * 3 + 0];
      const ny = normals[i * 3 + 1];
      const nz = normals[i * 3 + 2];

      const neighbors = adjacency[i];
      let sumWeights = 0.0;
      let totalOffset = 0.0;

      for (const n of neighbors) {
        const qx = vertices[n * 3 + 0];
        const qy = vertices[n * 3 + 1];
        const qz = vertices[n * 3 + 2];

        const dx = qx - vx;
        const dy = qy - vy;
        const dz = qz - vz;
        const dist2 = dx * dx + dy * dy + dz * dz;

        const h = dx * nx + dy * ny + dz * nz; // normal projection
        const wS = Math.exp(-dist2 / sigmaS2);
        const wN = Math.exp(-(h * h) / sigmaN2);
        const w = wS * wN;

        sumWeights += w;
        totalOffset += w * h;
      }

      if (sumWeights > 1e-8) {
        const d = totalOffset / sumWeights;
        nextVertices[i * 3 + 0] = vx + d * nx;
        nextVertices[i * 3 + 1] = vy + d * ny;
        nextVertices[i * 3 + 2] = vz + d * nz;
      } else {
        nextVertices[i * 3 + 0] = vx;
        nextVertices[i * 3 + 1] = vy;
        nextVertices[i * 3 + 2] = vz;
      }
    }

    return nextVertices;
  }

  /**
   * Compute smooth area-weighted vertex normals.
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {Float32Array}
   */
  static computeVertexNormals(vertices, indices) {
    const vCount = vertices.length / 3;
    const normals = new Float32Array(vertices.length);
    const fCount = indices.length / 3;

    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0];
      const i1 = indices[f * 3 + 1];
      const i2 = indices[f * 3 + 2];

      const x0 = vertices[i0 * 3 + 0], y0 = vertices[i0 * 3 + 1], z0 = vertices[i0 * 3 + 2];
      const x1 = vertices[i1 * 3 + 0], y1 = vertices[i1 * 3 + 1], z1 = vertices[i1 * 3 + 2];
      const x2 = vertices[i2 * 3 + 0], y2 = vertices[i2 * 3 + 1], z2 = vertices[i2 * 3 + 2];

      const e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
      const e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;

      const cx = e1y * e2z - e1z * e2y;
      const cy = e1z * e2x - e1x * e2z;
      const cz = e1x * e2y - e1y * e2x;

      normals[i0 * 3 + 0] += cx; normals[i0 * 3 + 1] += cy; normals[i0 * 3 + 2] += cz;
      normals[i1 * 3 + 0] += cx; normals[i1 * 3 + 1] += cy; normals[i1 * 3 + 2] += cz;
      normals[i2 * 3 + 0] += cx; normals[i2 * 3 + 1] += cy; normals[i2 * 3 + 2] += cz;
    }

    for (let i = 0; i < vCount; i++) {
      const idx = i * 3;
      const len = Math.hypot(normals[idx + 0], normals[idx + 1], normals[idx + 2]);
      if (len > 1e-8) {
        normals[idx + 0] /= len;
        normals[idx + 1] /= len;
        normals[idx + 2] /= len;
      } else {
        normals[idx + 0] = 0.0;
        normals[idx + 1] = 0.0;
        normals[idx + 2] = 1.0;
      }
    }

    return normals;
  }

  /**
   * Compute signed enclosed volume via divergence theorem.
   * @private
   */
  static _computeEnclosedVolume(vertices, indices) {
    let vol = 0.0;
    const fCount = indices.length / 3;
    for (let f = 0; f < fCount; f++) {
      const i0 = indices[f * 3 + 0] * 3;
      const i1 = indices[f * 3 + 1] * 3;
      const i2 = indices[f * 3 + 2] * 3;

      const x0 = vertices[i0], y0 = vertices[i0 + 1], z0 = vertices[i0 + 2];
      const x1 = vertices[i1], y1 = vertices[i1 + 1], z1 = vertices[i1 + 2];
      const x2 = vertices[i2], y2 = vertices[i2 + 1], z2 = vertices[i2 + 2];

      vol += (x0 * (y1 * z2 - z1 * y2) - y0 * (x1 * z2 - z1 * x2) + z0 * (x1 * y2 - y1 * x2)) / 6.0;
    }
    return Math.abs(vol);
  }

  /**
   * Compute vertex centroid [x, y, z].
   * @private
   */
  static _computeCentroid(vertices) {
    const vCount = vertices.length / 3;
    let sumX = 0, sumY = 0, sumZ = 0;
    for (let i = 0; i < vCount; i++) {
      sumX += vertices[i * 3 + 0];
      sumY += vertices[i * 3 + 1];
      sumZ += vertices[i * 3 + 2];
    }
    return [sumX / vCount, sumY / vCount, sumZ / vCount];
  }
}
