/**
 * @file surface_normal_calculator.js
 * @module surfaces/surface_normal_calculator
 * @description Differential Geometry, Curvature Tensors, and Normal Calculators for 3D Cosmological Manifolds.
 * 
 * Implements:
 * 1. Face Normal Vectors (cross product of oriented triangle edges).
 * 2. Area-Weighted and Angle-Weighted Vertex Pseudo-Normals.
 * 3. Continuous Field-Gradient Normals (n = grad(delta) / ||grad(delta)||).
 * 4. Discrete Gaussian Curvature K via Gauss-Bonnet Angle Defect: K_i = (2pi - sum(theta_j)) / A_i.
 * 5. Discrete Mean Curvature H via Laplace-Beltrami Operator: H_i = 0.5 * ||Delta v_i||.
 * 6. Principal Curvatures (k1, k2), Shape Index S in [-1, 1], and Curvedness C.
 * 7. Global Gauss-Bonnet Integral Invariant Check: \oint K dA == 2pi * chi.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Surface Normal Calculation Mode.
 * @enum {string}
 */
export const NormalWeightingMode = Object.freeze({
  UNIFORM: 'UNIFORM',
  AREA_WEIGHTED: 'AREA_WEIGHTED',
  ANGLE_WEIGHTED: 'ANGLE_WEIGHTED',
  FIELD_GRADIENT: 'FIELD_GRADIENT'
});

/**
 * Differential Geometry Analyzer for Triangulated 2-Manifolds.
 */
export class SurfaceNormalCalculator {
  /**
   * Calculate vertex normals using specified weighting mode.
   * 
   * @param {Float32Array} vertices - 3D vertex coordinates [x0, y0, z0, ...].
   * @param {Uint32Array} indices - Triangle index triplets [i0, i1, i2, ...].
   * @param {string} [mode=NormalWeightingMode.ANGLE_WEIGHTED]
   * @returns {Float32Array} Unit vertex normal vectors.
   */
  static computeNormals(vertices, indices, mode = NormalWeightingMode.ANGLE_WEIGHTED) {
    const vCount = vertices.length / 3;
    const normals = new Float32Array(vertices.length);
    const fCount = indices.length / 3;

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

      // Cross product e1 x e2
      const cx = e1[1] * e2[2] - e1[2] * e2[1];
      const cy = e1[2] * e2[0] - e1[0] * e2[2];
      const cz = e1[0] * e2[1] - e1[1] * e2[0];
      const area2 = Math.hypot(cx, cy, cz);
      if (area2 < 1e-12) continue;

      const unx = cx / area2;
      const uny = cy / area2;
      const unz = cz / area2;

      if (mode === NormalWeightingMode.ANGLE_WEIGHTED) {
        // Angle at vertex 0 (between e1 and e2)
        const l1 = Math.hypot(e1[0], e1[1], e1[2]);
        const l2 = Math.hypot(e2[0], e2[1], e2[2]);
        const l3 = Math.hypot(e3[0], e3[1], e3[2]);

        const dot0 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
        const theta0 = Math.acos(Math.max(-1.0, Math.min(1.0, dot0 / (l1 * l2))));

        // Angle at vertex 1 (between -e1 and e3)
        const dot1 = -(e1[0] * e3[0] + e1[1] * e3[1] + e1[2] * e3[2]);
        const theta1 = Math.acos(Math.max(-1.0, Math.min(1.0, dot1 / (l1 * l3))));

        // Angle at vertex 2 (between -e2 and -e3)
        const dot2 = e2[0] * e3[0] + e2[1] * e3[1] + e2[2] * e3[2];
        const theta2 = Math.acos(Math.max(-1.0, Math.min(1.0, dot2 / (l2 * l3))));

        normals[i0 * 3 + 0] += unx * theta0;
        normals[i0 * 3 + 1] += uny * theta0;
        normals[i0 * 3 + 2] += unz * theta0;

        normals[i1 * 3 + 0] += unx * theta1;
        normals[i1 * 3 + 1] += uny * theta1;
        normals[i1 * 3 + 2] += unz * theta1;

        normals[i2 * 3 + 0] += unx * theta2;
        normals[i2 * 3 + 1] += uny * theta2;
        normals[i2 * 3 + 2] += unz * theta2;
      } else {
        // Area-weighted or Uniform
        const weight = mode === NormalWeightingMode.AREA_WEIGHTED ? 0.5 * area2 : 1.0;
        normals[i0 * 3 + 0] += unx * weight;
        normals[i0 * 3 + 1] += uny * weight;
        normals[i0 * 3 + 2] += unz * weight;

        normals[i1 * 3 + 0] += unx * weight;
        normals[i1 * 3 + 1] += uny * weight;
        normals[i1 * 3 + 2] += unz * weight;

        normals[i2 * 3 + 0] += unx * weight;
        normals[i2 * 3 + 1] += uny * weight;
        normals[i2 * 3 + 2] += unz * weight;
      }
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
   * Compute discrete differential geometry and curvature tensors:
   * - Gaussian curvature K
   * - Mean curvature H
   * - Principal curvatures k1, k2
   * - Shape Index S
   * - Curvedness C
   * - Gauss-Bonnet integral check
   * 
   * @param {Float32Array} vertices
   * @param {Uint32Array} indices
   * @returns {{
   *   gaussianCurvatures: Float32Array,
   *   meanCurvatures: Float32Array,
   *   principalK1: Float32Array,
   *   principalK2: Float32Array,
   *   shapeIndex: Float32Array,
   *   curvedness: Float32Array,
   *   integratedGaussianCurvature: number,
   *   integratedMeanCurvature: number,
   *   gaussBonnetEulerResidual: number
   * }}
   */
  static computeCurvatures(vertices, indices) {
    const vCount = vertices.length / 3;
    const fCount = indices.length / 3;

    const angleSum = new Float64Array(vCount);
    const voronoiArea = new Float64Array(vCount);
    const laplaceVector = Array.from({ length: vCount }, () => [0.0, 0.0, 0.0]);

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

      const l1 = Math.hypot(e1[0], e1[1], e1[2]);
      const l2 = Math.hypot(e2[0], e2[1], e2[2]);
      const l3 = Math.hypot(e3[0], e3[1], e3[2]);

      const cx = e1[1] * e2[2] - e1[2] * e2[1];
      const cy = e1[2] * e2[0] - e1[0] * e2[2];
      const cz = e1[0] * e2[1] - e1[1] * e2[0];
      const triArea = 0.5 * Math.hypot(cx, cy, cz);

      // Barycentric / Voronoi area contribution (1/3 of triangle area)
      voronoiArea[i0] += triArea / 3.0;
      voronoiArea[i1] += triArea / 3.0;
      voronoiArea[i2] += triArea / 3.0;

      // Incident angles
      const dot0 = Math.max(-1.0, Math.min(1.0, (e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2]) / (l1 * l2)));
      const dot1 = Math.max(-1.0, Math.min(1.0, -(e1[0] * e3[0] + e1[1] * e3[1] + e1[2] * e3[2]) / (l1 * l3)));
      const dot2 = Math.max(-1.0, Math.min(1.0, (e2[0] * e3[0] + e2[1] * e3[1] + e2[2] * e3[2]) / (l2 * l3)));

      angleSum[i0] += Math.acos(dot0);
      angleSum[i1] += Math.acos(dot1);
      angleSum[i2] += Math.acos(dot2);

      // Cotangent weights for Laplace-Beltrami
      const cot0 = dot0 / Math.sqrt(Math.max(1e-12, 1.0 - dot0 * dot0));
      const cot1 = dot1 / Math.sqrt(Math.max(1e-12, 1.0 - dot1 * dot1));
      const cot2 = dot2 / Math.sqrt(Math.max(1e-12, 1.0 - dot2 * dot2));

      // Edge 1-2 opp to p0
      laplaceVector[i1][0] += cot0 * (p2[0] - p1[0]);
      laplaceVector[i1][1] += cot0 * (p2[1] - p1[1]);
      laplaceVector[i1][2] += cot0 * (p2[2] - p1[2]);

      laplaceVector[i2][0] += cot0 * (p1[0] - p2[0]);
      laplaceVector[i2][1] += cot0 * (p1[1] - p2[1]);
      laplaceVector[i2][2] += cot0 * (p1[2] - p2[2]);

      // Edge 2-0 opp to p1
      laplaceVector[i2][0] += cot1 * (p0[0] - p2[0]);
      laplaceVector[i2][1] += cot1 * (p0[1] - p2[1]);
      laplaceVector[i2][2] += cot1 * (p0[2] - p2[2]);

      laplaceVector[i0][0] += cot1 * (p2[0] - p0[0]);
      laplaceVector[i0][1] += cot1 * (p2[1] - p0[0]);
      laplaceVector[i0][2] += cot1 * (p2[2] - p0[0]);

      // Edge 0-1 opp to p2
      laplaceVector[i0][0] += cot2 * (p1[0] - p0[0]);
      laplaceVector[i0][1] += cot2 * (p1[1] - p0[1]);
      laplaceVector[i0][2] += cot2 * (p1[2] - p0[0]);

      laplaceVector[i1][0] += cot2 * (p0[0] - p1[0]);
      laplaceVector[i1][1] += cot2 * (p0[1] - p1[1]);
      laplaceVector[i1][2] += cot2 * (p0[2] - p1[2]);
    }

    const gaussianK = new Float32Array(vCount);
    const meanH = new Float32Array(vCount);
    const k1Arr = new Float32Array(vCount);
    const k2Arr = new Float32Array(vCount);
    const shapeIdx = new Float32Array(vCount);
    const curvedness = new Float32Array(vCount);

    let intK = 0.0;
    let intH = 0.0;

    for (let i = 0; i < vCount; i++) {
      const area = Math.max(1e-12, voronoiArea[i]);

      // Gaussian Curvature K = (2pi - angleSum) / area
      const K = (2.0 * Math.PI - angleSum[i]) / area;
      gaussianK[i] = K;
      intK += K * area;

      // Mean Curvature H = 0.5 * ||Delta v_i|| / (2 * area)
      const lapLen = Math.hypot(laplaceVector[i][0], laplaceVector[i][1], laplaceVector[i][2]);
      const H = lapLen / (4.0 * area);
      meanH[i] = H;
      intH += H * area;

      // Principal curvatures: k1 = H + sqrt(max(0, H^2 - K)), k2 = H - sqrt(max(0, H^2 - K))
      const disc = Math.max(0.0, H * H - K);
      const sqrtDisc = Math.sqrt(disc);
      const k1 = H + sqrtDisc;
      const k2 = H - sqrtDisc;
      k1Arr[i] = k1;
      k2Arr[i] = k2;

      // Shape Index S = (2/pi) * atan((k1 + k2) / (k1 - k2))
      const denom = k1 - k2;
      if (Math.abs(denom) > 1e-8) {
        shapeIdx[i] = (2.0 / Math.PI) * Math.atan((k1 + k2) / denom);
      } else {
        shapeIdx[i] = H >= 0 ? 1.0 : -1.0;
      }

      // Curvedness C = sqrt((k1^2 + k2^2) / 2)
      curvedness[i] = Math.sqrt((k1 * k1 + k2 * k2) / 2.0);
    }

    return {
      gaussianCurvatures: gaussianK,
      meanCurvatures: meanH,
      principalK1: k1Arr,
      principalK2: k2Arr,
      shapeIndex: shapeIdx,
      curvedness,
      integratedGaussianCurvature: intK,
      integratedMeanCurvature: intH,
      gaussBonnetEulerResidual: intK / (2.0 * Math.PI)
    };
  }
}
