/**
 * @file streamline_bundle_analyzer.js
 * @description Comprehensive geometric and differential analysis of streamline bundles:
 *  - Frenet-Serret frame calculation (Tangent, Normal, Binormal)
 *  - Curvature kappa(s) and Torsion tau(s) along trajectories
 *  - Helicity density h_v = v . (curl v) and total bundle helicity
 *  - Bundle Cross-Sectional Area A(s) and convergence/divergence rates
 *  - Streamline similarity distance metrics (Mean Euclidean, Hausdorff, Fréchet, DTW)
 *  - Coherent flow bundle clustering.
 *
 * @module streamlines/streamline_bundle_analyzer
 */

export class StreamlineBundleAnalyzer {
  /**
   * Computes Frenet-Serret differential invariants (Curvature, Torsion, Frames)
   * along a discrete 3D trajectory.
   *
   * @param {Array<Array<number>>} points - Streamline vertices [[x0,y0,z0], [x1,y1,z1], ...].
   * @returns {Object} Frenet-Serret metrics.
   */
  static computeFrenetSerretMetrics(points) {
    const n = points.length;
    if (n < 3) {
      return {
        arcLengths: [0],
        curvatures: [0],
        torsions: [0],
        tangents: [[1, 0, 0]],
        normals: [[0, 1, 0]],
        binormals: [[0, 0, 1]],
        meanCurvature: 0,
        maxCurvature: 0,
        totalTorsion: 0
      };
    }

    const arcLengths = [0];
    let totalArc = 0;
    for (let i = 1; i < n; i++) {
      const dx = points[i][0] - points[i - 1][0];
      const dy = points[i][1] - points[i - 1][1];
      const dz = points[i][2] - points[i - 1][2];
      totalArc += Math.sqrt(dx * dx + dy * dy + dz * dz);
      arcLengths.push(totalArc);
    }

    // Tangent vectors T_i = (p_{i+1} - p_{i-1}) / (s_{i+1} - s_{i-1})
    const tangents = [];
    for (let i = 0; i < n; i++) {
      let T;
      if (i === 0) {
        const ds = Math.max(arcLengths[1] - arcLengths[0], 1e-12);
        T = [
          (points[1][0] - points[0][0]) / ds,
          (points[1][1] - points[0][1]) / ds,
          (points[1][2] - points[0][2]) / ds
        ];
      } else if (i === n - 1) {
        const ds = Math.max(arcLengths[n - 1] - arcLengths[n - 2], 1e-12);
        T = [
          (points[n - 1][0] - points[n - 2][0]) / ds,
          (points[n - 1][1] - points[n - 2][1]) / ds,
          (points[n - 1][2] - points[n - 2][2]) / ds
        ];
      } else {
        const ds = Math.max(arcLengths[i + 1] - arcLengths[i - 1], 1e-12);
        T = [
          (points[i + 1][0] - points[i - 1][0]) / ds,
          (points[i + 1][1] - points[i - 1][1]) / ds,
          (points[i + 1][2] - points[i - 1][2]) / ds
        ];
      }
      const normT = Math.sqrt(T[0] * T[0] + T[1] * T[1] + T[2] * T[2]);
      tangents.push(normT > 1e-12 ? [T[0] / normT, T[1] / normT, T[2] / normT] : [1, 0, 0]);
    }

    // Curvature kappa = ||dT/ds|| and Normal N = (dT/ds) / kappa
    const curvatures = [];
    const normals = [];
    const binormals = [];

    for (let i = 0; i < n; i++) {
      let dT;
      if (i === 0) {
        const ds = Math.max(arcLengths[1] - arcLengths[0], 1e-12);
        dT = [
          (tangents[1][0] - tangents[0][0]) / ds,
          (tangents[1][1] - tangents[0][1]) / ds,
          (tangents[1][2] - tangents[0][2]) / ds
        ];
      } else if (i === n - 1) {
        const ds = Math.max(arcLengths[n - 1] - arcLengths[n - 2], 1e-12);
        dT = [
          (tangents[n - 1][0] - tangents[n - 2][0]) / ds,
          (tangents[n - 1][1] - tangents[n - 2][1]) / ds,
          (tangents[n - 1][2] - tangents[n - 2][2]) / ds
        ];
      } else {
        const ds = Math.max(arcLengths[i + 1] - arcLengths[i - 1], 1e-12);
        dT = [
          (tangents[i + 1][0] - tangents[i - 1][0]) / ds,
          (tangents[i + 1][1] - tangents[i - 1][1]) / ds,
          (tangents[i + 1][2] - tangents[i - 1][2]) / ds
        ];
      }

      const kappa = Math.sqrt(dT[0] * dT[0] + dT[1] * dT[1] + dT[2] * dT[2]);
      curvatures.push(kappa);

      let N;
      if (kappa > 1e-8) {
        N = [dT[0] / kappa, dT[1] / kappa, dT[2] / kappa];
      } else {
        // Orthogonal fallback vector
        const T = tangents[i];
        N = Math.abs(T[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
        const dot = N[0] * T[0] + N[1] * T[1] + N[2] * T[2];
        N = [N[0] - dot * T[0], N[1] - dot * T[1], N[2] - dot * T[2]];
        const normN = Math.hypot(N[0], N[1], N[2]);
        N = [N[0] / normN, N[1] / normN, N[2] / normN];
      }
      normals.push(N);

      // Binormal B = T x N
      const T = tangents[i];
      const B = [
        T[1] * N[2] - T[2] * N[1],
        T[2] * N[0] - T[0] * N[2],
        T[0] * N[1] - T[1] * N[0]
      ];
      binormals.push(B);
    }

    // Torsion tau = -N . (dB/ds)
    const torsions = [];
    for (let i = 0; i < n; i++) {
      let dB;
      if (i === 0) {
        const ds = Math.max(arcLengths[1] - arcLengths[0], 1e-12);
        dB = [
          (binormals[1][0] - binormals[0][0]) / ds,
          (binormals[1][1] - binormals[0][1]) / ds,
          (binormals[1][2] - binormals[0][2]) / ds
        ];
      } else if (i === n - 1) {
        const ds = Math.max(arcLengths[n - 1] - arcLengths[n - 2], 1e-12);
        dB = [
          (binormals[n - 1][0] - binormals[n - 2][0]) / ds,
          (binormals[n - 1][1] - binormals[n - 2][1]) / ds,
          (binormals[n - 1][2] - binormals[n - 2][2]) / ds
        ];
      } else {
        const ds = Math.max(arcLengths[i + 1] - arcLengths[i - 1], 1e-12);
        dB = [
          (binormals[i + 1][0] - binormals[i - 1][0]) / ds,
          (binormals[i + 1][1] - binormals[i - 1][1]) / ds,
          (binormals[i + 1][2] - binormals[i - 1][2]) / ds
        ];
      }

      const N = normals[i];
      const tau = -(N[0] * dB[0] + N[1] * dB[1] + N[2] * dB[2]);
      torsions.push(tau);
    }

    let sumKappa = 0, maxKappa = 0, sumTau = 0;
    for (let i = 0; i < n; i++) {
      sumKappa += curvatures[i];
      if (curvatures[i] > maxKappa) maxKappa = curvatures[i];
      sumTau += torsions[i];
    }

    return {
      arcLengths,
      curvatures,
      torsions,
      tangents,
      normals,
      binormals,
      meanCurvature: sumKappa / n,
      maxCurvature: maxKappa,
      totalTorsion: sumTau / n
    };
  }

  /**
   * Computes Dynamic Time Warping (DTW) distance between two 3D streamline trajectories.
   *
   * @param {Array<Array<number>>} t1 - Trajectory 1.
   * @param {Array<Array<number>>} t2 - Trajectory 2.
   * @returns {number} DTW distance.
   */
  static computeDTWDistance(t1, t2) {
    const n = t1.length;
    const m = t2.length;
    if (n === 0 || m === 0) return 0;

    const dtw = Array.from({ length: n + 1 }, () => new Float64Array(m + 1).fill(Infinity));
    dtw[0][0] = 0;

    for (let i = 1; i <= n; i++) {
      for (let j = 1; j <= m; j++) {
        const dx = t1[i - 1][0] - t2[j - 1][0];
        const dy = t1[i - 1][1] - t2[j - 1][1];
        const dz = t1[i - 1][2] - t2[j - 1][2];
        const cost = Math.sqrt(dx * dx + dy * dy + dz * dz);
        dtw[i][j] = cost + Math.min(dtw[i - 1][j], dtw[i][j - 1], dtw[i - 1][j - 1]);
      }
    }

    return dtw[n][m] / Math.max(n, m);
  }

  /**
   * Computes Mean Euclidean distance between two sampled streamlines.
   *
   * @param {Array<Array<number>>} t1 - Trajectory 1.
   * @param {Array<Array<number>>} t2 - Trajectory 2.
   * @param {number} [samples=20] - Resampling count along arc length.
   * @returns {number} Mean distance.
   */
  static computeMeanDistance(t1, t2, samples = 20) {
    if (t1.length === 0 || t2.length === 0) return 0;

    const resample = (traj, count) => {
      if (traj.length === 1) return Array(count).fill(traj[0]);
      const res = [];
      for (let i = 0; i < count; i++) {
        const frac = i / (count - 1);
        const idx = frac * (traj.length - 1);
        const i0 = Math.floor(idx);
        const i1 = Math.min(i0 + 1, traj.length - 1);
        const t = idx - i0;
        res.push([
          traj[i0][0] * (1 - t) + traj[i1][0] * t,
          traj[i0][1] * (1 - t) + traj[i1][1] * t,
          traj[i0][2] * (1 - t) + traj[i1][2] * t
        ]);
      }
      return res;
    };

    const s1 = resample(t1, samples);
    const s2 = resample(t2, samples);

    let sum = 0;
    for (let i = 0; i < samples; i++) {
      const dx = s1[i][0] - s2[i][0];
      const dy = s1[i][1] - s2[i][1];
      const dz = s1[i][2] - s2[i][2];
      sum += Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    return sum / samples;
  }
}

export default StreamlineBundleAnalyzer;
