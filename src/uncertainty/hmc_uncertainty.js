/**
 * @file hmc_uncertainty.js
 * @description CF4++ Hamiltonian Monte Carlo (HMC) Posterior Mean, RMS Grid Evaluator, and Uncertainty Propagator.
 * Implements:
 * 1. Streaming Welford algorithm for online posterior mean, RMS, and 3x3 covariance tensor calculation over 10,000 HMC steps
 * 2. Coordinate-to-velocity and velocity-to-coordinate uncertainty propagation via Jacobian tensors (Σ_v = J Σ_x J^T)
 * 3. Covariance error ellipsoid diagonalization, confidence scaling (68.3%, 95.4%, 99.7%), and fractional anisotropy
 */

import { jacobiDiagonalize3x3 } from '../topology/eigen_topology.js';
import { invert3x3 } from '../topology/root_finder.js';

/**
 * Streaming Welford accumulator for computing high-precision posterior mean, variance, RMS,
 * and 3x3 covariance across thousands of 3D grid realizations without high memory overhead.
 */
export class StreamingHMCAccumulator {
  /**
   * @param {number} totalVoxels - Total number of voxels N = nx * ny * nz
   */
  constructor(totalVoxels) {
    this.totalVoxels = totalVoxels;
    this.sampleCount = 0;

    // Running mean for Vx, Vy, Vz
    this.meanVx = new Float32Array(totalVoxels);
    this.meanVy = new Float32Array(totalVoxels);
    this.meanVz = new Float32Array(totalVoxels);

    // Sum of squared differences (Welford M2)
    this.M2_xx = new Float32Array(totalVoxels);
    this.M2_yy = new Float32Array(totalVoxels);
    this.M2_zz = new Float32Array(totalVoxels);
    this.M2_xy = new Float32Array(totalVoxels);
    this.M2_xz = new Float32Array(totalVoxels);
    this.M2_yz = new Float32Array(totalVoxels);
  }

  /**
   * Updates running statistics with a new HMC realization velocity grid buffer.
   * @param {Float32Array} realizationGrid - Velocity buffer [totalVoxels * 3]
   */
  update(realizationGrid) {
    this.sampleCount++;
    const n = this.sampleCount;
    const invN = 1.0 / n;

    for (let idx = 0; idx < this.totalVoxels; idx++) {
      const vx = realizationGrid[idx * 3];
      const vy = realizationGrid[idx * 3 + 1];
      const vz = realizationGrid[idx * 3 + 2];

      const dx1 = vx - this.meanVx[idx];
      const dy1 = vy - this.meanVy[idx];
      const dz1 = vz - this.meanVz[idx];

      this.meanVx[idx] += dx1 * invN;
      this.meanVy[idx] += dy1 * invN;
      this.meanVz[idx] += dz1 * invN;

      const dx2 = vx - this.meanVx[idx];
      const dy2 = vy - this.meanVy[idx];
      const dz2 = vz - this.meanVz[idx];

      this.M2_xx[idx] += dx1 * dx2;
      this.M2_yy[idx] += dy1 * dy2;
      this.M2_zz[idx] += dz1 * dz2;
      this.M2_xy[idx] += dx1 * dy2;
      this.M2_xz[idx] += dx1 * dz2;
      this.M2_yz[idx] += dy1 * dz2;
    }
  }

  /**
   * Finalizes and extracts posterior mean and RMS grids.
   * @returns {{
   *   meanGrid: Float32Array,
   *   rmsGrid: Float32Array,
   *   sampleCount: number
   * }}
   */
  finalize() {
    if (this.sampleCount < 2) {
      throw new Error('Need at least 2 HMC realization samples to finalize variance.');
    }

    const denom = this.sampleCount - 1;
    const meanGrid = new Float32Array(this.totalVoxels * 3);
    const rmsGrid = new Float32Array(this.totalVoxels * 3);

    for (let idx = 0; idx < this.totalVoxels; idx++) {
      meanGrid[idx * 3] = this.meanVx[idx];
      meanGrid[idx * 3 + 1] = this.meanVy[idx];
      meanGrid[idx * 3 + 2] = this.meanVz[idx];

      rmsGrid[idx * 3] = Math.sqrt(Math.max(0, this.M2_xx[idx] / denom));
      rmsGrid[idx * 3 + 1] = Math.sqrt(Math.max(0, this.M2_yy[idx] / denom));
      rmsGrid[idx * 3 + 2] = Math.sqrt(Math.max(0, this.M2_zz[idx] / denom));
    }

    return {
      meanGrid,
      rmsGrid,
      sampleCount: this.sampleCount
    };
  }

  /**
   * Gets the 3x3 covariance matrix at a specific voxel index.
   * @param {number} idx
   * @returns {number[][]}
   */
  getVoxelCovariance(idx) {
    if (this.sampleCount < 2) return [[0,0,0],[0,0,0],[0,0,0]];
    const denom = this.sampleCount - 1;

    const cxx = this.M2_xx[idx] / denom;
    const cyy = this.M2_yy[idx] / denom;
    const czz = this.M2_zz[idx] / denom;
    const cxy = this.M2_xy[idx] / denom;
    const cxz = this.M2_xz[idx] / denom;
    const cyz = this.M2_yz[idx] / denom;

    return [
      [cxx, cxy, cxz],
      [cxy, cyy, cyz],
      [cxz, cyz, czz]
    ];
  }
}

/**
 * Propagates spatial position uncertainty Σ_x into velocity dispersion Σ_v through the Jacobian tensor J:
 * Σ_v = J * Σ_x * J^T + Σ_intrinsic
 *
 * @param {number[][]} J - 3x3 Jacobian tensor dv_i / dx_j
 * @param {number[][]} sigmaX - 3x3 spatial covariance matrix in (Mpc/h)^2
 * @param {number[][]} [sigmaIntrinsic] - Optional intrinsic 3x3 velocity dispersion (km/s)^2
 * @returns {number[][]} 3x3 velocity covariance matrix Σ_v in (km/s)^2
 */
export function propagateSpatialToVelocityDispersion(J, sigmaX, sigmaIntrinsic = null) {
  // Compute M = J * sigmaX
  const M = [
    [
      J[0][0]*sigmaX[0][0] + J[0][1]*sigmaX[1][0] + J[0][2]*sigmaX[2][0],
      J[0][0]*sigmaX[0][1] + J[0][1]*sigmaX[1][1] + J[0][2]*sigmaX[2][1],
      J[0][0]*sigmaX[0][2] + J[0][1]*sigmaX[1][2] + J[0][2]*sigmaX[2][2]
    ],
    [
      J[1][0]*sigmaX[0][0] + J[1][1]*sigmaX[1][0] + J[1][2]*sigmaX[2][0],
      J[1][0]*sigmaX[0][1] + J[1][1]*sigmaX[1][1] + J[1][2]*sigmaX[2][1],
      J[1][0]*sigmaX[0][2] + J[1][1]*sigmaX[1][2] + J[1][2]*sigmaX[2][2]
    ],
    [
      J[2][0]*sigmaX[0][0] + J[2][1]*sigmaX[1][0] + J[2][2]*sigmaX[2][0],
      J[2][0]*sigmaX[0][1] + J[2][1]*sigmaX[1][1] + J[2][2]*sigmaX[2][1],
      J[2][0]*sigmaX[0][2] + J[2][1]*sigmaX[1][2] + J[2][2]*sigmaX[2][2]
    ]
  ];

  // Compute Σ_v = M * J^T
  const sigmaV = [
    [
      M[0][0]*J[0][0] + M[0][1]*J[0][1] + M[0][2]*J[0][2],
      M[0][0]*J[1][0] + M[0][1]*J[1][1] + M[0][2]*J[1][2],
      M[0][0]*J[2][0] + M[0][1]*J[2][1] + M[0][2]*J[2][2]
    ],
    [
      M[1][0]*J[0][0] + M[1][1]*J[0][1] + M[1][2]*J[0][2],
      M[1][0]*J[1][0] + M[1][1]*J[1][1] + M[1][2]*J[1][2],
      M[1][0]*J[2][0] + M[1][1]*J[2][1] + M[1][2]*J[2][2]
    ],
    [
      M[2][0]*J[0][0] + M[2][1]*J[0][1] + M[2][2]*J[0][2],
      M[2][0]*J[1][0] + M[2][1]*J[1][1] + M[2][2]*J[1][2],
      M[2][0]*J[2][0] + M[2][1]*J[2][1] + M[2][2]*J[2][2]
    ]
  ];

  if (sigmaIntrinsic) {
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        sigmaV[i][j] += sigmaIntrinsic[i][j];
      }
    }
  }

  return sigmaV;
}

/**
 * Propagates observational velocity uncertainty Σ_v into critical attractor position uncertainty Σ_x:
 * Σ_x = J^-1 * Σ_v * (J^-1)^T
 *
 * @param {number[][]} J - 3x3 Jacobian tensor dv_i / dx_j
 * @param {number[][]} sigmaV - 3x3 velocity covariance matrix in (km/s)^2
 * @returns {number[][]} 3x3 spatial covariance matrix Σ_x in (Mpc/h)^2
 */
export function propagateVelocityToSpatialDispersion(J, sigmaV) {
  const invJ = invert3x3(J);
  if (!invJ) {
    throw new Error('Singular Jacobian: cannot invert to propagate velocity uncertainty to spatial uncertainty.');
  }

  // M = invJ * sigmaV
  const M = [
    [
      invJ[0][0]*sigmaV[0][0] + invJ[0][1]*sigmaV[1][0] + invJ[0][2]*sigmaV[2][0],
      invJ[0][0]*sigmaV[0][1] + invJ[0][1]*sigmaV[1][1] + invJ[0][2]*sigmaV[2][1],
      invJ[0][0]*sigmaV[0][2] + invJ[0][1]*sigmaV[1][2] + invJ[0][2]*sigmaV[2][2]
    ],
    [
      invJ[1][0]*sigmaV[0][0] + invJ[1][1]*sigmaV[1][0] + invJ[1][2]*sigmaV[2][0],
      invJ[1][0]*sigmaV[0][1] + invJ[1][1]*sigmaV[1][1] + invJ[1][2]*sigmaV[2][1],
      invJ[1][0]*sigmaV[0][2] + invJ[1][1]*sigmaV[1][2] + invJ[1][2]*sigmaV[2][2]
    ],
    [
      invJ[2][0]*sigmaV[0][0] + invJ[2][1]*sigmaV[1][0] + invJ[2][2]*sigmaV[2][0],
      invJ[2][0]*sigmaV[0][1] + invJ[2][1]*sigmaV[1][1] + invJ[2][2]*sigmaV[2][1],
      invJ[2][0]*sigmaV[0][2] + invJ[2][1]*sigmaV[1][2] + invJ[2][2]*sigmaV[2][2]
    ]
  ];

  // Σ_x = M * invJ^T
  return [
    [
      M[0][0]*invJ[0][0] + M[0][1]*invJ[0][1] + M[0][2]*invJ[0][2],
      M[0][0]*invJ[1][0] + M[0][1]*invJ[1][1] + M[0][2]*invJ[1][2],
      M[0][0]*invJ[2][0] + M[0][1]*invJ[2][1] + M[0][2]*invJ[2][2]
    ],
    [
      M[1][0]*invJ[0][0] + M[1][1]*invJ[0][1] + M[1][2]*invJ[0][2],
      M[1][0]*invJ[1][0] + M[1][1]*invJ[1][1] + M[1][2]*invJ[1][2],
      M[1][0]*invJ[2][0] + M[1][1]*invJ[2][1] + M[1][2]*invJ[2][2]
    ],
    [
      M[2][0]*invJ[0][0] + M[2][1]*invJ[0][1] + M[2][2]*invJ[0][2],
      M[2][0]*invJ[1][0] + M[2][1]*invJ[1][1] + M[2][2]*invJ[1][2],
      M[2][0]*invJ[2][0] + M[2][1]*invJ[2][1] + M[2][2]*invJ[2][2]
    ]
  ];
}

/**
 * Diagonalizes a 3x3 covariance matrix and computes 3D error ellipsoid geometry.
 *
 * @param {number[][]} cov - 3x3 positive semi-definite covariance matrix
 * @param {'1sigma'|'2sigma'|'3sigma'} [confidenceLevel='1sigma']
 * @returns {{
 *   semiMajorAxis: number,
 *   semiIntermediateAxis: number,
 *   semiMinorAxis: number,
 *   eigenvectors: number[][],
 *   ellipsoidVolume: number,
 *   fractionalAnisotropy: number,
 *   confidenceScale: number
 * }}
 */
export function computeCovarianceEllipsoid(cov, confidenceLevel = '1sigma') {
  const { eigenvalues, eigenvectors } = jacobiDiagonalize3x3(cov);

  // Eigenvalues are sorted ascending: l1 <= l2 <= l3
  const varMin = Math.max(0, eigenvalues[0]);
  const varInt = Math.max(0, eigenvalues[1]);
  const varMaj = Math.max(0, eigenvalues[2]);

  // Chi-square critical values for 3 degrees of freedom:
  // 1-sigma (68.27%): sqrt(3.53) ≈ 1.8788
  // 2-sigma (95.45%): sqrt(8.02) ≈ 2.832
  // 3-sigma (99.73%): sqrt(14.16) ≈ 3.763
  let scale = 1.0;
  if (confidenceLevel === '1sigma') scale = Math.sqrt(3.53);
  else if (confidenceLevel === '2sigma') scale = Math.sqrt(8.02);
  else if (confidenceLevel === '3sigma') scale = Math.sqrt(14.16);

  const a = Math.sqrt(varMaj) * scale;
  const b = Math.sqrt(varInt) * scale;
  const c = Math.sqrt(varMin) * scale;

  const volume = (4.0 / 3.0) * Math.PI * a * b * c;

  // Fractional Anisotropy FA = sqrt(3/2) * sqrt((l1-l_mean)^2 + (l2-l_mean)^2 + (l3-l_mean)^2) / sqrt(l1^2 + l2^2 + l3^2)
  const lMean = (varMin + varInt + varMaj) / 3.0;
  const num = (varMin - lMean)**2 + (varInt - lMean)**2 + (varMaj - lMean)**2;
  const den = varMin**2 + varInt**2 + varMaj**2;
  const fa = den > 0 ? Math.sqrt(1.5) * (Math.sqrt(num) / Math.sqrt(den)) : 0;

  return {
    semiMajorAxis: a,
    semiIntermediateAxis: b,
    semiMinorAxis: c,
    eigenvectors,
    ellipsoidVolume: volume,
    fractionalAnisotropy: fa,
    confidenceScale: scale
  };
}
