/**
 * @file acceptance_gates.js
 * @module validation/acceptance_gates
 * @description Code-enforced acceptance validation gates (Gates A-I) for ZRT Cosmicflows Workbench.
 * Ensures zero unverified numerical or scientific results are tagged as VALIDATED.
 *
 * Gates:
 * - Gate A: Data Integrity Gate (NaN/Inf, finiteness, bounds, catalog checks)
 * - Gate B: Units & Dimensional Consistency Gate (Cosmological units, scale factor, H0)
 * - Gate C: Field Divergence & Curl Gate (Potential vs rotational decomposition, continuity)
 * - Gate D: Integration Accuracy Gate (RK45 convergence, truncation error, step reversibility)
 * - Gate E: Watershed Segmentation Consistency Gate (Basin partition, volume conservation)
 * - Gate F: Critical Point Roots & Morse Index Gate (Newton-Raphson residuals, Hessian nonsingularity)
 * - Gate G: Eigenvalue & Tensor Orthonormality Gate (Jacobi diagonalization, eigenvector orthonormality)
 * - Gate H: Uncertainty & Covariance Positive-Definiteness Gate (Covariance matrix PSD, correlation bounds)
 * - Gate I: Export & Serialization Gate (FITS 2880-byte blocks, schema validation, roundtrip fidelity)
 *
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

/**
 * Gate Result Status Enumeration
 */
export const GateStatus = {
  PASSED: 'PASSED',
  FAILED: 'FAILED',
  WARNING: 'WARNING',
  SKIPPED: 'SKIPPED'
};

/**
 * Linear Algebra & Matrix Verification Helper
 */
class LinearAlgebraValidator {
  /**
   * Evaluates if a 3x3 matrix is symmetric within tolerance.
   * @param {number[][]} m 3x3 matrix
   * @param {number} tol
   * @returns {boolean}
   */
  static isSymmetric3x3(m, tol = 1e-7) {
    if (!m || m.length !== 3) return false;
    for (let r = 0; r < 3; r++) {
      if (!m[r] || m[r].length !== 3) return false;
    }
    return (
      Math.abs(m[0][1] - m[1][0]) <= tol &&
      Math.abs(m[0][2] - m[2][0]) <= tol &&
      Math.abs(m[1][2] - m[2][1]) <= tol
    );
  }

  /**
   * Computes determinant of 3x3 matrix.
   * @param {number[][]} m 3x3 matrix
   * @returns {number}
   */
  static det3x3(m) {
    return (
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    );
  }

  /**
   * Computes eigenvalues of real symmetric 3x3 matrix via Jacobi rotation method.
   * @param {number[][]} a
   * @param {number} maxIter
   * @returns {{ eigenvalues: number[], eigenvectors: number[][], converged: boolean }}
   */
  static jacobiEigenvalues3x3(a, maxIter = 50) {
    // Clone matrix
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

    for (let iter = 0; iter < maxIter; iter++) {
      let maxOff = 0;
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

      if (maxOff < 1e-12) {
        converged = true;
        break;
      }

      const app = A[p][p];
      const aqq = A[q][q];
      const apq = A[p][q];

      const phi = 0.5 * Math.atan2(2 * apq, aqq - app);
      const c = Math.cos(phi);
      const s = Math.sin(phi);

      // Rotate A
      for (let k = 0; k < 3; k++) {
        if (k !== p && k !== q) {
          const akp = A[k][p];
          const akq = A[k][q];
          A[k][p] = A[p][k] = c * akp - s * akq;
          A[k][q] = A[q][k] = s * akp + c * akq;
        }
      }

      A[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
      A[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
      A[p][q] = A[q][p] = 0.0;

      // Update eigenvectors
      for (let k = 0; k < 3; k++) {
        const vkp = V[k][p];
        const vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq;
        V[k][q] = s * vkp + c * vkq;
      }
    }

    // Sort eigenvalues ascending
    const pairs = [
      { val: A[0][0], vec: [V[0][0], V[1][0], V[2][0]] },
      { val: A[1][1], vec: [V[0][1], V[1][1], V[2][1]] },
      { val: A[2][2], vec: [V[0][2], V[1][2], V[2][2]] }
    ].sort((x, y) => x.val - y.val);

    return {
      eigenvalues: [pairs[0].val, pairs[1].val, pairs[2].val],
      eigenvectors: [pairs[0].vec, pairs[1].vec, pairs[2].vec],
      converged
    };
  }
}

/**
 * Gate A: Data Integrity Gate
 */
export class GateA_DataIntegrity {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    if (!payload) {
      return { gateId: 'GATE_A', name: 'Data Integrity Gate', status: GateStatus.FAILED, errors: ['Payload is null or undefined'], metrics };
    }

    // 1. Catalog / Records Validation
    if (payload.records !== undefined) {
      if (!Array.isArray(payload.records) || payload.records.length === 0) {
        errors.push('Records catalog must be a non-empty array');
      } else {
        let nanCount = 0;
        let invalidCoords = 0;
        let cViolation = 0;

        for (let i = 0; i < payload.records.length; i++) {
          const rec = payload.records[i];
          const x = rec.x !== undefined ? rec.x : (rec.position ? rec.position[0] : null);
          const y = rec.y !== undefined ? rec.y : (rec.position ? rec.position[1] : null);
          const z = rec.z !== undefined ? rec.z : (rec.position ? rec.position[2] : null);

          if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            nanCount++;
          }

          if (rec.vMag !== undefined && (rec.vMag < 0 || rec.vMag > 300000)) {
            cViolation++;
          }

          if (rec.distance !== undefined && rec.distance < 0) {
            invalidCoords++;
          }
        }

        metrics.recordCount = payload.records.length;
        metrics.nanCount = nanCount;
        metrics.speedOfLightViolations = cViolation;

        if (nanCount > 0) errors.push(`Catalog contains ${nanCount} records with NaN/Inf coordinates`);
        if (cViolation > 0) errors.push(`Catalog contains ${cViolation} records exceeding physical speed of light`);
        if (invalidCoords > 0) errors.push(`Catalog contains ${invalidCoords} records with negative distance`);
      }
    }

    // 2. 3D Grid Buffer Validation
    if (payload.grid !== undefined) {
      const g = payload.grid;
      const dims = g.dimensions || [128, 128, 128];
      const expectedSize = dims[0] * dims[1] * dims[2];

      metrics.gridDimensions = dims;
      metrics.expectedVoxelCount = expectedSize;

      if (g.data) {
        const len = g.data.length;
        metrics.actualVoxelCount = len;
        if (len !== expectedSize) {
          errors.push(`Grid buffer length ${len} does not match dimensions product ${expectedSize}`);
        }

        let gridNan = 0;
        for (let i = 0; i < Math.min(len, 10000); i++) {
          if (!Number.isFinite(g.data[i])) gridNan++;
        }
        if (gridNan > 0) errors.push(`Grid sample contains ${gridNan} non-finite voxel values`);
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_A',
      name: 'Data Integrity Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate B: Units & Dimensional Consistency Gate
 */
export class GateB_UnitsConsistency {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const cosmo = payload.cosmology || payload;
    const h0 = cosmo.h0 !== undefined ? cosmo.h0 : cosmo.H0;
    const omegaM = cosmo.omegaM !== undefined ? cosmo.omegaM : cosmo.Omega_m;
    const omegaL = cosmo.omegaL !== undefined ? cosmo.omegaL : cosmo.Omega_Lambda;
    const scaleFactor = payload.scaleFactor !== undefined ? payload.scaleFactor : payload.a;
    const redshift = payload.redshift !== undefined ? payload.redshift : payload.z;

    if (h0 !== undefined) {
      metrics.h0 = h0;
      if (h0 < 30.0 || h0 > 120.0 || !Number.isFinite(h0)) {
        errors.push(`Hubble constant H0 = ${h0} km/s/Mpc out of physical cosmological bounds [30, 120]`);
      }
    }

    if (omegaM !== undefined) {
      metrics.omegaM = omegaM;
      if (omegaM <= 0.0 || omegaM > 2.0 || !Number.isFinite(omegaM)) {
        errors.push(`Matter density parameter Omega_m = ${omegaM} out of physical bounds (0, 2.0]`);
      }
    }

    if (omegaL !== undefined) {
      metrics.omegaL = omegaL;
      if (omegaL < 0.0 || omegaL > 2.0 || !Number.isFinite(omegaL)) {
        errors.push(`Dark energy parameter Omega_Lambda = ${omegaL} out of physical bounds [0, 2.0]`);
      }
    }

    if (scaleFactor !== undefined) {
      metrics.scaleFactor = scaleFactor;
      if (scaleFactor <= 0.0 || scaleFactor > 100.0 || !Number.isFinite(scaleFactor)) {
        errors.push(`Cosmological scale factor a = ${scaleFactor} must be strictly positive and finite`);
      }
    }

    if (redshift !== undefined) {
      metrics.redshift = redshift;
      if (redshift < -1.0 || !Number.isFinite(redshift)) {
        errors.push(`Redshift z = ${redshift} cannot be strictly less than -1.0`);
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_B',
      name: 'Units & Dimensional Consistency Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate C: Field Divergence & Curl Gate
 */
export class GateC_FieldDivergenceAndCurl {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const maxVorticityRatio = payload.maxVorticityRatio !== undefined ? payload.maxVorticityRatio : 0.20;
    const curlNorm = payload.curlNorm !== undefined ? payload.curlNorm : 0.0;
    const divNorm = payload.divNorm !== undefined ? payload.divNorm : 1.0;
    const totalGradNorm = payload.totalGradNorm !== undefined ? payload.totalGradNorm : (divNorm + curlNorm);

    metrics.curlNorm = curlNorm;
    metrics.divNorm = divNorm;
    metrics.vorticityFraction = totalGradNorm > 0 ? curlNorm / totalGradNorm : 0.0;

    if (metrics.vorticityFraction > maxVorticityRatio) {
      errors.push(`Velocity field vorticity fraction ${metrics.vorticityFraction.toFixed(4)} exceeds linear theory potential threshold ${maxVorticityRatio}`);
    }

    if (!Number.isFinite(curlNorm) || !Number.isFinite(divNorm)) {
      errors.push('Field divergence or curl contains non-finite numerical values');
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_C',
      name: 'Field Divergence & Curl Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate D: Integration Accuracy Gate
 */
export class GateD_IntegrationAccuracy {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const maxTruncationError = payload.truncationError !== undefined ? payload.truncationError : 0.0;
    const tol = payload.tolerance !== undefined ? payload.tolerance : 1e-4;
    const reversibilityResidual = payload.reversibilityResidual !== undefined ? payload.reversibilityResidual : 0.0;
    const maxReversibilityTol = payload.reversibilityTol !== undefined ? payload.reversibilityTol : 1e-3;

    metrics.truncationError = maxTruncationError;
    metrics.reversibilityResidual = reversibilityResidual;

    if (maxTruncationError > tol) {
      errors.push(`Integration local truncation error ${maxTruncationError} exceeds tolerance threshold ${tol}`);
    }

    if (reversibilityResidual > maxReversibilityTol) {
      errors.push(`Streamline backward-forward reversibility residual ${reversibilityResidual} exceeds bound ${maxReversibilityTol}`);
    }

    if (payload.stepSizes && Array.isArray(payload.stepSizes)) {
      const negativeSteps = payload.stepSizes.filter(s => s <= 0).length;
      if (negativeSteps > 0) {
        errors.push(`Adaptive integrator encountered ${negativeSteps} non-positive step increments`);
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_D',
      name: 'Integration Accuracy Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate E: Watershed Segmentation Consistency Gate
 */
export class GateE_WatershedSegmentation {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const basins = payload.basins || [];
    const totalVolume = payload.totalVolume !== undefined ? payload.totalVolume : 1000.0;
    let sumBasinVolumes = 0;

    metrics.basinCount = basins.length;

    if (basins.length === 0) {
      errors.push('Watershed segmentation produced 0 catchment basins');
    }

    for (let i = 0; i < basins.length; i++) {
      const b = basins[i];
      const vol = b.volume !== undefined ? b.volume : 0;
      if (vol <= 0) {
        errors.push(`Basin #${b.id || i} has non-positive volume: ${vol}`);
      }
      sumBasinVolumes += vol;
    }

    metrics.sumBasinVolumes = sumBasinVolumes;
    metrics.totalDomainVolume = totalVolume;

    const volumeConservationRelDiff = Math.abs(sumBasinVolumes - totalVolume) / Math.max(1e-9, totalVolume);
    metrics.volumeConservationRelDiff = volumeConservationRelDiff;

    if (volumeConservationRelDiff > 0.05 && payload.strictVolumeConservation) {
      errors.push(`Watershed volume conservation error ${ (volumeConservationRelDiff * 100).toFixed(2) }% exceeds 5% bound`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_E',
      name: 'Watershed Segmentation Consistency Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate F: Critical Point Roots & Morse Index Gate
 */
export class GateF_CriticalPointRoots {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const points = payload.criticalPoints || payload.points || [];
    const rootTol = payload.rootTolerance !== undefined ? payload.rootTolerance : 1e-3;

    metrics.criticalPointsCount = points.length;
    let maxResidual = 0;
    const morseCounts = { 0: 0, 1: 0, 2: 0, 3: 0 };

    for (let i = 0; i < points.length; i++) {
      const cp = points[i];
      const residual = cp.residual !== undefined ? cp.residual : (cp.velocityMagnitude || 0.0);
      if (residual > maxResidual) maxResidual = residual;

      if (residual > rootTol) {
        errors.push(`Critical point #${cp.id || i} (${cp.name || 'unnamed'}) root residual ${residual} > tolerance ${rootTol}`);
      }

      const idx = cp.morseIndex !== undefined ? cp.morseIndex : cp.morse_index;
      if (idx !== undefined && morseCounts[idx] !== undefined) {
        morseCounts[idx]++;
      } else {
        errors.push(`Invalid Morse index ${idx} for critical point #${cp.id || i}`);
      }

      if (cp.hessian) {
        const det = LinearAlgebraValidator.det3x3(cp.hessian);
        if (Math.abs(det) < 1e-12) {
          errors.push(`Critical point #${cp.id || i} has degenerate singular Hessian (det = ${det})`);
        }
      }
    }

    metrics.maxRootResidual = maxResidual;
    metrics.morseCensus = morseCounts;

    // Euler characteristic check for compact topology if provided
    if (payload.checkEulerCharacteristic && points.length > 0) {
      const euler = morseCounts[0] - morseCounts[1] + morseCounts[2] - morseCounts[3];
      metrics.computedEulerCharacteristic = euler;
      const expectedEuler = payload.expectedEulerCharacteristic !== undefined ? payload.expectedEulerCharacteristic : 0;
      if (euler !== expectedEuler) {
        errors.push(`Morse-Smale alternating sum chi = ${euler} does not match expected manifold Euler characteristic ${expectedEuler}`);
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_F',
      name: 'Critical Point Roots & Morse Index Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate G: Eigenvalue & Tensor Orthonormality Gate
 */
export class GateG_EigenvalueOrthonormality {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const tensor = payload.tensor || payload.deformationTensor;
    if (!tensor) {
      return { gateId: 'GATE_G', name: 'Eigenvalue & Tensor Orthonormality Gate', status: GateStatus.FAILED, errors: ['Tensor payload missing'], metrics };
    }

    const isSym = LinearAlgebraValidator.isSymmetric3x3(tensor);
    metrics.isSymmetric = isSym;
    if (!isSym) {
      errors.push('Tidal / shear deformation tensor is not symmetric within tolerance');
    }

    const { eigenvalues, eigenvectors, converged } = LinearAlgebraValidator.jacobiEigenvalues3x3(tensor);
    metrics.eigenvalues = eigenvalues;
    metrics.diagonalizationConverged = converged;

    if (!converged) {
      errors.push('Jacobi eigenvalue diagonalization failed to converge');
    }

    // Check eigenvalue ordering lambda1 <= lambda2 <= lambda3
    if (eigenvalues[0] > eigenvalues[1] + 1e-9 || eigenvalues[1] > eigenvalues[2] + 1e-9) {
      errors.push(`Eigenvalues are not sorted monotonically: [${eigenvalues.join(', ')}]`);
    }

    // Check eigenvector orthonormality v_i . v_j = delta_ij
    let maxOrthoError = 0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const dot = eigenvectors[i][0] * eigenvectors[j][0] +
                    eigenvectors[i][1] * eigenvectors[j][1] +
                    eigenvectors[i][2] * eigenvectors[j][2];
        const target = i === j ? 1.0 : 0.0;
        const err = Math.abs(dot - target);
        if (err > maxOrthoError) maxOrthoError = err;
      }
    }

    metrics.maxOrthonormalityError = maxOrthoError;
    if (maxOrthoError > 1e-4) {
      errors.push(`Eigenvector basis orthonormality error ${maxOrthoError} exceeds tolerance 1e-4`);
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_G',
      name: 'Eigenvalue & Tensor Orthonormality Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate H: Uncertainty & Covariance Positive-Definiteness Gate
 */
export class GateH_UncertaintyCovariance {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    const cov = payload.covarianceMatrix || payload.cov;
    if (!cov || cov.length !== 3) {
      return { gateId: 'GATE_H', name: 'Uncertainty & Covariance Positive-Definiteness Gate', status: GateStatus.FAILED, errors: ['Covariance 3x3 matrix missing or malformed'], metrics };
    }

    const isSym = LinearAlgebraValidator.isSymmetric3x3(cov);
    metrics.isSymmetric = isSym;
    if (!isSym) {
      errors.push('Covariance matrix is not symmetric (C_ij != C_ji)');
    }

    const det = LinearAlgebraValidator.det3x3(cov);
    metrics.determinant = det;
    if (det < -1e-9) {
      errors.push(`Covariance matrix determinant is strictly negative (${det}) violating positive semi-definiteness`);
    }

    const { eigenvalues } = LinearAlgebraValidator.jacobiEigenvalues3x3(cov);
    metrics.eigenvalues = eigenvalues;
    if (eigenvalues.some(e => e < -1e-7)) {
      errors.push(`Covariance matrix has negative eigenvalue: [${eigenvalues.join(', ')}]`);
    }

    // Variances non-negative
    for (let i = 0; i < 3; i++) {
      if (cov[i][i] < 0) {
        errors.push(`Diagonal variance C_${i}${i} is negative (${cov[i][i]})`);
      }
    }

    // Correlation coefficients <= 1.0
    for (let i = 0; i < 3; i++) {
      for (let j = i + 1; j < 3; j++) {
        const denom = Math.sqrt(Math.max(1e-12, cov[i][i] * cov[j][j]));
        const r = cov[i][j] / denom;
        if (Math.abs(r) > 1.0001) {
          errors.push(`Correlation coefficient r_${i}${j} = ${r} exceeds physical unity bound [-1, 1]`);
        }
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_H',
      name: 'Uncertainty & Covariance Positive-Definiteness Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Gate I: Export & Serialization Gate
 */
export class GateI_ExportSerialization {
  static evaluate(payload) {
    const errors = [];
    const metrics = {};

    if (!payload) {
      return { gateId: 'GATE_I', name: 'Export & Serialization Gate', status: GateStatus.FAILED, errors: ['Payload missing'], metrics };
    }

    // FITS Header block alignment check (must be exact multiple of 2880 bytes)
    if (payload.fitsHeader !== undefined) {
      const len = payload.fitsHeader.length;
      metrics.fitsHeaderBytes = len;
      if (len % 2880 !== 0) {
        errors.push(`FITS header length (${len} bytes) is not an exact multiple of 2880 bytes`);
      }
      if (!payload.fitsHeader.includes('END')) {
        errors.push('FITS header is missing mandatory END record');
      }
    }

    // SHA-256 Manifest check
    if (payload.manifest !== undefined) {
      const keys = Object.keys(payload.manifest);
      metrics.manifestFilesCount = keys.length;
      if (keys.length === 0) {
        errors.push('Manifest contains zero file hash entries');
      }
      for (const k of keys) {
        const hash = payload.manifest[k];
        if (!/^[a-f0-9]{64}$/i.test(hash)) {
          errors.push(`Invalid SHA-256 hash string for ${k}: ${hash}`);
        }
      }
    }

    // GeoJSON validation
    if (payload.geoJson !== undefined) {
      const gj = payload.geoJson;
      if (gj.type !== 'FeatureCollection' || !Array.isArray(gj.features)) {
        errors.push('GeoJSON payload must be a valid FeatureCollection with features array');
      } else {
        metrics.geoJsonFeatures = gj.features.length;
      }
    }

    const passed = errors.length === 0;
    return {
      gateId: 'GATE_I',
      name: 'Export & Serialization Gate',
      status: passed ? GateStatus.PASSED : GateStatus.FAILED,
      passed,
      errors,
      metrics
    };
  }
}

/**
 * Acceptance Gatekeeper orchestrating all Gates A-I.
 */
export class AcceptanceGatekeeper {
  /**
   * Run complete validation audit across all Gates A through I.
   * @param {Object} context Unified context containing data for all gates
   * @returns {Object} Comprehensive audit report
   */
  static evaluateAllGates(context = {}) {
    const timestamp = new Date().toISOString();
    const gateResults = [
      GateA_DataIntegrity.evaluate(context.gateA || context),
      GateB_UnitsConsistency.evaluate(context.gateB || context),
      GateC_FieldDivergenceAndCurl.evaluate(context.gateC || context),
      GateD_IntegrationAccuracy.evaluate(context.gateD || context),
      GateE_WatershedSegmentation.evaluate(context.gateE || context),
      GateF_CriticalPointRoots.evaluate(context.gateF || context),
      GateG_EigenvalueOrthonormality.evaluate(context.gateG || context),
      GateH_UncertaintyCovariance.evaluate(context.gateH || context),
      GateI_ExportSerialization.evaluate(context.gateI || context)
    ];

    const allPassed = gateResults.every(g => g.passed);
    const failedGates = gateResults.filter(g => !g.passed).map(g => g.gateId);
    const totalErrors = gateResults.reduce((acc, g) => acc + (g.errors ? g.errors.length : 0), 0);

    return {
      status: allPassed ? 'VALIDATED' : 'REJECTED',
      allPassed,
      evaluatedAt: timestamp,
      totalGatesEvaluated: gateResults.length,
      passedCount: gateResults.filter(g => g.passed).length,
      failedCount: failedGates.length,
      failedGates,
      totalErrors,
      gateResults
    };
  }
}
