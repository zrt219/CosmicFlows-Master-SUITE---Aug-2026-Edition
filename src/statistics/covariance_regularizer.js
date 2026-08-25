/**
 * @file covariance_regularizer.js
 * @description 3D Covariance Regularizer, Shrinkage Estimator, and Matrix Conditioner.
 *
 * Implements:
 * 1. Ledoit-Wolf optimal linear shrinkage covariance estimator: \Sigma^* = (1 - \lambda) S + \lambda \mu I
 *    with support for Constant Variance (Target A), Diagonal (Target B), Identity (Target C), and Constant Correlation (Target D).
 * 2. Oracle Approximating Shrinkage (OAS) and Rao-Blackwell Ledoit-Wolf (RBLW) estimators
 *    optimized for small sample sizes and high-dimensional cosmological velocity dispersion cubes.
 * 3. Matrix Conditioner, Eigendecomposition (Jacobi), Condition Number \kappa(C) = \lambda_{\max}/\lambda_{\min},
 *    and Tikhonov Regularizer C_{\epsilon} = C + \epsilon I (adaptive, relative, and GCV-optimized).
 * 4. Cholesky Factorization with guaranteed positive-definiteness (Modified Gill-Murray-Wright / Schnabel-Eskow),
 *    root-free LDL^T decomposition, matrix inversion, log-determinant, and Mahalanobis distance calculator.
 * 5. Spatial Covariance Tapering with Wendland (C^0, C^2, C^4, C^6), Spherical, and Gaussian compact support kernels,
 *    sparse representation, and Schur product theorem compliance.
 * 6. 3D Cosmological Velocity Dispersion Field & Cube Regularizer with triaxiality, anisotropy parameter \beta,
 *    and voxel-wise regularized dispersion tensors.
 *
 * @module statistics/covariance_regularizer
 */

// ============================================================================
// 1. FUNDAMENTAL MATRIX AND VECTOR UTILITIES (Float64Array based)
// ============================================================================

/**
 * Creates an N x M zero-initialized matrix stored in a 1D Float64Array (row-major).
 * @param {number} rows
 * @param {number} cols
 * @returns {Float64Array}
 */
export function createMatrix(rows, cols) {
  if (rows <= 0 || cols <= 0 || !Number.isInteger(rows) || !Number.isInteger(cols)) {
    throw new RangeError(`Matrix dimensions must be positive integers: received ${rows}x${cols}`);
  }
  return new Float64Array(rows * cols);
}

/**
 * Creates an N x N identity matrix.
 * @param {number} n
 * @returns {Float64Array}
 */
export function eye(n) {
  if (n <= 0 || !Number.isInteger(n)) {
    throw new RangeError(`Dimension n must be a positive integer: received ${n}`);
  }
  const out = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    out[i * n + i] = 1.0;
  }
  return out;
}

/**
 * Clones a matrix or vector.
 * @param {ArrayLike<number>} mat
 * @returns {Float64Array}
 */
export function cloneMatrix(mat) {
  return new Float64Array(mat);
}

/**
 * Sets matrix to 2D array or flat array.
 * @param {Array<Array<number>>|Array<number>|Float64Array} data
 * @param {number} [rows]
 * @param {number} [cols]
 * @returns {Float64Array}
 */
export function matFrom(data, rows, cols) {
  if (Array.isArray(data) && Array.isArray(data[0])) {
    const r = data.length;
    const c = data[0].length;
    const out = new Float64Array(r * c);
    for (let i = 0; i < r; i++) {
      for (let j = 0; j < c; j++) {
        out[i * c + j] = Number(data[i][j]);
      }
    }
    return out;
  }
  if (rows && cols) {
    if (data.length !== rows * cols) {
      throw new Error(`Data length (${data.length}) does not match specified dimensions (${rows}x${cols})`);
    }
    return new Float64Array(data);
  }
  return new Float64Array(data);
}

/**
 * Converts a flat matrix to a 2D JavaScript Array.
 * @param {Float64Array|Array<number>} mat
 * @param {number} rows
 * @param {number} cols
 * @returns {number[][]}
 */
export function matTo2D(mat, rows, cols) {
  const result = new Array(rows);
  for (let i = 0; i < rows; i++) {
    const row = new Array(cols);
    const offset = i * cols;
    for (let j = 0; j < cols; j++) {
      row[j] = mat[offset + j];
    }
    result[i] = row;
  }
  return result;
}

/**
 * Matrix multiplication C = A * B (rowsA x colsA by colsA x colsB).
 * @param {Float64Array|Array<number>} A
 * @param {Float64Array|Array<number>} B
 * @param {number} rowsA
 * @param {number} colsA
 * @param {number} colsB
 * @returns {Float64Array}
 */
export function matMul(A, B, rowsA, colsA, colsB) {
  const C = new Float64Array(rowsA * colsB);
  for (let i = 0; i < rowsA; i++) {
    const iOffsetA = i * colsA;
    const iOffsetC = i * colsB;
    for (let k = 0; k < colsA; k++) {
      const aVal = A[iOffsetA + k];
      if (aVal === 0.0) continue;
      const kOffsetB = k * colsB;
      for (let j = 0; j < colsB; j++) {
        C[iOffsetC + j] += aVal * B[kOffsetB + j];
      }
    }
  }
  return C;
}

/**
 * Matrix-vector multiplication y = A * x (rows x cols by cols).
 * @param {Float64Array|Array<number>} A
 * @param {Float64Array|Array<number>} x
 * @param {number} rows
 * @param {number} cols
 * @returns {Float64Array}
 */
export function matVecMul(A, x, rows, cols) {
  const y = new Float64Array(rows);
  for (let i = 0; i < rows; i++) {
    let sum = 0.0;
    const offset = i * cols;
    for (let j = 0; j < cols; j++) {
      sum += A[offset + j] * x[j];
    }
    y[i] = sum;
  }
  return y;
}

/**
 * Matrix transpose A^T.
 * @param {Float64Array|Array<number>} A
 * @param {number} rows
 * @param {number} cols
 * @returns {Float64Array}
 */
export function matTranspose(A, rows, cols) {
  const AT = new Float64Array(cols * rows);
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      AT[j * rows + i] = A[i * cols + j];
    }
  }
  return AT;
}

/**
 * Trace of a square matrix Tr(A).
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @returns {number}
 */
export function matTrace(A, n) {
  let tr = 0.0;
  for (let i = 0; i < n; i++) {
    tr += A[i * n + i];
  }
  return tr;
}

/**
 * Frobenius norm of a matrix ||A||_F = sqrt(sum_ij A_ij^2).
 * @param {Float64Array|Array<number>} A
 * @returns {number}
 */
export function matFrobeniusNorm(A) {
  let sumSq = 0.0;
  const len = A.length;
  for (let i = 0; i < len; i++) {
    const val = A[i];
    sumSq += val * val;
  }
  return Math.sqrt(sumSq);
}

/**
 * Squared Frobenius norm of a matrix ||A||_F^2 = sum_ij A_ij^2.
 * @param {Float64Array|Array<number>} A
 * @returns {number}
 */
export function matFrobeniusNormSq(A) {
  let sumSq = 0.0;
  const len = A.length;
  for (let i = 0; i < len; i++) {
    const val = A[i];
    sumSq += val * val;
  }
  return sumSq;
}

/**
 * Computes Tr(A^2) = Tr(A * A) for a square symmetric matrix A.
 * Note: For symmetric A, Tr(A^2) = ||A||_F^2.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @returns {number}
 */
export function matTraceSquare(A, n) {
  let tr = 0.0;
  for (let i = 0; i < n; i++) {
    const iOffset = i * n;
    for (let j = 0; j < n; j++) {
      tr += A[iOffset + j] * A[j * n + i];
    }
  }
  return tr;
}

/**
 * Schur (Hadamard / Element-wise) product C = A .* B.
 * @param {Float64Array|Array<number>} A
 * @param {Float64Array|Array<number>} B
 * @returns {Float64Array}
 */
export function matHadamard(A, B) {
  if (A.length !== B.length) {
    throw new Error(`Dimension mismatch: A (${A.length}) != B (${B.length})`);
  }
  const len = A.length;
  const C = new Float64Array(len);
  for (let i = 0; i < len; i++) {
    C[i] = A[i] * B[i];
  }
  return C;
}

/**
 * Matrix addition with scaling: C = alpha * A + beta * B.
 * @param {Float64Array|Array<number>} A
 * @param {Float64Array|Array<number>} B
 * @param {number} [alpha=1.0]
 * @param {number} [beta=1.0]
 * @returns {Float64Array}
 */
export function matAdd(A, B, alpha = 1.0, beta = 1.0) {
  if (A.length !== B.length) {
    throw new Error(`Dimension mismatch: A (${A.length}) != B (${B.length})`);
  }
  const len = A.length;
  const C = new Float64Array(len);
  for (let i = 0; i < len; i++) {
    C[i] = alpha * A[i] + beta * B[i];
  }
  return C;
}

/**
 * Symmetrizes a square matrix C = 0.5 * (A + A^T).
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @returns {Float64Array}
 */
export function matSymmetrize(A, n) {
  const C = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    C[i * n + i] = A[i * n + i];
    for (let j = i + 1; j < n; j++) {
      const avg = 0.5 * (A[i * n + j] + A[j * n + i]);
      C[i * n + j] = avg;
      C[j * n + i] = avg;
    }
  }
  return C;
}

/**
 * Checks if a matrix is symmetric within tolerance.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @param {number} [tol=1e-10]
 * @returns {boolean}
 */
export function isSymmetric(A, n, tol = 1e-10) {
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (Math.abs(A[i * n + j] - A[j * n + i]) > tol) {
        return false;
      }
    }
  }
  return true;
}

// ============================================================================
// 2. EIGENVALUE DECOMPOSITION & MATRIX FUNCTIONS
// ============================================================================

/**
 * Result of symmetric eigendecomposition.
 * @typedef {Object} EigenResult
 * @property {Float64Array} values - Eigenvalues in descending order (\lambda_1 \ge \lambda_2 \ge \dots)
 * @property {Float64Array} vectors - Column eigenvectors V, such that A = V * diag(values) * V^T
 * @property {number} iterations - Number of Jacobi sweeps
 */

/**
 * Computes exact eigendecomposition of a symmetric N x N real matrix using Cyclic Jacobi Algorithm with thresholding.
 *
 * Guaranteed convergence and orthogonal eigenvectors for any symmetric matrix.
 *
 * @param {Float64Array|Array<number>} A - Symmetric N x N matrix (will not be modified)
 * @param {number} n - Dimension of the matrix
 * @param {number} [maxSweeps=100] - Maximum sweeps
 * @param {number} [tol=1e-15] - Convergence tolerance
 * @returns {EigenResult}
 */
export function jacobiEigenvalues(A, n, maxSweeps = 100, tol = 1e-15) {
  if (n <= 0) throw new RangeError(`Dimension n must be > 0: ${n}`);
  
  // Copy A to working matrix
  const a = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) a[i] = A[i];

  // Initialize V = I
  const v = eye(n);
  const d = new Float64Array(n);
  for (let i = 0; i < n; i++) d[i] = a[i * n + i];

  let sweep = 0;
  while (sweep < maxSweeps) {
    // Sum of off-diagonal magnitudes
    let offDiagSum = 0.0;
    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        offDiagSum += Math.abs(a[i * n + j]);
      }
    }

    if (offDiagSum <= tol) break;

    // Threshold for first sweeps to accelerate convergence
    const threshold = sweep < 3 ? 0.2 * offDiagSum / (n * n) : 0.0;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = a[p * n + q];
        const g = 100.0 * Math.abs(apq);

        // Skip small off-diagonal elements
        if (sweep > 3 && (Math.abs(d[p]) + g === Math.abs(d[p])) && (Math.abs(d[q]) + g === Math.abs(d[q]))) {
          a[p * n + q] = 0.0;
          a[q * n + p] = 0.0;
          continue;
        }

        if (Math.abs(apq) <= threshold) continue;

        const h = d[q] - d[p];
        let t;
        if (Math.abs(h) + g === Math.abs(h)) {
          t = apq / h;
        } else {
          const theta = 0.5 * h / apq;
          t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta * theta));
          if (theta < 0.0) t = -t;
        }

        const c = 1.0 / Math.sqrt(1.0 + t * t);
        const s = t * c;
        const tau = s / (1.0 + c);
        const delta = t * apq;

        d[p] -= delta;
        d[q] += delta;
        a[p * n + q] = 0.0;
        a[q * n + p] = 0.0;

        for (let j = 0; j < p; j++) {
          const ajp = a[j * n + p];
          const ajq = a[j * n + q];
          a[j * n + p] = ajp - s * (ajq + ajp * tau);
          a[p * n + j] = a[j * n + p];
          a[j * n + q] = ajq + s * (ajp - ajq * tau);
          a[q * n + j] = a[j * n + q];
        }

        for (let j = p + 1; j < q; j++) {
          const apj = a[p * n + j];
          const ajq = a[j * n + q];
          a[p * n + j] = apj - s * (ajq + apj * tau);
          a[j * n + p] = a[p * n + j];
          a[j * n + q] = ajq + s * (apj - ajq * tau);
          a[q * n + j] = a[j * n + q];
        }

        for (let j = q + 1; j < n; j++) {
          const apj = a[p * n + j];
          const aqj = a[q * n + j];
          a[p * n + j] = apj - s * (aqj + apj * tau);
          a[j * n + p] = a[p * n + j];
          a[q * n + j] = aqj + s * (apj - aqj * tau);
          a[j * n + q] = a[q * n + j];
        }

        for (let j = 0; j < n; j++) {
          const vjp = v[j * n + p];
          const vjq = v[j * n + q];
          v[j * n + p] = vjp - s * (vjq + vjp * tau);
          v[j * n + q] = vjq + s * (vjp - vjq * tau);
        }
      }
    }
    sweep++;
  }

  // Sort eigenvalues and corresponding eigenvectors in descending order
  const indices = new Int32Array(n);
  for (let i = 0; i < n; i++) indices[i] = i;
  indices.sort((i, j) => d[j] - d[i]);

  const sortedValues = new Float64Array(n);
  const sortedVectors = new Float64Array(n * n);

  for (let i = 0; i < n; i++) {
    const idx = indices[i];
    sortedValues[i] = d[idx];
    for (let r = 0; r < n; r++) {
      sortedVectors[r * n + i] = v[r * n + idx];
    }
  }

  return {
    values: sortedValues,
    vectors: sortedVectors,
    iterations: sweep
  };
}

/**
 * Computes analytical eigenvalues and eigenvectors for 3x3 symmetric matrix (Velocity dispersion tensor).
 * @param {Float64Array|Array<number>} A - 3x3 symmetric matrix
 * @returns {EigenResult}
 */
export function eigen3x3Symmetric(A) {
  return jacobiEigenvalues(A, 3, 50, 1e-15);
}

/**
 * Projects a symmetric matrix onto the cone of Positive Semi-Definite (PSD) matrices (Higham 1988).
 * Any negative eigenvalue is clamped to minEig (default 0 or small positive epsilon).
 *
 * @param {Float64Array|Array<number>} A - Symmetric matrix
 * @param {number} n - Dimension
 * @param {number} [minEig=1e-12] - Minimum allowed eigenvalue
 * @returns {Float64Array}
 */
export function nearestPositiveSemiDefinite(A, n, minEig = 1e-12) {
  const eig = jacobiEigenvalues(A, n);
  const clampedVals = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    clampedVals[i] = Math.max(minEig, eig.values[i]);
  }

  const out = new Float64Array(n * n);
  const V = eig.vectors;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += V[i * n + k] * clampedVals[k] * V[j * n + k];
      }
      out[i * n + j] = sum;
    }
  }
  return matSymmetrize(out, n);
}

/**
 * Computes matrix power A^p for a symmetric positive semi-definite matrix.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @param {number} power
 * @returns {Float64Array}
 */
export function matPowerSymmetric(A, n, power) {
  const eig = jacobiEigenvalues(A, n);
  const poweredVals = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (eig.values[i] < 0 && power % 1 !== 0) {
      throw new Error(`Cannot compute fractional power ${power} of negative eigenvalue ${eig.values[i]}`);
    }
    poweredVals[i] = Math.pow(Math.max(0, eig.values[i]), power);
  }

  const out = new Float64Array(n * n);
  const V = eig.vectors;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += V[i * n + k] * poweredVals[k] * V[j * n + k];
      }
      out[i * n + j] = sum;
    }
  }
  return matSymmetrize(out, n);
}

/**
 * Matrix square root A^{1/2} for symmetric positive semi-definite matrix.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @returns {Float64Array}
 */
export function matSqrtSymmetric(A, n) {
  return matPowerSymmetric(A, n, 0.5);
}

/**
 * Matrix inverse square root A^{-1/2} for symmetric positive definite matrix.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @param {number} [eps=1e-12]
 * @returns {Float64Array}
 */
export function matInvSqrtSymmetric(A, n, eps = 1e-12) {
  const eig = jacobiEigenvalues(A, n);
  const invSqrtVals = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const val = Math.max(eps, eig.values[i]);
    invSqrtVals[i] = 1.0 / Math.sqrt(val);
  }

  const out = new Float64Array(n * n);
  const V = eig.vectors;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += V[i * n + k] * invSqrtVals[k] * V[j * n + k];
      }
      out[i * n + j] = sum;
    }
  }
  return matSymmetrize(out, n);
}

/**
 * Matrix logarithm ln(A) for symmetric positive definite matrix.
 * @param {Float64Array|Array<number>} A
 * @param {number} n
 * @param {number} [eps=1e-15]
 * @returns {Float64Array}
 */
export function matLogSymmetric(A, n, eps = 1e-15) {
  const eig = jacobiEigenvalues(A, n);
  const logVals = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const val = Math.max(eps, eig.values[i]);
    logVals[i] = Math.log(val);
  }

  const out = new Float64Array(n * n);
  const V = eig.vectors;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += V[i * n + k] * logVals[k] * V[j * n + k];
      }
      out[i * n + j] = sum;
    }
  }
  return matSymmetrize(out, n);
}

// ============================================================================
// 3. MATRIX CONDITIONER AND TIKHONOV REGULARIZATION
// ============================================================================

/**
 * Complete Condition Number Analysis Result.
 * @typedef {Object} ConditionReport
 * @property {number} kappa2 - 2-norm condition number \kappa_2(C) = \lambda_{\max} / \lambda_{\min}
 * @property {number} rcond - Reciprocal condition number 1 / \kappa_2(C)
 * @property {number} lambdaMax - Largest eigenvalue
 * @property {number} lambdaMin - Smallest eigenvalue
 * @property {boolean} isPositiveDefinite - Whether \lambda_{\min} > 0
 * @property {boolean} isIllConditioned - Whether \kappa_2 > 1e12 or \lambda_{\min} <= 0
 * @property {string} classification - 'well-conditioned', 'moderately-conditioned', 'ill-conditioned', 'singular', or 'indefinite'
 */

/**
 * Evaluates the condition number and spectral health of a symmetric matrix.
 * @param {Float64Array|Array<number>} C - Symmetric matrix
 * @param {number} n - Dimension
 * @returns {ConditionReport}
 */
export function computeConditionReport(C, n) {
  const eig = jacobiEigenvalues(C, n);
  const lambdaMax = eig.values[0];
  const lambdaMin = eig.values[n - 1];

  let kappa2;
  let rcond;
  let isPositiveDefinite = lambdaMin > 0;
  let classification;

  if (lambdaMin <= 0) {
    kappa2 = Infinity;
    rcond = 0.0;
    classification = lambdaMin < -1e-12 ? 'indefinite' : 'singular';
  } else {
    kappa2 = lambdaMax / lambdaMin;
    rcond = lambdaMin / lambdaMax;
    if (kappa2 < 1e3) {
      classification = 'well-conditioned';
    } else if (kappa2 < 1e7) {
      classification = 'moderately-conditioned';
    } else {
      classification = 'ill-conditioned';
    }
  }

  const isIllConditioned = kappa2 > 1e12 || !isPositiveDefinite;

  return {
    kappa2,
    rcond,
    lambdaMax,
    lambdaMin,
    isPositiveDefinite,
    isIllConditioned,
    classification
  };
}

/**
 * Computes standard 2-norm condition number \kappa_2(C) = \lambda_{\max} / \lambda_{\min}.
 * @param {Float64Array|Array<number>} C
 * @param {number} n
 * @returns {number}
 */
export function conditionNumber(C, n) {
  return computeConditionReport(C, n).kappa2;
}

/**
 * Tikhonov Matrix Regularization: C_\epsilon = C + \epsilon I.
 *
 * @param {Float64Array|Array<number>} C - Square matrix
 * @param {number} n - Dimension
 * @param {number} epsilon - Ridge / Tikhonov parameter
 * @returns {Float64Array}
 */
export function tikhonovRegularize(C, n, epsilon) {
  if (epsilon < 0) {
    throw new RangeError(`Tikhonov parameter epsilon must be >= 0: received ${epsilon}`);
  }
  const reg = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) {
    reg[i] = C[i];
  }
  for (let i = 0; i < n; i++) {
    reg[i * n + i] += epsilon;
  }
  return reg;
}

/**
 * Relative Tikhonov Regularization: C_{\epsilon} = C + (\alpha \cdot \frac{\text{Tr}(C)}{n}) I.
 *
 * @param {Float64Array|Array<number>} C
 * @param {number} n
 * @param {number} [alpha=1e-4]
 * @returns {Float64Array}
 */
export function relativeTikhonovRegularize(C, n, alpha = 1e-4) {
  const tr = matTrace(C, n);
  const avgVar = Math.max(1e-12, tr / n);
  const eps = alpha * avgVar;
  return tikhonovRegularize(C, n, eps);
}

/**
 * Adaptive Condition-Capped Regularizer.
 * Dynamically computes minimal \epsilon \ge 0 such that \kappa_2(C + \epsilon I) \le \kappa_{\max}.
 *
 * Formula: If \kappa(C) > \kappa_{\max}, \epsilon = \frac{\lambda_{\max} - \kappa_{\max} \lambda_{\min}}{\kappa_{\max} - 1}.
 *
 * @param {Float64Array|Array<number>} C - Symmetric matrix
 * @param {number} n - Dimension
 * @param {number} [targetKappa=1e4] - Target maximum condition number
 * @returns {{ regularized: Float64Array, epsilon: number, originalKappa: number, finalKappa: number }}
 */
export function adaptiveConditionCapping(C, n, targetKappa = 1e4) {
  if (targetKappa <= 1.0) {
    throw new RangeError(`targetKappa must be > 1.0: received ${targetKappa}`);
  }
  const report = computeConditionReport(C, n);
  const lambdaMax = report.lambdaMax;
  const lambdaMin = report.lambdaMin;

  let epsilon = 0.0;
  if (report.kappa2 > targetKappa || lambdaMin <= 0) {
    const safeMin = Math.max(0.0, lambdaMin);
    epsilon = (lambdaMax - targetKappa * safeMin) / (targetKappa - 1.0);
    if (epsilon <= 0 || !Number.isFinite(epsilon)) {
      epsilon = 1e-6 * Math.max(1.0, Math.abs(lambdaMax));
    }
  }

  const regularized = tikhonovRegularize(C, n, epsilon);
  const finalReport = computeConditionReport(regularized, n);

  return {
    regularized,
    epsilon,
    originalKappa: report.kappa2,
    finalKappa: finalReport.kappa2
  };
}

/**
 * Generalized Cross-Validation (GCV) optimal Tikhonov parameter selector for linear inversion.
 *
 * @param {Float64Array|Array<number>} A - Forward operator matrix (M x N)
 * @param {Float64Array|Array<number>} y - Observation vector (length M)
 * @param {number} m - Observations
 * @param {number} n - State dimension
 * @param {number[]} [gridEps] - Candidate epsilon values to evaluate
 * @returns {{ bestEpsilon: number, bestGcvScore: number, evaluations: Array<{epsilon: number, gcv: number}> }}
 */
export function optimalTikhonovGCV(A, y, m, n, gridEps) {
  const candidates = gridEps || [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0];
  const AT = matTranspose(A, m, n);
  const ATA = matMul(AT, A, n, m, n);
  const ATy = matVecMul(AT, y, n, m);

  let bestEps = candidates[0];
  let bestScore = Infinity;
  const evaluations = [];

  for (const eps of candidates) {
    const regATA = tikhonovRegularize(ATA, n, eps);
    const ch = choleskyFactorization(regATA, n);
    if (!ch.isPositiveDefinite) continue;

    const x = choleskySolve(ch.L, ATy, n);
    const Ax = matVecMul(A, x, m, n);
    let resNormSq = 0.0;
    for (let i = 0; i < m; i++) {
      const diff = y[i] - Ax[i];
      resNormSq += diff * diff;
    }

    const eig = jacobiEigenvalues(ATA, n);
    let dof = 0.0;
    for (let i = 0; i < n; i++) {
      const l = Math.max(0, eig.values[i]);
      dof += l / (l + eps);
    }

    const trDenom = m - dof;
    const gcvScore = trDenom > 0 ? (resNormSq / (trDenom * trDenom)) : Infinity;

    evaluations.push({ epsilon: eps, gcv: gcvScore });
    if (gcvScore < bestScore) {
      bestScore = gcvScore;
      bestEps = eps;
    }
  }

  return {
    bestEpsilon: bestEps,
    bestGcvScore: bestScore,
    evaluations
  };
}

// ============================================================================
// 4. CHOLESKY FACTORIZATION, LDL^T DECOMPOSITION, AND SOLVERS
// ============================================================================

/**
 * Result of Cholesky Factorization.
 * @typedef {Object} CholeskyResult
 * @property {Float64Array} L - Lower triangular matrix L such that C = L * L^T
 * @property {boolean} isPositiveDefinite - Whether original matrix was strictly positive definite
 * @property {Float64Array} diagonalAdded - Perturbation \Delta D added to guarantee positive definiteness
 * @property {number} maxPerturbation - Maximum diagonal perturbation added
 */

/**
 * Standard Cholesky Factorization: C = L * L^T.
 *
 * @param {Float64Array|Array<number>} C - Symmetric N x N matrix
 * @param {number} n - Dimension
 * @param {number} [tol=1e-12] - Small positive tolerance
 * @returns {CholeskyResult}
 */
export function choleskyFactorization(C, n, tol = 1e-12) {
  const L = new Float64Array(n * n);
  let isPositiveDefinite = true;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0.0;
      const iOffset = i * n;
      const jOffset = j * n;
      for (let k = 0; k < j; k++) {
        sum += L[iOffset + k] * L[jOffset + k];
      }

      if (i === j) {
        const val = C[iOffset + i] - sum;
        if (val <= tol) {
          isPositiveDefinite = false;
          L[iOffset + i] = 0.0;
        } else {
          L[iOffset + i] = Math.sqrt(val);
        }
      } else {
        const diagJ = L[jOffset + j];
        if (diagJ <= tol) {
          L[iOffset + j] = 0.0;
        } else {
          L[iOffset + j] = (C[iOffset + j] - sum) / diagJ;
        }
      }
    }
  }

  return {
    L,
    isPositiveDefinite,
    diagonalAdded: new Float64Array(n),
    maxPerturbation: 0.0
  };
}

/**
 * Modified Gill-Murray-Wright / Schnabel-Eskow Cholesky Factorization.
 *
 * Guarantees that (C + E) = L * L^T where E is a non-negative diagonal matrix
 * with minimal ||E||_\infty, ensuring positive definiteness for any symmetric matrix.
 *
 * @param {Float64Array|Array<number>} C - Symmetric N x N matrix
 * @param {number} n - Dimension
 * @param {number} [delta=1e-6] - Minimum eigenvalue guarantee parameter
 * @returns {CholeskyResult}
 */
export function modifiedCholesky(C, n, delta = 1e-6) {
  let gamma = 0.0;
  let xi = 0.0;
  for (let i = 0; i < n; i++) {
    gamma = Math.max(gamma, Math.abs(C[i * n + i]));
    for (let j = 0; j < i; j++) {
      xi = Math.max(xi, Math.abs(C[i * n + j]));
    }
  }

  const nu = Math.max(1.0, Math.sqrt(n * n - 1));
  const beta2 = Math.max(gamma, xi / nu, 1e-12);
  const beta = Math.sqrt(beta2);

  const L = new Float64Array(n * n);
  const D = new Float64Array(n);
  const diagAdded = new Float64Array(n);
  let maxPerturb = 0.0;

  const A = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) A[i] = C[i];

  for (let j = 0; j < n; j++) {
    let thetaJ = 0.0;
    for (let i = j + 1; i < n; i++) {
      thetaJ = Math.max(thetaJ, Math.abs(A[i * n + j]));
    }

    const thetaJ_beta = thetaJ / beta;
    const requiredDiag = Math.max(Math.abs(A[j * n + j]), (thetaJ_beta * thetaJ_beta), delta);
    const added = requiredDiag - A[j * n + j];
    const e_j = added > 0 ? added : 0.0;

    diagAdded[j] = e_j;
    maxPerturb = Math.max(maxPerturb, e_j);
    D[j] = A[j * n + j] + e_j;

    L[j * n + j] = Math.sqrt(D[j]);

    for (let i = j + 1; i < n; i++) {
      L[i * n + j] = A[i * n + j] / L[j * n + j];
      for (let k = j + 1; k <= i; k++) {
        A[i * n + k] -= L[i * n + j] * L[k * n + j];
        A[k * n + i] = A[i * n + k];
      }
    }
  }

  return {
    L,
    isPositiveDefinite: maxPerturb === 0.0,
    diagonalAdded: diagAdded,
    maxPerturbation: maxPerturb
  };
}

/**
 * Root-Free LDL^T Decomposition: C = L * D * L^T.
 * L is unit lower triangular (diagonal = 1), D is diagonal.
 *
 * @param {Float64Array|Array<number>} C - Symmetric matrix
 * @param {number} n - Dimension
 * @returns {{ L: Float64Array, D: Float64Array, isPositiveDefinite: boolean }}
 */
export function ldltDecomposition(C, n) {
  const L = eye(n);
  const D = new Float64Array(n);
  let isPositiveDefinite = true;

  for (let j = 0; j < n; j++) {
    let sumD = 0.0;
    for (let k = 0; k < j; k++) {
      sumD += L[j * n + k] * L[j * n + k] * D[k];
    }
    D[j] = C[j * n + j] - sumD;
    if (D[j] <= 0) {
      isPositiveDefinite = false;
    }

    const dJ = D[j];
    for (let i = j + 1; i < n; i++) {
      let sumL = 0.0;
      for (let k = 0; k < j; k++) {
        sumL += L[i * n + k] * L[j * n + k] * D[k];
      }
      L[i * n + j] = Math.abs(dJ) > 1e-15 ? (C[i * n + j] - sumL) / dJ : 0.0;
    }
  }

  return { L, D, isPositiveDefinite };
}

/**
 * Forward substitution: Solves L * y = b for lower triangular matrix L.
 * @param {Float64Array|Array<number>} L - Lower triangular matrix (n x n)
 * @param {Float64Array|Array<number>} b - RHS vector (length n)
 * @param {number} n - Dimension
 * @returns {Float64Array}
 */
export function forwardSolve(L, b, n) {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let sum = 0.0;
    const offset = i * n;
    for (let j = 0; j < i; j++) {
      sum += L[offset + j] * y[j];
    }
    const diag = L[offset + i];
    if (Math.abs(diag) < 1e-15) {
      throw new Error(`Zero pivot encountered on diagonal at index ${i} in forwardSolve`);
    }
    y[i] = (b[i] - sum) / diag;
  }
  return y;
}

/**
 * Backward substitution: Solves L^T * x = y for upper triangular matrix L^T.
 * @param {Float64Array|Array<number>} L - Lower triangular matrix (n x n)
 * @param {Float64Array|Array<number>} y - RHS vector (length n)
 * @param {number} n - Dimension
 * @returns {Float64Array}
 */
export function backwardSolve(L, y, n) {
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let sum = 0.0;
    for (let j = i + 1; j < n; j++) {
      sum += L[j * n + i] * x[j];
    }
    const diag = L[i * n + i];
    if (Math.abs(diag) < 1e-15) {
      throw new Error(`Zero pivot encountered on diagonal at index ${i} in backwardSolve`);
    }
    x[i] = (y[i] - sum) / diag;
  }
  return x;
}

/**
 * Solves C * x = b using Cholesky factor L (C = L * L^T).
 * @param {Float64Array|Array<number>} L - Lower triangular factor
 * @param {Float64Array|Array<number>} b - RHS vector
 * @param {number} n - Dimension
 * @returns {Float64Array}
 */
export function choleskySolve(L, b, n) {
  const y = forwardSolve(L, b, n);
  return backwardSolve(L, y, n);
}

/**
 * Computes exact matrix inverse C^{-1} via Cholesky decomposition.
 * @param {Float64Array|Array<number>} C - Symmetric positive definite matrix
 * @param {number} n - Dimension
 * @returns {Float64Array}
 */
export function choleskyInvert(C, n) {
  const ch = choleskyFactorization(C, n);
  let L = ch.L;
  if (!ch.isPositiveDefinite) {
    const mod = modifiedCholesky(C, n);
    L = mod.L;
  }

  const Linv = new Float64Array(n * n);
  const col = new Float64Array(n);

  for (let j = 0; j < n; j++) {
    col.fill(0.0);
    col[j] = 1.0;
    const y = forwardSolve(L, col, n);
    for (let i = 0; i < n; i++) {
      Linv[i * n + j] = y[i];
    }
  }

  const Cinv = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += Linv[k * n + i] * Linv[k * n + j];
      }
      Cinv[i * n + j] = sum;
      Cinv[j * n + i] = sum;
    }
  }

  return Cinv;
}

/**
 * Computes numerically stable log-determinant \ln \det(C) = 2 \sum_{i=1}^n \ln(L_{ii}).
 * @param {Float64Array|Array<number>} C - Symmetric positive definite matrix
 * @param {number} n - Dimension
 * @returns {number}
 */
export function logDeterminant(C, n) {
  const ch = choleskyFactorization(C, n);
  if (!ch.isPositiveDefinite) {
    const eig = jacobiEigenvalues(C, n);
    let logDet = 0.0;
    for (let i = 0; i < n; i++) {
      if (eig.values[i] <= 0) return -Infinity;
      logDet += Math.log(eig.values[i]);
    }
    return logDet;
  }

  let logDet = 0.0;
  for (let i = 0; i < n; i++) {
    logDet += 2.0 * Math.log(ch.L[i * n + i]);
  }
  return logDet;
}

/**
 * Computes exact determinant \det(C) = \prod L_{ii}^2.
 * @param {Float64Array|Array<number>} C - Symmetric matrix
 * @param {number} n - Dimension
 * @returns {number}
 */
export function determinant(C, n) {
  return Math.exp(logDeterminant(C, n));
}

/**
 * Computes Mahalanobis distance d_M(x, \mu) = \sqrt{(x - \mu)^T C^{-1} (x - \mu)}.
 *
 * @param {Float64Array|Array<number>} x - Vector 1
 * @param {Float64Array|Array<number>} mu - Center vector (or Vector 2)
 * @param {Float64Array|Array<number>} C - Covariance matrix
 * @param {number} n - Dimension
 * @returns {number}
 */
export function mahalanobisDistance(x, mu, C, n) {
  const diff = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    diff[i] = x[i] - mu[i];
  }

  const ch = choleskyFactorization(C, n);
  let L = ch.L;
  if (!ch.isPositiveDefinite) {
    L = modifiedCholesky(C, n).L;
  }

  const z = forwardSolve(L, diff, n);

  let sumSq = 0.0;
  for (let i = 0; i < n; i++) {
    sumSq += z[i] * z[i];
  }
  return Math.sqrt(sumSq);
}

// ============================================================================
// 5. LEDOIT-WOLF OPTIMAL LINEAR SHRINKAGE ESTIMATOR
// ============================================================================

/**
 * Target matrix types for shrinkage estimation.
 * @readonly
 * @enum {string}
 */
export const ShrinkageTarget = {
  SCALED_IDENTITY: 'scaled_identity', // F = \mu I (Target A, Ledoit-Wolf 2004)
  DIAGONAL: 'diagonal',               // F = \text{diag}(S) (Target B, unequal variances)
  IDENTITY: 'identity',               // F = I (Target C, isotropic unit variance)
  CONSTANT_CORRELATION: 'constant_correlation' // Target D, equal off-diagonal correlation
};

/**
 * Complete Ledoit-Wolf Shrinkage Result.
 * @typedef {Object} LedoitWolfResult
 * @property {Float64Array} regularizedCovariance - Shrunk covariance \Sigma^* = (1-\lambda) S + \lambda F
 * @property {Float64Array} sampleCovariance - Raw empirical sample covariance S
 * @property {Float64Array} targetMatrix - Target matrix F
 * @property {number} shrinkageIntensity - Optimal shrinkage coefficient \lambda^* \in [0, 1]
 * @property {number} deltaSq - ||S - F||_F^2
 * @property {number} bSq - Asymptotic estimation variance b^2
 * @property {number} mu - Mean variance parameter \mu = Tr(S)/p
 * @property {number} p - Dimension
 * @property {number} n - Sample count
 */

/**
 * Computes the empirical sample covariance matrix from an N_samples x P_dim data matrix.
 *
 * @param {Float64Array|Array<number>|Array<Array<number>>} X - Data matrix (n x p)
 * @param {number} n - Number of observation samples
 * @param {number} p - Number of variables / dimensions
 * @param {boolean} [biasCorrected=true] - If true, divide by (n-1), else n
 * @returns {{ S: Float64Array, means: Float64Array, centeredX: Float64Array }}
 */
export function computeSampleCovariance(X, n, p, biasCorrected = true) {
  if (n < 2) {
    throw new RangeError(`Sample size n must be >= 2 for covariance estimation: received ${n}`);
  }
  const flatX = matFrom(X, n, p);
  const means = new Float64Array(p);

  for (let i = 0; i < n; i++) {
    const rowOffset = i * p;
    for (let j = 0; j < p; j++) {
      means[j] += flatX[rowOffset + j];
    }
  }
  for (let j = 0; j < p; j++) {
    means[j] /= n;
  }

  const centeredX = new Float64Array(n * p);
  for (let i = 0; i < n; i++) {
    const rowOffset = i * p;
    for (let j = 0; j < p; j++) {
      centeredX[rowOffset + j] = flatX[rowOffset + j] - means[j];
    }
  }

  const denom = biasCorrected ? (n - 1) : n;
  const S = new Float64Array(p * p);

  for (let i = 0; i < p; i++) {
    for (let j = i; j < p; j++) {
      let sum = 0.0;
      for (let k = 0; k < n; k++) {
        sum += centeredX[k * p + i] * centeredX[k * p + j];
      }
      const val = sum / denom;
      S[i * p + j] = val;
      S[j * p + i] = val;
    }
  }

  return { S, means, centeredX };
}

/**
 * Computes Ledoit-Wolf (2004) Optimal Linear Shrinkage Covariance Matrix.
 *
 * Analytical Formulation:
 * - \Sigma^* = (1 - \lambda^*) S + \lambda^* F
 * - Target A: F = \mu I where \mu = \text{Tr}(S) / p
 * - Asymptotically optimal MSE minimization without Gaussianity assumption:
 *   \delta^2 = ||S - F||_F^2
 *   \bar{b}^2 = \frac{1}{n^2} \sum_{k=1}^n || y_k - S ||_F^2 where y_k = x_k x_k^T (centered)
 *   b^2 = \min(\bar{b}^2, \delta^2)
 *   \lambda^* = b^2 / \delta^2
 *
 * @param {Float64Array|Array<number>|Array<Array<number>>} X - Data matrix (n x p)
 * @param {number} n - Number of samples
 * @param {number} p - Dimension
 * @param {ShrinkageTarget} [targetType=ShrinkageTarget.SCALED_IDENTITY] - Target prior
 * @returns {LedoitWolfResult}
 */
export function ledoitWolfShrinkage(X, n, p, targetType = ShrinkageTarget.SCALED_IDENTITY) {
  const { S, centeredX } = computeSampleCovariance(X, n, p, false);

  const trS = matTrace(S, p);
  const mu = trS / p;

  const F = new Float64Array(p * p);
  if (targetType === ShrinkageTarget.SCALED_IDENTITY) {
    for (let i = 0; i < p; i++) F[i * p + i] = mu;
  } else if (targetType === ShrinkageTarget.DIAGONAL) {
    for (let i = 0; i < p; i++) F[i * p + i] = S[i * p + i];
  } else if (targetType === ShrinkageTarget.IDENTITY) {
    for (let i = 0; i < p; i++) F[i * p + i] = 1.0;
  } else if (targetType === ShrinkageTarget.CONSTANT_CORRELATION) {
    let sumCorr = 0.0;
    let countCorr = 0;
    for (let i = 0; i < p; i++) {
      const varI = Math.max(1e-15, S[i * p + i]);
      for (let j = i + 1; j < p; j++) {
        const varJ = Math.max(1e-15, S[j * p + j]);
        const r_ij = S[i * p + j] / Math.sqrt(varI * varJ);
        sumCorr += r_ij;
        countCorr++;
      }
    }
    const rBar = countCorr > 0 ? sumCorr / countCorr : 0.0;
    for (let i = 0; i < p; i++) {
      const stdI = Math.sqrt(Math.max(1e-15, S[i * p + i]));
      for (let j = 0; j < p; j++) {
        const stdJ = Math.sqrt(Math.max(1e-15, S[j * p + j]));
        F[i * p + j] = i === j ? S[i * p + i] : rBar * stdI * stdJ;
      }
    }
  }

  let deltaSq = 0.0;
  for (let i = 0; i < p * p; i++) {
    const diff = S[i] - F[i];
    deltaSq += diff * diff;
  }

  let sumNormSq = 0.0;
  for (let k = 0; k < n; k++) {
    const kOffset = k * p;
    for (let i = 0; i < p; i++) {
      const xki = centeredX[kOffset + i];
      for (let j = 0; j < p; j++) {
        const xkj = centeredX[kOffset + j];
        const y_kij = xki * xkj;
        const diff = y_kij - S[i * p + j];
        sumNormSq += diff * diff;
      }
    }
  }

  const bBarSq = sumNormSq / (n * n);
  const bSq = Math.min(bBarSq, deltaSq);

  let lambda = deltaSq > 1e-15 ? (bSq / deltaSq) : 0.0;
  lambda = Math.max(0.0, Math.min(1.0, lambda));

  const regularized = new Float64Array(p * p);
  for (let i = 0; i < p * p; i++) {
    regularized[i] = (1.0 - lambda) * S[i] + lambda * F[i];
  }

  return {
    regularizedCovariance: regularized,
    sampleCovariance: S,
    targetMatrix: F,
    shrinkageIntensity: lambda,
    deltaSq,
    bSq,
    mu,
    p,
    n
  };
}

// ============================================================================
// 6. ORACLE APPROXIMATING SHRINKAGE (OAS) & RBLW ESTIMATOR
// ============================================================================

/**
 * Result of Oracle Approximating Shrinkage.
 * @typedef {Object} OASResult
 * @property {Float64Array} regularizedCovariance - OAS Shrunk covariance matrix
 * @property {Float64Array} sampleCovariance - Sample covariance matrix S
 * @property {number} shrinkageIntensity - OAS optimal shrinkage intensity \rho_{OAS} \in [0, 1]
 * @property {number} mu - Trace mean \text{Tr}(S)/p
 * @property {number} trS - Trace of S
 * @property {number} trS2 - Trace of S^2
 * @property {number} p - Dimension
 * @property {number} n - Number of samples
 */

/**
 * Oracle Approximating Shrinkage (OAS) Covariance Estimator (Chen et al. 2010).
 *
 * Formula:
 * \mu = \text{Tr}(S) / p
 * \hat{\rho}_{OAS} = \frac{ (1 - 2/p) \text{Tr}(S^2) + \text{Tr}^2(S) }{ (n + 1 - 2/p) (\text{Tr}(S^2) - \frac{1}{p}\text{Tr}^2(S)) }
 * \rho_{OAS}^* = \min\left(1, \max\left(0, \hat{\rho}_{OAS}\right)\right)
 * \Sigma_{OAS} = (1 - \rho_{OAS}^*) S + \rho_{OAS}^* \mu I
 *
 * @param {Float64Array|Array<number>|Array<Array<number>>} X - Data matrix (n x p)
 * @param {number} n - Samples
 * @param {number} p - Dimension
 * @returns {OASResult}
 */
export function oracleApproximatingShrinkage(X, n, p) {
  const { S } = computeSampleCovariance(X, n, p, false);

  const trS = matTrace(S, p);
  const mu = trS / p;
  const trS2 = matTraceSquare(S, p);
  const tr2S = trS * trS;

  const num = (1.0 - 2.0 / p) * trS2 + tr2S;
  const denom = (n + 1.0 - 2.0 / p) * (trS2 - tr2S / p);

  let rho = 0.0;
  if (Math.abs(denom) > 1e-15) {
    rho = num / denom;
  }
  rho = Math.max(0.0, Math.min(1.0, rho));

  const regularized = new Float64Array(p * p);
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      const sVal = S[i * p + j];
      const targetVal = i === j ? mu : 0.0;
      regularized[i * p + j] = (1.0 - rho) * sVal + rho * targetVal;
    }
  }

  return {
    regularizedCovariance: regularized,
    sampleCovariance: S,
    shrinkageIntensity: rho,
    mu,
    trS,
    trS2,
    p,
    n
  };
}

/**
 * Rao-Blackwell Ledoit-Wolf (RBLW) Estimator (Chen et al. 2010).
 *
 * Formula:
 * \hat{\rho}_{RBLW} = \frac{ \frac{n-2}{n} \text{Tr}(S^2) + \text{Tr}^2(S) }{ (n + 2) (\text{Tr}(S^2) - \frac{1}{p}\text{Tr}^2(S)) }
 *
 * @param {Float64Array|Array<number>|Array<Array<number>>} X - Data matrix (n x p)
 * @param {number} n - Samples
 * @param {number} p - Dimension
 * @returns {OASResult}
 */
export function raoBlackwellLedoitWolf(X, n, p) {
  const { S } = computeSampleCovariance(X, n, p, false);

  const trS = matTrace(S, p);
  const mu = trS / p;
  const trS2 = matTraceSquare(S, p);
  const tr2S = trS * trS;

  const num = ((n - 2.0) / n) * trS2 + tr2S;
  const denom = (n + 2.0) * (trS2 - tr2S / p);

  let rho = 0.0;
  if (Math.abs(denom) > 1e-15) {
    rho = num / denom;
  }
  rho = Math.max(0.0, Math.min(1.0, rho));

  const regularized = new Float64Array(p * p);
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      const sVal = S[i * p + j];
      const targetVal = i === j ? mu : 0.0;
      regularized[i * p + j] = (1.0 - rho) * sVal + rho * targetVal;
    }
  }

  return {
    regularizedCovariance: regularized,
    sampleCovariance: S,
    shrinkageIntensity: rho,
    mu,
    trS,
    trS2,
    p,
    n
  };
}

// ============================================================================
// 7. SPATIAL COVARIANCE TAPERING & COMPACT SUPPORT KERNELS
// ============================================================================

/**
 * Compact Support Kernel Functions for Spatial Covariance Tapering (Furrer et al. 2006).
 */
export const TaperKernels = {
  /**
   * Wendland C^0 Kernel in R^3:
   * K(r) = (1 - r)_+^2
   * @param {number} r - Scaled distance d / theta
   * @returns {number}
   */
  wendlandC0(r) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    const t = 1.0 - r;
    return t * t;
  },

  /**
   * Wendland C^2 Kernel in R^3 (Compactly supported, positive definite on R^3):
   * K(r) = (1 - r)_+^4 (1 + 4r)
   * @param {number} r - Scaled distance d / theta
   * @returns {number}
   */
  wendlandC2(r) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    const t = 1.0 - r;
    const t2 = t * t;
    return t2 * t2 * (1.0 + 4.0 * r);
  },

  /**
   * Wendland C^4 Kernel in R^3:
   * K(r) = (1 - r)_+^6 (1 + 6r + (35/3) r^2)
   * @param {number} r - Scaled distance d / theta
   * @returns {number}
   */
  wendlandC4(r) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    const t = 1.0 - r;
    const t3 = t * t * t;
    return t3 * t3 * (1.0 + 6.0 * r + (35.0 / 3.0) * r * r);
  },

  /**
   * Wendland C^6 Kernel in R^3:
   * K(r) = (1 - r)_+^8 (1 + 8r + 25r^2 + 32r^3)
   * @param {number} r - Scaled distance d / theta
   * @returns {number}
   */
  wendlandC6(r) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    const t = 1.0 - r;
    const t4 = (t * t) * (t * t);
    return t4 * t4 * (1.0 + 8.0 * r + 25.0 * r * r + 32.0 * r * r * r);
  },

  /**
   * Spherical Tapering Kernel:
   * K(r) = (1 - 1.5r + 0.5r^3) for r <= 1, 0 otherwise
   * @param {number} r - Scaled distance d / theta
   * @returns {number}
   */
  spherical(r) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    return 1.0 - 1.5 * r + 0.5 * r * r * r;
  },

  /**
   * Gaussian Compactly-Tapered Kernel:
   * K(r) = \exp(-0.5 (r/\sigma)^2) \cdot (1 - r^2)_+^2
   * @param {number} r - Scaled distance d / theta
   * @param {number} [sigma=0.35] - Gaussian width
   * @returns {number}
   */
  gaussianTapered(r, sigma = 0.35) {
    if (r >= 1.0 || r < 0.0) return 0.0;
    const g = Math.exp(-0.5 * (r / sigma) * (r / sigma));
    const t = 1.0 - r * r;
    return g * t * t;
  }
};

/**
 * Result of Spatial Covariance Tapering.
 * @typedef {Object} TaperedCovarianceResult
 * @property {Float64Array} taperedMatrix - Tapered covariance C_{tap} = C \circ K(D/\theta)
 * @property {Float64Array} taperKernelMatrix - Kernel correlation matrix K(D/\theta)
 * @property {number} sparsity - Percentage of zero entries (0.0 to 1.0)
 * @property {number} nonZeroCount - Number of non-zero entries
 * @property {number} theta - Taper range / support radius
 * @property {string} kernelName - Name of compact kernel used
 */

/**
 * Computes spatial tapered covariance matrix C_{tap} = C \circ K(D/\theta).
 *
 * Preserves positive definiteness under the Schur Product Theorem:
 * If C is PSD and K is a positive definite kernel on R^d, then C \circ K is strictly PSD.
 *
 * @param {Float64Array|Array<number>} C - Base covariance matrix (N x N)
 * @param {Array<[number, number, number]>|Float64Array} coords - 3D Coordinates (N points, [x, y, z])
 * @param {number} n - Number of spatial points N
 * @param {number} theta - Taper support radius (cutoff distance)
 * @param {string} [kernelType='wendlandC2'] - Kernel type: 'wendlandC0', 'wendlandC2', 'wendlandC4', 'wendlandC6', 'spherical', 'gaussianTapered'
 * @returns {TaperedCovarianceResult}
 */
export function spatialCovarianceTapering(C, coords, n, theta, kernelType = 'wendlandC2') {
  if (theta <= 0) {
    throw new RangeError(`Taper range theta must be positive: received ${theta}`);
  }

  const kernelFn = TaperKernels[kernelType] || TaperKernels.wendlandC2;
  const tapered = new Float64Array(n * n);
  const kernelMat = new Float64Array(n * n);

  let flatCoords;
  if (coords instanceof Float64Array || Array.isArray(coords[0]) === false) {
    flatCoords = coords;
  } else {
    flatCoords = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) {
      flatCoords[i * 3 + 0] = coords[i][0];
      flatCoords[i * 3 + 1] = coords[i][1];
      flatCoords[i * 3 + 2] = coords[i][2];
    }
  }

  let nonZeroCount = 0;

  for (let i = 0; i < n; i++) {
    const x1 = flatCoords[i * 3 + 0];
    const y1 = flatCoords[i * 3 + 1];
    const z1 = flatCoords[i * 3 + 2];

    tapered[i * n + i] = C[i * n + i];
    kernelMat[i * n + i] = 1.0;
    nonZeroCount++;

    for (let j = i + 1; j < n; j++) {
      const dx = x1 - flatCoords[j * 3 + 0];
      const dy = y1 - flatCoords[j * 3 + 1];
      const dz = z1 - flatCoords[j * 3 + 2];
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const r = dist / theta;

      const kVal = kernelFn(r);
      kernelMat[i * n + j] = kVal;
      kernelMat[j * n + i] = kVal;

      if (kVal > 0.0) {
        const cVal = C[i * n + j] * kVal;
        tapered[i * n + j] = cVal;
        tapered[j * n + i] = cVal;
        nonZeroCount += 2;
      } else {
        tapered[i * n + j] = 0.0;
        tapered[j * n + i] = 0.0;
      }
    }
  }

  const totalEntries = n * n;
  const sparsity = (totalEntries - nonZeroCount) / totalEntries;

  return {
    taperedMatrix: tapered,
    taperKernelMatrix: kernelMat,
    sparsity,
    nonZeroCount,
    theta,
    kernelName: kernelType
  };
}

/**
 * Compressed Sparse Row (CSR) Sparse Matrix Representation for Large Tapered Covariances.
 */
export class CSRSparseMatrix {
  /**
   * @param {number} rows
   * @param {number} cols
   * @param {Float64Array} values
   * @param {Int32Array} colIndices
   * @param {Int32Array} rowPointers
   */
  constructor(rows, cols, values, colIndices, rowPointers) {
    this.rows = rows;
    this.cols = cols;
    this.values = values;
    this.colIndices = colIndices;
    this.rowPointers = rowPointers;
  }

  /**
   * Creates CSR matrix from dense matrix with given threshold.
   * @param {Float64Array|Array<number>} dense
   * @param {number} rows
   * @param {number} cols
   * @param {number} [threshold=1e-15]
   * @returns {CSRSparseMatrix}
   */
  static fromDense(dense, rows, cols, threshold = 1e-15) {
    const valList = [];
    const colList = [];
    const rowPtrList = new Int32Array(rows + 1);

    let count = 0;
    for (let i = 0; i < rows; i++) {
      rowPtrList[i] = count;
      const offset = i * cols;
      for (let j = 0; j < cols; j++) {
        const val = dense[offset + j];
        if (Math.abs(val) > threshold) {
          valList.push(val);
          colList.push(j);
          count++;
        }
      }
    }
    rowPtrList[rows] = count;

    return new CSRSparseMatrix(
      rows,
      cols,
      new Float64Array(valList),
      new Int32Array(colList),
      rowPtrList
    );
  }

  /**
   * Sparse Matrix-Vector product y = A * x.
   * @param {Float64Array|Array<number>} x
   * @returns {Float64Array}
   */
  multiplyVector(x) {
    const y = new Float64Array(this.rows);
    for (let i = 0; i < this.rows; i++) {
      let sum = 0.0;
      const start = this.rowPointers[i];
      const end = this.rowPointers[i + 1];
      for (let idx = start; idx < end; idx++) {
        sum += this.values[idx] * x[this.colIndices[idx]];
      }
      y[i] = sum;
    }
    return y;
  }
}

// ============================================================================
// 8. 3D COSMOLOGICAL VELOCITY DISPERSION FIELD & CUBE REGULARIZER
// ============================================================================

/**
 * 3D Velocity Dispersion Tensor Diagnostics.
 * @typedef {Object} DispersionDiagnostics
 * @property {Float64Array} tensor - 3x3 symmetric dispersion tensor \sigma_{ij}
 * @property {number} trace - Trace \sigma_{xx} + \sigma_{yy} + \sigma_{zz}
 * @property {number} isotropicDispersion - Isotropic 1D dispersion \sigma_{1D} = \sqrt{\text{Tr}(\sigma)/3}
 * @property {number} sigma3D - Total 3D velocity dispersion \sqrt{\text{Tr}(\sigma)}
 * @property {number[]} principalDispersions - Principal axes dispersions [\sigma_1, \sigma_2, \sigma_3] in descending order
 * @property {Float64Array} principalAxes - 3x3 Matrix of principal eigenvectors
 * @property {number} anisotropyParameter - Anisotropy \beta = 1 - (\sigma_\theta^2 + \sigma_\phi^2) / (2 \sigma_r^2)
 * @property {number} triaxiality - Triaxiality parameter T = (\sigma_1^2 - \sigma_2^2) / (\sigma_1^2 - \sigma_3^2)
 * @property {number} ellipticity - Ellipticity e = (\sigma_1 - \sigma_3) / (2 (\sigma_1 + \sigma_2 + \sigma_3))
 * @property {number} prolateness - Prolateness p = (\sigma_1 - 2\sigma_2 + \sigma_3) / (2 (\sigma_1 + \sigma_2 + \sigma_3))
 */

/**
 * Computes 3D Velocity Dispersion Tensor \sigma_{ij} = \langle v_i v_j \rangle - \langle v_i \rangle \langle v_j \rangle
 * and cosmological shape/anisotropy parameters.
 *
 * @param {Array<[number, number, number]>|Float64Array} velocities - Array of 3D velocity vectors [vx, vy, vz]
 * @param {number[]} [weights] - Optional particle weights / kernel values
 * @param {[number, number, number]} [radialDirection] - Unit radial vector \hat{r} for anisotropy \beta computation
 * @returns {DispersionDiagnostics}
 */
export function computeVelocityDispersionTensor(velocities, weights = null, radialDirection = null) {
  const n = Array.isArray(velocities[0]) ? velocities.length : (velocities.length / 3);
  if (n < 2) {
    throw new RangeError(`At least 2 velocity samples required: received ${n}`);
  }

  let sumW = 0.0;
  const meanV = new Float64Array(3);
  const secondMoments = new Float64Array(9);

  for (let k = 0; k < n; k++) {
    let vx, vy, vz;
    if (Array.isArray(velocities[0])) {
      vx = velocities[k][0];
      vy = velocities[k][1];
      vz = velocities[k][2];
    } else {
      vx = velocities[k * 3 + 0];
      vy = velocities[k * 3 + 1];
      vz = velocities[k * 3 + 2];
    }

    const w = weights ? weights[k] : 1.0;
    sumW += w;

    meanV[0] += w * vx;
    meanV[1] += w * vy;
    meanV[2] += w * vz;

    secondMoments[0] += w * vx * vx;
    secondMoments[1] += w * vx * vy;
    secondMoments[2] += w * vx * vz;
    secondMoments[3] += w * vy * vx;
    secondMoments[4] += w * vy * vy;
    secondMoments[5] += w * vy * vz;
    secondMoments[6] += w * vz * vx;
    secondMoments[7] += w * vz * vy;
    secondMoments[8] += w * vz * vz;
  }

  if (sumW <= 0) {
    throw new Error('Sum of weights must be strictly positive');
  }

  meanV[0] /= sumW;
  meanV[1] /= sumW;
  meanV[2] /= sumW;

  const tensor = new Float64Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      tensor[i * 3 + j] = (secondMoments[i * 3 + j] / sumW) - (meanV[i] * meanV[j]);
    }
  }
  const symTensor = matSymmetrize(tensor, 3);

  const eig = eigen3x3Symmetric(symTensor);
  const s1_sq = Math.max(0.0, eig.values[0]);
  const s2_sq = Math.max(0.0, eig.values[1]);
  const s3_sq = Math.max(0.0, eig.values[2]);

  const s1 = Math.sqrt(s1_sq);
  const s2 = Math.sqrt(s2_sq);
  const s3 = Math.sqrt(s3_sq);

  const tr = symTensor[0] + symTensor[4] + symTensor[8];
  const isoDisp = Math.sqrt(Math.max(0.0, tr / 3.0));
  const sigma3D = Math.sqrt(Math.max(0.0, tr));

  const denomTriax = s1_sq - s3_sq;
  const triaxiality = denomTriax > 1e-12 ? (s1_sq - s2_sq) / denomTriax : 0.5;

  const sumS = s1 + s2 + s3;
  const ellipticity = sumS > 1e-12 ? (s1 - s3) / (2.0 * sumS) : 0.0;
  const prolateness = sumS > 1e-12 ? (s1 - 2.0 * s2 + s3) / (2.0 * sumS) : 0.0;

  let beta = 0.0;
  if (radialDirection) {
    const rx = radialDirection[0];
    const ry = radialDirection[1];
    const rz = radialDirection[2];
    const rNorm = Math.sqrt(rx * rx + ry * ry + rz * rz);
    if (rNorm > 1e-12) {
      const ur = [rx / rNorm, ry / rNorm, rz / rNorm];
      const sigmaR_sq = ur[0] * (symTensor[0] * ur[0] + symTensor[1] * ur[1] + symTensor[2] * ur[2]) +
                        ur[1] * (symTensor[3] * ur[0] + symTensor[4] * ur[1] + symTensor[5] * ur[2]) +
                        ur[2] * (symTensor[6] * ur[0] + symTensor[7] * ur[1] + symTensor[8] * ur[2]);

      const sigmaTan_sq = Math.max(0.0, tr - sigmaR_sq);
      if (sigmaR_sq > 1e-12) {
        beta = 1.0 - sigmaTan_sq / (2.0 * sigmaR_sq);
      }
    }
  }

  return {
    tensor: symTensor,
    trace: tr,
    isotropicDispersion: isoDisp,
    sigma3D,
    principalDispersions: [s1, s2, s3],
    principalAxes: eig.vectors,
    anisotropyParameter: beta,
    triaxiality,
    ellipticity,
    prolateness
  };
}

/**
 * 3D Velocity Dispersion Regularizer for Cosmological Voxel Grids.
 */
export class VelocityDispersionCubeRegularizer {
  /**
   * @param {Object} options
   * @param {number} [options.globalIsotropicDispersion=300.0] - Prior isotropic dispersion \sigma_0 (km/s)
   * @param {number} [options.minSamplesForDirect=10] - Sample threshold below which shrinkage is active
   * @param {string} [options.shrinkageMethod='OAS'] - 'OAS', 'LedoitWolf', or 'Tikhonov'
   * @param {number} [options.tikhonovEpsilon=100.0] - Fixed ridge term for Tikhonov mode
   */
  constructor(options = {}) {
    this.sigma0 = options.globalIsotropicDispersion || 300.0;
    this.minSamples = options.minSamplesForDirect || 10;
    this.shrinkageMethod = options.shrinkageMethod || 'OAS';
    this.tikhonovEpsilon = options.tikhonovEpsilon || 100.0;
  }

  /**
   * Regularizes a single voxel dispersion tensor given its sample particle velocities.
   *
   * @param {Array<[number, number, number]>|Float64Array} sampleVelocities - Velocities in voxel
   * @returns {DispersionDiagnostics}
   */
  regularizeVoxel(sampleVelocities) {
    const count = Array.isArray(sampleVelocities[0]) ? sampleVelocities.length : (sampleVelocities.length / 3);

    // If voxel is completely empty, return isotropic prior
    if (count === 0) {
      const priorTensor = new Float64Array(9);
      const var0 = this.sigma0 * this.sigma0;
      priorTensor[0] = var0;
      priorTensor[4] = var0;
      priorTensor[8] = var0;
      return {
        tensor: priorTensor,
        trace: 3.0 * var0,
        isotropicDispersion: this.sigma0,
        sigma3D: Math.sqrt(3.0) * this.sigma0,
        principalDispersions: [this.sigma0, this.sigma0, this.sigma0],
        principalAxes: eye(3),
        anisotropyParameter: 0.0,
        triaxiality: 0.5,
        ellipticity: 0.0,
        prolateness: 0.0
      };
    }

    if (count === 1) {
      const priorTensor = new Float64Array(9);
      const var0 = this.sigma0 * this.sigma0;
      priorTensor[0] = var0;
      priorTensor[4] = var0;
      priorTensor[8] = var0;
      return {
        tensor: priorTensor,
        trace: 3.0 * var0,
        isotropicDispersion: this.sigma0,
        sigma3D: Math.sqrt(3.0) * this.sigma0,
        principalDispersions: [this.sigma0, this.sigma0, this.sigma0],
        principalAxes: eye(3),
        anisotropyParameter: 0.0,
        triaxiality: 0.5,
        ellipticity: 0.0,
        prolateness: 0.0
      };
    }

    const rawDiag = computeVelocityDispersionTensor(sampleVelocities);

    let regTensor;
    if (this.shrinkageMethod === 'OAS') {
      const oas = oracleApproximatingShrinkage(sampleVelocities, count, 3);
      regTensor = nearestPositiveSemiDefinite(oas.regularizedCovariance, 3, 1e-4);
    } else if (this.shrinkageMethod === 'LedoitWolf') {
      const lw = ledoitWolfShrinkage(sampleVelocities, count, 3, ShrinkageTarget.SCALED_IDENTITY);
      regTensor = nearestPositiveSemiDefinite(lw.regularizedCovariance, 3, 1e-4);
    } else {
      regTensor = tikhonovRegularize(rawDiag.tensor, 3, this.tikhonovEpsilon);
    }

    const eig = eigen3x3Symmetric(regTensor);
    const s1_sq = Math.max(0.0, eig.values[0]);
    const s2_sq = Math.max(0.0, eig.values[1]);
    const s3_sq = Math.max(0.0, eig.values[2]);

    const s1 = Math.sqrt(s1_sq);
    const s2 = Math.sqrt(s2_sq);
    const s3 = Math.sqrt(s3_sq);

    const tr = regTensor[0] + regTensor[4] + regTensor[8];
    const isoDisp = Math.sqrt(Math.max(0.0, tr / 3.0));
    const sigma3D = Math.sqrt(Math.max(0.0, tr));

    const denomTriax = s1_sq - s3_sq;
    const triaxiality = denomTriax > 1e-12 ? (s1_sq - s2_sq) / denomTriax : 0.5;

    const sumS = s1 + s2 + s3;
    const ellipticity = sumS > 1e-12 ? (s1 - s3) / (2.0 * sumS) : 0.0;
    const prolateness = sumS > 1e-12 ? (s1 - 2.0 * s2 + s3) / (2.0 * sumS) : 0.0;

    return {
      tensor: regTensor,
      trace: tr,
      isotropicDispersion: isoDisp,
      sigma3D,
      principalDispersions: [s1, s2, s3],
      principalAxes: eig.vectors,
      anisotropyParameter: rawDiag.anisotropyParameter,
      triaxiality,
      ellipticity,
      prolateness
    };
  }
}

// ============================================================================
// 9. RIEMANNIAN GEOMETRY AND MATRIX DISTANCE METRICS
// ============================================================================

/**
 * Computes the Affine-Invariant Riemannian Metric (AIRM) distance between two PSD matrices:
 * \delta_R(A, B) = || \ln(A^{-1/2} B A^{-1/2}) ||_F = \sqrt{\sum_{i=1}^n \ln^2(\lambda_i)}
 * where \lambda_i are generalized eigenvalues of (B, A).
 *
 * @param {Float64Array|Array<number>} A - Symmetric PSD matrix (n x n)
 * @param {Float64Array|Array<number>} B - Symmetric PSD matrix (n x n)
 * @param {number} n - Dimension
 * @returns {number}
 */
export function affineInvariantRiemannianDistance(A, B, n) {
  const A_inv_sqrt = matInvSqrtSymmetric(A, n);
  const mid = matMul(A_inv_sqrt, B, n, n, n);
  const symMid = matMul(mid, A_inv_sqrt, n, n, n);
  const logMid = matLogSymmetric(symMid, n);
  return matFrobeniusNorm(logMid);
}

/**
 * Computes the Bures-Wasserstein distance between two zero-mean Gaussians:
 * d_{BW}^2(\Sigma_1, \Sigma_2) = \text{Tr}(\Sigma_1) + \text{Tr}(\Sigma_2) - 2 \text{Tr}\left((\Sigma_1^{1/2} \Sigma_2 \Sigma_1^{1/2})^{1/2}\right)
 *
 * @param {Float64Array|Array<number>} S1 - Covariance 1 (n x n)
 * @param {Float64Array|Array<number>} S2 - Covariance 2 (n x n)
 * @param {number} n - Dimension
 * @returns {number}
 */
export function buresWassersteinDistance(S1, S2, n) {
  const S1_sqrt = matSqrtSymmetric(S1, n);
  const mid = matMul(S1_sqrt, S2, n, n, n);
  const M = matMul(mid, S1_sqrt, n, n, n);
  const M_sqrt = matSqrtSymmetric(M, n);

  const tr1 = matTrace(S1, n);
  const tr2 = matTrace(S2, n);
  const trM = matTrace(M_sqrt, n);

  const distSq = Math.max(0.0, tr1 + tr2 - 2.0 * trM);
  return Math.sqrt(distSq);
}

/**
 * Computes Kullback-Leibler (KL) Divergence D_{KL}(\mathcal{N}_0 || \mathcal{N}_1) between zero-mean multivariate Gaussians:
 * D_{KL} = \frac{1}{2} \left[ \text{Tr}(\Sigma_1^{-1} \Sigma_0) - n + \ln \frac{\det \Sigma_1}{\det \Sigma_0} \right]
 *
 * @param {Float64Array|Array<number>} S0 - Covariance 0 (n x n)
 * @param {Float64Array|Array<number>} S1 - Covariance 1 (n x n)
 * @param {number} n - Dimension
 * @returns {number}
 */
export function gaussianKLDivergence(S0, S1, n) {
  const S1_inv = choleskyInvert(S1, n);
  const S1inv_S0 = matMul(S1_inv, S0, n, n, n);
  const trTerm = matTrace(S1inv_S0, n);

  const logDet0 = logDeterminant(S0, n);
  const logDet1 = logDeterminant(S1, n);

  const kl = 0.5 * (trTerm - n + (logDet1 - logDet0));
  return Math.max(0.0, kl);
}
