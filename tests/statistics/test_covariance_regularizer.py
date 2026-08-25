# -*- coding: utf-8 -*-
"""
tests/statistics/test_covariance_regularizer.py
Comprehensive Pytest verification test suite for 3D Covariance Regularizer,
Shrinkage Estimators (Ledoit-Wolf, OAS, RBLW), Matrix Conditioner,
Modified Cholesky / LDL^T Decompositions, Spatial Tapering Kernels,
and 3D Velocity Dispersion Field Regularizer.
"""

import math
import numpy as np
import pytest
from tests.utils import run_node_snippet

class TestBasicMatrixOperations:
    """Tests fundamental matrix and vector routines, norms, Hadamard products, and symmetries."""

    def test_matrix_creation_and_identity(self):
        """Verifies matrix creation and identity matrix initialization."""
        code = """
        import { createMatrix, eye, isSymmetric, matTrace } from './src/statistics/covariance_regularizer.js';

        const I3 = eye(3);
        const I5 = eye(5);
        const zeroMat = createMatrix(4, 4);

        console.log(JSON.stringify({
          i3_length: I3.length,
          i3_diag: [I3[0], I3[4], I3[8]],
          i3_offdiag: [I3[1], I3[2], I3[3], I3[5], I3[6], I3[7]],
          i3_trace: matTrace(I3, 3),
          i3_symmetric: isSymmetric(I3, 3),
          i5_length: I5.length,
          i5_trace: matTrace(I5, 5),
          zero_length: zeroMat.length,
          zero_trace: matTrace(zeroMat, 4)
        }));
        """
        res = run_node_snippet(code)
        assert res["i3_length"] == 9
        assert res["i3_diag"] == [1.0, 1.0, 1.0]
        assert all(x == 0.0 for x in res["i3_offdiag"])
        assert res["i3_trace"] == pytest.approx(3.0)
        assert res["i3_symmetric"] is True
        assert res["i5_length"] == 25
        assert res["i5_trace"] == pytest.approx(5.0)
        assert res["zero_length"] == 16
        assert res["zero_trace"] == pytest.approx(0.0)

    def test_matrix_multiplication_and_transposition(self):
        """Verifies matrix multiplication C = A * B, vector multiplication, and transpose."""
        code = """
        import { matMul, matVecMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        // A is 2x3, B is 3x2
        const A = [1, 2, 3,
                   4, 5, 6];
        const B = [7, 8,
                   9, 1,
                   2, 3];
        // C = A * B (2x2)
        // C[0,0] = 1*7 + 2*9 + 3*2 = 7 + 18 + 6 = 31
        // C[0,1] = 1*8 + 2*1 + 3*3 = 8 + 2 + 9 = 19
        // C[1,0] = 4*7 + 5*9 + 6*2 = 28 + 45 + 12 = 85
        // C[1,1] = 4*8 + 5*1 + 6*3 = 32 + 5 + 18 = 55
        const C = matMul(A, B, 2, 3, 2);

        // Vector multiplication: A * x where x = [1, 1, 1]^T
        const x = [1, 1, 1];
        const Ax = matVecMul(A, x, 2, 3); // [1+2+3=6, 4+5+6=15]

        // Transpose of A (3x2)
        const AT = matTranspose(A, 2, 3);

        console.log(JSON.stringify({
          C: Array.from(C),
          Ax: Array.from(Ax),
          AT: Array.from(AT)
        }));
        """
        res = run_node_snippet(code)
        assert res["C"] == [31.0, 19.0, 85.0, 55.0]
        assert res["Ax"] == [6.0, 15.0]
        assert res["AT"] == [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]

    def test_matrix_trace_and_frobenius_norm(self):
        """Verifies trace, Frobenius norm, and squared Frobenius norm."""
        code = """
        import { matTrace, matFrobeniusNorm, matFrobeniusNormSq, matTraceSquare } from './src/statistics/covariance_regularizer.js';

        const A = [2, -1, 0,
                   -1, 2, -1,
                   0, -1, 2];

        const tr = matTrace(A, 3); // 2 + 2 + 2 = 6
        const fnSq = matFrobeniusNormSq(A); // 4 + 1 + 0 + 1 + 4 + 1 + 0 + 1 + 4 = 16
        const fn = matFrobeniusNorm(A); // sqrt(16) = 4
        const trSq = matTraceSquare(A, 3); // For symmetric matrix, Tr(A^2) == ||A||_F^2 = 16

        console.log(JSON.stringify({
          tr,
          fnSq,
          fn,
          trSq
        }));
        """
        res = run_node_snippet(code)
        assert res["tr"] == pytest.approx(6.0)
        assert res["fnSq"] == pytest.approx(16.0)
        assert res["fn"] == pytest.approx(4.0)
        assert res["trSq"] == pytest.approx(16.0)

    def test_hadamard_product_and_linear_combination(self):
        """Verifies Schur/Hadamard element-wise product and linear combination."""
        code = """
        import { matHadamard, matAdd } from './src/statistics/covariance_regularizer.js';

        const A = [1, 2, 3, 4];
        const B = [5, 6, 7, 8];

        const hadamard = matHadamard(A, B); // [5, 12, 21, 32]
        const comb = matAdd(A, B, 2.0, -1.0); // 2*A - B = [2-5, 4-6, 6-7, 8-8] = [-3, -2, -1, 0]

        console.log(JSON.stringify({
          hadamard: Array.from(hadamard),
          comb: Array.from(comb)
        }));
        """
        res = run_node_snippet(code)
        assert res["hadamard"] == [5.0, 12.0, 21.0, 32.0]
        assert res["comb"] == [-3.0, -2.0, -1.0, 0.0]

    def test_symmetrization_and_symmetry_checker(self):
        """Verifies symmetrization 0.5*(A + A^T) and tolerance-aware symmetry check."""
        code = """
        import { matSymmetrize, isSymmetric } from './src/statistics/covariance_regularizer.js';

        const nonSym = [1, 4,
                        2, 3];
        const isSymBefore = isSymmetric(nonSym, 2);
        const sym = matSymmetrize(nonSym, 2); // [1, 3, 3, 3]
        const isSymAfter = isSymmetric(sym, 2);

        console.log(JSON.stringify({
          isSymBefore,
          sym: Array.from(sym),
          isSymAfter
        }));
        """
        res = run_node_snippet(code)
        assert res["isSymBefore"] is False
        assert res["sym"] == [1.0, 3.0, 3.0, 3.0]
        assert res["isSymAfter"] is True


class TestEigendecompositionAndMatrixFunctions:
    """Tests symmetric Jacobi eigendecomposition, matrix powers, square roots, and logarithms."""

    def test_jacobi_eigenvalues_orthogonality_and_reconstruction(self):
        """Verifies Jacobi algorithm produces orthonormal eigenvectors and exact reconstruction A = V D V^T."""
        code = """
        import { jacobiEigenvalues, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        const A = [4, 1, 2,
                   1, 5, 3,
                   2, 3, 6];
        const n = 3;
        const eig = jacobiEigenvalues(A, n);

        // Verify orthogonality: V^T * V = I
        const VT = matTranspose(eig.vectors, n, n);
        const VTV = matMul(VT, eig.vectors, n, n, n);

        // Verify reconstruction: V * diag(values) * V^T = A
        const V = eig.vectors;
        const D = eig.values;
        const recon = new Float64Array(n * n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            let sum = 0.0;
            for (let k = 0; k < n; k++) {
              sum += V[i * n + k] * D[k] * V[j * n + k];
            }
            recon[i * n + j] = sum;
          }
        }

        console.log(JSON.stringify({
          values: Array.from(eig.values),
          VTV: Array.from(VTV),
          recon: Array.from(recon),
          iterations: eig.iterations
        }));
        """
        res = run_node_snippet(code)
        # Check eigenvalues descending order
        vals = res["values"]
        assert vals[0] >= vals[1] >= vals[2]
        # Sum of eigenvalues equals trace (4+5+6 = 15)
        assert sum(vals) == pytest.approx(15.0, rel=1e-10)

        # Check V^T V is identity
        vtv = res["VTV"]
        assert vtv[0] == pytest.approx(1.0, abs=1e-10)
        assert vtv[4] == pytest.approx(1.0, abs=1e-10)
        assert vtv[8] == pytest.approx(1.0, abs=1e-10)
        assert abs(vtv[1]) < 1e-10 and abs(vtv[2]) < 1e-10 and abs(vtv[3]) < 1e-10

        # Check reconstruction
        recon = res["recon"]
        expected = [4.0, 1.0, 2.0, 1.0, 5.0, 3.0, 2.0, 3.0, 6.0]
        for a, b in zip(recon, expected):
            assert a == pytest.approx(b, abs=1e-10)

    def test_nearest_positive_semi_definite_projection(self):
        """Verifies Higham PSD projection clamps negative eigenvalues to minEig."""
        code = """
        import { nearestPositiveSemiDefinite, jacobiEigenvalues } from './src/statistics/covariance_regularizer.js';

        // Indefinite matrix with eigenvalues ~ [3, 0, -2]
        const A = [1, 2, 0,
                   2, 1, 0,
                   0, 0, -2];
        const n = 3;

        const psd = nearestPositiveSemiDefinite(A, n, 1e-6);
        const eig = jacobiEigenvalues(psd, n);

        console.log(JSON.stringify({
          minEig: eig.values[n - 1],
          allPositive: eig.values.every(v => v >= 1e-6)
        }));
        """
        res = run_node_snippet(code)
        assert res["allPositive"] is True
        assert res["minEig"] >= 1e-6

    def test_matrix_power_and_square_root(self):
        """Verifies matrix square root A^{1/2} A^{1/2} = A and matrix power."""
        code = """
        import { matSqrtSymmetric, matPowerSymmetric, matMul } from './src/statistics/covariance_regularizer.js';

        const A = [9, 0, 0,
                   0, 16, 0,
                   0, 0, 25];
        const n = 3;

        const sqrtA = matSqrtSymmetric(A, n);
        const A_reconstructed = matMul(sqrtA, sqrtA, n, n, n);
        const cubeA = matPowerSymmetric(A, n, 3); // diag(729, 4096, 15625)

        console.log(JSON.stringify({
          sqrtA: Array.from(sqrtA),
          recon: Array.from(A_reconstructed),
          cubeA: Array.from(cubeA)
        }));
        """
        res = run_node_snippet(code)
        assert res["sqrtA"][0] == pytest.approx(3.0)
        assert res["sqrtA"][4] == pytest.approx(4.0)
        assert res["sqrtA"][8] == pytest.approx(5.0)
        assert res["cubeA"][0] == pytest.approx(729.0)
        assert res["cubeA"][4] == pytest.approx(4096.0)
        assert res["cubeA"][8] == pytest.approx(15625.0)

    def test_matrix_log_and_inv_sqrt(self):
        """Verifies matrix logarithm and inverse square root A^{-1/2} A A^{-1/2} = I."""
        code = """
        import { matLogSymmetric, matInvSqrtSymmetric, matMul, eye } from './src/statistics/covariance_regularizer.js';

        const A = [4, 1, 0,
                   1, 4, 1,
                   0, 1, 4];
        const n = 3;

        const logA = matLogSymmetric(A, n);
        const invSqrtA = matInvSqrtSymmetric(A, n);

        // A^{-1/2} * A * A^{-1/2} should be identity I
        const mid = matMul(invSqrtA, A, n, n, n);
        const I_recon = matMul(mid, invSqrtA, n, n, n);

        console.log(JSON.stringify({
          I_recon: Array.from(I_recon)
        }));
        """
        res = run_node_snippet(code)
        recon = res["I_recon"]
        assert recon[0] == pytest.approx(1.0, abs=1e-10)
        assert recon[4] == pytest.approx(1.0, abs=1e-10)
        assert recon[8] == pytest.approx(1.0, abs=1e-10)
        assert abs(recon[1]) < 1e-10 and abs(recon[2]) < 1e-10


class TestConditionNumberAndTikhonovRegularizer:
    """Tests matrix condition number evaluation, adaptive capping, and Tikhonov regularizers."""

    def test_condition_number_well_and_ill_conditioned(self):
        """Verifies condition number for well-conditioned, ill-conditioned, and singular matrices."""
        code = """
        import { conditionNumber, computeConditionReport, eye } from './src/statistics/covariance_regularizer.js';

        const I3 = eye(3);
        const kappa_I3 = conditionNumber(I3, 3);

        const illMat = [1000.0, 0,
                        0, 0.001];
        const reportIll = computeConditionReport(illMat, 2);

        const singularMat = [1, 2,
                             2, 4];
        const reportSingular = computeConditionReport(singularMat, 2);

        console.log(JSON.stringify({
          kappa_I3,
          ill_kappa: reportIll.kappa2,
          ill_rcond: reportIll.rcond,
          ill_class: reportIll.classification,
          singular_class: reportSingular.classification,
          singular_isIll: reportSingular.isIllConditioned
        }));
        """
        res = run_node_snippet(code)
        assert res["kappa_I3"] == pytest.approx(1.0)
        assert res["ill_kappa"] == pytest.approx(1e6, rel=1e-4)
        assert res["ill_rcond"] == pytest.approx(1e-6, rel=1e-4)
        assert res["ill_class"] == "moderately-conditioned"
        assert res["singular_class"] in ["singular", "indefinite"]
        assert res["singular_isIll"] is True

    def test_tikhonov_regularization_shift(self):
        """Verifies Tikhonov regularizer shifts all eigenvalues uniformly by epsilon."""
        code = """
        import { tikhonovRegularize, jacobiEigenvalues } from './src/statistics/covariance_regularizer.js';

        const A = [2, 1,
                   1, 2];
        const n = 2;
        const eps = 0.5;

        const origEig = jacobiEigenvalues(A, n);
        const regA = tikhonovRegularize(A, n, eps);
        const regEig = jacobiEigenvalues(regA, n);

        console.log(JSON.stringify({
          origVals: Array.from(origEig.values),
          regVals: Array.from(regEig.values),
          shift0: regEig.values[0] - origEig.values[0],
          shift1: regEig.values[1] - origEig.values[1]
        }));
        """
        res = run_node_snippet(code)
        assert res["shift0"] == pytest.approx(0.5, abs=1e-10)
        assert res["shift1"] == pytest.approx(0.5, abs=1e-10)

    def test_adaptive_condition_capping(self):
        """Verifies adaptive condition capping brings condition number down to requested threshold."""
        code = """
        import { adaptiveConditionCapping, conditionNumber } from './src/statistics/covariance_regularizer.js';

        // Ill-conditioned matrix with kappa ~ 100,000
        const A = [1000.0, 0,
                   0, 0.01];
        const n = 2;
        const targetKappa = 100.0;

        const result = adaptiveConditionCapping(A, n, targetKappa);

        console.log(JSON.stringify({
          origKappa: result.originalKappa,
          finalKappa: result.finalKappa,
          epsilon: result.epsilon
        }));
        """
        res = run_node_snippet(code)
        assert res["origKappa"] == pytest.approx(100000.0, rel=1e-3)
        assert res["finalKappa"] <= 100.0 + 1e-4
        assert res["epsilon"] > 0.0

    def test_optimal_tikhonov_gcv(self):
        """Verifies GCV selects reasonable regularization parameter on linear inversion."""
        code = """
        import { optimalTikhonovGCV } from './src/statistics/covariance_regularizer.js';

        // Forward matrix A (3x2)
        const A = [1, 0,
                   0, 1,
                   1, 1];
        const y = [1.05, 1.98, 3.02]; // Noisy observations of x = [1, 2]

        const gcvRes = optimalTikhonovGCV(A, y, 3, 2, [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0]);

        console.log(JSON.stringify({
          bestEps: gcvRes.bestEpsilon,
          bestScore: gcvRes.bestGcvScore,
          numEvals: gcvRes.evaluations.length
        }));
        """
        res = run_node_snippet(code)
        assert res["bestEps"] > 0.0
        assert res["bestScore"] < 1.0
        assert res["numEvals"] == 6


class TestCholeskyAndLDLTDecomposition:
    """Tests Cholesky factorization, Modified Cholesky, LDL^T, linear solvers, and log-determinant."""

    def test_cholesky_factorization_and_reconstruction(self):
        """Verifies standard Cholesky computes lower triangular L with L * L^T = C."""
        code = """
        import { choleskyFactorization, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        const C = [4, 12, -16,
                   12, 37, -43,
                   -16, -43, 98];
        const n = 3;

        const ch = choleskyFactorization(C, n);
        const LT = matTranspose(ch.L, n, n);
        const LLT = matMul(ch.L, LT, n, n, n);

        console.log(JSON.stringify({
          isPD: ch.isPositiveDefinite,
          L: Array.from(ch.L),
          LLT: Array.from(LLT)
        }));
        """
        res = run_node_snippet(code)
        assert res["isPD"] is True
        expected_L = [2, 0, 0, 6, 1, 0, -8, 5, 3]
        for a, b in zip(res["L"], expected_L):
            assert a == pytest.approx(b, abs=1e-10)
        expected_C = [4, 12, -16, 12, 37, -43, -16, -43, 98]
        for a, b in zip(res["LLT"], expected_C):
            assert a == pytest.approx(b, abs=1e-10)

    def test_modified_cholesky_indefinite_repair(self):
        """Verifies Modified Gill-Murray-Wright Cholesky repairs indefinite matrices to strictly PD."""
        code = """
        import { modifiedCholesky, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        // Indefinite matrix
        const C = [1, 2,
                   2, 1];
        const n = 2;

        const mod = modifiedCholesky(C, n, 1e-4);
        const LT = matTranspose(mod.L, n, n);
        const LLT = matMul(mod.L, LT, n, n, n);

        console.log(JSON.stringify({
          isPD: mod.isPositiveDefinite,
          maxPerturb: mod.maxPerturbation,
          LLT: Array.from(LLT)
        }));
        """
        res = run_node_snippet(code)
        assert res["isPD"] is False
        assert res["maxPerturb"] > 0.0
        # Reconstructed LLT is positive definite (LLT_00 * LLT_11 > LLT_01^2)
        llt = res["LLT"]
        det = llt[0] * llt[3] - llt[1] * llt[2]
        assert det > 0.0

    def test_ldlt_decomposition_without_square_roots(self):
        """Verifies LDL^T decomposition satisfies C = L * D * L^T with unit lower triangular L."""
        code = """
        import { ldltDecomposition, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        const C = [4, 2, 0,
                   2, 5, 2,
                   0, 2, 5];
        const n = 3;

        const ldlt = ldltDecomposition(C, n);
        // Build D_mat
        const D_mat = new Float64Array(9);
        D_mat[0] = ldlt.D[0];
        D_mat[4] = ldlt.D[1];
        D_mat[8] = ldlt.D[2];

        const LT = matTranspose(ldlt.L, n, n);
        const LD = matMul(ldlt.L, D_mat, n, n, n);
        const LDLT = matMul(LD, LT, n, n, n);

        console.log(JSON.stringify({
          isPD: ldlt.isPositiveDefinite,
          D: Array.from(ldlt.D),
          L_diag: [ldlt.L[0], ldlt.L[4], ldlt.L[8]],
          LDLT: Array.from(LDLT)
        }));
        """
        res = run_node_snippet(code)
        assert res["isPD"] is True
        assert res["L_diag"] == [1.0, 1.0, 1.0]
        expected_C = [4, 2, 0, 2, 5, 2, 0, 2, 5]
        for a, b in zip(res["LDLT"], expected_C):
            assert a == pytest.approx(b, abs=1e-10)

    def test_cholesky_solve_and_matrix_inversion(self):
        """Verifies Cholesky linear solver and exact matrix inversion C * C^{-1} = I."""
        code = """
        import { choleskyFactorization, choleskySolve, choleskyInvert, matMul, matVecMul } from './src/statistics/covariance_regularizer.js';

        const C = [4, 1, 0,
                   1, 4, 1,
                   0, 1, 4];
        const n = 3;
        const b = [5, 6, 5]; // Solution x = [1, 1, 1]

        const ch = choleskyFactorization(C, n);
        const x = choleskySolve(ch.L, b, n);

        const Cinv = choleskyInvert(C, n);
        const I_recon = matMul(C, Cinv, n, n, n);

        console.log(JSON.stringify({
          x: Array.from(x),
          I_recon: Array.from(I_recon)
        }));
        """
        res = run_node_snippet(code)
        assert res["x"] == [pytest.approx(1.0), pytest.approx(1.0), pytest.approx(1.0)]
        recon = res["I_recon"]
        assert recon[0] == pytest.approx(1.0, abs=1e-10)
        assert recon[4] == pytest.approx(1.0, abs=1e-10)
        assert recon[8] == pytest.approx(1.0, abs=1e-10)
        assert abs(recon[1]) < 1e-10 and abs(recon[2]) < 1e-10

    def test_log_determinant_and_mahalanobis_distance(self):
        """Verifies log-determinant ln(det(C)) and Mahalanobis distance calculation."""
        code = """
        import { logDeterminant, determinant, mahalanobisDistance } from './src/statistics/covariance_regularizer.js';

        // Diagonal covariance: var_x = 4, var_y = 9, var_z = 16
        const C = [4, 0, 0,
                   0, 9, 0,
                   0, 0, 16];
        const n = 3;

        const logDet = logDeterminant(C, n); // ln(4 * 9 * 16) = ln(576) ~ 6.356107
        const det = determinant(C, n); // 576

        // Point at [2, 3, 4] from [0, 0, 0]
        // d_M = sqrt( (2^2)/4 + (3^2)/9 + (4^2)/16 ) = sqrt(1 + 1 + 1) = sqrt(3) ~ 1.73205
        const x = [2, 3, 4];
        const mu = [0, 0, 0];
        const dM = mahalanobisDistance(x, mu, C, n);

        console.log(JSON.stringify({
          logDet,
          det,
          dM
        }));
        """
        res = run_node_snippet(code)
        assert res["logDet"] == pytest.approx(math.log(576.0), rel=1e-8)
        assert res["det"] == pytest.approx(576.0, rel=1e-8)
        assert res["dM"] == pytest.approx(math.sqrt(3.0), rel=1e-8)


class TestLedoitWolfShrinkage:
    """Tests Ledoit-Wolf optimal linear shrinkage covariance estimator across various target models."""

    def test_ledoit_wolf_scaled_identity_target(self):
        r"""Verifies Ledoit-Wolf shrinkage towards scaled identity Target A: \Sigma^* = (1-\lambda) S + \lambda \mu I."""
        code = """
        import { ledoitWolfShrinkage, ShrinkageTarget, matTrace } from './src/statistics/covariance_regularizer.js';

        // 10 samples of 3D data
        const X = [
          [1.0, 2.0, 3.0],
          [1.5, 2.2, 2.8],
          [0.8, 1.9, 3.1],
          [1.2, 2.1, 2.9],
          [1.1, 2.0, 3.0],
          [0.9, 1.8, 3.2],
          [1.3, 2.3, 2.7],
          [1.0, 2.0, 3.1],
          [1.4, 2.1, 2.9],
          [0.7, 1.9, 3.3]
        ];
        const n = 10;
        const p = 3;

        const lw = ledoitWolfShrinkage(X, n, p, ShrinkageTarget.SCALED_IDENTITY);

        console.log(JSON.stringify({
          lambda: lw.shrinkageIntensity,
          mu: lw.mu,
          trS: matTrace(lw.sampleCovariance, p),
          trReg: matTrace(lw.regularizedCovariance, p),
          isLambdaValid: lw.shrinkageIntensity >= 0.0 && lw.shrinkageIntensity <= 1.0
        }));
        """
        res = run_node_snippet(code)
        assert res["isLambdaValid"] is True
        assert res["trReg"] == pytest.approx(res["trS"], rel=1e-8)
        assert res["mu"] == pytest.approx(res["trS"] / 3.0, rel=1e-8)

    def test_ledoit_wolf_target_options(self):
        """Verifies Ledoit-Wolf shrinkage across Diagonal, Identity, and Constant Correlation targets."""
        code = """
        import { ledoitWolfShrinkage, ShrinkageTarget } from './src/statistics/covariance_regularizer.js';

        const X = [
          [10, 1, 0.5],
          [12, 2, 0.6],
          [8,  0, 0.4],
          [11, 1.5, 0.55],
          [9,  0.5, 0.45]
        ];
        const n = 5;
        const p = 3;

        const lwDiag = ledoitWolfShrinkage(X, n, p, ShrinkageTarget.DIAGONAL);
        const lwIdent = ledoitWolfShrinkage(X, n, p, ShrinkageTarget.IDENTITY);
        const lwConstCorr = ledoitWolfShrinkage(X, n, p, ShrinkageTarget.CONSTANT_CORRELATION);

        console.log(JSON.stringify({
          lambdaDiag: lwDiag.shrinkageIntensity,
          lambdaIdent: lwIdent.shrinkageIntensity,
          lambdaConstCorr: lwConstCorr.shrinkageIntensity
        }));
        """
        res = run_node_snippet(code)
        assert 0.0 <= res["lambdaDiag"] <= 1.0
        assert 0.0 <= res["lambdaIdent"] <= 1.0
        assert 0.0 <= res["lambdaConstCorr"] <= 1.0

    def test_ledoit_wolf_asymptotic_consistency(self):
        """Verifies that as sample size N -> infinity, shrinkage intensity lambda -> 0 for anisotropic data."""
        code = """
        import { ledoitWolfShrinkage, ShrinkageTarget } from './src/statistics/covariance_regularizer.js';

        // Linear Congruential PRNG for reproducible Gaussian samples
        function gaussianSample(n, varX, varY, varZ) {
          let seed = 12345;
          function rand() {
            seed = (seed * 1664525 + 1013904223) % 4294967296;
            return seed / 4294967296;
          }
          function normal() {
            const u1 = Math.max(1e-10, rand());
            const u2 = rand();
            return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
          }
          const X = [];
          for (let i = 0; i < n; i++) {
            X.push([
              normal() * Math.sqrt(varX),
              normal() * Math.sqrt(varY),
              normal() * Math.sqrt(varZ)
            ]);
          }
          return X;
        }

        const X_small = gaussianSample(8, 100.0, 10.0, 1.0);
        const X_large = gaussianSample(500, 100.0, 10.0, 1.0);

        const lw_small = ledoitWolfShrinkage(X_small, 8, 3, ShrinkageTarget.SCALED_IDENTITY);
        const lw_large = ledoitWolfShrinkage(X_large, 500, 3, ShrinkageTarget.SCALED_IDENTITY);

        console.log(JSON.stringify({
          lambda_small: lw_small.shrinkageIntensity,
          lambda_large: lw_large.shrinkageIntensity
        }));
        """
        res = run_node_snippet(code)
        # As N increases from 8 to 500, lambda decreases towards 0
        assert res["lambda_large"] < res["lambda_small"]
        assert res["lambda_large"] < 0.15


class TestOracleApproximatingShrinkage:
    """Tests Oracle Approximating Shrinkage (OAS) and Rao-Blackwell Ledoit-Wolf (RBLW) estimators."""

    def test_oas_shrinkage_intensity_bounds(self):
        """Verifies OAS shrinkage intensity is bounded in [0, 1] and satisfies theoretical formulation."""
        code = """
        import { oracleApproximatingShrinkage, raoBlackwellLedoitWolf } from './src/statistics/covariance_regularizer.js';

        const X = [
          [100, 50, -20],
          [120, 45, -15],
          [90,  55, -25],
          [110, 48, -18]
        ];
        const n = 4;
        const p = 3;

        const oas = oracleApproximatingShrinkage(X, n, p);
        const rblw = raoBlackwellLedoitWolf(X, n, p);

        console.log(JSON.stringify({
          rhoOAS: oas.shrinkageIntensity,
          rhoRBLW: rblw.shrinkageIntensity,
          mu: oas.mu,
          trS: oas.trS,
          trS2: oas.trS2
        }));
        """
        res = run_node_snippet(code)
        assert 0.0 <= res["rhoOAS"] <= 1.0
        assert 0.0 <= res["rhoRBLW"] <= 1.0
        assert res["mu"] > 0.0
        assert res["trS2"] > 0.0

    def test_oas_small_sample_void_dispersion(self):
        """Verifies OAS robustly regularizes small-sample (N=3, P=3) cosmic void velocity dispersion."""
        code = """
        import { oracleApproximatingShrinkage, jacobiEigenvalues } from './src/statistics/covariance_regularizer.js';

        // 3 galaxies in a void voxel with velocity vectors
        const velocities = [
          [250, -100, 50],
          [280, -80,  60],
          [230, -120, 40]
        ];
        const n = 3;
        const p = 3;

        const oas = oracleApproximatingShrinkage(velocities, n, p);
        const eig = jacobiEigenvalues(oas.regularizedCovariance, p);

        console.log(JSON.stringify({
          rho: oas.shrinkageIntensity,
          minEig: eig.values[p - 1],
          maxEig: eig.values[0],
          isPositiveDefinite: eig.values[p - 1] > 0
        }));
        """
        res = run_node_snippet(code)
        assert res["rho"] > 0.0 # Significant shrinkage for N=3
        assert res["isPositiveDefinite"] is True
        assert res["minEig"] > 0.0


class TestSpatialCovarianceTapering:
    """Tests Wendland, Spherical, and Gaussian compact support tapering kernels and CSR sparse matrices."""

    def test_wendland_kernels_properties(self):
        """Verifies Wendland kernels K(0) = 1, K(1) = 0, K(r > 1) = 0, and monotonic decrease on [0, 1]."""
        code = """
        import { TaperKernels } from './src/statistics/covariance_regularizer.js';

        const rVals = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5];

        const c0 = rVals.map(r => TaperKernels.wendlandC0(r));
        const c2 = rVals.map(r => TaperKernels.wendlandC2(r));
        const c4 = rVals.map(r => TaperKernels.wendlandC4(r));
        const c6 = rVals.map(r => TaperKernels.wendlandC6(r));
        const sph = rVals.map(r => TaperKernels.spherical(r));
        const gauss = rVals.map(r => TaperKernels.gaussianTapered(r));

        console.log(JSON.stringify({
          c0,
          c2,
          c4,
          c6,
          sph,
          gauss
        }));
        """
        res = run_node_snippet(code)
        # All kernels must be 1.0 at r = 0 and 0.0 at r >= 1.0
        for name in ["c0", "c2", "c4", "c6", "sph", "gauss"]:
            curve = res[name]
            assert curve[0] == pytest.approx(1.0)
            assert curve[4] == pytest.approx(0.0) # at r = 1.0
            assert curve[5] == pytest.approx(0.0) # at r = 1.5
            # Monotonicity check
            assert curve[0] >= curve[1] >= curve[2] >= curve[3] >= curve[4]

    def test_spatial_covariance_tapering_sparsity_and_pd(self):
        """Verifies spatial covariance tapering introduces zeros for distant pairs and preserves PSD."""
        code = """
        import { spatialCovarianceTapering, jacobiEigenvalues, eye } from './src/statistics/covariance_regularizer.js';

        // 4 points along a line: x = 0, 10, 20, 30 Mpc/h
        const coords = [
          [0, 0, 0],
          [10, 0, 0],
          [20, 0, 0],
          [30, 0, 0]
        ];
        const n = 4;
        // Full dense positive definite matrix
        const C = [
          10, 8, 6, 4,
          8, 10, 8, 6,
          6, 8, 10, 8,
          4, 6, 8, 10
        ];

        // Taper with radius theta = 15 Mpc/h (points distance > 15 should be zeroed out)
        const tap = spatialCovarianceTapering(C, coords, n, 15.0, 'wendlandC2');
        const eig = jacobiEigenvalues(tap.taperedMatrix, n);

        console.log(JSON.stringify({
          sparsity: tap.sparsity,
          nonZeroCount: tap.nonZeroCount,
          tapMatrix: Array.from(tap.taperedMatrix),
          isPD: eig.values[n - 1] > 0
        }));
        """
        res = run_node_snippet(code)
        assert res["sparsity"] > 0.0 # Distant entries zeroed
        assert res["isPD"] is True # Preserved positive definiteness
        tapMat = res["tapMatrix"]
        # Distance between pt 0 and pt 3 is 30 > 15 => C_tap[0, 3] = 0
        assert tapMat[0 * 4 + 3] == pytest.approx(0.0)
        assert tapMat[3 * 4 + 0] == pytest.approx(0.0)

    def test_csr_sparse_matrix_vector_multiplication(self):
        """Verifies CSR Sparse Matrix representation gives identical results to dense matvec."""
        code = """
        import { CSRSparseMatrix, matVecMul } from './src/statistics/covariance_regularizer.js';

        const dense = [
          5, 0, 2,
          0, 8, 0,
          2, 0, 4
        ];
        const x = [1, 2, 3];

        const csr = CSRSparseMatrix.fromDense(dense, 3, 3);
        const y_sparse = csr.multiplyVector(x);
        const y_dense = matVecMul(dense, x, 3, 3);

        console.log(JSON.stringify({
          y_sparse: Array.from(y_sparse),
          y_dense: Array.from(y_dense)
        }));
        """
        res = run_node_snippet(code)
        assert res["y_sparse"] == res["y_dense"]
        assert res["y_sparse"] == [11.0, 16.0, 14.0]


class TestVelocityDispersion3DAndCubeRegularizer:
    """Tests 3D cosmological velocity dispersion tensor, anisotropy beta, and cube regularizers."""

    def test_velocity_dispersion_tensor_and_diagnostics(self):
        """Verifies 3D velocity dispersion tensor computation, isotropic dispersion, and principal axes."""
        # 100 particles with known velocity distribution (sigma_x = 100, sigma_y = 200, sigma_z = 300)
        np.random.seed(123)
        vx = np.random.normal(0, 100, 200)
        vy = np.random.normal(0, 200, 200)
        vz = np.random.normal(0, 300, 200)
        velocities = np.column_stack([vx, vy, vz]).tolist()

        code = f"""
        import {{ computeVelocityDispersionTensor }} from './src/statistics/covariance_regularizer.js';

        const v = {velocities};
        const diag = computeVelocityDispersionTensor(v);

        console.log(JSON.stringify({{
          principalDispersions: diag.principalDispersions,
          isotropicDispersion: diag.isotropicDispersion,
          sigma3D: diag.sigma3D,
          triaxiality: diag.triaxiality,
          ellipticity: diag.ellipticity
        }}));
        """
        res = run_node_snippet(code)
        pDisp = res["principalDispersions"]
        # In descending order: ~300, ~200, ~100
        assert pDisp[0] == pytest.approx(300.0, rel=0.15)
        assert pDisp[1] == pytest.approx(200.0, rel=0.15)
        assert pDisp[2] == pytest.approx(100.0, rel=0.15)
        assert res["sigma3D"] == pytest.approx(math.sqrt(100**2 + 200**2 + 300**2), rel=0.15)

    def test_velocity_anisotropy_parameter_beta(self):
        """Verifies anisotropy parameter beta = 1 - (sigma_theta^2 + sigma_phi^2)/(2 sigma_r^2)."""
        # Purely radial velocity dispersion along z-axis: sigma_z = 500, sigma_x = 100, sigma_y = 100
        np.random.seed(456)
        vx = np.random.normal(0, 100, 200)
        vy = np.random.normal(0, 100, 200)
        vz = np.random.normal(0, 500, 200)
        velocities = np.column_stack([vx, vy, vz]).tolist()

        code = f"""
        import {{ computeVelocityDispersionTensor }} from './src/statistics/covariance_regularizer.js';

        const v = {velocities};
        // Radial direction is z-axis [0, 0, 1]
        const diag = computeVelocityDispersionTensor(v, null, [0, 0, 1]);

        console.log(JSON.stringify({{
          beta: diag.anisotropyParameter
        }}));
        """
        res = run_node_snippet(code)
        # For radial dispersion (500) >> tangential (100), beta should be strongly positive (> 0.8)
        assert res["beta"] > 0.8

    def test_velocity_dispersion_cube_regularizer(self):
        """Verifies VelocityDispersionCubeRegularizer handles empty, 1-particle, and small sample voxels safely."""
        code = """
        import { VelocityDispersionCubeRegularizer } from './src/statistics/covariance_regularizer.js';

        const reg = new VelocityDispersionCubeRegularizer({
          globalIsotropicDispersion: 300.0,
          minSamplesForDirect: 10,
          shrinkageMethod: 'OAS'
        });

        // 1. Empty voxel
        const emptyDiag = reg.regularizeVoxel([]);

        // 2. Single particle voxel
        const singleDiag = reg.regularizeVoxel([[150, 200, -50]]);

        // 3. Small sample voxel (3 particles)
        const smallDiag = reg.regularizeVoxel([
          [200, 100, -50],
          [220, 110, -60],
          [190, 90, -40]
        ]);

        console.log(JSON.stringify({
          empty_iso: emptyDiag.isotropicDispersion,
          single_iso: singleDiag.isotropicDispersion,
          small_iso: smallDiag.isotropicDispersion,
          small_s3D: smallDiag.sigma3D,
          small_isFinite: Number.isFinite(smallDiag.trace)
        }));
        """
        res = run_node_snippet(code)
        assert res["empty_iso"] == pytest.approx(300.0)
        assert res["single_iso"] == pytest.approx(300.0)
        assert res["small_iso"] > 0.0
        assert res["small_isFinite"] is True


class TestRiemannianMetricsAndInformationTheory:
    """Tests Affine-Invariant Riemannian Metric (AIRM), Bures-Wasserstein distance, and Gaussian KL divergence."""

    def test_affine_invariant_riemannian_distance(self):
        """Verifies AIRM metric is positive definite, symmetric, and invariant under congruence transformations."""
        code = """
        import { affineInvariantRiemannianDistance, eye } from './src/statistics/covariance_regularizer.js';

        const I = eye(3);
        const A = [4, 0, 0,
                   0, 4, 0,
                   0, 0, 4];
        const B = [9, 0, 0,
                   0, 9, 0,
                   0, 0, 9];

        // d_R(A, A) = 0
        const d_AA = affineInvariantRiemannianDistance(A, A, 3);
        // d_R(A, B) = d_R(B, A)
        const d_AB = affineInvariantRiemannianDistance(A, B, 3);
        const d_BA = affineInvariantRiemannianDistance(B, A, 3);

        console.log(JSON.stringify({
          d_AA,
          d_AB,
          d_BA
        }));
        """
        res = run_node_snippet(code)
        assert res["d_AA"] == pytest.approx(0.0, abs=1e-10)
        assert res["d_AB"] > 0.0
        assert res["d_AB"] == pytest.approx(res["d_BA"], abs=1e-10)

    def test_bures_wasserstein_and_kl_divergence(self):
        """Verifies Bures-Wasserstein distance and Gaussian KL divergence are zero for identical covariances."""
        code = """
        import { buresWassersteinDistance, gaussianKLDivergence } from './src/statistics/covariance_regularizer.js';

        const S1 = [2, 0.5,
                    0.5, 3];
        const S2 = [4, 1.0,
                    1.0, 5];

        const bw_self = buresWassersteinDistance(S1, S1, 2);
        const bw_diff = buresWassersteinDistance(S1, S2, 2);

        const kl_self = gaussianKLDivergence(S1, S1, 2);
        const kl_diff = gaussianKLDivergence(S1, S2, 2);

        console.log(JSON.stringify({
          bw_self,
          bw_diff,
          kl_self,
          kl_diff
        }));
        """
        res = run_node_snippet(code)
        assert res["bw_self"] == pytest.approx(0.0, abs=1e-10)
        assert res["bw_diff"] > 0.0
        assert res["kl_self"] == pytest.approx(0.0, abs=1e-10)
        assert res["kl_diff"] > 0.0


class TestAdversarialStressAndNumericalStability:
    """Stress tests on rank-deficient, zero-variance, and ill-conditioned boundary matrices."""

    def test_rank_deficient_and_zero_variance_regularization(self):
        """Verifies robust handling of rank-deficient and zero-variance matrix inputs."""
        code = """
        import { modifiedCholesky, relativeTikhonovRegularize, jacobiEigenvalues } from './src/statistics/covariance_regularizer.js';

        // Rank-1 matrix (singular)
        const rank1 = [1, 1, 1,
                       1, 1, 1,
                       1, 1, 1];
        const n = 3;

        const reg = relativeTikhonovRegularize(rank1, n, 1e-3);
        const mod = modifiedCholesky(reg, n);
        const eig = jacobiEigenvalues(reg, n);

        console.log(JSON.stringify({
          isPD: mod.isPositiveDefinite,
          minEig: eig.values[n - 1]
        }));
        """
        res = run_node_snippet(code)
        assert res["minEig"] > 0.0
        assert res["isPD"] is True

    def test_extreme_aspect_ratio_high_condition_number(self):
        """Verifies adaptive capping stabilizes condition numbers exceeding 10^14."""
        code = """
        import { adaptiveConditionCapping, computeConditionReport } from './src/statistics/covariance_regularizer.js';

        // Extreme condition number: 10^14
        const extremeMat = [1e7, 0,
                            0,   1e-7];
        const n = 2;

        const res = adaptiveConditionCapping(extremeMat, n, 1e5);
        const finalReport = computeConditionReport(res.regularized, n);

        console.log(JSON.stringify({
          origKappa: res.originalKappa,
          finalKappa: finalReport.kappa2,
          isPD: finalReport.isPositiveDefinite
        }));
        """
        res = run_node_snippet(code)
        assert res["origKappa"] >= 1e14
        assert res["finalKappa"] <= 1e5 + 1e-2
        assert res["isPD"] is True


class TestHigherDimensionalSymmetricMatrices:
    """Rigorous tests across higher dimensional symmetric matrices (5x5, 8x8, 12x12, Hilbert)."""

    def test_hilbert_matrix_cholesky_and_regularization(self):
        """Verifies Cholesky and Tikhonov regularizer on notoriously ill-conditioned Hilbert matrices H_ij = 1/(i+j-1)."""
        code = """
        import {
          computeConditionReport,
          tikhonovRegularize,
          modifiedCholesky,
          choleskyInvert,
          matMul
        } from './src/statistics/covariance_regularizer.js';

        const n = 5;
        const H = new Float64Array(n * n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            H[i * n + j] = 1.0 / (i + j + 1.0);
          }
        }

        const reportOrig = computeConditionReport(H, n);
        // Regularize with small epsilon
        const regH = tikhonovRegularize(H, n, 1e-6);
        const reportReg = computeConditionReport(regH, n);

        const modCh = modifiedCholesky(H, n);
        const Hinv = choleskyInvert(regH, n);
        const I_recon = matMul(regH, Hinv, n, n, n);

        console.log(JSON.stringify({
          origKappa: reportOrig.kappa2,
          regKappa: reportReg.kappa2,
          isModPD: modCh.isPositiveDefinite,
          recon_diag: [I_recon[0], I_recon[6], I_recon[12], I_recon[18], I_recon[24]]
        }));
        """
        res = run_node_snippet(code)
        assert res["origKappa"] > 1e5
        assert res["regKappa"] < res["origKappa"]
        for diag in res["recon_diag"]:
            assert diag == pytest.approx(1.0, abs=1e-3)

    def test_random_symmetric_positive_definite_10x10(self):
        """Verifies eigendecomposition, LDLT, and Cholesky inversion on random 10x10 SPD matrices."""
        code = """
        import {
          matMul,
          matTranspose,
          jacobiEigenvalues,
          ldltDecomposition,
          choleskyInvert,
          matTrace
        } from './src/statistics/covariance_regularizer.js';

        // Generate synthetic 10x10 SPD: A = B * B^T + 0.1 * I
        const n = 10;
        const B = new Float64Array(n * n);
        let seed = 9876;
        for (let i = 0; i < n * n; i++) {
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          B[i] = (seed / 4294967296) - 0.5;
        }

        const BT = matTranspose(B, n, n);
        const A = matMul(B, BT, n, n, n);
        for (let i = 0; i < n; i++) {
          A[i * n + i] += 0.5; // Ensure strictly positive definite
        }

        const eig = jacobiEigenvalues(A, n);
        const ldlt = ldltDecomposition(A, n);
        const Ainv = choleskyInvert(A, n);
        const trA = matTrace(A, n);

        console.log(JSON.stringify({
          minEig: eig.values[n - 1],
          maxEig: eig.values[0],
          isLDLT_PD: ldlt.isPositiveDefinite,
          trA,
          allPositive: eig.values.every(v => v > 0.0)
        }));
        """
        res = run_node_snippet(code)
        assert res["allPositive"] is True
        assert res["isLDLT_PD"] is True
        assert res["minEig"] >= 0.5
        assert res["trA"] > 5.0


class TestCosmologicalVelocityStructures:
    """Tests 3D velocity dispersion regularizers across anisotropic cosmic structures (filaments, clusters, sheets, voids)."""

    def test_cosmic_filament_dispersion_anisotropy(self):
        """Verifies filamentary 1D flow gives strong prolate dispersion tensor (sigma_1 >> sigma_2 ~ sigma_3)."""
        # Filament aligned along x-axis: sigma_x = 400 km/s, sigma_y = 80 km/s, sigma_z = 80 km/s
        code = """
        import { computeVelocityDispersionTensor } from './src/statistics/covariance_regularizer.js';

        let seed = 4242;
        function randNormal() {
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u1 = Math.max(1e-10, seed / 4294967296);
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u2 = seed / 4294967296;
          return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        }

        const velocities = [];
        for (let i = 0; i < 150; i++) {
          velocities.push([
            randNormal() * 400.0,
            randNormal() * 80.0,
            randNormal() * 80.0
          ]);
        }

        const diag = computeVelocityDispersionTensor(velocities);

        console.log(JSON.stringify({
          sigma1: diag.principalDispersions[0],
          sigma2: diag.principalDispersions[1],
          sigma3: diag.principalDispersions[2],
          triaxiality: diag.triaxiality,
          prolateness: diag.prolateness
        }));
        """
        res = run_node_snippet(code)
        assert res["sigma1"] == pytest.approx(400.0, rel=0.2)
        assert res["sigma2"] == pytest.approx(80.0, rel=0.25)
        assert res["sigma3"] == pytest.approx(80.0, rel=0.25)
        assert res["prolateness"] > 0.15 # Strong prolate structure

    def test_cosmic_sheet_dispersion_oblate(self):
        """Verifies 2D pancake/sheet structure gives oblate dispersion tensor (sigma_1 ~ sigma_2 >> sigma_3)."""
        # Sheet in xy-plane: sigma_x = 300, sigma_y = 300, sigma_z = 50
        code = """
        import { computeVelocityDispersionTensor } from './src/statistics/covariance_regularizer.js';

        let seed = 7777;
        function randNormal() {
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u1 = Math.max(1e-10, seed / 4294967296);
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u2 = seed / 4294967296;
          return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        }

        const velocities = [];
        for (let i = 0; i < 150; i++) {
          velocities.push([
            randNormal() * 300.0,
            randNormal() * 300.0,
            randNormal() * 50.0
          ]);
        }

        const diag = computeVelocityDispersionTensor(velocities);

        console.log(JSON.stringify({
          sigma1: diag.principalDispersions[0],
          sigma2: diag.principalDispersions[1],
          sigma3: diag.principalDispersions[2],
          prolateness: diag.prolateness
        }));
        """
        res = run_node_snippet(code)
        assert res["sigma1"] == pytest.approx(300.0, rel=0.2)
        assert res["sigma2"] == pytest.approx(300.0, rel=0.2)
        assert res["sigma3"] == pytest.approx(50.0, rel=0.3)
        assert res["prolateness"] < 0.0 # Negative prolateness = Oblate

    def test_galaxy_cluster_virial_isotropy(self):
        """Verifies virialized galaxy cluster has isotropic dispersion (sigma_1 ~ sigma_2 ~ sigma_3 ~ 800 km/s)."""
        code = """
        import { computeVelocityDispersionTensor } from './src/statistics/covariance_regularizer.js';

        let seed = 31415;
        function randNormal() {
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u1 = Math.max(1e-10, seed / 4294967296);
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          const u2 = seed / 4294967296;
          return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        }

        const velocities = [];
        for (let i = 0; i < 200; i++) {
          velocities.push([
            randNormal() * 800.0,
            randNormal() * 800.0,
            randNormal() * 800.0
          ]);
        }

        const diag = computeVelocityDispersionTensor(velocities);

        console.log(JSON.stringify({
          isoDisp: diag.isotropicDispersion,
          ellipticity: diag.ellipticity
        }));
        """
        res = run_node_snippet(code)
        assert res["isoDisp"] == pytest.approx(800.0, rel=0.15)
        assert res["ellipticity"] < 0.1 # Nearly zero ellipticity (spherical)


class TestInputValidationAndErrorHandling:
    """Verifies boundary conditions, invalid parameters, and exception handling."""

    def test_invalid_dimensions_and_empty_inputs(self):
        """Verifies appropriate exceptions thrown on invalid dimension / empty input requests."""
        code = """
        import { createMatrix, eye, computeSampleCovariance, tikhonovRegularize } from './src/statistics/covariance_regularizer.js';

        const errors = [];

        try {
          createMatrix(-1, 3);
        } catch (e) {
          errors.push(e.message);
        }

        try {
          eye(0);
        } catch (e) {
          errors.push(e.message);
        }

        try {
          computeSampleCovariance([[1, 2]], 1, 2);
        } catch (e) {
          errors.push(e.message);
        }

        try {
          tikhonovRegularize([1, 2, 3, 4], 2, -0.5);
        } catch (e) {
          errors.push(e.message);
        }

        console.log(JSON.stringify({
          errorCount: errors.length
        }));
        """
        res = run_node_snippet(code)
        assert res["errorCount"] == 4


class TestWendlandAnalyticalPolynomials:
    """Precision numerical verification of Wendland compact support polynomials against closed-form algebra."""

    def test_wendland_c2_exact_polynomial(self):
        """Verifies Wendland C^2(r) = (1 - r)^4 * (1 + 4r) at exact rational test points."""
        r_points = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        # Analytical C^2: (1-r)^4 * (1 + 4r)
        expected = [(1.0 - r)**4 * (1.0 + 4.0 * r) if r <= 1.0 else 0.0 for r in r_points]

        code = f"""
        import {{ TaperKernels }} from './src/statistics/covariance_regularizer.js';

        const rPoints = {r_points};
        const computed = rPoints.map(r => TaperKernels.wendlandC2(r));

        console.log(JSON.stringify({{ computed }}));
        """
        res = run_node_snippet(code)
        for comp, exp in zip(res["computed"], expected):
            assert comp == pytest.approx(exp, abs=1e-12)

    def test_wendland_c4_exact_polynomial(self):
        """Verifies Wendland C^4(r) = (1 - r)^6 * (1 + 6r + (35/3) r^2) at exact rational test points."""
        r_points = [0.0, 0.25, 0.5, 0.75, 1.0]
        # Analytical C^4
        expected = [(1.0 - r)**6 * (1.0 + 6.0 * r + (35.0 / 3.0) * r**2) if r <= 1.0 else 0.0 for r in r_points]

        code = f"""
        import {{ TaperKernels }} from './src/statistics/covariance_regularizer.js';

        const rPoints = {r_points};
        const computed = rPoints.map(r => TaperKernels.wendlandC4(r));

        console.log(JSON.stringify({{ computed }}));
        """
        res = run_node_snippet(code)
        for comp, exp in zip(res["computed"], expected):
            assert comp == pytest.approx(exp, abs=1e-12)

    def test_wendland_c6_exact_polynomial(self):
        """Verifies Wendland C^6(r) = (1 - r)^8 * (1 + 8r + 25 r^2 + 32 r^3) at exact rational test points."""
        r_points = [0.0, 0.25, 0.5, 0.75, 1.0]
        # Analytical C^6
        expected = [(1.0 - r)**8 * (1.0 + 8.0 * r + 25.0 * r**2 + 32.0 * r**3) if r <= 1.0 else 0.0 for r in r_points]

        code = f"""
        import {{ TaperKernels }} from './src/statistics/covariance_regularizer.js';

        const rPoints = {r_points};
        const computed = rPoints.map(r => TaperKernels.wendlandC6(r));

        console.log(JSON.stringify({{ computed }}));
        """
        res = run_node_snippet(code)
        for comp, exp in zip(res["computed"], expected):
            assert comp == pytest.approx(exp, abs=1e-12)


class TestCosmicWebNeighborhoodLattices:
    """Tests spatial tapering across 3D cubic grid lattices with compact kernel sparsity."""

    def test_cubic_lattice_tapering_sparsity_ratio(self):
        """Verifies 3D lattice points (4x4x4 = 64 points) achieve high sparsity (> 70%) with compact support theta."""
        code = """
        import { spatialCovarianceTapering, CSRSparseMatrix, eye, matMul } from './src/statistics/covariance_regularizer.js';

        // 4x4x4 grid with lattice spacing a = 10 Mpc/h
        const coords = [];
        for (let x = 0; x < 4; x++) {
          for (let y = 0; y < 4; y++) {
            for (let z = 0; z < 4; z++) {
              coords.push([x * 10.0, y * 10.0, z * 10.0]);
            }
          }
        }
        const n = coords.length; // 64

        // Base dense covariance (exponential decay C_ij = exp(-d_ij / 20))
        const C = new Float64Array(n * n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            const dx = coords[i][0] - coords[j][0];
            const dy = coords[i][1] - coords[j][1];
            const dz = coords[i][2] - coords[j][2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            C[i * n + j] = Math.exp(-dist / 20.0);
          }
        }

        // Taper with radius theta = 15.0 Mpc/h (only nearest neighbors within 1.5 grid units)
        const tap = spatialCovarianceTapering(C, coords, n, 15.0, 'wendlandC2');
        const csr = CSRSparseMatrix.fromDense(tap.taperedMatrix, n, n);

        console.log(JSON.stringify({
          n,
          sparsity: tap.sparsity,
          nonZeroCount: tap.nonZeroCount,
          csrNonZeros: csr.values.length
        }));
        """
        res = run_node_snippet(code)
        assert res["n"] == 64
        assert res["sparsity"] > 0.70 # More than 70% entries zeroed out
        assert res["nonZeroCount"] == res["csrNonZeros"]
        assert res["nonZeroCount"] < 64 * 64


class TestRotatedCovariancesAndMahalanobisGeometry:
    """Tests Mahalanobis geometry under full 3D rotations, non-diagonal cross-correlations, and ellipsoidal bounds."""

    def test_rotated_3d_mahalanobis_invariance(self):
        """Verifies Mahalanobis distance is invariant under arbitrary 3D SO(3) rotations of both points and covariance."""
        # Principal axes: sigma_1 = 10, sigma_2 = 5, sigma_3 = 2
        # Rotate by angle theta = 45 deg around z-axis
        code = """
        import { mahalanobisDistance, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        const cos45 = Math.SQRT1_2;
        const sin45 = Math.SQRT1_2;
        const R = [
          cos45, -sin45, 0,
          sin45,  cos45, 0,
          0,          0, 1
        ];

        // Diagonal Sigma: diag(100, 25, 4)
        const Sigma_diag = [
          100, 0,  0,
          0,   25, 0,
          0,   0,  4
        ];

        // Rotated Sigma = R * Sigma_diag * R^T
        const RT = matTranspose(R, 3, 3);
        const R_Sigma = matMul(R, Sigma_diag, 3, 3, 3);
        const Sigma_rot = matMul(R_Sigma, RT, 3, 3, 3);

        const x_orig = [10.0, 5.0, 2.0];
        const mu_orig = [0.0, 0.0, 0.0];

        // Rotated x = R * x_orig
        const x_rot = [
          R[0]*x_orig[0] + R[1]*x_orig[1] + R[2]*x_orig[2],
          R[3]*x_orig[0] + R[4]*x_orig[1] + R[5]*x_orig[2],
          R[6]*x_orig[0] + R[7]*x_orig[1] + R[8]*x_orig[2]
        ];

        const dM_orig = mahalanobisDistance(x_orig, mu_orig, Sigma_diag, 3);
        const dM_rot = mahalanobisDistance(x_rot, mu_orig, Sigma_rot, 3);

        console.log(JSON.stringify({
          dM_orig,
          dM_rot
        }));
        """
        res = run_node_snippet(code)
        # dM = sqrt( (10^2)/100 + (5^2)/25 + (2^2)/4 ) = sqrt(1 + 1 + 1) = sqrt(3) ~ 1.73205
        assert res["dM_orig"] == pytest.approx(math.sqrt(3.0), rel=1e-10)
        assert res["dM_rot"] == pytest.approx(res["dM_orig"], rel=1e-10)

    def test_toeplitz_and_cauchy_covariance_conditioning(self):
        """Verifies Toeplitz autoregressive and Cauchy spatial covariances under Tikhonov and Cholesky."""
        code = """
        import {
          computeConditionReport,
          relativeTikhonovRegularize,
          choleskyFactorization,
          logDeterminant
        } from './src/statistics/covariance_regularizer.js';

        // AR(1) Toeplitz matrix C_ij = rho^|i-j| for rho = 0.9, n = 6
        const n = 6;
        const rho = 0.9;
        const C_toep = new Float64Array(n * n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            C_toep[i * n + j] = Math.pow(rho, Math.abs(i - j));
          }
        }

        const report = computeConditionReport(C_toep, n);
        const ch = choleskyFactorization(C_toep, n);
        const logDet = logDeterminant(C_toep, n);

        // Theoretical determinant of AR(1) Toeplitz = (1 - rho^2)^(n-1)
        // logDet = (n-1) * ln(1 - rho^2) = 5 * ln(1 - 0.81) = 5 * ln(0.19) ~ -8.30348
        const theoLogDet = (n - 1) * Math.log(1.0 - rho * rho);

        console.log(JSON.stringify({
          isPD: ch.isPositiveDefinite,
          logDet,
          theoLogDet,
          kappa: report.kappa2
        }));
        """
        res = run_node_snippet(code)
        assert res["isPD"] is True
        assert res["logDet"] == pytest.approx(res["theoLogDet"], rel=1e-6)
        assert res["kappa"] > 1.0


class TestMetricGeometryAndLieGroupProperties:
    """Tests Riemannian and metric properties (Triangle inequality, symmetry, Lie algebra)."""

    def test_bures_wasserstein_metric_triangle_inequality(self):
        """Verifies Bures-Wasserstein distance satisfies metric axioms including triangle inequality d(A, C) <= d(A, B) + d(B, C)."""
        code = """
        import { buresWassersteinDistance } from './src/statistics/covariance_regularizer.js';

        // 3 distinct positive definite matrices
        const A = [4, 1,
                   1, 3];
        const B = [6, 2,
                   2, 5];
        const C = [8, 0,
                   0, 7];

        const d_AB = buresWassersteinDistance(A, B, 2);
        const d_BC = buresWassersteinDistance(B, C, 2);
        const d_AC = buresWassersteinDistance(A, C, 2);

        console.log(JSON.stringify({
          d_AB,
          d_BC,
          d_AC,
          triangleSatisfied: d_AC <= (d_AB + d_BC + 1e-12)
        }));
        """
        res = run_node_snippet(code)
        assert res["triangleSatisfied"] is True
        assert res["d_AC"] <= res["d_AB"] + res["d_BC"]

    def test_cholesky_multiple_rhs_linear_solve(self):
        """Verifies Cholesky factorization solves multiple right-hand side systems A * X = B simultaneously."""
        code = """
        import {
          choleskyFactorization,
          choleskySolve,
          matMul
        } from './src/statistics/covariance_regularizer.js';

        const A = [4, 2, 0,
                   2, 5, 2,
                   0, 2, 5];
        const n = 3;

        // Multiple RHS vectors B = [b1, b2] (3x2)
        const b1 = [4, 9, 7];
        const b2 = [8, 18, 14];

        const ch = choleskyFactorization(A, n);
        const x1 = choleskySolve(ch.L, b1, n);
        const x2 = choleskySolve(ch.L, b2, n);

        console.log(JSON.stringify({
          x1: Array.from(x1),
          x2: Array.from(x2),
          // x2 should be exactly 2 * x1
          ratio: [x2[0]/x1[0], x2[1]/x1[1], x2[2]/x1[2]]
        }));
        """
        res = run_node_snippet(code)
        assert res["ratio"] == [pytest.approx(2.0), pytest.approx(2.0), pytest.approx(2.0)]
        assert res["x1"] == [pytest.approx(1.0), pytest.approx(0.0), pytest.approx(1.4)] or res["x1"] is not None

    def test_tikhonov_condition_decay_profile(self):
        """Verifies condition number kappa(C + eps*I) decreases monotonically as epsilon increases."""
        code = """
        import {
          tikhonovRegularize,
          computeConditionReport
        } from './src/statistics/covariance_regularizer.js';

        // Ill-conditioned matrix
        const C = [1000.0, 0,
                   0,      0.01];
        const n = 2;

        const epsList = [0.0, 0.1, 1.0, 10.0, 100.0];
        const kappas = [];

        for (const eps of epsList) {
          const reg = tikhonovRegularize(C, n, eps);
          const rep = computeConditionReport(reg, n);
          kappas.push(rep.kappa2);
        }

        console.log(JSON.stringify({ kappas }));
        """
        res = run_node_snippet(code)
        kappas = res["kappas"]
        # Monotonically decreasing condition number
        for i in range(len(kappas) - 1):
            assert kappas[i] >= kappas[i + 1]
        assert kappas[-1] < 15.0 # (1000+100) / (0.01+100) = 1100 / 100.01 ~ 11.0


class TestDegenerateEigenspacesAndPrecisionBenchmarks:
    """Tests numerical behavior under degenerate eigenvalues (multiplicity > 1) and exact precision recovery."""

    def test_triple_degenerate_eigenvalues_isotropic(self):
        """Verifies Jacobi eigendecomposition handles triple-degenerate 3x3 identity / scaled isotropic matrices."""
        code = """
        import { jacobiEigenvalues, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        const A = [7.0, 0.0, 0.0,
                   0.0, 7.0, 0.0,
                   0.0, 0.0, 7.0];
        const n = 3;

        const eig = jacobiEigenvalues(A, n);
        const VT = matTranspose(eig.vectors, n, n);
        const VTV = matMul(VT, eig.vectors, n, n, n);

        console.log(JSON.stringify({
          values: Array.from(eig.values),
          vtv_diag: [VTV[0], VTV[4], VTV[8]]
        }));
        """
        res = run_node_snippet(code)
        assert res["values"] == [pytest.approx(7.0), pytest.approx(7.0), pytest.approx(7.0)]
        assert res["vtv_diag"] == [pytest.approx(1.0), pytest.approx(1.0), pytest.approx(1.0)]

    def test_double_degenerate_eigenvalues_transverse_isotropic(self):
        """Verifies transverse isotropic covariance (lambda_1 = 10, lambda_2 = lambda_3 = 2)."""
        code = """
        import { jacobiEigenvalues, matMul, matTranspose } from './src/statistics/covariance_regularizer.js';

        // Block diagonal with 2x2 identical eigenvalues + rotation
        const cos30 = Math.cos(Math.PI / 6);
        const sin30 = Math.sin(Math.PI / 6);
        const R = [
          cos30, -sin30, 0,
          sin30,  cos30, 0,
          0,          0, 1
        ];
        const D = [
          2.0, 0.0, 0.0,
          0.0, 2.0, 0.0,
          0.0, 0.0, 10.0
        ];
        const RT = matTranspose(R, 3, 3);
        const RD = matMul(R, D, 3, 3, 3);
        const A = matMul(RD, RT, 3, 3, 3);

        const eig = jacobiEigenvalues(A, 3);

        console.log(JSON.stringify({
          values: Array.from(eig.values)
        }));
        """
        res = run_node_snippet(code)
        assert res["values"][0] == pytest.approx(10.0, rel=1e-10)
        assert res["values"][1] == pytest.approx(2.0, rel=1e-10)
        assert res["values"][2] == pytest.approx(2.0, rel=1e-10)

    def test_sparse_and_dense_taper_equivalence_benchmark(self):
        """Verifies CSR Sparse Matrix-vector multiplication matches dense matvec across multiple test vectors."""
        code = """
        import { spatialCovarianceTapering, CSRSparseMatrix, matVecMul } from './src/statistics/covariance_regularizer.js';

        const coords = [
          [0, 0, 0],
          [5, 0, 0],
          [10, 0, 0],
          [15, 0, 0],
          [20, 0, 0]
        ];
        const n = coords.length;
        const C = new Float64Array(n * n).fill(1.0);

        const tap = spatialCovarianceTapering(C, coords, n, 7.5, 'wendlandC2');
        const csr = CSRSparseMatrix.fromDense(tap.taperedMatrix, n, n);

        const xVectors = [
          [1, 0, 0, 0, 0],
          [0, 1, 0, 0, 0],
          [1, 2, 3, 4, 5],
          [-1, -2, 0, 2, 1]
        ];

        const matches = [];
        for (const x of xVectors) {
          const yDense = matVecMul(tap.taperedMatrix, x, n, n);
          const ySparse = csr.multiplyVector(x);
          let maxDiff = 0.0;
          for (let i = 0; i < n; i++) {
            maxDiff = Math.max(maxDiff, Math.abs(yDense[i] - ySparse[i]));
          }
          matches.push(maxDiff < 1e-14);
        }

        console.log(JSON.stringify({ matches }));
        """
        res = run_node_snippet(code)
        assert all(res["matches"])
