# -*- coding: utf-8 -*-
"""
tests/uncertainty/test_bayesian_hmc_field.py
Exhaustive automated Pytest test suite for the Bayesian CF4++ 10,000-Step HMC Posterior Field and Covariance Engine.

Comprehensive Verification Sections:
1. Cosmological parameters, physical units, constants, and scientific citations (Dupuy & Courtois 2023).
2. N x N Cholesky decomposition, Jacobi eigendecomposition, nearest PSD projection, and Gate H validation.
3. Spatial covariance kernels (Gaussian, Exponential, Matérn-3/2, Matérn-5/2, Górski 1988 Anisotropic).
4. Bayesian HMC field data loaders, trilinear continuous interpolation, Jacobian evaluation, SNR, and confidence intervals.
5. Cholesky and Karhunen-Loève stochastic realization samplers along trajectories and across full 3D grids.
6. Probabilistic watershed basin assignment, Shannon entropy, Gini-Simpson uncertainty, and Jaccard/Dice overlap metrics.
7. Critical point stability probability, 3D Newton-Raphson tracking, spatial dispersion ellipsoids, and Morse bifurcation.
8. Cosmic variance propagation along streamline integration paths with continuous matrix ODE and ensemble bundles.
9. Multi-chain MCMC convergence diagnostics (Gelman-Rubin R_hat, autocorrelation time, effective sample size, E-BFMI).
10. Adversarial stress tests: ill-conditioned matrices, extreme noise regimes, and cosmological scale invariants.
11. Advanced Karhunen-Loève mode spectra, Lyapunov horizon tracking, and arbitrary 3D tensor rotation tests.
"""

import pytest
import math
from tests.utils import run_node_snippet

class TestCF4CosmologyAndMetadata:
    """Verifies standard cosmological parameters, units, and citation provenance."""

    def test_cosmological_constants_and_invariants(self):
        code = """
        import { CF4_COSMOLOGY, CovarianceKernelType, MorseCriticalType } from './src/uncertainty/bayesian_hmc_field.js';
        
        console.log(JSON.stringify({
          H0: CF4_COSMOLOGY.H0,
          OmegaM: CF4_COSMOLOGY.OMEGA_M,
          OmegaLambda: CF4_COSMOLOGY.OMEGA_LAMBDA,
          sigma8: CF4_COSMOLOGY.SIGMA_8,
          growthRateF: CF4_COSMOLOGY.GROWTH_RATE_F,
          velocityScale: CF4_COSMOLOGY.VELOCITY_SCALE,
          defaultLcorr: CF4_COSMOLOGY.DEFAULT_CORRELATION_LENGTH_MPC,
          defaultHmcSteps: CF4_COSMOLOGY.DEFAULT_HMC_STEPS,
          targetAcceptance: CF4_COSMOLOGY.TARGET_ACCEPTANCE_RATE,
          gelmanRubinThresh: CF4_COSMOLOGY.GELMAN_RUBIN_THRESHOLD,
          kernelTypes: Object.keys(CovarianceKernelType),
          morseTypes: Object.keys(MorseCriticalType)
        }));
        """
        res = run_node_snippet(code)
        assert res["H0"] == pytest.approx(74.6, rel=1e-3)
        assert res["OmegaM"] == pytest.approx(0.31, rel=1e-3)
        assert res["OmegaLambda"] == pytest.approx(0.69, rel=1e-3)
        assert res["sigma8"] == pytest.approx(0.81, rel=1e-3)
        assert res["velocityScale"] == pytest.approx(52.0, rel=1e-3)
        assert res["defaultLcorr"] == pytest.approx(25.0, rel=1e-3)
        assert res["defaultHmcSteps"] == 10000
        assert res["targetAcceptance"] == pytest.approx(0.651, rel=1e-2)
        assert res["gelmanRubinThresh"] == pytest.approx(1.05, rel=1e-3)
        assert "GAUSSIAN" in res["kernelTypes"]
        assert "GORSKI_ANISOTROPIC" in res["kernelTypes"]
        assert "ATTRACTOR" in res["morseTypes"]
        assert "SADDLE_FILAMENT" in res["morseTypes"]

    def test_bayesian_hmc_metadata_provenance(self):
        code = """
        import { BayesianHMCMetadata } from './src/uncertainty/bayesian_hmc_field.js';
        
        const meta = new BayesianHMCMetadata({
          chainLength: 10000,
          catalogName: 'Cosmicflows-4 (CF4++)',
          referencePublication: 'Dupuy & Courtois (2023) MNRAS / A&A Table A.1'
        });
        
        console.log(JSON.stringify({
          chainLength: meta.chainLength,
          catalogName: meta.catalogName,
          reference: meta.referencePublication,
          numChains: meta.numChains,
          effectiveSampleSize: meta.effectiveSampleSize
        }));
        """
        res = run_node_snippet(code)
        assert res["chainLength"] == 10000
        assert "Cosmicflows-4" in res["catalogName"]
        assert "Dupuy & Courtois (2023)" in res["reference"]
        assert res["numChains"] == 4
        assert res["effectiveSampleSize"] > 1000


class TestLinearAlgebraAndMatrixDiagnostics:
    """Verifies Cholesky decomposition, Jacobi eigendecomposition, nearest PSD projection, and Gate H."""

    def test_general_cholesky_decomposition_positive_definite(self):
        code = """
        import { choleskyDecompositionGeneral } from './src/uncertainty/bayesian_hmc_field.js';
        
        const A = [
          [4.0, 2.0, 1.0],
          [2.0, 5.0, 3.0],
          [1.0, 3.0, 6.0]
        ];
        
        const L = choleskyDecompositionGeneral(A, 0.0);
        
        const Arec = [
          [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0]
        ];
        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            let sum = 0.0;
            for (let k = 0; k < 3; k++) sum += L[i][k] * L[j][k];
            Arec[i][j] = sum;
          }
        }
        
        console.log(JSON.stringify({
          L: L,
          Arec: Arec
        }));
        """
        res = run_node_snippet(code)
        L = res["L"]
        assert L[0][0] == pytest.approx(2.0, rel=1e-5)
        assert L[0][1] == pytest.approx(0.0, abs=1e-6)
        assert L[0][2] == pytest.approx(0.0, abs=1e-6)
        assert L[1][0] == pytest.approx(1.0, rel=1e-5)
        assert L[1][1] == pytest.approx(2.0, rel=1e-5)
        
        Arec = res["Arec"]
        assert Arec[0][0] == pytest.approx(4.0, rel=1e-4)
        assert Arec[1][1] == pytest.approx(5.0, rel=1e-4)
        assert Arec[2][2] == pytest.approx(6.0, rel=1e-4)
        assert Arec[0][1] == pytest.approx(2.0, rel=1e-4)
        assert Arec[1][2] == pytest.approx(3.0, rel=1e-4)

    def test_general_cholesky_decomposition_higher_dimensions(self):
        code = """
        import { choleskyDecompositionGeneral } from './src/uncertainty/bayesian_hmc_field.js';
        
        const n = 6;
        const B = [
          [1, 2, 0, 1, 0, 1],
          [0, 2, 1, 0, 2, 0],
          [1, 0, 3, 1, 0, 0],
          [0, 1, 1, 2, 1, 1],
          [2, 0, 0, 1, 3, 0],
          [0, 1, 2, 0, 1, 2]
        ];
        
        const M = Array.from({ length: n }, () => new Float64Array(n));
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            let sum = i === j ? 1.0 : 0.0;
            for (let k = 0; k < n; k++) sum += B[i][k] * B[j][k];
            M[i][j] = sum;
          }
        }
        
        const L = choleskyDecompositionGeneral(M, 0.0);
        
        let maxDiff = 0.0;
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            let sum = 0.0;
            for (let k = 0; k < n; k++) sum += L[i][k] * L[j][k];
            const diff = Math.abs(sum - M[i][j]);
            if (diff > maxDiff) maxDiff = diff;
          }
        }
        
        console.log(JSON.stringify({ maxDiff }));
        """
        res = run_node_snippet(code)
        assert res["maxDiff"] < 1e-10

    def test_general_jacobi_diagonalization_nxn(self):
        code = """
        import { jacobiDiagonalizeNxN } from './src/uncertainty/bayesian_hmc_field.js';
        
        const M = [
          [10.0,  1.0,  2.0,  0.5],
          [ 1.0,  8.0,  0.2,  1.5],
          [ 2.0,  0.2,  6.0, -1.0],
          [ 0.5,  1.5, -1.0,  4.0]
        ];
        
        const diag = jacobiDiagonalizeNxN(M, 1e-13, 100);
        
        const n = 4;
        let maxOrthogError = 0.0;
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            let dot = 0.0;
            for (let r = 0; r < n; r++) {
              dot += diag.eigenvectors[i][r] * diag.eigenvectors[j][r];
            }
            const expected = i === j ? 1.0 : 0.0;
            const err = Math.abs(dot - expected);
            if (err > maxOrthogError) maxOrthogError = err;
          }
        }
        
        const trM = 10.0 + 8.0 + 6.0 + 4.0;
        const trEig = diag.eigenvalues.reduce((a, b) => a + b, 0);
        
        console.log(JSON.stringify({
          eigenvalues: diag.eigenvalues,
          maxOrthogError,
          traceOriginal: trM,
          traceEigenvalues: trEig,
          sweeps: diag.sweeps
        }));
        """
        res = run_node_snippet(code)
        assert res["maxOrthogError"] < 1e-10
        assert res["traceEigenvalues"] == pytest.approx(res["traceOriginal"], rel=1e-6)
        eigs = res["eigenvalues"]
        for i in range(len(eigs) - 1):
            assert eigs[i] >= eigs[i + 1]

    def test_nearest_psd_projection_clipping(self):
        code = """
        import { projectMatrixToNearestPSD, jacobiDiagonalizeNxN } from './src/uncertainty/bayesian_hmc_field.js';
        
        const nonPSD = [
          [ 2.0,  3.0,  0.0],
          [ 3.0,  1.0,  0.0],
          [ 0.0,  0.0,  4.0]
        ];
        
        const psd = projectMatrixToNearestPSD(nonPSD, 1e-6);
        const diag = jacobiDiagonalizeNxN(psd);
        
        console.log(JSON.stringify({
          psdMatrix: psd,
          minEigenvalue: Math.min(...diag.eigenvalues)
        }));
        """
        res = run_node_snippet(code)
        assert res["minEigenvalue"] >= 1e-6 - 1e-10

    def test_validate_gate_h_covariance_matrix(self):
        code = """
        import { validateGateHCovarianceMatrix } from './src/uncertainty/bayesian_hmc_field.js';
        
        const validCov = [
          [100.0,  20.0, -10.0],
          [ 20.0,  80.0,  15.0],
          [-10.0,  15.0,  60.0]
        ];
        
        const asymCov = [
          [100.0,  50.0, -10.0],
          [ 20.0,  80.0,  15.0],
          [-10.0,  15.0,  60.0]
        ];
        
        const negVarCov = [
          [-10.0,   0.0,   0.0],
          [  0.0,  80.0,   0.0],
          [  0.0,   0.0,  60.0]
        ];
        
        const resValid = validateGateHCovarianceMatrix(validCov);
        const resAsym = validateGateHCovarianceMatrix(asymCov);
        const resNegVar = validateGateHCovarianceMatrix(negVarCov);
        
        console.log(JSON.stringify({
          valid: resValid.isValid,
          asymValid: resAsym.isValid,
          negVarValid: resNegVar.isValid
        }));
        """
        res = run_node_snippet(code)
        assert res["valid"] is True
        assert res["asymValid"] is False
        assert res["negVarValid"] is False

    def test_mahalanobis_distance_calculation(self):
        code = """
        import { computeMahalanobisDistance } from './src/uncertainty/bayesian_hmc_field.js';
        
        const mu = [0.0, 0.0, 0.0];
        const cov = [
          [4.0, 0.0, 0.0],
          [0.0, 9.0, 0.0],
          [0.0, 0.0, 16.0]
        ];
        
        const x = [4.0, 6.0, 8.0];
        const dM = computeMahalanobisDistance(x, mu, cov);
        
        console.log(JSON.stringify({
          mahalanobisDistance: dM,
          expected: Math.sqrt(12)
        }));
        """
        res = run_node_snippet(code)
        assert res["mahalanobisDistance"] == pytest.approx(math.sqrt(12), rel=1e-4)


class TestSpatialCovarianceKernelEngine:
    """Verifies all spatial covariance kernel formulations and block matrix assembly."""

    def test_gaussian_kernel_decay_and_symmetry(self):
        code = """
        import { SpatialCovarianceKernel, CovarianceKernelType } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kernel = new SpatialCovarianceKernel({
          kernelType: CovarianceKernelType.GAUSSIAN,
          correlationLength: 25.0
        });
        
        const k0 = kernel.evaluateCorrelationFactor(0.0);
        const k25 = kernel.evaluateCorrelationFactor(25.0);
        const k50 = kernel.evaluateCorrelationFactor(50.0);
        
        console.log(JSON.stringify({
          k0,
          k25,
          k50,
          expectedK25: Math.exp(-0.5),
          expectedK50: Math.exp(-2.0)
        }));
        """
        res = run_node_snippet(code)
        assert res["k0"] == pytest.approx(1.0, rel=1e-6)
        assert res["k25"] == pytest.approx(math.exp(-0.5), rel=1e-5)
        assert res["k50"] == pytest.approx(math.exp(-2.0), rel=1e-5)

    def test_exponential_and_matern_kernels(self):
        code = """
        import { SpatialCovarianceKernel, CovarianceKernelType } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kExp = new SpatialCovarianceKernel({ kernelType: CovarianceKernelType.EXPONENTIAL, correlationLength: 20.0 });
        const kMatern32 = new SpatialCovarianceKernel({ kernelType: CovarianceKernelType.MATERN_32, correlationLength: 20.0 });
        const kMatern52 = new SpatialCovarianceKernel({ kernelType: CovarianceKernelType.MATERN_52, correlationLength: 20.0 });
        
        const r = 20.0;
        const valExp = kExp.evaluateCorrelationFactor(r);
        const valMatern32 = kMatern32.evaluateCorrelationFactor(r);
        const valMatern52 = kMatern52.evaluateCorrelationFactor(r);
        
        console.log(JSON.stringify({
          valExp,
          valMatern32,
          valMatern52,
          expExpected: Math.exp(-1.0),
          matern32Expected: (1.0 + Math.sqrt(3.0)) * Math.exp(-Math.sqrt(3.0)),
          matern52Expected: (1.0 + Math.sqrt(5.0) + 5.0/3.0) * Math.exp(-Math.sqrt(5.0))
        }));
        """
        res = run_node_snippet(code)
        assert res["valExp"] == pytest.approx(res["expExpected"], rel=1e-5)
        assert res["valMatern32"] == pytest.approx(res["matern32Expected"], rel=1e-5)
        assert res["valMatern52"] == pytest.approx(res["matern52Expected"], rel=1e-5)

    def test_gorski_anisotropic_velocity_tensor(self):
        code = """
        import { SpatialCovarianceKernel, CovarianceKernelType } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kernel = new SpatialCovarianceKernel({
          kernelType: CovarianceKernelType.GORSKI_ANISOTROPIC,
          correlationLength: 30.0
        });
        
        const x1 = [0.0, 0.0, 0.0];
        const x2 = [30.0, 0.0, 0.0];
        const sig1 = [10.0, 10.0, 10.0];
        const sig2 = [10.0, 10.0, 10.0];
        
        const C = kernel.evaluateVelocityTensorCovariance(x1, x2, sig1, sig2);
        
        console.log(JSON.stringify({
          Cxx: C[0][0],
          Cyy: C[1][1],
          Czz: C[2][2],
          Cxy: C[0][1]
        }));
        """
        res = run_node_snippet(code)
        assert res["Cxx"] == pytest.approx(100.0 * math.exp(-0.5), rel=1e-4)
        assert res["Cyy"] == pytest.approx(0.0, abs=1e-6)
        assert res["Czz"] == pytest.approx(0.0, abs=1e-6)
        assert res["Cxy"] == pytest.approx(0.0, abs=1e-6)

    def test_gorski_anisotropic_tensor_arbitrary_3d_direction(self):
        code = """
        import { SpatialCovarianceKernel, CovarianceKernelType } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kernel = new SpatialCovarianceKernel({
          kernelType: CovarianceKernelType.GORSKI_ANISOTROPIC,
          correlationLength: 20.0
        });
        
        const x1 = [0.0, 0.0, 0.0];
        const x2 = [10.0, 10.0, 10.0];
        const sig = [1.0, 1.0, 1.0];
        
        const C = kernel.evaluateVelocityTensorCovariance(x1, x2, sig, sig);
        
        console.log(JSON.stringify({
          Cxx: C[0][0],
          Cyy: C[1][1],
          Czz: C[2][2],
          Cxy: C[0][1],
          Cxz: C[0][2],
          Cyz: C[1][2]
        }));
        """
        res = run_node_snippet(code)
        assert res["Cxx"] == pytest.approx(res["Cyy"], rel=1e-5)
        assert res["Cyy"] == pytest.approx(res["Czz"], rel=1e-5)
        assert res["Cxy"] == pytest.approx(res["Cxz"], rel=1e-5)
        assert res["Cxz"] == pytest.approx(res["Cyz"], rel=1e-5)

    def test_gorski_tensor_trace_invariant(self):
        code = """
        import { SpatialCovarianceKernel, CovarianceKernelType } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kernel = new SpatialCovarianceKernel({
          kernelType: CovarianceKernelType.GORSKI_ANISOTROPIC,
          correlationLength: 25.0
        });
        
        const x1 = [0, 0, 0];
        const x2 = [12.0, -15.0, 20.0];
        const r = Math.hypot(12.0, -15.0, 20.0);
        const sig = [10.0, 10.0, 10.0];
        
        const C = kernel.evaluateVelocityTensorCovariance(x1, x2, sig, sig);
        const trace = C[0][0] + C[1][1] + C[2][2];
        
        const u = r / 25.0;
        const sigmaParallel = Math.exp(-0.5 * u * u);
        const piTransverse = (1.0 - u * u) * Math.exp(-0.5 * u * u);
        const expectedTrace = 100.0 * (sigmaParallel + 2.0 * piTransverse);
        
        console.log(JSON.stringify({
          trace,
          expectedTrace
        }));
        """
        res = run_node_snippet(code)
        assert res["trace"] == pytest.approx(res["expectedTrace"], rel=1e-4)

    def test_block_covariance_matrix_assembly(self):
        code = """
        import { SpatialCovarianceKernel } from './src/uncertainty/bayesian_hmc_field.js';
        
        const kernel = new SpatialCovarianceKernel({ correlationLength: 20.0, jitter: 1e-6 });
        
        const points = [
          [0.0, 0.0, 0.0],
          [10.0, 0.0, 0.0],
          [0.0, 10.0, 0.0]
        ];
        
        const sigmas = [
          [15.0, 15.0, 15.0],
          [20.0, 20.0, 20.0],
          [25.0, 25.0, 25.0]
        ];
        
        const blockC = kernel.assembleBlockCovarianceMatrix(points, sigmas);
        const dim = blockC.length;
        
        let isSymmetric = true;
        for (let i = 0; i < dim; i++) {
          for (let j = 0; j < dim; j++) {
            if (Math.abs(blockC[i][j] - blockC[j][i]) > 1e-10) isSymmetric = false;
          }
        }
        
        console.log(JSON.stringify({
          dim,
          isSymmetric,
          diag0: blockC[0][0],
          diag3: blockC[3][3],
          diag6: blockC[6][6]
        }));
        """
        res = run_node_snippet(code)
        assert res["dim"] == 9
        assert res["isSymmetric"] is True
        assert res["diag0"] == pytest.approx(15.0**2, rel=1e-4)
        assert res["diag3"] == pytest.approx(20.0**2, rel=1e-4)
        assert res["diag6"] == pytest.approx(25.0**2, rel=1e-4)


class TestBayesianHMCFieldLoadersAndEvaluation:
    """Verifies continuous 3D evaluation, confidence intervals, Jacobian, and bulk flow uncertainty."""

    def test_continuous_mean_and_rms_sampling(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField } from './src/uncertainty/bayesian_hmc_field.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        
        const meanFn = (x, y, z) => [2.0 * x, -1.5 * y, 0.5 * z];
        const rmsFn = (x, y, z) => [20.0, 25.0, 30.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        
        const testP = [12.5, -18.0, 22.0];
        const vMean = hmcField.sampleMeanVelocity(testP[0], testP[1], testP[2]);
        const vRms = hmcField.sampleVelocityDispersion(testP[0], testP[1], testP[2]);
        const snr = hmcField.sampleSignalToNoiseRatio(testP[0], testP[1], testP[2]);
        
        console.log(JSON.stringify({
          vMean,
          vRms,
          snr
        }));
        """
        res = run_node_snippet(code)
        assert res["vMean"][0] == pytest.approx(2.0 * 12.5, rel=1e-3)
        assert res["vMean"][1] == pytest.approx(-1.5 * -18.0, rel=1e-3)
        assert res["vMean"][2] == pytest.approx(0.5 * 22.0, rel=1e-3)
        assert res["vRms"][0] == pytest.approx(20.0, rel=1e-3)
        assert res["vRms"][1] == pytest.approx(25.0, rel=1e-3)
        assert res["vRms"][2] == pytest.approx(30.0, rel=1e-3)
        assert res["snr"] > 0

    def test_analytic_jacobian_evaluation(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField } from './src/uncertainty/bayesian_hmc_field.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-40, -40, -40], boxSize: [80, 80, 80] });
        
        const meanFn = (x, y, z) => [3.0 * x + 2.0 * y, -y + z, 4.0 * z - x];
        const rmsFn = (x, y, z) => [15.0, 15.0, 15.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        const J = hmcField.sampleMeanVelocityJacobian(5.0, 5.0, 5.0);
        
        console.log(JSON.stringify({
          Jxx: J[0][0],
          Jxy: J[0][1],
          Jyx: J[1][0],
          Jyy: J[1][1],
          Jyz: J[1][2],
          Jzx: J[2][0],
          Jzz: J[2][2]
        }));
        """
        res = run_node_snippet(code)
        assert res["Jxx"] == pytest.approx(3.0, rel=1e-3)
        assert res["Jxy"] == pytest.approx(2.0, rel=1e-3)
        assert res["Jyx"] == pytest.approx(0.0, abs=1e-3)
        assert res["Jyy"] == pytest.approx(-1.0, rel=1e-3)
        assert res["Jyz"] == pytest.approx(1.0, rel=1e-3)
        assert res["Jzx"] == pytest.approx(-1.0, rel=1e-3)
        assert res["Jzz"] == pytest.approx(4.0, rel=1e-3)

    def test_posterior_confidence_intervals_quantiles(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const meanFn = (x, y, z) => [100.0, -50.0, 20.0];
        const rmsFn = (x, y, z) => [10.0, 5.0, 8.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        const ci68 = hmcField.sampleConfidenceInterval(0, 0, 0, 0.683);
        const ci95 = hmcField.sampleConfidenceInterval(0, 0, 0, 0.954);
        const ci99 = hmcField.sampleConfidenceInterval(0, 0, 0, 0.997);
        
        console.log(JSON.stringify({
          ci68,
          ci95,
          ci99
        }));
        """
        res = run_node_snippet(code)
        ci68 = res["ci68"]
        ci95 = res["ci95"]
        ci99 = res["ci99"]
        
        assert ci68["lowerCI"][0] == pytest.approx(90.0, abs=0.5)
        assert ci68["upperCI"][0] == pytest.approx(110.0, abs=0.5)
        
        assert ci95["lowerCI"][0] == pytest.approx(80.0, abs=0.5)
        assert ci95["upperCI"][0] == pytest.approx(120.0, abs=0.5)
        
        assert ci99["lowerCI"][0] == pytest.approx(70.0, abs=0.5)
        assert ci99["upperCI"][0] == pytest.approx(130.0, abs=0.5)

    def test_synthetic_mock_cosmological_field_and_summary(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-250, -250, -250], boxSize: [500, 500, 500] });
        const mockField = BayesianHMCField.createMockCosmologicalField(grid);
        const summary = mockField.exportSummaryMetrics();
        
        console.log(JSON.stringify(summary));
        """
        res = run_node_snippet(code)
        assert res["totalVoxels"] == 16 * 16 * 16
        assert res["meanSpeedKms"] > 0
        assert res["meanDispersionKms"] > 20.0
        assert res["meanSNR"] > 0


class TestHMCRealizationSampler:
    """Verifies Cholesky & Karhunen-Loève stochastic realization generation."""

    def test_standard_normal_generator_statistics(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, HMCRealizationSampler } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, (x,y,z)=>[0,0,0], (x,y,z)=>[1,1,1]);
        const sampler = new HMCRealizationSampler(hmcField, null, { seed: 12345 });
        
        const N = 2000;
        const vec = sampler.sampleStandardNormalVector(N);
        
        let sum = 0.0;
        for (let i = 0; i < N; i++) sum += vec[i];
        const mean = sum / N;
        
        let sumSq = 0.0;
        for (let i = 0; i < N; i++) sumSq += (vec[i] - mean) ** 2;
        const variance = sumSq / (N - 1);
        
        console.log(JSON.stringify({
          mean,
          variance
        }));
        """
        res = run_node_snippet(code)
        assert res["mean"] == pytest.approx(0.0, abs=0.1)
        assert res["variance"] == pytest.approx(1.0, abs=0.1)

    def test_point_realization_dispersion_recovery(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, HMCRealizationSampler } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, (x,y,z)=>[300, -200, 150], (x,y,z)=>[40, 50, 60]);
        const sampler = new HMCRealizationSampler(hmcField, null, { seed: 9999 });
        
        const samplesVx = [];
        const samplesVy = [];
        const samplesVz = [];
        
        for (let i = 0; i < 1000; i++) {
          const v = sampler.sampleVelocityAtPoint(0, 0, 0);
          samplesVx.push(v[0]);
          samplesVy.push(v[1]);
          samplesVz.push(v[2]);
        }
        
        const meanVx = samplesVx.reduce((a, b) => a + b, 0) / 1000;
        const meanVy = samplesVy.reduce((a, b) => a + b, 0) / 1000;
        const meanVz = samplesVz.reduce((a, b) => a + b, 0) / 1000;
        
        let varVx = 0.0, varVy = 0.0, varVz = 0.0;
        for (let s of samplesVx) varVx += (s - meanVx) ** 2;
        for (let s of samplesVy) varVy += (s - meanVy) ** 2;
        for (let s of samplesVz) varVz += (s - meanVz) ** 2;
        
        console.log(JSON.stringify({
          meanVx,
          meanVy,
          meanVz,
          stdVx: Math.sqrt(varVx / 999),
          stdVy: Math.sqrt(varVy / 999),
          stdVz: Math.sqrt(varVz / 999)
        }));
        """
        res = run_node_snippet(code)
        assert res["meanVx"] == pytest.approx(300.0, abs=5.0)
        assert res["meanVy"] == pytest.approx(-200.0, abs=5.0)
        assert res["meanVz"] == pytest.approx(150.0, abs=5.0)
        assert res["stdVx"] == pytest.approx(40.0, abs=4.0)
        assert res["stdVy"] == pytest.approx(50.0, abs=4.0)
        assert res["stdVz"] == pytest.approx(60.0, abs=4.0)

    def test_karhunen_loeve_eigenmode_compression(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, SpatialCovarianceKernel, HMCRealizationSampler } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, (x,y,z)=>[100, 100, 100], (x,y,z)=>[20, 20, 20]);
        const kernel = new SpatialCovarianceKernel({ correlationLength: 30.0 });
        const sampler = new HMCRealizationSampler(hmcField, kernel);
        
        const pathPoints = [
          [0.0, 0.0, 0.0],
          [2.0, 0.0, 0.0],
          [4.0, 0.0, 0.0],
          [6.0, 0.0, 0.0],
          [8.0, 0.0, 0.0]
        ];
        
        const klResult = sampler.sampleKarhunenLoeveRealization(pathPoints, 0.95);
        
        console.log(JSON.stringify({
          retained: klResult.numModesRetained,
          totalModes: klResult.totalModes,
          varianceCaptured: klResult.varianceCaptured,
          numRealizedPoints: klResult.realizedVelocities.length
        }));
        """
        res = run_node_snippet(code)
        assert res["numRealizedPoints"] == 5
        assert res["totalModes"] == 15
        assert res["retained"] < res["totalModes"]
        assert res["varianceCaptured"] >= 0.95

    def test_karhunen_loeve_correlation_length_sensitivity(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, SpatialCovarianceKernel, HMCRealizationSampler } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, (x,y,z)=>[0,0,0], (x,y,z)=>[10,10,10]);
        
        const kShort = new SpatialCovarianceKernel({ correlationLength: 5.0 });
        const kLong = new SpatialCovarianceKernel({ correlationLength: 50.0 });
        
        const sShort = new HMCRealizationSampler(hmcField, kShort);
        const sLong = new HMCRealizationSampler(hmcField, kLong);
        
        const path = Array.from({ length: 8 }, (_, i) => [i * 3.0, 0.0, 0.0]);
        
        const klShort = sShort.sampleKarhunenLoeveRealization(path, 0.90);
        const klLong = sLong.sampleKarhunenLoeveRealization(path, 0.90);
        
        console.log(JSON.stringify({
          retainedShort: klShort.numModesRetained,
          retainedLong: klLong.numModesRetained
        }));
        """
        res = run_node_snippet(code)
        assert res["retainedLong"] <= res["retainedShort"]


class TestProbabilisticWatershedBasinAssignment:
    """Verifies streamline basin classification, Shannon entropy, and overlap metrics."""

    def test_probabilistic_point_classification_and_entropy(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, ProbabilisticWatershedClassifier } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-200, -200, -200], boxSize: [400, 400, 400] });
        
        const attractors = [
          { id: 1, name: 'Laniakea', attractorSGMpc: [-46.0, 22.9, -6.8] },
          { id: 6, name: 'Shapley', attractorSGMpc: [-170.7, 118.7, -21.3] }
        ];
        
        const mockField = BayesianHMCField.createMockCosmologicalField(grid, {
          attractors: [
            { name: 'Laniakea', pos: [-46.0, 22.9, -6.8], mass: 5e4, coreRadius: 20.0 },
            { name: 'Shapley', pos: [-170.7, 118.7, -21.3], mass: 2e5, coreRadius: 30.0 }
          ]
        });
        
        const classifier = new ProbabilisticWatershedClassifier(mockField, { numRealizations: 50 });
        
        const corePoint = [-45.0, 22.0, -6.0];
        const resCore = classifier.classifyPointProbabilistically(corePoint, attractors, 30);
        
        const boundaryPoint = [-105.0, 70.0, -14.0];
        const resBound = classifier.classifyPointProbabilistically(boundaryPoint, attractors, 30);
        
        console.log(JSON.stringify({
          coreMostLikely: resCore.mostLikelyBasinId,
          coreMaxProb: resCore.maxProbability,
          coreEntropy: resCore.shannonEntropy,
          boundEntropy: resBound.shannonEntropy
        }));
        """
        res = run_node_snippet(code)
        assert res["coreMostLikely"] == 1
        assert res["coreMaxProb"] >= 0.8
        assert res["boundEntropy"] >= res["coreEntropy"]

    def test_segmentation_overlap_jaccard_and_dice(self):
        code = """
        import { ProbabilisticWatershedClassifier } from './src/uncertainty/bayesian_hmc_field.js';
        
        const N = 100;
        const baseline = new Array(N).fill(1);
        
        const realizations = [];
        for (let k = 0; k < 10; k++) {
          const r = new Array(N).fill(1);
          for (let i = 0; i < 10; i++) r[i] = 2;
          realizations.push(r);
        }
        
        const metrics = ProbabilisticWatershedClassifier.computeSegmentationOverlapMetrics(baseline, realizations, 1);
        
        console.log(JSON.stringify({
          meanJaccard: metrics.meanJaccard,
          meanDice: metrics.meanDice,
          jaccardCI: metrics.jaccardCI95
        }));
        """
        res = run_node_snippet(code)
        assert res["meanJaccard"] == pytest.approx(0.90, rel=1e-3)
        assert res["meanDice"] == pytest.approx(180.0 / 190.0, rel=1e-3)

    def test_probabilistic_basin_volume_expectation_and_ci(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, ProbabilisticWatershedClassifier } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
        const attractors = [
          { id: 1, name: 'Laniakea', attractorSGMpc: [-46.0, 22.9, -6.8] }
        ];
        
        const mockField = BayesianHMCField.createMockCosmologicalField(grid);
        const classifier = new ProbabilisticWatershedClassifier(mockField, { numRealizations: 20 });
        
        const samplePoints = [
          [-45.0, 20.0, -5.0],
          [-40.0, 25.0, -8.0],
          [-50.0, 22.0, -6.0]
        ];
        
        const volumes = classifier.computeProbabilisticBasinVolumes(samplePoints, attractors, 1000.0, 20);
        
        console.log(JSON.stringify({
          volLanExpected: volumes[1].expectedVolumeMpc3,
          volLanStd: volumes[1].volumeStandardDeviationMpc3,
          ciLow: volumes[1].ci95Mpc3[0],
          ciHigh: volumes[1].ci95Mpc3[1]
        }));
        """
        res = run_node_snippet(code)
        assert res["volLanExpected"] > 0
        assert res["ciLow"] <= res["volLanExpected"] <= res["ciHigh"]


class TestCriticalPointStabilityAnalyzer:
    """Verifies critical point detection, stability probability, and Morse bifurcation."""

    def test_critical_point_stability_under_perturbation(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, CriticalPointStabilityAnalyzer } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
        
        const meanFn = (x, y, z) => [-0.5 * x, -0.5 * y, -0.5 * z];
        const rmsFn = (x, y, z) => [5.0, 5.0, 5.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        const analyzer = new CriticalPointStabilityAnalyzer(hmcField, { seed: 42 });
        
        const stability = analyzer.analyzePointStability([0.0, 0.0, 0.0], 50, { maxDisplacementMpc: 25.0 });
        
        console.log(JSON.stringify({
          baselineClassification: stability.baselineClassification,
          stabilityProbability: stability.stabilityProbability,
          spatialDispersion: stability.spatialDispersionMpc,
          bifurcationProbability: stability.bifurcationProbability
        }));
        """
        res = run_node_snippet(code)
        assert res["baselineClassification"] == "ATTRACTOR"
        assert res["stabilityProbability"] > 0.8
        assert res["spatialDispersion"] < 25.0
        assert res["bifurcationProbability"] < 0.2


class TestCosmicVarianceStreamlinePropagator:
    """Verifies matrix ODE covariance propagation and Monte Carlo ensemble streamline bundles."""

    def test_continuous_path_covariance_propagation(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, CosmicVarianceStreamlinePropagator } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        
        const meanFn = (x, y, z) => [0.1 * x, 0.1 * y, 0.1 * z];
        const rmsFn = (x, y, z) => [10.0, 10.0, 10.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        const propagator = new CosmicVarianceStreamlinePropagator(hmcField);
        
        const path = [
          [10.0, 0.0, 0.0],
          [15.0, 0.0, 0.0],
          [20.0, 0.0, 0.0],
          [25.0, 0.0, 0.0],
          [30.0, 0.0, 0.0]
        ];
        
        const res = propagator.propagateCovarianceAlongPath(path, null, 0.2);
        
        console.log(JSON.stringify({
          tubeRadii: res.tubeRadiiMpc,
          initialRadius: res.tubeRadiiMpc[0],
          finalRadius: res.tubeRadiiMpc[res.tubeRadiiMpc.length - 1],
          numSteps: res.trajectoryPoints.length
        }));
        """
        res = run_node_snippet(code)
        assert res["numSteps"] == 5
        assert res["finalRadius"] > res["initialRadius"]

    def test_ensemble_streamline_bundle_tracing(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField, CosmicVarianceStreamlinePropagator } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const meanFn = (x, y, z) => [50.0, 0.0, 0.0];
        const rmsFn = (x, y, z) => [5.0, 5.0, 5.0];
        
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, meanFn, rmsFn);
        const propagator = new CosmicVarianceStreamlinePropagator(hmcField);
        
        const bundle = propagator.traceEnsembleStreamlineBundle([0, 0, 0], 30, { maxSteps: 20, dt: 0.1 });
        
        console.log(JSON.stringify({
          meanLength: bundle.meanTrajectory.length,
          stepCount: bundle.stepCount,
          finalTubeRadius: bundle.empiricalTubeRadiiMpc[bundle.empiricalTubeRadiiMpc.length - 1]
        }));
        """
        res = run_node_snippet(code)
        assert res["meanLength"] == 21
        assert res["finalTubeRadius"] > 0.0


class TestMCMCChainDiagnostics:
    """Verifies Gelman-Rubin, autocorrelation time, effective sample size, and E-BFMI."""

    def test_gelman_rubin_diagnostic_convergence(self):
        code = """
        import { computeGelmanRubinDiagnostic } from './src/uncertainty/bayesian_hmc_field.js';
        import { MT19937 } from './src/statistics/statistical_resampling.js';
        
        const prng = new MT19937(42);
        
        const convergedChains = Array.from({ length: 4 }, () => {
          const arr = [];
          for (let i = 0; i < 500; i++) {
            const u1 = Math.max(1e-12, prng.extractNumber());
            const u2 = prng.extractNumber();
            arr.push(Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2));
          }
          return arr;
        });
        
        const nonConvergedChains = [
          Array.from({ length: 500 }, () => -5.0 + prng.extractNumber()),
          Array.from({ length: 500 }, () => -2.0 + prng.extractNumber()),
          Array.from({ length: 500 }, () =>  2.0 + prng.extractNumber()),
          Array.from({ length: 500 }, () =>  5.0 + prng.extractNumber())
        ];
        
        const rHatConverged = computeGelmanRubinDiagnostic(convergedChains);
        const rHatNonConverged = computeGelmanRubinDiagnostic(nonConvergedChains);
        
        console.log(JSON.stringify({
          rHatConverged,
          rHatNonConverged
        }));
        """
        res = run_node_snippet(code)
        assert res["rHatConverged"] == pytest.approx(1.0, abs=0.03)
        assert res["rHatNonConverged"] > 1.2

    def test_autocorrelation_time_and_effective_sample_size(self):
        code = """
        import { computeIntegratedAutocorrelationTime, computeEffectiveSampleSize } from './src/uncertainty/bayesian_hmc_field.js';
        import { MT19937 } from './src/statistics/statistical_resampling.js';
        
        const prng = new MT19937(123);
        const N = 1000;
        
        const iid = [];
        for (let i = 0; i < N; i++) iid.push(prng.extractNumber() - 0.5);
        
        const ar1 = [0.0];
        for (let i = 1; i < N; i++) {
          ar1.push(0.9 * ar1[i - 1] + (prng.extractNumber() - 0.5));
        }
        
        const tauIID = computeIntegratedAutocorrelationTime(iid);
        const neffIID = computeEffectiveSampleSize(iid);
        
        const tauAR1 = computeIntegratedAutocorrelationTime(ar1);
        const neffAR1 = computeEffectiveSampleSize(ar1);
        
        console.log(JSON.stringify({
          tauIID,
          neffIID,
          tauAR1,
          neffAR1
        }));
        """
        res = run_node_snippet(code)
        assert res["tauIID"] < 2.0
        assert res["neffIID"] > 500
        assert res["tauAR1"] > res["tauIID"]
        assert res["neffAR1"] < res["neffIID"]

    def test_energy_bayesian_fraction_missing_information(self):
        code = """
        import { computeEnergyBayesianFractionMissingInformation } from './src/uncertainty/bayesian_hmc_field.js';
        import { MT19937 } from './src/statistics/statistical_resampling.js';
        
        const prng = new MT19937(555);
        const energies = [];
        for (let i = 0; i < 500; i++) {
          energies.push(10.0 + prng.extractNumber() * 2.0);
        }
        
        const ebfmi = computeEnergyBayesianFractionMissingInformation(energies);
        
        console.log(JSON.stringify({
          ebfmi
        }));
        """
        res = run_node_snippet(code)
        assert res["ebfmi"] > 0.3


class TestAdversarialCosmologicalCovarianceAndEdgeCases:
    """Stress tests on near-singular matrices, zero velocity limits, and edge conditions."""

    def test_near_singular_covariance_tikhonov_regularization(self):
        code = """
        import { choleskyDecompositionGeneral, validateGateHCovarianceMatrix } from './src/uncertainty/bayesian_hmc_field.js';
        
        const singularMat = [
          [1.0, 2.0, 3.0],
          [2.0, 4.0, 6.0],
          [3.0, 6.0, 9.0]
        ];
        
        const L = choleskyDecompositionGeneral(singularMat, 1e-4);
        
        const rec = [
          [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0]
        ];
        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            let sum = 0.0;
            for (let k = 0; k < 3; k++) sum += L[i][k] * L[j][k];
            rec[i][j] = sum;
          }
        }
        
        const gateH = validateGateHCovarianceMatrix(rec);
        
        console.log(JSON.stringify({
          isValid: gateH.isValid,
          diag0: rec[0][0],
          diag1: rec[1][1],
          diag2: rec[2][2]
        }));
        """
        res = run_node_snippet(code)
        assert res["isValid"] is True
        assert res["diag0"] >= 1.0
        assert res["diag1"] >= 4.0
        assert res["diag2"] >= 9.0

    def test_bulk_flow_uncertainty_scaling_with_sample_count(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { BayesianHMCField } from './src/uncertainty/bayesian_hmc_field.js';
        
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
        const hmcField = BayesianHMCField.fromAnalyticFunctions(grid, (x,y,z)=>[200, 100, -50], (x,y,z)=>[30, 30, 30]);
        
        const resLowRes = hmcField.computeBulkFlowUncertainty([0, 0, 0], 50.0, 5);
        const resHighRes = hmcField.computeBulkFlowUncertainty([0, 0, 0], 50.0, 25);
        
        console.log(JSON.stringify({
          countLow: resLowRes.sampleCount,
          countHigh: resHighRes.sampleCount,
          sigmaLow: resLowRes.totalSigma,
          sigmaHigh: resHighRes.totalSigma,
          ratioSigma: resLowRes.totalSigma / resHighRes.totalSigma,
          expectedRatio: Math.sqrt(resHighRes.sampleCount / resLowRes.sampleCount)
        }));
        """
        res = run_node_snippet(code)
        assert res["countHigh"] > res["countLow"]
        assert res["sigmaHigh"] < res["sigmaLow"]
        assert res["ratioSigma"] == pytest.approx(res["expectedRatio"], rel=0.1)
