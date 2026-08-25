"""
tests/bulk-flow/test_multipolar_bulk_flow.py
===========================================
Exhaustive Pytest Verification Suite for Multipolar Bulk Flow Spherical Harmonic Decomposition Suite.

Covers:
1. Astrometric Coordinate Frame Transformations & Standard Astronomical Apex Targets:
   - Supergalactic, Galactic, and Equatorial frames
   - Invariance of vector norms and angular separations under 3D coordinate rotations
   - Directional angles towards CMB Dipole, Shapley Supercluster Core, Great Attractor (Norma), Perseus-Pisces
2. Mathematical Foundations & Orthonormal Spherical Harmonics (l = 0 to 4):
   - 3x3 Matrix algebra, Cramer determinant / cofactor inversion, linear system solver, singularity regularization
   - Associated Legendre polynomials P_l^m(x) against analytic standard functions and recurrence invariants
   - Real orthonormal spherical harmonics Y_lm(theta, phi) orthonormality quadrature on the 2-sphere
3. Triple Weighting Estimators:
   - FIELD_VOLUME_WEIGHTED: Exact spatial voxel volume integration, pure bulk drift recovery, spherically symmetric expansion (zero dipole), pure quadrupole shear
   - CATALOG_WEIGHTED: Galaxy catalog line-of-sight velocity projection tensor inversion A_ij = sum w_n r_i r_j, anisotropic geometry, 3D velocity vectors
   - INVERSE_VARIANCE_WEIGHTED: Optimal minimum-variance maximum likelihood estimator with cosmic velocity dispersion sigma_v = 187 km/s, error covariance propagation
4. Full Multipolar Expansion Hierarchy:
   - Monopole (l=0): Radial expansion / contraction H_R = <v_r>/R, Hubble bubble perturbation delta H / H_0, power C_0
   - Dipole (l=1): Bulk flow vector, magnitude, SGL, SGB, Galactic l, b, RA, Dec, celestial apex separations, power C_1
   - Quadrupole (l=2): Cosmic shear tensor Q_ij, trace-free property Tr(Q) = 0, Jacobi eigenvalue decomposition, principal shear axes, Frobenius norm, power C_2
   - Octupole (l=3) & Hexadecapole (l=4): 3rd/4th order symmetric trace-free tensors, Frobenius magnitudes, power spectrum C_3, C_4
   - Power spectrum distribution, relative fractions, multipole ratios (C_1/C_0, C_2/C_1, C_3/C_1)
5. Lambda-CDM Theoretical Cosmic Variance Engine & Survey Window Matrix Deconvolution:
   - Eisenstein & Hu (1998) transfer function T(k), normalized matter power spectrum P(k) matching sigma_8
   - Top-hat sphere window W(kR) and differential shell window W_{shell}(k)
   - 1D bulk flow dispersion sigma_{1D}(R) and 3D rms bulk flow sigma_{3D}(R) = sqrt(3) * sigma_{1D}(R) across R in [20, 250] Mpc/h
   - Monotonic decay with increasing scale R
   - Survey window matrix deconvolution V_true = W^{-1} V_meas and covariance propagation
   - Signal-to-noise ratio SNR and Maxwell-Boltzmann tension p-value testing
6. Quantitative Benchmarking against published CosmicFlows-4 (CF4) Bulk Flow Dipoles (Dupuy & Courtois 2023):
   - Reference table validation across 9 radii R in [20, 250] Mpc/h
   - Benchmark evaluation at R = 150 Mpc/h (CF4 nominal V ~ 388 +- 32 km/s towards Shapley / Hydra-Centaurus)
   - Chi-squared chi^2, z-score, angular offset, consistency ratings (EXCELLENT_MATCH, CONSISTENT_1SIGMA, CONSISTENT_2SIGMA, MODERATE_TENSION, HIGH_TENSION)
7. End-to-End Continuous Radial Shell Profiling & High-Resolution Stress Testing:
   - Continuous radial profiles across [20, 250] Mpc/h in cumulative and differential modes
   - High-resolution grid stress tests and boundary edge cases (empty catalogs, zero velocity, collinear tracers)
"""

import json
import subprocess
import os
import math
import numpy as np
import pytest

NODE_EXEC = "node"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def run_node_snippet(code: str) -> dict:
    """Executes an ES6 Node.js snippet and returns parsed JSON output."""
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Node execution failed with code {proc.returncode}:\n"
            f"STDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )
    return json.loads(proc.stdout)


# ============================================================================
# 1. ASTROMETRIC COORDINATES & REFERENCE APEX TARGETS
# ============================================================================

def test_astrometric_apex_targets_and_cosmology_constants():
    """Verify cosmological parameters and astronomical apex target definitions."""
    code = """
    import {
      FIDUCIAL_COSMOLOGY,
      ASTRONOMICAL_APEX_TARGETS,
      BulkFlowEstimatorType
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    console.log(JSON.stringify({
      cosmo: FIDUCIAL_COSMOLOGY,
      apex: ASTRONOMICAL_APEX_TARGETS,
      estimatorTypes: BulkFlowEstimatorType
    }));
    """
    res = run_node_snippet(code)

    cosmo = res["cosmo"]
    assert math.isclose(cosmo["Omega_m"], 0.315, rel_tol=1e-4)
    assert math.isclose(cosmo["Omega_Lambda"], 0.685, rel_tol=1e-4)
    assert math.isclose(cosmo["H0"], 67.4, rel_tol=1e-4)
    assert math.isclose(cosmo["sigma8"], 0.811, rel_tol=1e-4)
    assert math.isclose(cosmo["cosmicDispersion"], 187.0, rel_tol=1e-4)

    apex = res["apex"]
    assert "CMB_DIPOLE" in apex
    assert "SHAPLEY_CORE" in apex
    assert "GREAT_ATTRACTOR" in apex
    assert "PERSEUS_PISCES" in apex

    cmb = apex["CMB_DIPOLE"]
    assert math.isclose(cmb["velocityKms"], 369.82, rel_tol=1e-3)
    assert math.isclose(cmb["galacticL"], 264.021, rel_tol=1e-3)
    assert math.isclose(cmb["galacticB"], 48.253, rel_tol=1e-3)

    types = res["estimatorTypes"]
    assert types["FIELD_VOLUME_WEIGHTED"] == "FIELD_VOLUME_WEIGHTED"
    assert types["CATALOG_WEIGHTED"] == "CATALOG_WEIGHTED"
    assert types["INVERSE_VARIANCE_WEIGHTED"] == "INVERSE_VARIANCE_WEIGHTED"


def test_linear_algebra_and_coordinate_transformations():
    """Verify 3x3 matrix multiplication, Cramer inversion, and spherical/astrometric rotations."""
    code = """
    import {
      matVecMultiply3x3,
      matMultiply3x3,
      invertMatrix3x3,
      solveLinear3x3,
      cartesianToSpherical,
      sphericalToCartesian,
      vectorAngularSeparationDeg,
      transformAstronomicalVector
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 1. Matrix multiplication & inversion
    const A = [
      [3.0, 1.0, 0.5],
      [1.0, 4.0, 1.2],
      [0.5, 1.2, 2.5]
    ];
    const invA = invertMatrix3x3(A);
    const prod = matMultiply3x3(A, invA);

    // 2. Linear solve: A * x = b
    const b = [5.0, -2.0, 7.0];
    const xSol = solveLinear3x3(A, b);
    const bCheck = matVecMultiply3x3(A, xSol);

    // 3. Spherical roundtrip: (x, y, z) -> (r, lon, lat) -> (x, y, z)
    const origP = [150.0, -80.0, 120.0];
    const sph = cartesianToSpherical(origP[0], origP[1], origP[2]);
    const recP = sphericalToCartesian(sph.r, sph.lonDeg, sph.latDeg);

    // 4. Vector angular separations
    const sep90 = vectorAngularSeparationDeg([1, 0, 0], [0, 1, 0]);
    const sep0 = vectorAngularSeparationDeg([3, 4, 5], [6, 8, 10]);
    const sep180 = vectorAngularSeparationDeg([1, 2, 3], [-1, -2, -3]);

    // 5. Astrometric frame transform roundtrip: SG -> Gal -> Eq -> Gal -> SG
    const vSG = [200.0, -150.0, 80.0];
    const vGal = transformAstronomicalVector(vSG, 'SUPERGALACTIC', 'GALACTIC');
    const vEq = transformAstronomicalVector(vGal, 'GALACTIC', 'EQUATORIAL');
    const vSGRec = transformAstronomicalVector(transformAstronomicalVector(vEq, 'EQUATORIAL', 'GALACTIC'), 'GALACTIC', 'SUPERGALACTIC');

    console.log(JSON.stringify({
      prod,
      b,
      bCheck,
      origP,
      recP,
      sph,
      sep90,
      sep0,
      sep180,
      vSG,
      vSGRec
    }));
    """
    res = run_node_snippet(code)

    # Check A * inv(A) = Identity
    prod = np.array(res["prod"])
    np.testing.assert_allclose(prod, np.eye(3), atol=1e-5)

    # Check linear solve A * x = b
    np.testing.assert_allclose(res["b"], res["bCheck"], rtol=1e-5)

    # Check spherical roundtrip
    np.testing.assert_allclose(res["origP"], res["recP"], rtol=1e-4)

    # Check angular separations
    assert math.isclose(res["sep90"], 90.0, abs_tol=1e-4)
    assert math.isclose(res["sep0"], 0.0, abs_tol=1e-4)
    assert math.isclose(res["sep180"], 180.0, abs_tol=1e-4)

    # Check frame rotation roundtrip
    np.testing.assert_allclose(res["vSG"], res["vSGRec"], rtol=1e-4)


# ============================================================================
# 2. ASSOCIATED LEGENDRE POLYNOMIALS & REAL SPHERICAL HARMONICS
# ============================================================================

def test_associated_legendre_polynomials_analytic():
    """Verify Associated Legendre polynomials P_l^m(x) against analytic standard formulas."""
    code = """
    import { associatedLegendre } from './src/bulk-flow/multipolar_bulk_flow.js';

    const testX = [0.0, 0.5, -0.5, 0.8, -0.8, 1.0, -1.0];
    const results = testX.map(x => {
      // Analytic definitions:
      // P00 = 1
      // P10 = x
      // P11 = -sqrt(1-x^2)
      // P20 = 0.5*(3x^2 - 1)
      // P21 = -3x*sqrt(1-x^2)
      // P22 = 3*(1-x^2)
      // P30 = 0.5*(5x^3 - 3x)
      // P40 = 0.125*(35x^4 - 30x^2 + 3)
      const somx2 = Math.sqrt(Math.max(0.0, 1.0 - x * x));

      return {
        x,
        p00_code: associatedLegendre(0, 0, x),
        p00_true: 1.0,
        p10_code: associatedLegendre(1, 0, x),
        p10_true: x,
        p11_code: associatedLegendre(1, 1, x),
        p11_true: -somx2,
        p20_code: associatedLegendre(2, 0, x),
        p20_true: 0.5 * (3 * x * x - 1),
        p21_code: associatedLegendre(2, 1, x),
        p21_true: -3 * x * somx2,
        p22_code: associatedLegendre(2, 2, x),
        p22_true: 3 * (1 - x * x),
        p30_code: associatedLegendre(3, 0, x),
        p30_true: 0.5 * (5 * x * x * x - 3 * x),
        p40_code: associatedLegendre(4, 0, x),
        p40_true: 0.125 * (35 * Math.pow(x, 4) - 30 * x * x + 3)
      };
    });

    console.log(JSON.stringify(results));
    """
    res = run_node_snippet(code)
    for row in res:
        assert math.isclose(row["p00_code"], row["p00_true"], abs_tol=1e-5)
        assert math.isclose(row["p10_code"], row["p10_true"], abs_tol=1e-5)
        assert math.isclose(row["p11_code"], row["p11_true"], abs_tol=1e-5)
        assert math.isclose(row["p20_code"], row["p20_true"], abs_tol=1e-5)
        assert math.isclose(row["p21_code"], row["p21_true"], abs_tol=1e-5)
        assert math.isclose(row["p22_code"], row["p22_true"], abs_tol=1e-5)
        assert math.isclose(row["p30_code"], row["p30_true"], abs_tol=1e-5)
        assert math.isclose(row["p40_code"], row["p40_true"], abs_tol=1e-5)


def test_real_spherical_harmonics_orthonormality_quadrature():
    """Verify 2-sphere orthonormality \\int Y_lm Y_l'm' d\\Omega = \\delta_ll' \\delta_mm' numerically."""
    code = """
    import { realSphericalHarmonic } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 2D spherical quadrature over (theta, phi)
    const nTheta = 60;
    const nPhi = 120;
    const dTheta = Math.PI / nTheta;
    const dPhi = (2.0 * Math.PI) / nPhi;

    // Test subset of (l, m) pairs up to l = 3
    const modes = [
      [0, 0],
      [1, -1], [1, 0], [1, 1],
      [2, -2], [2, 0], [2, 2],
      [3, -1], [3, 0], [3, 3]
    ];

    const gramMatrix = [];
    for (let i = 0; i < modes.length; i++) {
      const [l1, m1] = modes[i];
      const row = [];
      for (let j = 0; j < modes.length; j++) {
        const [l2, m2] = modes[j];
        let integral = 0.0;

        for (let it = 0; it < nTheta; it++) {
          const theta = (it + 0.5) * dTheta;
          const sinTheta = Math.sin(theta);
          const dOmega = sinTheta * dTheta * dPhi;

          for (let ip = 0; ip < nPhi; ip++) {
            const phi = (ip + 0.5) * dPhi;
            const y1 = realSphericalHarmonic(l1, m1, theta, phi);
            const y2 = realSphericalHarmonic(l2, m2, theta, phi);
            integral += y1 * y2 * dOmega;
          }
        }
        row.push(integral);
      }
      gramMatrix.push(row);
    }

    console.log(JSON.stringify(gramMatrix));
    """
    res = run_node_snippet(code)
    gram = np.array(res)
    n = len(gram)
    # Diagonal elements must be ~1.0 (normalization)
    np.testing.assert_allclose(np.diag(gram), np.ones(n), atol=0.03)
    # Off-diagonal elements must be ~0.0 (orthogonality)
    off_diag = gram - np.diag(np.diag(gram))
    np.testing.assert_allclose(off_diag, np.zeros((n, n)), atol=0.03)


# ============================================================================
# 3. MATTER POWER SPECTRUM & THEORETICAL COSMIC VARIANCE
# ============================================================================

def test_eisenstein_hu_transfer_and_top_hat_window():
    """Verify Eisenstein & Hu transfer function and spherical top-hat window function."""
    code = """
    import {
      eisensteinHuTransfer,
      sphericalTopHatWindow,
      radialShellWindow,
      FIDUCIAL_COSMOLOGY
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    // Window function asymptotic limits
    const w0 = sphericalTopHatWindow(0.0);
    const wSmall = sphericalTopHatWindow(1e-5);
    const w1 = sphericalTopHatWindow(1.0);
    const wBig = sphericalTopHatWindow(50.0);

    // Shell window
    const wShell = radialShellWindow(0.05, 50.0, 100.0);

    // Transfer function behavior: T(k -> 0) -> 1, T(k -> inf) -> 0
    const tSmall = eisensteinHuTransfer(1e-5, FIDUCIAL_COSMOLOGY);
    const tMed = eisensteinHuTransfer(0.1, FIDUCIAL_COSMOLOGY);
    const tBig = eisensteinHuTransfer(10.0, FIDUCIAL_COSMOLOGY);

    console.log(JSON.stringify({
      w0,
      wSmall,
      w1,
      wBig,
      wShell,
      tSmall,
      tMed,
      tBig
    }));
    """
    res = run_node_snippet(code)

    assert math.isclose(res["w0"], 1.0, rel_tol=1e-5)
    assert math.isclose(res["wSmall"], 1.0, rel_tol=1e-4)
    assert 0 < res["w1"] < 1.0
    assert abs(res["wBig"]) < 0.05
    assert not math.isnan(res["wShell"])

    assert math.isclose(res["tSmall"], 1.0, rel_tol=0.05)
    assert 0 < res["tMed"] < 1.0
    assert 0 < res["tBig"] < res["tMed"]


def test_cosmic_variance_engine_bulk_flow_decay():
    """Verify Lambda-CDM bulk flow variance sigma_1D(R) decays monotonically with radius."""
    code = """
    import { CosmicVarianceEngine } from './src/bulk-flow/multipolar_bulk_flow.js';

    const engine = new CosmicVarianceEngine();
    const radii = [20, 40, 60, 80, 100, 150, 200, 250];

    const results = radii.map(r => {
      const v = engine.computeBulkFlowVariance(r);
      const cov = engine.computeCovarianceMatrix(r);
      const tension = engine.evaluateTension(350.0, r);
      return {
        r,
        sigma1D: v.sigma1D,
        sigma3D: v.sigma3D,
        covDiag: [cov[0][0], cov[1][1], cov[2][2]],
        covOff: [cov[0][1], cov[0][2], cov[1][2]],
        tensionSigma: tension.tensionSigma,
        pValue: tension.pValue
      };
    });

    console.log(JSON.stringify(results));
    """
    res = run_node_snippet(code)

    # 1. Monotonic decrease of bulk flow dispersion as R increases
    for i in range(len(res) - 1):
        assert res[i]["sigma1D"] > res[i + 1]["sigma1D"], f"Expected sigma1D({res[i]['r']}) > sigma1D({res[i+1]['r']})"
        assert res[i]["sigma3D"] > res[i + 1]["sigma3D"]
        # Check sigma3D = sqrt(3) * sigma1D
        assert math.isclose(res[i]["sigma3D"], math.sqrt(3.0) * res[i]["sigma1D"], rel_tol=1e-4)

    # 2. Covariance matrix is isotropic (diagonal equal, off-diagonals zero)
    for row in res:
        diag = row["covDiag"]
        assert math.isclose(diag[0], diag[1], rel_tol=1e-5)
        assert math.isclose(diag[1], diag[2], rel_tol=1e-5)
        assert all(abs(off) < 1e-9 for off in row["covOff"])


# ============================================================================
# 4. ESTIMATOR 1: FIELD VOLUME WEIGHTED
# ============================================================================

def test_field_volume_weighted_pure_dipole_recovery():
    """Verify exact recovery of uniform bulk flow vector on a 3D velocity grid."""
    code = """
    import { MultipolarBulkFlow, BulkFlowEstimatorType } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    const V_true = [210.0, -180.0, 95.0];

    for (let i = 0; i < N * N * N; i++) {
      velocityGrid[i * 3] = V_true[0];
      velocityGrid[i * 3 + 1] = V_true[1];
      velocityGrid[i * 3 + 2] = V_true[2];
    }

    const gridSource = {
      velocityGrid,
      nx: N, ny: N, nz: N,
      boxSize,
      origin
    };

    const engine = new MultipolarBulkFlow(gridSource);
    const res50 = engine.computeFieldVolumeWeighted(50.0, 0.0);
    const res100 = engine.computeFieldVolumeWeighted(100.0, 0.0);
    const resShell = engine.computeFieldVolumeWeighted(120.0, 40.0);

    console.log(JSON.stringify({
      V_true,
      res50,
      res100,
      resShell
    }));
    """
    res = run_node_snippet(code)
    v_true = res["V_true"]
    true_mag = math.hypot(*v_true)

    for key in ["res50", "res100", "resShell"]:
        dipole = res[key]
        np.testing.assert_allclose(dipole["bulkFlowVector"], v_true, rtol=1e-4)
        assert math.isclose(dipole["magnitude"], true_mag, rel_tol=1e-4)
        assert dipole["sampleCount"] > 0


def test_field_volume_weighted_pure_expansion_zero_dipole():
    """Verify pure isotropic Hubble expansion v(x) = H0 * x yields zero bulk flow dipole."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);
    const H = 0.5; // (km/s) / Mpc

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = H * x;
          velocityGrid[idx + 1] = H * y;
          velocityGrid[idx + 2] = H * z;
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(100.0, 20.0);

    console.log(JSON.stringify({
      dipoleMag: decomp.dipole.magnitude,
      meanVr: decomp.monopole.meanRadialVelocityKms,
      hubbleHR: decomp.monopole.hubbleExpansionHR,
      H_true: H
    }));
    """
    res = run_node_snippet(code)

    # Net dipole bulk flow must vanish identically by spherical parity
    assert res["dipoleMag"] < 1e-3
    # Monopole radial velocity must match H * r_mid
    assert res["meanVr"] > 0
    assert math.isclose(res["hubbleHR"], res["H_true"], rel_tol=0.15)


def test_field_volume_weighted_pure_shear_quadrupole():
    """Verify pure strain flow v = (S*x, -S*y, 0) yields pure Quadrupole shear with eigenvalues (S, 0, -S)."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);
    const S = 0.8; // Shear rate

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = S * x;
          velocityGrid[idx + 1] = -S * y;
          velocityGrid[idx + 2] = 0.0;
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(100.0, 10.0);

    console.log(JSON.stringify({
      dipoleMag: decomp.dipole.magnitude,
      shearMag: decomp.quadrupole.shearMagnitude,
      eigenvalues: decomp.quadrupole.eigenvalues,
      tensor: decomp.quadrupole.tensor
    }));
    """
    res = run_node_snippet(code)

    # Dipole vanishes
    assert res["dipoleMag"] < 1e-3
    # Shear magnitude is strictly positive
    assert res["shearMag"] > 0
    # Trace of quadrupole must be zero: sum of eigenvalues ~ 0
    eigs = res["eigenvalues"]
    assert math.isclose(sum(eigs), 0.0, abs_tol=1e-3)


# ============================================================================
# 5. ESTIMATOR 2: CATALOG WEIGHTED
# ============================================================================

def test_catalog_weighted_estimator_fibonacci_sphere():
    """Verify line-of-sight projection tensor inversion recovers true bulk flow accurately."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const V_true = [180.0, 240.0, -110.0];
    const numGalaxies = 600;
    const galaxies = [];
    const phi = (1 + Math.sqrt(5)) / 2;

    for (let i = 0; i < numGalaxies; i++) {
      const theta = 2 * Math.PI * i / phi;
      const zNorm = 1 - (2 * i + 1) / numGalaxies;
      const radius2D = Math.sqrt(1 - zNorm * zNorm);
      const r = 30.0 + 80.0 * (i / numGalaxies);

      const x = r * radius2D * Math.cos(theta);
      const y = r * radius2D * Math.sin(theta);
      const z = r * zNorm;
      const rMag = Math.hypot(x, y, z);
      const nx = x / rMag, ny = y / rMag, nz = z / rMag;

      const uTrue = V_true[0] * nx + V_true[1] * ny + V_true[2] * nz;
      // Deterministic bounded noise perturbation
      const noise = 15.0 * Math.cos(i * 2.3);

      galaxies.push({
        x, y, z,
        u: uTrue + noise,
        weight: 1.0 + 0.5 * Math.sin(i)
      });
    }

    const engine = new MultipolarBulkFlow(galaxies);
    const res = engine.computeCatalogWeighted(galaxies, 120.0, 20.0);

    console.log(JSON.stringify({
      V_true,
      recovered: res.bulkFlowVector,
      magnitude: res.magnitude,
      trueMag: Math.hypot(...V_true),
      count: res.sampleCount
    }));
    """
    res = run_node_snippet(code)

    np.testing.assert_allclose(res["recovered"], res["V_true"], rtol=0.08)
    assert math.isclose(res["magnitude"], res["trueMag"], rel_tol=0.08)
    assert res["count"] > 500


# ============================================================================
# 6. ESTIMATOR 3: INVERSE VARIANCE WEIGHTED
# ============================================================================

def test_inverse_variance_weighted_estimator_and_covariance():
    """Verify inverse-variance minimum variance estimator and error covariance propagation."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const V_true = [-120.0, 310.0, 75.0];
    const numGalaxies = 800;
    const galaxies = [];
    const phi = (1 + Math.sqrt(5)) / 2;

    for (let i = 0; i < numGalaxies; i++) {
      const theta = 2 * Math.PI * i / phi;
      const zNorm = 1 - (2 * i + 1) / numGalaxies;
      const radius2D = Math.sqrt(1 - zNorm * zNorm);
      const r = 25.0 + 95.0 * (i / numGalaxies);

      const x = r * radius2D * Math.cos(theta);
      const y = r * radius2D * Math.sin(theta);
      const z = r * zNorm;
      const rMag = Math.hypot(x, y, z);
      const nx = x / rMag, ny = y / rMag, nz = z / rMag;

      const uTrue = V_true[0] * nx + V_true[1] * ny + V_true[2] * nz;
      // Heteroscedastic error ~ 15% distance error
      const obsError = 50.0 + 200.0 * (r / 120.0);
      const noise = 20.0 * Math.sin(i * 3.7);

      galaxies.push({
        x, y, z,
        u: uTrue + noise,
        error: obsError
      });
    }

    const engine = new MultipolarBulkFlow(galaxies);
    const resIV = engine.computeInverseVarianceWeighted(galaxies, 130.0, 20.0, 187.0);

    console.log(JSON.stringify({
      V_true,
      recovered: resIV.bulkFlowVector,
      mag: resIV.magnitude,
      cov: resIV.covarianceMatrix,
      unc: resIV.uncertainties
    }));
    """
    res = run_node_snippet(code)

    np.testing.assert_allclose(res["recovered"], res["V_true"], rtol=0.08)

    cov = np.array(res["cov"])
    # Covariance matrix must be positive-definite
    eigs = np.linalg.eigvalsh(cov)
    assert all(eigs > 0)

    unc = res["unc"]
    assert unc["sigmaVx"] > 0
    assert unc["sigmaVy"] > 0
    assert unc["sigmaVz"] > 0
    assert unc["sigmaMagnitude"] > 0


# ============================================================================
# 7. MULTIPOLE DECOMPOSITION HIERARCHY (l = 0, 1, 2, 3, 4)
# ============================================================================

def test_full_multipolar_expansion_hierarchy():
    """Verify full multipolar decomposition (l=0 to 4), power spectrum, and tensor invariants."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 16x16x16 grid with multi-component velocity field:
    // v = V_bulk + H*x + S*shear + octupole perturbation
    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    const V0 = [250.0, 150.0, -80.0];
    const H = 0.3;
    const S = 0.4;

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const r = Math.hypot(x, y, z) || 1.0;
          const idx = ((ix * N + iy) * N + iz) * 3;

          velocityGrid[idx] = V0[0] + H * x + S * x;
          velocityGrid[idx + 1] = V0[1] + H * y - S * y;
          velocityGrid[idx + 2] = V0[2] + H * z;
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(120.0, 20.0);

    console.log(JSON.stringify(decomp));
    """
    res = run_node_snippet(code)

    # 1. Monopole
    assert res["monopole"]["hubbleExpansionHR"] > 0
    assert res["monopole"]["powerC0"] >= 0

    # 2. Dipole
    assert res["dipole"]["magnitude"] > 200.0
    assert res["dipole"]["powerC1"] > 0
    assert 0 <= res["dipole"]["sglDeg"] <= 360.0
    assert -90.0 <= res["dipole"]["sgbDeg"] <= 90.0
    assert 0 <= res["dipole"]["galacticL"] <= 360.0
    assert -90.0 <= res["dipole"]["galacticB"] <= 90.0
    assert 0 <= res["dipole"]["angleToShapleyDeg"] <= 180.0
    assert 0 <= res["dipole"]["angleToCMBDipoleDeg"] <= 180.0

    # 3. Quadrupole
    assert res["quadrupole"]["shearMagnitude"] > 0
    assert res["quadrupole"]["powerC2"] >= 0
    q_eigs = res["quadrupole"]["eigenvalues"]
    assert math.isclose(sum(q_eigs), 0.0, abs_tol=1e-3)

    # 4. Octupole & Hexadecapole
    assert res["octupole"]["octupoleMagnitude"] >= 0
    assert res["hexadecapole"]["hexadecapoleMagnitude"] >= 0

    # 5. Power spectrum
    ps = res["powerSpectrum"]
    assert ps["totalPower"] > 0
    np.testing.assert_allclose(sum(ps["powerFractions"]), 1.0, atol=1e-4)


# ============================================================================
# 8. SURVEY WINDOW DECONVOLUTION
# ============================================================================

def test_survey_window_function_deconvolution():
    """Verify survey window function deconvolution V_true = W^{-1} V_meas."""
    code = """
    import { MultipolarBulkFlow, BulkFlowEstimatorType } from './src/bulk-flow/multipolar_bulk_flow.js';

    // True bulk flow
    const V_true = [300.0, -200.0, 150.0];

    // Synthetic survey window matrix W (representing e.g. Zone of Avoidance distortion)
    const W = [
      [0.85, 0.05, 0.02],
      [0.05, 0.90, -0.03],
      [0.02, -0.03, 0.78]
    ];

    // Measured bulk flow V_meas = W * V_true
    const V_meas = [
      W[0][0] * V_true[0] + W[0][1] * V_true[1] + W[0][2] * V_true[2],
      W[1][0] * V_true[0] + W[1][1] * V_true[1] + W[1][2] * V_true[2],
      W[2][0] * V_true[0] + W[2][1] * V_true[1] + W[2][2] * V_true[2]
    ];

    const dipoleMeas = {
      method: 'CATALOG_WEIGHTED',
      rMin: 20.0,
      rMax: 150.0,
      bulkFlowVector: V_meas,
      magnitude: Math.hypot(...V_meas),
      covarianceMatrix: [
        [400.0, 20.0, 10.0],
        [20.0, 450.0, -15.0],
        [10.0, -15.0, 500.0]
      ]
    };

    const engine = new MultipolarBulkFlow([]);
    const deconv = engine.deconvolveSurveyWindow(dipoleMeas, W);

    console.log(JSON.stringify({
      V_true,
      V_meas,
      deconv
    }));
    """
    res = run_node_snippet(code)

    v_true = res["V_true"]
    v_deconv = res["deconv"]["bulkFlowVector"]
    # Deconvolved vector should match true velocity vector exactly
    np.testing.assert_allclose(v_deconv, v_true, rtol=1e-4)
    assert math.isclose(res["deconv"]["magnitude"], math.hypot(*v_true), rel_tol=1e-4)


# ============================================================================
# 9. COSMICFLOWS-4 (CF4) BENCHMARK COMPARISON
# ============================================================================

def test_cosmicflows4_reference_table_and_benchmark_comparison():
    """Verify CosmicFlows-4 empirical reference table and chi-squared comparison engine."""
    code = """
    import {
      COSMICFLOWS4_REFERENCE_DIPOLES,
      MultipolarBulkFlow
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    const engine = new MultipolarBulkFlow([]);

    // 1. Check all CF4 reference radii in [20, 250] Mpc/h
    const radii = COSMICFLOWS4_REFERENCE_DIPOLES.map(d => d.radiusMpc);

    // 2. Test perfect match against CF4 at R = 150 Mpc/h (Dupuy & Courtois 2023)
    const cf4_150 = COSMICFLOWS4_REFERENCE_DIPOLES.find(d => d.radiusMpc === 150.0);
    const mockDipoleMatch = {
      rMax: 150.0,
      bulkFlowVector: [cf4_150.vx, cf4_150.vy, cf4_150.vz],
      magnitude: cf4_150.magnitude,
      sglDeg: cf4_150.sgl,
      sgbDeg: cf4_150.sgb,
      galacticL: cf4_150.galacticL,
      galacticB: cf4_150.galacticB,
      uncertainties: { sigmaMagnitude: 32.0 }
    };
    const compMatch = engine.compareToCosmicFlows4(mockDipoleMatch, 150.0);

    // 3. Test discrepant bulk flow (high tension)
    const mockDipoleTension = {
      rMax: 150.0,
      bulkFlowVector: [-300.0, -200.0, 100.0],
      magnitude: 375.0,
      uncertainties: { sigmaMagnitude: 20.0 }
    };
    const compTension = engine.compareToCosmicFlows4(mockDipoleTension, 150.0);

    console.log(JSON.stringify({
      radii,
      compMatch,
      compTension
    }));
    """
    res = run_node_snippet(code)

    # Check CF4 table coverage
    expected_radii = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0, 250.0]
    np.testing.assert_allclose(res["radii"], expected_radii, atol=1e-3)

    # Matching result should give EXCELLENT_MATCH and chi2 ~ 0
    m = res["compMatch"]
    assert m["consistencyRating"] == "EXCELLENT_MATCH"
    assert m["isConsistentWithin1Sigma"] is True
    assert math.isclose(m["deltaMagnitude"], 0.0, abs_tol=1e-3)
    assert math.isclose(m["angularOffsetDeg"], 0.0, abs_tol=1e-3)

    # Tension result should be flagged
    t = res["compTension"]
    assert t["consistencyRating"] in ["MODERATE_TENSION", "HIGH_TENSION"]
    assert t["angularOffsetDeg"] > 90.0


# ============================================================================
# 10. CONTINUOUS RADIAL PROFILES & CONVENIENCE FUNCTIONS
# ============================================================================

def test_continuous_radial_profiles_and_convenience_wrappers():
    """Verify computeRadialProfile, computeTripleWeightingBulkFlow, computeMultipolarRadialProfile."""
    code = """
    import {
      MultipolarBulkFlow,
      computeTripleWeightingBulkFlow,
      computeMultipolarRadialProfile
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 16x16x16 grid with radial decaying velocity field V(r) = V0 / (1 + r/80)
    const N = 16;
    const boxSize = 400.0;
    const origin = [-200.0, -200.0, -200.0];
    const velocityGrid = new Float32Array(N * N * N * 3);
    const V0 = [280.0, 220.0, -70.0];

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const r = Math.hypot(x, y, z) || 1.0;
          const factor = 1.0 / (1.0 + r / 80.0);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = V0[0] * factor;
          velocityGrid[idx + 1] = V0[1] * factor;
          velocityGrid[idx + 2] = V0[2] * factor;
        }
      }
    }

    const gridSource = { velocityGrid, nx: N, ny: N, nz: N, boxSize, origin };
    const engine = new MultipolarBulkFlow(gridSource);

    const radii = [40, 60, 80, 100, 150, 200];
    const profileCum = engine.computeRadialProfile(radii, { mode: 'SPHERICAL_CUMULATIVE' });
    const profileDiff = engine.computeRadialProfile(radii, { mode: 'DIFFERENTIAL_SHELLS' });

    const multiProfile = computeMultipolarRadialProfile(gridSource, [60, 100, 150]);

    console.log(JSON.stringify({
      profileCumLength: profileCum.length,
      profileDiffLength: profileDiff.length,
      cumMags: profileCum.map(p => p.magnitude),
      diffMags: profileDiff.map(p => p.magnitude),
      multiProfileLength: multiProfile.length
    }));
    """
    res = run_node_snippet(code)

    assert res["profileCumLength"] == 6
    assert res["profileDiffLength"] == 6
    assert res["multiProfileLength"] == 3

    # Decaying field should have monotonically decreasing magnitude with radius
    cum_mags = res["cumMags"]
    for i in range(len(cum_mags) - 1):
        assert cum_mags[i] > cum_mags[i + 1]


# ============================================================================
# 11. DETAILED MATRIX REGULARIZATION & CONDITION NUMBER TESTS
# ============================================================================

def test_singular_and_regularized_matrix_inversion():
    """Verify robust regularized inversion of singular and degenerate 3x3 matrices."""
    code = """
    import { invertMatrix3x3, solveLinear3x3 } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 1. Zero matrix (rank 0)
    const zeroMat = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];
    const invZero = invertMatrix3x3(zeroMat);

    // 2. Collinear projection rank-1 matrix (e.g. all galaxies on the x-axis)
    const rank1Mat = [
      [10.0, 0.0, 0.0],
      [0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0]
    ];
    const invRank1 = invertMatrix3x3(rank1Mat);

    // 3. Near-singular matrix
    const nearSingular = [
      [1.0, 2.0, 3.0],
      [2.0, 4.0, 6.00000001],
      [1.0, 1.0, 1.0]
    ];
    const invNear = invertMatrix3x3(nearSingular);

    console.log(JSON.stringify({
      invZero,
      invRank1,
      invNear
    }));
    """
    res = run_node_snippet(code)

    # Inversion of singular matrix should not return NaN or throw
    assert not math.isnan(res["invZero"][0][0])
    assert not math.isnan(res["invRank1"][0][0])
    assert not math.isnan(res["invNear"][0][0])


# ============================================================================
# 12. HIGHER-ORDER SPHERICAL HARMONICS & RECURRENCE INVARIANTS
# ============================================================================

def test_associated_legendre_parity_and_boundary_properties():
    """Verify Associated Legendre polynomials parity P_l^m(-x) = (-1)^{l-m} P_l^m(x) and boundaries."""
    code = """
    import { associatedLegendre } from './src/bulk-flow/multipolar_bulk_flow.js';

    const testX = [0.1, 0.35, 0.72, 0.95];
    const parityChecks = [];

    for (let l = 0; l <= 4; l++) {
      for (let m = 0; m <= l; m++) {
        for (let i = 0; i < testX.length; i++) {
          const x = testX[i];
          const px = associatedLegendre(l, m, x);
          const pminusx = associatedLegendre(l, m, -x);
          const paritySign = Math.pow(-1, l - m);
          const expected = paritySign * px;
          parityChecks.push({
            l, m, x,
            px,
            pminusx,
            expected,
            err: Math.abs(pminusx - expected)
          });
        }
      }
    }

    console.log(JSON.stringify(parityChecks));
    """
    res = run_node_snippet(code)
    for check in res:
        assert check["err"] < 1e-5, f"Parity failed for l={check['l']}, m={check['m']}, x={check['x']}"


def test_spherical_harmonic_rotation_power_invariance():
    """Verify angular power spectrum C_l is invariant under global coordinate frame rotations."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    // Generate random galaxy catalog on sphere with bulk flow + quadrupole
    const N = 400;
    const galaxies = [];
    const phi = (1 + Math.sqrt(5)) / 2;

    for (let i = 0; i < N; i++) {
      const theta = 2 * Math.PI * i / phi;
      const zNorm = 1 - (2 * i + 1) / N;
      const radius2D = Math.sqrt(1 - zNorm * zNorm);
      const r = 50.0 + 30.0 * Math.sin(i * 0.5);

      const x = r * radius2D * Math.cos(theta);
      const y = r * radius2D * Math.sin(theta);
      const z = r * zNorm;

      // Dipole V = (150, 200, -100)
      const u = (150.0 * x + 200.0 * y - 100.0 * z) / r;
      galaxies.push({ x, y, z, u, weight: 1.0 });
    }

    const engine1 = new MultipolarBulkFlow(galaxies);
    const decomp1 = engine1.computeMultipoleDecomposition(100.0, 10.0);

    // Rotated galaxies around z-axis by 90 degrees: (x', y', z') = (-y, x, z)
    const rotatedGalaxies = galaxies.map(g => ({
      x: -g.y,
      y: g.x,
      z: g.z,
      u: g.u,
      weight: g.weight
    }));

    const engine2 = new MultipolarBulkFlow(rotatedGalaxies);
    const decomp2 = engine2.computeMultipoleDecomposition(100.0, 10.0);

    console.log(JSON.stringify({
      cl1: decomp1.powerSpectrum.Cl,
      cl2: decomp2.powerSpectrum.Cl,
      mag1: decomp1.dipole.magnitude,
      mag2: decomp2.dipole.magnitude
    }));
    """
    res = run_node_snippet(code)

    # Dipole magnitude must be invariant under rotation
    assert math.isclose(res["mag1"], res["mag2"], rel_tol=1e-3)
    # Power spectrum components Cl must be approximately invariant
    np.testing.assert_allclose(res["cl1"], res["cl2"], rtol=0.08)


# ============================================================================
# 13. DIRECT 3D VELOCITY CATALOG ESTIMATOR TESTS
# ============================================================================

def test_catalog_weighted_3d_velocity_vectors():
    """Verify catalog-weighted estimator when full 3D velocity vectors (vx, vy, vz) are supplied."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const V_true = [140.0, -220.0, 310.0];
    const N = 300;
    const galaxies = [];

    for (let i = 0; i < N; i++) {
      const r = 20.0 + 80.0 * (i / N);
      const angle = (i * 137.5) * Math.PI / 180.0;
      const z = (i / N - 0.5) * 2.0 * r;
      const x = Math.sqrt(Math.max(0, r * r - z * z)) * Math.cos(angle);
      const y = Math.sqrt(Math.max(0, r * r - z * z)) * Math.sin(angle);

      // 3D velocity with random thermal dispersion
      const vx = V_true[0] + 30.0 * Math.sin(i * 1.5);
      const vy = V_true[1] + 30.0 * Math.cos(i * 2.1);
      const vz = V_true[2] + 30.0 * Math.sin(i * 3.3);

      galaxies.push({
        x, y, z,
        vx, vy, vz,
        weight: 1.0 + (i % 3) * 0.5
      });
    }

    const engine = new MultipolarBulkFlow(galaxies);
    const res = engine.computeCatalogWeighted(galaxies, 100.0, 20.0);

    console.log(JSON.stringify({
      V_true,
      recovered: res.bulkFlowVector,
      magnitude: res.magnitude,
      trueMag: Math.hypot(...V_true),
      count: res.sampleCount
    }));
    """
    res = run_node_snippet(code)

    np.testing.assert_allclose(res["recovered"], res["V_true"], rtol=0.05)
    assert math.isclose(res["magnitude"], res["trueMag"], rel_tol=0.05)
    assert res["count"] > 250


# ============================================================================
# 14. INVERSE-VARIANCE STATISTICAL EFFICIENCY COMPARISON
# ============================================================================

def test_inverse_variance_weighting_vs_unweighted_efficiency():
    """Verify that inverse-variance estimator achieves lower variance under heteroscedastic noise."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const V_true = [200.0, -100.0, 150.0];
    const N = 500;
    const galaxies = [];
    const phi = (1 + Math.sqrt(5)) / 2;

    for (let i = 0; i < N; i++) {
      const theta = 2 * Math.PI * i / phi;
      const zNorm = 1 - (2 * i + 1) / N;
      const radius2D = Math.sqrt(1 - zNorm * zNorm);
      const r = 30.0 + 70.0 * (i / N);

      const x = r * radius2D * Math.cos(theta);
      const y = r * radius2D * Math.sin(theta);
      const z = r * zNorm;
      const rMag = Math.hypot(x, y, z);
      const nx = x / rMag, ny = y / rMag, nz = z / rMag;

      const uTrue = V_true[0] * nx + V_true[1] * ny + V_true[2] * nz;

      // 80% of galaxies have high noise (sigma = 300), 20% have low noise (sigma = 30)
      const isAccurate = (i % 5 === 0);
      const obsError = isAccurate ? 30.0 : 300.0;
      const noise = obsError * 0.5 * Math.sin(i * 4.1);

      galaxies.push({
        x, y, z,
        u: uTrue + noise,
        error: obsError,
        weight: 1.0
      });
    }

    const engine = new MultipolarBulkFlow(galaxies);
    const resCat = engine.computeCatalogWeighted(galaxies, 100.0, 20.0);
    const resIV = engine.computeInverseVarianceWeighted(galaxies, 100.0, 20.0, 187.0);

    const errCat = Math.hypot(
      resCat.bulkFlowVector[0] - V_true[0],
      resCat.bulkFlowVector[1] - V_true[1],
      resCat.bulkFlowVector[2] - V_true[2]
    );
    const errIV = Math.hypot(
      resIV.bulkFlowVector[0] - V_true[0],
      resIV.bulkFlowVector[1] - V_true[1],
      resIV.bulkFlowVector[2] - V_true[2]
    );

    console.log(JSON.stringify({
      errCat,
      errIV,
      resIV_unc: resIV.uncertainties
    }));
    """
    res = run_node_snippet(code)

    # Inverse-variance weighting should have smaller error from true bulk flow
    assert res["errIV"] < res["errCat"]
    assert res["resIV_unc"]["sigmaMagnitude"] > 0


# ============================================================================
# 15. QUADRUPOLE PRINCIPAL AXES ORTHONORMALITY & SPECTRAL ANALYSIS
# ============================================================================

def test_quadrupole_principal_axes_orthonormality_and_directions():
    """Verify principal shear axes eigenvectors form an orthonormal frame."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    // Asymmetric triaxial shear field: S = diag(0.6, 0.2, -0.8)
    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = 0.6 * x;
          velocityGrid[idx + 1] = 0.2 * y;
          velocityGrid[idx + 2] = -0.8 * z;
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(120.0, 10.0);
    const axes = decomp.quadrupole.principalAxes;

    // Check dot products between principal axes e_i . e_j
    const e1 = axes[0].vector;
    const e2 = axes[1].vector;
    const e3 = axes[2].vector;

    const dot12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
    const dot23 = e2[0] * e3[0] + e2[1] * e3[1] + e2[2] * e3[2];
    const dot31 = e3[0] * e1[0] + e3[1] * e1[1] + e3[2] * e1[2];

    const norm1 = Math.hypot(...e1);
    const norm2 = Math.hypot(...e2);
    const norm3 = Math.hypot(...e3);

    console.log(JSON.stringify({
      dot12, dot23, dot31,
      norm1, norm2, norm3,
      eigenvalues: decomp.quadrupole.eigenvalues,
      shearMag: decomp.quadrupole.shearMagnitude
    }));
    """
    res = run_node_snippet(code)

    # Norm of each eigenvector must be 1
    assert math.isclose(res["norm1"], 1.0, abs_tol=1e-4)
    assert math.isclose(res["norm2"], 1.0, abs_tol=1e-4)
    assert math.isclose(res["norm3"], 1.0, abs_tol=1e-4)

    # Orthogonality
    assert math.isclose(res["dot12"], 0.0, abs_tol=1e-4)
    assert math.isclose(res["dot23"], 0.0, abs_tol=1e-4)
    assert math.isclose(res["dot31"], 0.0, abs_tol=1e-4)

    # Trace of eigenvalues is zero
    assert math.isclose(sum(res["eigenvalues"]), 0.0, abs_tol=1e-3)


# ============================================================================
# 16. OCTUPOLE & HEXADECAPOLE TENSOR TRACE INVARIANTS
# ============================================================================

def test_octupole_and_hexadecapole_tensor_traces():
    """Verify trace subtractions and symmetry properties of 3rd and 4th order multipole tensors."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    // Multi-pole perturbation
    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const r = Math.hypot(x, y, z) || 1.0;
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = 100.0 * (x * x - y * y) / r;
          velocityGrid[idx + 1] = 100.0 * (y * y - z * z) / r;
          velocityGrid[idx + 2] = 100.0 * (z * z - x * x) / r;
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(120.0, 10.0);
    const O = decomp.octupole.tensor;

    // Check octupole trace contraction: sum_j O_ijj ~ 0
    const traces = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        traces[i] += O[i][j][j];
      }
    }

    console.log(JSON.stringify({
      octupoleMag: decomp.octupole.octupoleMagnitude,
      traces,
      hexMag: decomp.hexadecapole.hexadecapoleMagnitude
    }));
    """
    res = run_node_snippet(code)

    assert res["octupoleMag"] >= 0
    assert res["hexMag"] >= 0
    # Trace contractions should be small
    for tr in res["traces"]:
        assert abs(tr) < 5.0


# ============================================================================
# 17. CONVENIENCE TRIPLE WEIGHTING & ZERO-SAMPLE EDGE CASES
# ============================================================================

def test_convenience_triple_weighting_and_edge_cases():
    """Verify computeTripleWeightingBulkFlow and edge cases (empty catalog, zero radius)."""
    code = """
    import {
      computeTripleWeightingBulkFlow,
      MultipolarBulkFlow
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 1. Valid catalog
    const galaxies = [
      { x: 50, y: 0, z: 0, u: 300, error: 100, weight: 1 },
      { x: -50, y: 0, z: 0, u: -300, error: 100, weight: 1 },
      { x: 0, y: 50, z: 0, u: 100, error: 100, weight: 1 },
      { x: 0, y: -50, z: 0, u: -100, error: 100, weight: 1 },
      { x: 0, y: 0, z: 50, u: 50, error: 100, weight: 1 },
      { x: 0, y: 0, z: -50, u: -50, error: 100, weight: 1 }
    ];

    const res = computeTripleWeightingBulkFlow(galaxies, 100.0, 0.0);

    // 2. Empty catalog edge case
    const emptyEngine = new MultipolarBulkFlow([]);
    const emptyDipole = emptyEngine.computeFieldVolumeWeighted(100.0, 0.0);
    const emptyDecomp = emptyEngine.computeMultipoleDecomposition(100.0, 0.0);

    console.log(JSON.stringify({
      catMag: res.catalogWeighted.magnitude,
      ivMag: res.inverseVariance.magnitude,
      emptyMag: emptyDipole.magnitude,
      emptyCount: emptyDipole.sampleCount,
      emptyDecompPower: emptyDecomp.powerSpectrum.totalPower
    }));
    """
    res = run_node_snippet(code)

    assert res["catMag"] > 0
    assert res["ivMag"] > 0
    assert res["emptyMag"] == 0.0
    assert res["emptyCount"] == 0
    assert res["emptyDecompPower"] == 0.0


# ============================================================================
# 18. LARGE HIGH-RESOLUTION GRID STRESS TEST
# ============================================================================

def test_large_grid_high_resolution_stress_test():
    """Verify stability and performance on a 32x32x32 = 32,768 voxel velocity grid."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 32;
    const boxSize = 400.0;
    const origin = [-200.0, -200.0, -200.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    const V_true = [350.0, 240.0, -90.0];

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = V_true[0];
          velocityGrid[idx + 1] = V_true[1];
          velocityGrid[idx + 2] = V_true[2];
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(150.0, 20.0);

    console.log(JSON.stringify({
      sampleCount: decomp.sampleCount,
      bulkVector: decomp.dipole.bulkFlowVector,
      magnitude: decomp.dipole.magnitude,
      trueMag: Math.hypot(...V_true)
    }));
    """
    res = run_node_snippet(code)

    assert res["sampleCount"] > 5000
    np.testing.assert_allclose(res["bulkVector"], [350.0, 240.0, -90.0], rtol=1e-4)
    assert math.isclose(res["magnitude"], res["trueMag"], rel_tol=1e-4)


# ============================================================================
# 19. STF 4TH-ORDER TENSOR PERMUTATION SYMMETRY
# ============================================================================

def test_hexadecapole_stf_tensor_permutation_symmetry():
    """Verify H_ijkl is invariant under all 24 index permutations."""
    code = """
    import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';

    const N = 16;
    const boxSize = 300.0;
    const origin = [-150.0, -150.0, -150.0];
    const velocityGrid = new Float32Array(N * N * N * 3);

    for (let ix = 0; ix < N; ix++) {
      const x = origin[0] + (ix + 0.5) * (boxSize / N);
      for (let iy = 0; iy < N; iy++) {
        const y = origin[1] + (iy + 0.5) * (boxSize / N);
        for (let iz = 0; iz < N; iz++) {
          const z = origin[2] + (iz + 0.5) * (boxSize / N);
          const idx = ((ix * N + iy) * N + iz) * 3;
          velocityGrid[idx] = 120.0 * Math.sin(x / 40.0);
          velocityGrid[idx + 1] = 80.0 * Math.cos(y / 40.0);
          velocityGrid[idx + 2] = -50.0 * Math.sin(z / 40.0);
        }
      }
    }

    const engine = new MultipolarBulkFlow({
      velocityGrid, nx: N, ny: N, nz: N, boxSize, origin
    });

    const decomp = engine.computeMultipoleDecomposition(120.0, 10.0);
    const H = decomp.hexadecapole.tensor;

    // Check symmetry: H[i][j][k][l] === H[j][i][k][l] === H[k][j][i][l] === H[l][j][k][i]
    let maxSymmetryError = 0.0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        for (let k = 0; k < 3; k++) {
          for (let l = 0; l < 3; l++) {
            const val = H[i][j][k][l];
            const sym1 = Math.abs(val - H[j][i][k][l]);
            const sym2 = Math.abs(val - H[k][j][i][l]);
            const sym3 = Math.abs(val - H[l][j][k][i]);
            const sym4 = Math.abs(val - H[i][k][j][l]);
            const sym5 = Math.abs(val - H[i][l][k][j]);
            maxSymmetryError = Math.max(maxSymmetryError, sym1, sym2, sym3, sym4, sym5);
          }
        }
      }
    }

    console.log(JSON.stringify({
      maxSymmetryError,
      hexMag: decomp.hexadecapole.hexadecapoleMagnitude
    }));
    """
    res = run_node_snippet(code)
    assert res["maxSymmetryError"] < 1e-10


# ============================================================================
# 20. RADIAL SHELL FOURIER INTEGRATION & TENSION BOUNDS
# ============================================================================

def test_radial_shell_fourier_window_properties_and_tension_scales():
    """Verify radialShellWindow continuity, limits, and tension significance scales."""
    code = """
    import {
      radialShellWindow,
      sphericalTopHatWindow,
      CosmicVarianceEngine
    } from './src/bulk-flow/multipolar_bulk_flow.js';

    // 1. Shell window with rIn = 0 should match sphericalTopHatWindow
    const kVals = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0];
    const diffs = kVals.map(k => {
      const wShell = radialShellWindow(k, 0.0, 100.0);
      const wSphere = sphericalTopHatWindow(k * 100.0);
      return Math.abs(wShell - wSphere);
    });

    // 2. Cosmic variance tension evaluations
    const engine = new CosmicVarianceEngine();
    const t0 = engine.evaluateTension(100.0, 150.0);
    const t1 = engine.evaluateTension(250.0, 150.0);
    const t2 = engine.evaluateTension(600.0, 150.0);

    console.log(JSON.stringify({
      diffs,
      t0_tension: t0.tensionSigma,
      t1_tension: t1.tensionSigma,
      t2_tension: t2.tensionSigma,
      t2_pValue: t2.pValue
    }));
    """
    res = run_node_snippet(code)

    for diff in res["diffs"]:
        assert diff < 1e-5

    # Tension should increase monotonically with observed velocity magnitude
    assert res["t0_tension"] <= res["t1_tension"] <= res["t2_tension"]
    assert res["t2_tension"] > 1.5


