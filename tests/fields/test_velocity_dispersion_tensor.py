"""
tests/fields/test_velocity_dispersion_tensor.py
================================================
Exhaustive automated Pytest test suite for the 3D Velocity Dispersion Tensor,
Cosmic Thermal and Kinetic Pressure Fields, Spherical Velocity Anisotropy Parameter beta(r),
Jeans Equation Cosmological Mass Estimator, and Supercluster Multiscale Virial Diagnostics.

Test Categories:
----------------
1. Mathematical Invariants, Tensor Properties & Jacobi Eigensolver (Exact 3x3 Diagonalization, Cardano roots, Orthonormality, Trace, Triaxiality, Invariants)
2. Spatial Smoothing Kernels & Convolution Mechanics (Gaussian, Top-Hat, SPH Cubic Spline, Wendland C^4, Laminar flow, Linear shear)
3. Discrete Particle Catalog Dispersion Estimator & Weighting
4. Spherical Halo Coordinate Transformation & Anisotropy Parameter beta(x) (Isotropic beta=0, Pure Radial beta=1, Circular beta->-inf, Symmetrized gamma_beta, Osipkov-Merritt)
5. Radial Anisotropy Profile Shell Binning & Regime Classification
6. Cosmic Thermal and Kinetic Pressure Tensors & Sound Speed (P_ij = rho * sigma_ij^2, Isotropic P, Sound Speed c_s, Gas Temperature Proxy, Mach Number)
7. Spherical Jeans Equation Cosmological Mass Estimator (Logarithmic Slopes, Singular Isothermal Sphere SIS, Hernquist Profile, NFW Profile, Anisotropy Bias, Monotonicity)
8. Supercluster Multiscale Virial Ratio 2K / |W| (Center of Mass Bulk Velocity, Internal Dispersion Kinetic Energy, Gravitational Potential, Virial Equilibrium Classification)
9. W3C PROV-O Cryptographic Lineage & NIST SHA-256 Audit Manifest
10. Robustness, Singularity Handling & Boundary Edge Cases

Author: Scientific Computational Cosmology Engineer
"""

import json
import math
import os
import subprocess
import numpy as np
import pytest

NODE_EXEC = "node"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

# Cosmological constants
G_COSMO = 4.3009172706e-9  # (km/s)^2 * Mpc / M_sun
RHO_CRIT_0 = 2.77536627e11  # (M_sun/h) / (Mpc/h)^3 for h=1.0
OMEGA_M = 0.3111


def run_node_snippet(code: str) -> dict:
    """Executes an inline ES6 Node.js snippet and returns parsed JSON output."""
    wrapped_code = f"""
    {code}
    """
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", wrapped_code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Node execution failed with code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as err:
        raise RuntimeError(f"Failed to parse JSON from node output:\n{proc.stdout}\nError: {err}")


# ============================================================================
# 1. MATHEMATICAL INVARIANTS, TENSOR PROPERTIES & JACOBI EIGENSOLVER
# ============================================================================

class TestMathematicalInvariantsAndJacobiEigensolver:
    """Verify 3x3 symmetric matrix diagonalization, Cardano eigensolver, invariants, and tensor properties."""

    def test_exact_jacobi_diagonalization_diagonal_matrix(self):
        """Diagonal matrix should yield exact diagonal entries as eigenvalues and identity as eigenvectors."""
        code = """
        import { diagonalizeSymmetric3x3 } from './src/fields/velocity_dispersion_tensor.js';
        const mat = new Float64Array([
          500.0,   0.0,   0.0,
            0.0, 300.0,   0.0,
            0.0,   0.0, 100.0
        ]);
        const res = diagonalizeSymmetric3x3(mat);
        console.log(JSON.stringify({
          eigenvalues: Array.from(res.eigenvalues),
          eigenvectors: Array.from(res.eigenvectors),
          sweeps: res.sweepCount
        }));
        """
        out = run_node_snippet(code)
        evals = out["eigenvalues"]
        evecs = np.array(out["eigenvectors"]).reshape(3, 3)

        assert pytest.approx(evals[0], rel=1e-12) == 500.0
        assert pytest.approx(evals[1], rel=1e-12) == 300.0
        assert pytest.approx(evals[2], rel=1e-12) == 100.0

        # Eigenvectors should be standard basis vectors
        assert pytest.approx(abs(evecs[0, 0]), rel=1e-12) == 1.0
        assert pytest.approx(abs(evecs[1, 1]), rel=1e-12) == 1.0
        assert pytest.approx(abs(evecs[2, 2]), rel=1e-12) == 1.0

    def test_jacobi_vs_numpy_eigh_general_symmetric(self):
        """General full 3x3 symmetric dispersion tensor: Jacobi matches NumPy eigh with machine precision."""
        np_mat = np.array([
          [450.0,  120.0,  -65.0],
          [120.0,  320.0,   85.0],
          [-65.0,   85.0,  210.0]
        ], dtype=np.float64)

        # Expected eigenvalues in descending order
        np_evals, np_evecs = np.linalg.eigh(np_mat)
        np_evals_sorted = np_evals[::-1]

        code = f"""
        import {{ diagonalizeSymmetric3x3 }} from './src/fields/velocity_dispersion_tensor.js';
        const mat = new Float64Array({json.dumps(np_mat.flatten().tolist())});
        const res = diagonalizeSymmetric3x3(mat);
        console.log(JSON.stringify({{
          eigenvalues: Array.from(res.eigenvalues),
          eigenvectors: Array.from(res.eigenvectors)
        }}));
        """
        out = run_node_snippet(code)
        evals = np.array(out["eigenvalues"])
        evecs = np.array(out["eigenvectors"]).reshape(3, 3)

        np.testing.assert_allclose(evals, np_evals_sorted, rtol=1e-12)

        # Verify orthonormal eigenvectors: V^T V = I
        identity = np.eye(3)
        np.testing.assert_allclose(evecs.T @ evecs, identity, atol=1e-12)

        # Verify right-handed determinant det(V) = +1
        assert pytest.approx(np.linalg.det(evecs), rel=1e-12) == 1.0

        # Verify spectral decomposition: V * Lambda * V^T = A
        recon = evecs @ np.diag(evals) @ evecs.T
        np.testing.assert_allclose(recon, np_mat, rtol=1e-12)

    def test_cardano_eigensolver_consistency(self):
        """Cardano analytical eigensolver produces identical eigenvalues to Jacobi method."""
        code = """
        import { cardanoEigenvaluesSymmetric3x3, diagonalizeSymmetric3x3 } from './src/fields/velocity_dispersion_tensor.js';
        const mat = new Float64Array([
          620.0, -150.0,  80.0,
         -150.0,  410.0, -90.0,
           80.0,  -90.0, 250.0
        ]);
        const cardanoEvals = cardanoEigenvaluesSymmetric3x3(mat);
        const jacobiRes = diagonalizeSymmetric3x3(mat);
        console.log(JSON.stringify({
          cardano: Array.from(cardanoEvals),
          jacobi: Array.from(jacobiRes.eigenvalues)
        }));
        """
        out = run_node_snippet(code)
        np.testing.assert_allclose(out["cardano"], out["jacobi"], rtol=1e-11)

    def test_dispersion_tensor_point_invariants(self):
        """Verify trace, 1D/3D velocity dispersions, determinant, second invariant, and triaxiality."""
        sxx, syy, szz = 900.0, 400.0, 100.0
        sxy, sxz, syz = 50.0, -30.0, 20.0

        code = f"""
        import {{ DispersionTensorPoint }} from './src/fields/velocity_dispersion_tensor.js';
        const pt = new DispersionTensorPoint({sxx}, {sxy}, {sxz}, {syy}, {syz}, {szz}, [10.0, -5.0, 15.0]);
        console.log(JSON.stringify({{
          trace: pt.trace,
          sigma1D: pt.sigma1D,
          sigma3D: pt.sigma3D,
          determinant: pt.determinant,
          secondInvariant: pt.secondInvariant,
          triaxiality: pt.triaxiality,
          ellipticity: pt.ellipticity,
          prolateness: pt.prolateness,
          conditionNumber: pt.conditionNumber,
          isPSD: pt.isPositiveSemiDefinite()
        }}));
        """
        out = run_node_snippet(code)

        expected_trace = sxx + syy + szz
        expected_sig1d = math.sqrt(expected_trace / 3.0)
        expected_sig3d = math.sqrt(expected_trace)

        assert pytest.approx(out["trace"], rel=1e-12) == expected_trace
        assert pytest.approx(out["sigma1D"], rel=1e-12) == expected_sig1d
        assert pytest.approx(out["sigma3D"], rel=1e-12) == expected_sig3d
        assert out["isPSD"] is True
        assert 0.0 <= out["triaxiality"] <= 1.0

        # Determinant check with numpy
        mat = np.array([
            [sxx, sxy, sxz],
            [sxy, syy, syz],
            [sxz, syz, szz]
        ])
        assert pytest.approx(out["determinant"], rel=1e-11) == np.linalg.det(mat)

    def test_triaxiality_prolate_oblate_limits(self):
        """Triaxiality T=1 for prolate (lambda_1 > lambda_2 = lambda_3), T=0 for oblate (lambda_1 = lambda_2 > lambda_3)."""
        code = """
        import { DispersionTensorPoint } from './src/fields/velocity_dispersion_tensor.js';
        // Prolate: lambda_1 = 400, lambda_2 = 100, lambda_3 = 100 -> T = (400-100)/(400-100) = 1.0
        const prolate = new DispersionTensorPoint(400.0, 0.0, 0.0, 100.0, 0.0, 100.0);
        // Oblate: lambda_1 = 400, lambda_2 = 400, lambda_3 = 100 -> T = (400-400)/(400-100) = 0.0
        const oblate = new DispersionTensorPoint(400.0, 0.0, 0.0, 400.0, 0.0, 100.0);
        // Spherical: lambda_1 = lambda_2 = lambda_3 = 250 -> T = 0.5
        const sphere = new DispersionTensorPoint(250.0, 0.0, 0.0, 250.0, 0.0, 250.0);

        console.log(JSON.stringify({
          tProlate: prolate.triaxiality,
          tOblate: oblate.triaxiality,
          tSphere: sphere.triaxiality
        }));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["tProlate"], abs=1e-7) == 1.0
        assert pytest.approx(out["tOblate"], abs=1e-7) == 0.0
        assert pytest.approx(out["tSphere"], abs=1e-7) == 0.5


# ============================================================================
# 2. SPATIAL SMOOTHING KERNELS & CONVOLUTION MECHANICS
# ============================================================================

class TestSpatialSmoothingKernelsAndMechanics:
    """Verify kernel weights, normalization, support, and spatial convolution dispersion estimation."""

    @pytest.mark.parametrize("kernel_type", ["gaussian", "tophat", "cubic_spline", "wendland_c4"])
    def test_kernel_normalization_3d_numerical_integration(self, kernel_type):
        """Spherical volume integral 4*pi*int_0^R_cut r^2 W(r, R) dr must equal 1.0."""
        R = 5.0  # Mpc/h
        r_max = 5.0 * R if kernel_type == "gaussian" else 2.0 * R
        N_pts = 5000

        code = f"""
        import {{ evaluateSmoothingKernel }} from './src/fields/velocity_dispersion_tensor.js';
        const R = {R};
        const rMax = {r_max};
        const N = {N_pts};
        const dr = rMax / N;

        let integral = 0.0;
        for (let i = 0; i < N; i++) {{
          const rMid = (i + 0.5) * dr;
          const w = evaluateSmoothingKernel(rMid, R, '{kernel_type}');
          integral += 4.0 * Math.PI * rMid * rMid * w * dr;
        }}

        console.log(JSON.stringify({{ integral }}));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["integral"], rel=1e-3) == 1.0

    def test_constant_velocity_field_yields_zero_dispersion(self):
        """Laminar constant flow v(x) = (u0, v0, w0) must have identically zero velocity dispersion."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { VelocityDispersionTensor3D } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 12, ny: 12, nz: 12,
          origin: [-30.0, -30.0, -30.0],
          boxSize: [60.0, 60.0, 60.0]
        });

        const vField = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [350.0, -200.0, 150.0]);
        const disp = VelocityDispersionTensor3D.computeFromVelocityField(vField, 5.0);

        const sig1dField = disp.getScalarSigma1DField();
        const stats = sig1dField.getStatistics();

        console.log(JSON.stringify({ maxSig1D: stats.max }));
        """
        out = run_node_snippet(code)
        assert out["maxSig1D"] < 1e-4

    def test_linear_shear_flow_analytical_dispersion(self):
        """
        Linear shear flow vx(y) = S0 * y, vy = vz = 0.
        With Gaussian kernel smoothing of scale R:
        <vx> = S0 * y, <vx^2> = S0^2 * <y'^2> = S0^2 * (y^2 + R^2)
        => sigma_xx^2 = <vx^2> - <vx>^2 = S0^2 * R^2.
        All other dispersion components = 0.
        """
        S0 = 25.0  # km/s / (Mpc/h)
        R = 4.0   # Mpc/h
        expected_sxx = (S0 * R) ** 2  # 100^2 = 10000 (km/s)^2

        code = f"""
        import {{ GridIndexer }} from './src/fields/grid_indexer.js';
        import {{ VelocityField }} from './src/fields/velocity_field.js';
        import {{ VelocityDispersionTensor3D }} from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({{
          nx: 25, ny: 25, nz: 25,
          origin: [-50.0, -50.0, -50.0],
          boxSize: [100.0, 100.0, 100.0]
        }});

        const S0 = {S0};
        const R = {R};

        const vField = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [S0 * y, 0.0, 0.0]);
        const disp = VelocityDispersionTensor3D.computeFromVelocityField(vField, R, {{ kernelType: 'gaussian' }});

        // Sample at interior grid center (ix=12, iy=12, iz=12) where x=0, y=0, z=0
        const pt = disp.getTensorAt(12, 12, 12);

        console.log(JSON.stringify({{
          sxx: pt.sxx,
          syy: pt.syy,
          szz: pt.szz,
          sxy: pt.sxy,
          expectedSxx: {expected_sxx}
        }}));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sxx"], rel=0.05) == expected_sxx
        assert out["syy"] < 1e-5
        assert out["szz"] < 1e-5
        assert abs(out["sxy"]) < 1e-5

    def test_continuous_trilinear_interpolation(self):
        """Continuous evaluation at non-grid point matches trilinear interpolated tensor."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityDispersionTensor3D } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 8, ny: 8, nz: 8,
          origin: [0.0, 0.0, 0.0],
          boxSize: [70.0, 70.0, 70.0] // dx=10
        });

        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          100.0 + x * 2.0, 10.0, 0.0,
          200.0 + y * 3.0, 0.0,
          300.0 + z * 4.0
        ]);

        // Evaluate at continuous position (15.0, 25.0, 35.0)
        const pt = disp.evaluateAt(15.0, 25.0, 35.0);

        console.log(JSON.stringify({
          sxx: pt.sxx,
          syy: pt.syy,
          szz: pt.szz,
          sxy: pt.sxy
        }));
        """
        out = run_node_snippet(code)
        # Linear function is reconstructed exactly by trilinear interpolation
        assert pytest.approx(out["sxx"], rel=1e-10) == 100.0 + 15.0 * 2.0
        assert pytest.approx(out["syy"], rel=1e-10) == 200.0 + 25.0 * 3.0
        assert pytest.approx(out["szz"], rel=1e-10) == 300.0 + 35.0 * 4.0
        assert pytest.approx(out["sxy"], rel=1e-10) == 10.0


# ============================================================================
# 3. DISCRETE PARTICLE CATALOG DISPERSION ESTIMATOR
# ============================================================================

class TestParticleCatalogDispersionEstimator:
    """Verify particle-to-grid dispersion tensor computation from discrete tracer catalogs."""

    def test_two_stream_particle_distribution(self):
        """
        Two equal-mass counter-streaming populations at the same location:
        Half moving with +V0, half moving with -V0 along x-axis.
        <vx> = 0, <vx^2> = V0^2 => sigma_xx^2 = V0^2.
        """
        V0 = 400.0  # km/s
        code = f"""
        import {{ GridIndexer }} from './src/fields/grid_indexer.js';
        import {{ VelocityDispersionTensor3D }} from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({{
          nx: 11, ny: 11, nz: 11,
          origin: [-20.0, -20.0, -20.0],
          boxSize: [40.0, 40.0, 40.0]
        }});

        const positions = [];
        const velocities = [];
        const V0 = {V0};

        // Create 200 particles distributed symmetrically near origin (0, 0, 0)
        for (let i = 0; i < 100; i++) {{
          const dx = (Math.random() - 0.5) * 2.0;
          const dy = (Math.random() - 0.5) * 2.0;
          const dz = (Math.random() - 0.5) * 2.0;
          positions.push([dx, dy, dz]);
          velocities.push([V0, 0.0, 0.0]); // Stream +V0
        }}
        for (let i = 0; i < 100; i++) {{
          const dx = (Math.random() - 0.5) * 2.0;
          const dy = (Math.random() - 0.5) * 2.0;
          const dz = (Math.random() - 0.5) * 2.0;
          positions.push([dx, dy, dz]);
          velocities.push([-V0, 0.0, 0.0]); // Stream -V0
        }}

        const disp = VelocityDispersionTensor3D.computeFromParticleCatalog(
          positions, velocities, grid, 6.0
        );

        // Sample at grid node (5, 5, 5) which is exactly at (0, 0, 0)
        const pt = disp.getTensorAt(5, 5, 5);

        console.log(JSON.stringify({{
          sxx: pt.sxx,
          syy: pt.syy,
          szz: pt.szz,
          expectedSxx: {V0 * V0}
        }}));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sxx"], rel=0.01) == V0 * V0
        assert out["syy"] < 1e-5
        assert out["szz"] < 1e-5


# ============================================================================
# 4. SPHERICAL HALO COORDINATE TRANSFORMATION & ANISOTROPY PARAMETER BETA
# ============================================================================

class TestSphericalHaloAnisotropyParameter:
    """Verify coordinate transformations from Cartesian to Spherical frames and Binney beta evaluation."""

    def test_isotropic_dispersion_tensor_yields_beta_zero(self):
        """Spherically isotropic dispersion Sigma = sigma_0^2 * I yields sigma_r^2 = sigma_theta^2 = sigma_phi^2 and beta = 0."""
        sigma0 = 350.0
        code = f"""
        import {{ DispersionTensorPoint, SphericalAnisotropyAnalyzer }} from './src/fields/velocity_dispersion_tensor.js';
        const s0sq = {sigma0 * sigma0};
        const pt = new DispersionTensorPoint(s0sq, 0.0, 0.0, s0sq, 0.0, s0sq, [12.0, -8.0, 15.0]);
        const sph = SphericalAnisotropyAnalyzer.transformToSpherical(pt, [0.0, 0.0, 0.0]);

        console.log(JSON.stringify({{
          sigmaR2: sph.sigmaR2,
          sigmaTheta2: sph.sigmaTheta2,
          sigmaPhi2: sph.sigmaPhi2,
          beta: sph.beta,
          symmetrizedBeta: sph.symmetrizedBeta,
          regime: sph.classifyRegime()
        }}));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sigmaR2"], rel=1e-12) == sigma0 ** 2
        assert pytest.approx(out["sigmaTheta2"], rel=1e-12) == sigma0 ** 2
        assert pytest.approx(out["sigmaPhi2"], rel=1e-12) == sigma0 ** 2
        assert pytest.approx(out["beta"], abs=1e-12) == 0.0
        assert pytest.approx(out["symmetrizedBeta"], abs=1e-12) == 0.0
        assert out["regime"] == "isotropic"

    def test_pure_radial_orbit_anisotropy_beta_one(self):
        """Pure radial dispersion (e.g. position along x-axis with non-zero sxx only) yields beta = 1.0."""
        code = """
        import { DispersionTensorPoint, SphericalAnisotropyAnalyzer } from './src/fields/velocity_dispersion_tensor.js';
        // Point on x-axis: pos = [10.0, 0.0, 0.0], halo center = [0.0, 0.0, 0.0]
        // Radial direction is x-axis. sxx = 500^2, syy = szz = 0
        const pt = new DispersionTensorPoint(250000.0, 0.0, 0.0, 0.0, 0.0, 0.0, [10.0, 0.0, 0.0]);
        const sph = SphericalAnisotropyAnalyzer.transformToSpherical(pt, [0.0, 0.0, 0.0]);

        console.log(JSON.stringify({
          sigmaR2: sph.sigmaR2,
          sigmaTheta2: sph.sigmaTheta2,
          sigmaPhi2: sph.sigmaPhi2,
          sigmaTan2: sph.sigmaTan2,
          beta: sph.beta,
          symmetrizedBeta: sph.symmetrizedBeta,
          regime: sph.classifyRegime()
        }));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sigmaR2"], rel=1e-12) == 250000.0
        assert pytest.approx(out["sigmaTan2"], abs=1e-12) == 0.0
        assert pytest.approx(out["beta"], abs=1e-12) == 1.0
        assert pytest.approx(out["symmetrizedBeta"], abs=1e-12) == 1.0
        assert out["regime"] == "pure_radial"

    def test_pure_tangential_orbit_anisotropy(self):
        """Pure tangential motion (sigma_r^2 = 0, sigma_t^2 > 0) yields beta -> -inf and symmetrized beta = -1.0."""
        code = """
        import { DispersionTensorPoint, SphericalAnisotropyAnalyzer } from './src/fields/velocity_dispersion_tensor.js';
        // Point on x-axis: pos = [10.0, 0.0, 0.0] => radial is x.
        // Set syy = szz = 100000 (tangential), sxx = 0 (radial = 0)
        const pt = new DispersionTensorPoint(0.0, 0.0, 0.0, 100000.0, 0.0, 100000.0, [10.0, 0.0, 0.0]);
        const sph = SphericalAnisotropyAnalyzer.transformToSpherical(pt, [0.0, 0.0, 0.0]);

        console.log(JSON.stringify({
          sigmaR2: sph.sigmaR2,
          sigmaTan2: sph.sigmaTan2,
          isNegativeInfinity: !isFinite(sph.beta) && sph.beta < 0,
          symmetrizedBeta: sph.symmetrizedBeta,
          regime: sph.classifyRegime()
        }));
        """
        out = run_node_snippet(code)
        assert out["isNegativeInfinity"] is True
        assert pytest.approx(out["symmetrizedBeta"], abs=1e-12) == -1.0
        assert out["regime"] == "circular_tangential"

    def test_osipkov_merritt_anisotropy_profile(self):
        """Osipkov-Merritt beta(r) = r^2 / (r^2 + ra^2) matches analytical curve."""
        code = """
        import { SphericalAnisotropyAnalyzer } from './src/fields/velocity_dispersion_tensor.js';
        const ra = 5.0; // Mpc/h
        const radii = [0.0, 2.5, 5.0, 10.0, 20.0];
        const betas = radii.map(r => SphericalAnisotropyAnalyzer.osipkovMerrittBeta(r, ra));

        console.log(JSON.stringify({ betas }));
        """
        out = run_node_snippet(code)
        betas = out["betas"]
        # r=0 -> beta=0
        assert pytest.approx(betas[0], abs=1e-12) == 0.0
        # r=2.5 -> beta = 2.5^2 / (2.5^2 + 5^2) = 6.25 / 31.25 = 0.2
        assert pytest.approx(betas[1], rel=1e-12) == 0.2
        # r=5.0 -> beta = 5^2 / (5^2 + 5^2) = 0.5
        assert pytest.approx(betas[2], rel=1e-12) == 0.5
        # r=10.0 -> beta = 100 / 125 = 0.8
        assert pytest.approx(betas[3], rel=1e-12) == 0.8


# ============================================================================
# 5. RADIAL ANISOTROPY PROFILE SHELL BINNING
# ============================================================================

class TestRadialAnisotropyProfileShellBinning:
    """Verify radial shell binning and spherical dispersion profile reconstruction."""

    def test_radial_profile_binning_linear_spacing(self):
        """Radial shell binning over 3D isotropic dispersion field reconstructs flat beta=0 profile."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityDispersionTensor3D, SphericalAnisotropyAnalyzer } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 16, ny: 16, nz: 16,
          origin: [-20.0, -20.0, -20.0],
          boxSize: [40.0, 40.0, 40.0]
        });

        // Constant isotropic 300 km/s dispersion
        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          90000.0, 0.0, 0.0,
          90000.0, 0.0,
          90000.0
        ]);

        const profile = SphericalAnisotropyAnalyzer.computeRadialAnisotropyProfile(
          disp, [0.0, 0.0, 0.0], { rMin: 2.0, rMax: 16.0, numBins: 7, logSpacing: false }
        );

        console.log(JSON.stringify({ profile }));
        """
        out = run_node_snippet(code)
        profile = out["profile"]
        assert len(profile) == 7
        for bin_res in profile:
            if bin_res["count"] > 0:
                assert pytest.approx(bin_res["sigmaR"], rel=1e-5) == 300.0
                assert pytest.approx(bin_res["sigmaTan"], rel=1e-5) == 300.0
                assert pytest.approx(bin_res["beta"], abs=1e-5) == 0.0
                assert bin_res["regime"] == "isotropic"


# ============================================================================
# 6. COSMIC THERMAL & KINETIC PRESSURE FIELDS & SOUND SPEED
# ============================================================================

class TestCosmicThermalAndKineticPressureFields:
    """Verify pressure tensor P_ij = rho * sigma_ij^2, scalar pressure P, sound speed c_s, and Mach number."""

    def test_pressure_tensor_and_effective_sound_speed(self):
        """
        Homogeneous medium with density rho0, isotropic dispersion sigma0.
        P = rho0 * sigma0^2.
        Sound speed c_s = sqrt(gamma * P / rho0) = sqrt(gamma) * sigma0.
        With gamma = 5/3: c_s = sqrt(5/3) * sigma0.
        """
        sigma0 = 300.0  # km/s
        gamma = 5.0 / 3.0
        expected_cs = math.sqrt(gamma) * sigma0  # ~387.298 km/s

        code = f"""
        import {{ GridIndexer }} from './src/fields/grid_indexer.js';
        import {{ VelocityDispersionTensor3D, CosmicPressureField }} from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({{
          nx: 8, ny: 8, nz: 8,
          origin: [-10.0, -10.0, -10.0],
          boxSize: [20.0, 20.0, 20.0]
        }});

        const s0sq = {sigma0 * sigma0};
        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          s0sq, 0.0, 0.0,
          s0sq, 0.0,
          s0sq
        ]);

        const pressure = new CosmicPressureField(disp, null, {{
          gammaAdiabatic: {gamma}
        }});

        const pField = pressure.getIsotropicPressureField();
        const csField = pressure.getSoundSpeedField();
        const tempField = pressure.getTemperatureField();

        console.log(JSON.stringify({{
          sampleP: pField.get(4, 4, 4),
          sampleCs: csField.get(4, 4, 4),
          sampleTemp: tempField.get(4, 4, 4),
          expectedCs: {expected_cs}
        }}));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sampleCs"], rel=1e-10) == expected_cs
        assert out["sampleP"] > 0.0
        assert out["sampleTemp"] > 0.0

    def test_cosmic_mach_number_field(self):
        """Supersonic vs subsonic flow regime classification via Mach number M = |v| / c_s."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { VelocityDispersionTensor3D, CosmicPressureField } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 8, ny: 8, nz: 8,
          origin: [-10.0, -10.0, -10.0],
          boxSize: [20.0, 20.0, 20.0]
        });

        // Dispersion sigma1D = 200 km/s -> c_s = sqrt(5/3)*200 = 258.199 km/s
        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          40000.0, 0.0, 0.0,
          40000.0, 0.0,
          40000.0
        ]);

        // Peculiar velocity: vx = 500 km/s, vy = vz = 0 -> supersonic M ~ 1.9365
        const vField = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [500.0, 0.0, 0.0]);

        const pressure = new CosmicPressureField(disp, null);
        const machField = pressure.computeMachNumberField(vField);

        const sampleMach = machField.get(4, 4, 4);

        console.log(JSON.stringify({
          sampleMach,
          expectedMach: 500.0 / Math.sqrt((5.0 / 3.0) * 40000.0)
        }));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["sampleMach"], rel=1e-10) == out["expectedMach"]
        assert out["sampleMach"] > 1.0  # Supersonic


# ============================================================================
# 7. SPHERICAL JEANS EQUATION COSMOLOGICAL MASS ESTIMATOR
# ============================================================================

class TestSphericalJeansEquationMassEstimator:
    """Verify Jeans equation dynamical mass estimator against exact analytical profiles."""

    def test_singular_isothermal_sphere_exact_solution(self):
        """
        Singular Isothermal Sphere (SIS):
          rho(r) = sigma^2 / (2 pi G r^2) => d ln(rho) / d ln(r) = -2
          sigma_r^2(r) = sigma^2 = const  => d ln(sigma_r^2) / d ln(r) = 0
          beta(r) = 0
          Exact enclosed mass: M(<r) = 2 * sigma^2 * r / G.
        """
        sigma = 220.0  # km/s (typical Milky Way / spiral halo)
        radii = np.logspace(np.log10(0.01), np.log10(0.3), 30)  # Mpc/h

        # Analytical density profile
        density = (sigma ** 2) / (2.0 * np.pi * G_COSMO * (radii ** 2))
        sigma_r2 = np.full_like(radii, sigma ** 2)
        beta = np.zeros_like(radii)

        code = f"""
        import {{ JeansMassEstimator, G_COSMO_MPC_MSUN }} from './src/fields/velocity_dispersion_tensor.js';
        const radii = {json.dumps(radii.tolist())};
        const density = {json.dumps(density.tolist())};
        const sigmaR2 = {json.dumps(sigma_r2.tolist())};
        const beta = {json.dumps(beta.tolist())};

        const estimated = JeansMassEstimator.estimateEnclosedMassProfile(
          radii, density, sigmaR2, beta, {{ G: G_COSMO_MPC_MSUN }}
        );

        const exact = radii.map(r => JeansMassEstimator.singularIsothermalSphereExactMass(r, {sigma}, G_COSMO_MPC_MSUN));

        console.log(JSON.stringify({{
          estimatedMasses: estimated.map(p => p.mass),
          exactMasses: exact
        }}));
        """
        out = run_node_snippet(code)
        est = np.array(out["estimatedMasses"])
        exact = np.array(out["exactMasses"])

        # Interior points should match exact SIS solution within < 0.1%
        np.testing.assert_allclose(est[2:-2], exact[2:-2], rtol=1e-3)

    def test_hernquist_exact_mass_profile(self):
        """Verify Hernquist profile analytical mass M(<r) = M_tot * r^2 / (r + a)^2."""
        M_tot = 1.0e12  # M_sun
        a = 0.02        # Mpc/h (20 kpc/h)
        radii = [0.01, 0.02, 0.04, 0.1]

        code = f"""
        import {{ JeansMassEstimator }} from './src/fields/velocity_dispersion_tensor.js';
        const radii = {json.dumps(radii)};
        const masses = radii.map(r => JeansMassEstimator.hernquistExactMass(r, {M_tot}, {a}));
        console.log(JSON.stringify({{ masses }}));
        """
        out = run_node_snippet(code)
        masses = out["masses"]

        # r = a -> M(<a) = M_tot * 1/4 = 0.25e12
        assert pytest.approx(masses[1], rel=1e-10) == 0.25 * M_tot
        # r = 2a (0.04) -> M(<2a) = M_tot * (2/3)^2 = 4/9 * M_tot
        assert pytest.approx(masses[2], rel=1e-10) == (4.0 / 9.0) * M_tot

    def test_nfw_exact_mass_profile(self):
        """Verify NFW profile analytical mass formula."""
        rho0 = 1.0e14  # M_sun / (Mpc/h)^3
        rs = 0.05      # Mpc/h (50 kpc/h)
        r = 0.1        # x = r/rs = 2.0

        code = f"""
        import {{ JeansMassEstimator }} from './src/fields/velocity_dispersion_tensor.js';
        const mass = JeansMassEstimator.nfwExactMass({r}, {rho0}, {rs});
        console.log(JSON.stringify({{ mass }}));
        """
        out = run_node_snippet(code)
        # Expected: 4 * pi * rho0 * rs^3 * [ln(3) - 2/3]
        expected = 4.0 * math.pi * rho0 * (rs ** 3) * (math.log(3.0) - (2.0 / 3.0))
        assert pytest.approx(out["mass"], rel=1e-10) == expected

    def test_anisotropy_parameter_bias_on_mass_estimation(self):
        """
        Jeans Mass: M(<r) = (r * sigma_r^2 / G) * (gamma_rho + gamma_sigma - 2*beta).
        Increasing beta (radial bias) decreases the inferred mass;
        Decreasing beta (tangential bias) increases the inferred mass.
        """
        code = """
        import { JeansMassEstimator, G_COSMO_MPC_MSUN } from './src/fields/velocity_dispersion_tensor.js';
        const radii = [0.05, 0.10, 0.15, 0.20];
        const density = [1e15, 2e14, 5e13, 1e13];
        const sigmaR2 = [90000.0, 85000.0, 80000.0, 75000.0];

        const betaIsotropic = [0.0, 0.0, 0.0, 0.0];
        const betaRadial = [0.5, 0.5, 0.5, 0.5];
        const betaTangential = [-0.5, -0.5, -0.5, -0.5];

        const mIso = JeansMassEstimator.estimateEnclosedMassProfile(radii, density, sigmaR2, betaIsotropic);
        const mRad = JeansMassEstimator.estimateEnclosedMassProfile(radii, density, sigmaR2, betaRadial);
        const mTan = JeansMassEstimator.estimateEnclosedMassProfile(radii, density, sigmaR2, betaTangential);

        console.log(JSON.stringify({
          mIso: mIso[2].mass,
          mRad: mRad[2].mass,
          mTan: mTan[2].mass
        }));
        """
        out = run_node_snippet(code)
        assert out["mRad"] < out["mIso"] < out["mTan"]


# ============================================================================
# 8. SUPERCLUSTER MULTISCALE VIRIAL RATIO (2K / |W|)
# ============================================================================

class TestSuperclusterMultiscaleVirialRatio:
    """Verify kinetic energy integration, gravitational potential, and virial ratio evaluation."""

    def test_center_of_mass_velocity_subtraction(self):
        """Adding a uniform boost v_boost to the velocity field does not change relative bulk kinetic energy."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { VelocityDispersionTensor3D, SuperclusterVirialAnalyzer } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 8, ny: 8, nz: 8,
          origin: [-10.0, -10.0, -10.0],
          boxSize: [20.0, 20.0, 20.0]
        });

        // Base shearing velocity field
        const vField1 = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [x * 10.0, y * 5.0, -z * 15.0]);
        // Boosted velocity field
        const vField2 = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [x * 10.0 + 500.0, y * 5.0 - 300.0, -z * 15.0 + 200.0]);

        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          10000.0, 0.0, 0.0,
          10000.0, 0.0,
          10000.0
        ]);

        const kRes1 = SuperclusterVirialAnalyzer.computeKineticEnergy(grid, vField1, disp);
        const kRes2 = SuperclusterVirialAnalyzer.computeKineticEnergy(grid, vField2, disp);

        console.log(JSON.stringify({
          kBulk1: kRes1.kBulk,
          kBulk2: kRes2.kBulk,
          kDisp1: kRes1.kDisp,
          kDisp2: kRes2.kDisp,
          vcm1: Array.from(kRes1.centerOfMassVelocity),
          vcm2: Array.from(kRes2.centerOfMassVelocity)
        }));
        """
        out = run_node_snippet(code)
        assert pytest.approx(out["kBulk1"], rel=1e-10) == out["kBulk2"]
        assert pytest.approx(out["kDisp1"], rel=1e-10) == out["kDisp2"]
        assert pytest.approx(500.0, rel=1e-10) == (out["vcm2"][0] - out["vcm1"][0])
        assert pytest.approx(-300.0, rel=1e-10) == (out["vcm2"][1] - out["vcm1"][1])
        assert pytest.approx(200.0, rel=1e-10) == (out["vcm2"][2] - out["vcm1"][2])

    def test_virial_state_classification(self):
        """Verify 2K / |W| calculation and classification into equilibrium, collapsing, and unbound regimes."""
        code = """
        import { SuperclusterVirialAnalyzer, VirialState } from './src/fields/velocity_dispersion_tensor.js';

        // 1. Virial equilibrium: 2K = |W| -> 2K/|W| = 1.0
        const resEquil = SuperclusterVirialAnalyzer.evaluateVirialRatio(500.0, -1000.0);

        // 2. Collapsing (super-virial): 2K < |W| -> 2K/|W| = 0.5 < 0.85
        const resCollapse = SuperclusterVirialAnalyzer.evaluateVirialRatio(250.0, -1000.0);

        // 3. Unbound / expanding (sub-virial): 2K > |W| -> 2K/|W| = 1.5 > 1.15
        const resUnbound = SuperclusterVirialAnalyzer.evaluateVirialRatio(750.0, -1000.0);

        console.log(JSON.stringify({
          equil: resEquil,
          collapse: resCollapse,
          unbound: resUnbound
        }));
        """
        out = run_node_snippet(code)

        assert pytest.approx(out["equil"]["virialRatio"], rel=1e-12) == 1.0
        assert out["equil"]["isVirialized"] is True
        assert out["equil"]["virialState"] == "virialized_equilibrium"
        assert pytest.approx(out["equil"]["virialDeficit"], abs=1e-12) == 0.0

        assert pytest.approx(out["collapse"]["virialRatio"], rel=1e-12) == 0.5
        assert out["collapse"]["isVirialized"] is False
        assert out["collapse"]["virialState"] == "collapsing_super_virial"

        assert pytest.approx(out["unbound"]["virialRatio"], rel=1e-12) == 1.5
        assert out["unbound"]["isVirialized"] is False
        assert out["unbound"]["virialState"] == "unbound_sub_virial"


# ============================================================================
# 9. W3C PROV-O CRYPTOGRAPHIC LINEAGE & AUDIT MANIFEST
# ============================================================================

class TestAuditReportAndProvenance:
    """Verify cryptographic audit report and W3C PROV-O manifest generation."""

    def test_dispersion_audit_report_generation(self):
        """Audit report generates valid JSON and NIST SHA-256 integrity digest."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityDispersionTensor3D, CosmicPressureField, generateDispersionAuditReport } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({
          nx: 6, ny: 6, nz: 6,
          origin: [-10.0, -10.0, -10.0],
          boxSize: [20.0, 20.0, 20.0]
        });

        const disp = VelocityDispersionTensor3D.fromAnalyticFunction(grid, (x, y, z) => [
          40000.0, 1000.0, 0.0,
          25000.0, 0.0,
          10000.0
        ]);

        const pressure = new CosmicPressureField(disp, null);
        const report = generateDispersionAuditReport(disp, pressure, { targetSupercluster: 'Laniakea' });

        console.log(JSON.stringify(report));
        """
        out = run_node_snippet(code)

        assert out["activity"] == "CosmicFlows-4:VelocityDispersionAndCosmicPressureAnalysis"
        assert "provenance" in out
        assert len(out["provenance"]["integrityDigest"]) == 64  # SHA-256 hex length
        assert out["dispersionField"]["sigma1D"]["mean"] > 0.0
        assert out["customMetadata"]["targetSupercluster"] == "Laniakea"


# ============================================================================
# 10. ROBUSTNESS, SINGULARITY HANDLING & BOUNDARY EDGE CASES
# ============================================================================

class TestRobustnessAndEdgeCases:
    """Verify parameter bounds, error handling, and singularity protections."""

    def test_smoothing_scale_must_be_strictly_positive(self):
        """Non-positive smoothing radius must throw RangeError."""
        code = """
        import { evaluateSmoothingKernel } from './src/fields/velocity_dispersion_tensor.js';
        try {
          evaluateSmoothingKernel(1.0, 0.0);
          console.log(JSON.stringify({ error: null }));
        } catch (err) {
          console.log(JSON.stringify({ error: err.name }));
        }
        """
        out = run_node_snippet(code)
        assert out["error"] == "RangeError"

    def test_null_buffer_construction_throws_error(self):
        """Constructing VelocityDispersionTensor3D with missing buffers throws Error."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityDispersionTensor3D } from './src/fields/velocity_dispersion_tensor.js';

        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4 });
        try {
          new VelocityDispersionTensor3D(grid, null, null, null, null, null, null);
          console.log(JSON.stringify({ error: null }));
        } catch (err) {
          console.log(JSON.stringify({ error: err.name }));
        }
        """
        out = run_node_snippet(code)
        assert out["error"] == "Error"
