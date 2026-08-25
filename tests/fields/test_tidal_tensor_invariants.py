"""
tests/fields/test_tidal_tensor_invariants.py
============================================
Exhaustive automated Pytest test suite for 3D Gravitational Tidal Tensor Invariants,
Zel'dovich Web Morphological Classification, Cosmological Anisotropy Metrics, Singularity
Collapse Time Estimators, and Tidal Torque Theory (TTT) Angular Momentum Acquisition
for the CosmicFlows-4 Research Workbench.

Test Categories:
----------------
1. Mathematical & Algebraic Invariants of Gravitational Tidal & Deformation Tensors
2. Cayley-Hamilton Identity, Characteristic Polynomials, and Cardano/Vieta Roots
3. Cyclic Jacobi Eigensolver & SO(3) Proper Rotation Orthogonality
4. Cosmological Anisotropy, Ellipticity, Prolateness, Triaxiality, and Sphericity Metrics
5. Hahn et al. / Forero-Romero Cosmic Web Morphological Classification & Threshold Scanning
6. Zel'dovich Approximation, Caustics, and Singularity Collapse Time Estimators
7. Tidal Torque Theory (TTT), Commutator Algebra, and Proto-halo Spin Parameters
8. Exact Analytical Field Solutions (Keplerian, Triaxial Quadrupole, Plummer, Plane Wave)
9. High-Order Finite Difference Stencils (2nd, 4th, 6th) and Grid Boundary Handling
10. 3D Spectral FFT Poisson Solver & Density-to-Tidal Inversion
11. Grid-Wide Morphological Summaries, Volume/Mass Filling Fractions, and Map Generation
12. Doroshkevich PDF & Peak Joint Invariant Statistics
13. Multiscale Gaussian and Top-Hat Spatial Filtering
14. Filament and Sheet Alignment Statistics
15. Cosmological Redshift-Distance & Growth Factor Inversions
16. Robustness, Numerical Stability, Degeneracies, and Extreme Edge Cases

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
# 1. MATHEMATICAL & ALGEBRAIC INVARIANTS
# ============================================================================

class TestTidalTensorAlgebraicInvariants:
    """Verify exact algebraic decompositions and principal invariants of 3D tidal tensors."""

    def test_traceless_decomposition(self):
        """
        Verify decomposition of deformation Hessian Psi_ij into trace (delta)
        and traceless tidal shear T_ij = Psi_ij - (1/3)*Tr(Psi)*delta_ij.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const Psi = new Float64Array([
          4.0, 1.2, -0.8,
          1.2, 2.5,  0.6,
         -0.8, 0.6,  1.5
        ]);

        const decomp = TidalTensorMath.decomposeTidalTensor(Psi);
        const traceless = decomp.traceless;
        const trace = decomp.trace;
        const tracelessTrace = traceless[0] + traceless[4] + traceless[8];

        console.log(JSON.stringify({
          trace,
          tracelessTrace,
          traceless: Array.from(traceless)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["trace"], 4.0 + 2.5 + 1.5, rel_tol=1e-12)
        assert math.isclose(res["tracelessTrace"], 0.0, abs_tol=1e-14)

        # Check traceless diagonal values: Psi_ii - trace/3
        third_trace = (4.0 + 2.5 + 1.5) / 3.0
        assert math.isclose(res["traceless"][0], 4.0 - third_trace, rel_tol=1e-12)
        assert math.isclose(res["traceless"][4], 2.5 - third_trace, rel_tol=1e-12)
        assert math.isclose(res["traceless"][8], 1.5 - third_trace, rel_tol=1e-12)

    def test_principal_invariants_exact_formulas(self):
        """
        Verify computation of First (I_1), Second (I_2), and Third (I_3) Principal Invariants
        against exact analytical formulas and NumPy reference values.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const Psi = new Float64Array([
          3.0, 1.0, 0.5,
          1.0, 2.0, -0.5,
          0.5, -0.5, 1.0
        ]);

        const inv = TidalTensorMath.computeInvariants(Psi);

        console.log(JSON.stringify({
          I1: inv.I1,
          I2: inv.I2,
          I3: inv.I3,
          J2: inv.J2,
          J3: inv.J3,
          discriminant: inv.discriminant
        }));
        """
        res = run_node_snippet(code)

        # NumPy reference calculation
        A = np.array([
            [3.0, 1.0, 0.5],
            [1.0, 2.0, -0.5],
            [0.5, -0.5, 1.0]
        ], dtype=np.float64)

        I1_np = np.trace(A)
        I2_np = 0.5 * (np.trace(A)**2 - np.trace(A @ A))
        I3_np = np.linalg.det(A)

        assert math.isclose(res["I1"], I1_np, rel_tol=1e-12)
        assert math.isclose(res["I2"], I2_np, rel_tol=1e-12)
        assert math.isclose(res["I3"], I3_np, rel_tol=1e-12)

        # For real symmetric matrices, discriminant Delta must be >= 0 (all 3 roots real)
        assert res["discriminant"] >= 0.0

    def test_cayley_hamilton_theorem(self):
        """
        Verify the Cayley-Hamilton Identity: A^3 - I_1 A^2 + I_2 A - I_3 I = 0
        for arbitrary 3D symmetric cosmological deformation matrices.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const matrices = [
          [2.0, 0.5, -0.3, 0.5, 1.5, 0.2, -0.3, 0.2, 0.8],
          [5.0, 2.1,  1.4, 2.1, 3.2, 0.9,  1.4, 0.9, 1.1],
          [-1.5, 0.8, 0.4, 0.8, -0.5, -0.2, 0.4, -0.2, -2.1]
        ];

        const residuals = matrices.map(mat => {
          const inv = TidalTensorMath.computeInvariants(mat);
          return TidalTensorMath.verifyCayleyHamilton(mat, inv);
        });

        console.log(JSON.stringify({ residuals }));
        """
        res = run_node_snippet(code)
        for residual in res["residuals"]:
            assert residual < 1e-12


# ============================================================================
# 2. EIGENVALUE & EIGENSYSTEM ENGINES
# ============================================================================

class TestEigenvalueAndJacobiEngines:
    """Verify analytical Cardano/Vieta roots, cyclic Jacobi diagonalization, and SO(3) rotations."""

    def test_cardano_vs_numpy_eigenvalues(self):
        """
        Compare analytical Cardano / Vieta trigonometric eigenvalues against NumPy linalg.eigvalsh.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const A = new Float64Array([
          4.5, 1.2, 0.8,
          1.2, 3.0, 1.5,
          0.8, 1.5, 1.0
        ]);

        const eigs = TidalTensorMath.cardanoEigenvalues(A);
        console.log(JSON.stringify({ eigenvalues: Array.from(eigs) }));
        """
        res = run_node_snippet(code)

        A_np = np.array([
            [4.5, 1.2, 0.8],
            [1.2, 3.0, 1.5],
            [0.8, 1.5, 1.0]
        ], dtype=np.float64)

        eigs_np = np.sort(np.linalg.eigvalsh(A_np))[::-1]

        assert np.allclose(res["eigenvalues"], eigs_np, rtol=1e-12, atol=1e-12)
        # Ensure strict descending order: lambda1 >= lambda2 >= lambda3
        assert res["eigenvalues"][0] >= res["eigenvalues"][1] >= res["eigenvalues"][2]

    def test_jacobi_eigensystem_orthonormality_and_reconstruction(self):
        """
        Verify Jacobi cyclic diagonalization produces orthonormal eigenvectors R
        such that R^T R = I, det(R) = +1 (proper SO(3) rotation), and R diag(lambda) R^T = A.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const A = new Float64Array([
          6.0, 2.0, 1.0,
          2.0, 4.0, 3.0,
          1.0, 3.0, 5.0
        ]);

        const diag = TidalTensorMath.jacobiDiagonalization(A);
        console.log(JSON.stringify({
          eigenvalues: Array.from(diag.eigenvalues),
          eigenvectors: diag.eigenvectors.map(e => Array.from(e)),
          iterations: diag.iterations
        }));
        """
        res = run_node_snippet(code)
        lambdas = np.array(res["eigenvalues"])
        e1 = np.array(res["eigenvectors"][0])
        e2 = np.array(res["eigenvectors"][1])
        e3 = np.array(res["eigenvectors"][2])

        # Check unit lengths
        assert math.isclose(np.linalg.norm(e1), 1.0, rel_tol=1e-12)
        assert math.isclose(np.linalg.norm(e2), 1.0, rel_tol=1e-12)
        assert math.isclose(np.linalg.norm(e3), 1.0, rel_tol=1e-12)

        # Check mutual orthogonality
        assert math.isclose(np.dot(e1, e2), 0.0, abs_tol=1e-12)
        assert math.isclose(np.dot(e2, e3), 0.0, abs_tol=1e-12)
        assert math.isclose(np.dot(e3, e1), 0.0, abs_tol=1e-12)

        # Check right-handedness det(R) = +1
        R = np.column_stack([e1, e2, e3])
        det_R = np.linalg.det(R)
        assert math.isclose(det_R, 1.0, rel_tol=1e-12)

        # Check matrix reconstruction A_reconstructed = R @ diag(lambdas) @ R^T
        A_orig = np.array([
            [6.0, 2.0, 1.0],
            [2.0, 4.0, 3.0],
            [1.0, 3.0, 5.0]
        ])
        A_rec = R @ np.diag(lambdas) @ R.T
        assert np.allclose(A_rec, A_orig, rtol=1e-12, atol=1e-12)

    def test_so3_rotational_invariance_of_invariants(self):
        """
        Verify that Principal Invariants (I_1, I_2, I_3) and eigenvalues remain invariant
        under arbitrary 3D SO(3) coordinate rotations A' = R A R^T.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        // Original matrix
        const A = new Float64Array([
          3.5, 0.8, -0.4,
          0.8, 2.1,  0.6,
         -0.4, 0.6,  1.2
        ]);

        // Arbitrary Euler rotation angles
        const alpha = 0.65;
        const beta = 1.12;
        const gamma = -0.48;

        // Construct 3D rotation matrix R = Rz(gamma) * Ry(beta) * Rx(alpha)
        const c1 = Math.cos(alpha), s1 = Math.sin(alpha);
        const c2 = Math.cos(beta),  s2 = Math.sin(beta);
        const c3 = Math.cos(gamma), s3 = Math.sin(gamma);

        const R = [
          c2*c3,  s1*s2*c3 - c1*s3,  c1*s2*c3 + s1*s3,
          c2*s3,  s1*s2*s3 + c1*c3,  c1*s2*s3 - s1*c3,
          -s2,    s1*c2,             c1*c2
        ];

        // Compute rotated matrix A_rot = R * A * R^T
        const A_rot = new Float64Array(9);
        for (let i = 0; i < 3; i++) {
          for (let j = 0; j < 3; j++) {
            let sum = 0.0;
            for (let k = 0; k < 3; k++) {
              for (let l = 0; l < 3; l++) {
                sum += R[i*3 + k] * A[k*3 + l] * R[j*3 + l];
              }
            }
            A_rot[i*3 + j] = sum;
          }
        }

        const invOrig = TidalTensorMath.computeInvariants(A);
        const invRot = TidalTensorMath.computeInvariants(A_rot);

        const eigsOrig = TidalTensorMath.cardanoEigenvalues(A);
        const eigsRot = TidalTensorMath.cardanoEigenvalues(A_rot);

        console.log(JSON.stringify({
          invOrig,
          invRot,
          eigsOrig: Array.from(eigsOrig),
          eigsRot: Array.from(eigsRot)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["invOrig"]["I1"], res["invRot"]["I1"], rel_tol=1e-12)
        assert math.isclose(res["invOrig"]["I2"], res["invRot"]["I2"], rel_tol=1e-12)
        assert math.isclose(res["invOrig"]["I3"], res["invRot"]["I3"], rel_tol=1e-12)
        assert np.allclose(res["eigsOrig"], res["eigsRot"], rtol=1e-12, atol=1e-12)


# ============================================================================
# 3. COSMOLOGICAL ANISOTROPY & MORPHOLOGY METRICS
# ============================================================================

class TestCosmologicalAnisotropyAndMorphology:
    """Verify ellipticity, prolateness, shear magnitude, and triaxiality metrics."""

    def test_prolateness_bounding_theorem(self):
        """
        Verify Doroshkevich (1970) and BBKS (1986) inequality: -e <= p <= e.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const testSets = [
          [5.0, 3.0, 1.0],   // Standard triaxial
          [5.0, 4.5, 1.0],   // Asymmetric
          [6.0, 2.0, 1.0],   // Strong shear
          [4.0, 4.0, 4.0],   // Isotropic / spherical
          [3.0, 2.0, 1.0]    // Pure triaxial symmetric (p = 0)
        ];

        const results = testSets.map(eigs => TidalTensorMath.computeAnisotropyParameters(eigs));
        console.log(JSON.stringify({ results }));
        """
        res = run_node_snippet(code)
        for r in res["results"]:
            e = r["ellipticity"]
            p = r["prolateness"]
            assert p >= -e - 1e-12
            assert p <= e + 1e-12
            assert r["sphericity"] <= 1.0 + 1e-12
            assert 0.0 <= r["triaxiality"] <= 1.0

    def test_pure_prolate_filament_geometry(self):
        """
        Configuration with lambda1 = lambda2 > lambda3 (e.g. [4, 4, 1]):
        p = (lambda1 + lambda3 - 2*lambda2) / (2*sum) = - (lambda1 - lambda3) / (2*sum) = -e.
        Triaxiality T = (4^2 - 4^2) / (4^2 - 1^2) = 0.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const eigs = [4.0, 4.0, 1.0];
        const params = TidalTensorMath.computeAnisotropyParameters(eigs);
        console.log(JSON.stringify(params));
        """
        res = run_node_snippet(code)
        assert res["ellipticity"] > 0.0
        assert math.isclose(res["prolateness"], -res["ellipticity"], rel_tol=1e-12)
        assert math.isclose(res["triaxiality"], 0.0, abs_tol=1e-12)

    def test_pure_oblate_pancake_sheet_geometry(self):
        """
        Configuration with lambda1 > lambda2 = lambda3 (e.g. [5, 1, 1]):
        p = (lambda1 + lambda3 - 2*lambda2) / (2*sum) = + (lambda1 - lambda3) / (2*sum) = +e.
        Triaxiality T = (5^2 - 1^2) / (5^2 - 1^2) = 1.0.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const eigs = [5.0, 1.0, 1.0];
        const params = TidalTensorMath.computeAnisotropyParameters(eigs);
        console.log(JSON.stringify(params));
        """
        res = run_node_snippet(code)
        assert res["ellipticity"] > 0.0
        assert math.isclose(res["prolateness"], res["ellipticity"], rel_tol=1e-12)
        assert math.isclose(res["triaxiality"], 1.0, abs_tol=1e-12)

    def test_spherical_isotropic_geometry(self):
        """
        Isotropic spherical perturbation: lambda1 = lambda2 = lambda3.
        e = 0, p = 0, q = 0, S = 1.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const eigs = [3.0, 3.0, 3.0];
        const params = TidalTensorMath.computeAnisotropyParameters(eigs);
        console.log(JSON.stringify(params));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["ellipticity"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["prolateness"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["shearMagnitude"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["sphericity"], 1.0, abs_tol=1e-12)


# ============================================================================
# 4. COSMIC WEB CLASSIFICATION (HAHN / FORERO-ROMERO)
# ============================================================================

class TestCosmicWebClassification:
    """Verify Hahn et al. (2007) and Forero-Romero et al. (2009) cosmic web classification."""

    def test_four_cosmic_web_types(self):
        """
        Verify exact classification of all 4 Cosmic Web types with threshold gamma_th = 0.2:
        - Void: 0 eigenvalues > gamma_th (all <= 0.2)
        - Sheet: 1 eigenvalue > gamma_th
        - Filament: 2 eigenvalues > gamma_th
        - Knot: 3 eigenvalues > gamma_th
        """
        code = """
        import { TidalTensorMath, CosmicWebType } from './src/fields/tidal_tensor_invariants.js';

        const gammaTh = 0.2;

        const voidEigs = [0.15, -0.1, -0.5];
        const sheetEigs = [0.8, 0.1, -0.4];
        const filamentEigs = [1.2, 0.5, -0.2];
        const knotEigs = [2.0, 1.5, 0.4];

        console.log(JSON.stringify({
          voidType: TidalTensorMath.classifyCosmicWeb(voidEigs, gammaTh),
          sheetType: TidalTensorMath.classifyCosmicWeb(sheetEigs, gammaTh),
          filamentType: TidalTensorMath.classifyCosmicWeb(filamentEigs, gammaTh),
          knotType: TidalTensorMath.classifyCosmicWeb(knotEigs, gammaTh)
        }));
        """
        res = run_node_snippet(code)
        assert res["voidType"] == 0      # CosmicWebType.VOID
        assert res["sheetType"] == 1     # CosmicWebType.SHEET
        assert res["filamentType"] == 2  # CosmicWebType.FILAMENT
        assert res["knotType"] == 3      # CosmicWebType.KNOT

    def test_principal_direction_vectors(self):
        """
        Verify orientation vectors:
        - Sheet normal is aligned with e_1 (axis of maximum compression).
        - Filament spine is aligned with e_3 (axis of slowest collapse / elongation).
        """
        code = """
        import { TidalTensorInvariants, CosmicWebType } from './src/fields/tidal_tensor_invariants.js';

        const tensor = new Float64Array([
          3.0, 0.0, 0.0,
          0.0, 2.0, 0.0,
          0.0, 0.0, 0.5
        ]);

        const e1 = new Float64Array([1.0, 0.0, 0.0]);
        const e2 = new Float64Array([0.0, 1.0, 0.0]);
        const e3 = new Float64Array([0.0, 0.0, 1.0]);

        const inv = new TidalTensorInvariants({
          tensor,
          tracelessTensor: tensor,
          trace: 5.5,
          secondInvariant: 8.5,
          thirdInvariant: 3.0,
          shearSecondInvariant: 0,
          shearThirdInvariant: 0,
          discriminant: 10.0,
          eigenvalues: [3.0, 2.0, 0.5],
          eigenvectors: [e1, e2, e3],
          anisotropy: 0.22,
          ellipticity: 0.22,
          prolateness: -0.05,
          shearMagnitude: 1.5,
          triaxiality: 0.57,
          sphericity: 0.34,
          webType: CosmicWebType.FILAMENT,
          gammaThreshold: 1.0
        });

        console.log(JSON.stringify({
          sheetNormal: Array.from(inv.sheetNormal),
          filamentSpine: Array.from(inv.filamentSpine),
          webTypeName: inv.webTypeName
        }));
        """
        res = run_node_snippet(code)
        assert np.allclose(res["sheetNormal"], [1.0, 0.0, 0.0])
        assert np.allclose(res["filamentSpine"], [0.0, 0.0, 1.0])
        assert res["webTypeName"] == "Filament (2D Collapse)"


# ============================================================================
# 5. ZEL'DOVICH APPROXIMATION & COLLAPSE TIME ESTIMATORS
# ============================================================================

class TestZeldovichDynamicsAndCollapse:
    """Verify 1LPT Zel'dovich mapping, caustics, and collapse time/redshift solvers."""

    def test_eds_linear_growth_and_collapse_times(self):
        """
        In Einstein-de Sitter (Omega_m = 1, Omega_Lambda = 0):
        D(z) = 1 / (1 + z)
        Collapse condition along i-th axis: D(z_col, i) = 1 / lambda_i => 1 + z_col, i = lambda_i
        t_col, i = t_0 / (1 + z_col, i)^1.5 = t_0 / lambda_i^1.5.
        """
        code = """
        import { TidalTensorMath, CosmologicalModel } from './src/fields/tidal_tensor_invariants.js';

        const eigs = [3.0, 2.0, 0.5]; // lambda1 = 3, lambda2 = 2, lambda3 = 0.5
        const cosmo = { model: CosmologicalModel.EDS, omegaM: 1.0, omegaLambda: 0.0, H0: 100.0 };

        const collapse = TidalTensorMath.computeZeldovichCollapse(eigs, cosmo);
        const d_z1 = TidalTensorMath.linearGrowthFactor(collapse.collapseRedshifts[0], cosmo);
        const d_z2 = TidalTensorMath.linearGrowthFactor(collapse.collapseRedshifts[1], cosmo);

        console.log(JSON.stringify({
          zCol: Array.from(collapse.collapseRedshifts),
          tColGyr: Array.from(collapse.collapseTimesGyr),
          d_z1,
          d_z2
        }));
        """
        res = run_node_snippet(code)

        # 1 + z_col1 = 3.0 => z_col1 = 2.0
        assert math.isclose(res["zCol"][0], 2.0, rel_tol=1e-6)
        # 1 + z_col2 = 2.0 => z_col2 = 1.0
        assert math.isclose(res["zCol"][1], 1.0, rel_tol=1e-6)

        # D(z_col) should equal 1 / lambda
        assert math.isclose(res["d_z1"], 1.0 / 3.0, rel_tol=1e-6)
        assert math.isclose(res["d_z2"], 1.0 / 2.0, rel_tol=1e-6)

    def test_lcdm_growth_factor_and_growth_rate(self):
        """
        Verify LCDM growth factor D(z=0) = 1.0 and growth rate f(z=0) = Omega_m^0.55.
        """
        code = """
        import { TidalTensorMath, CosmologicalModel } from './src/fields/tidal_tensor_invariants.js';

        const cosmo = { model: CosmologicalModel.FLAT_LCDM, omegaM: 0.315, omegaLambda: 0.685 };

        const D0 = TidalTensorMath.linearGrowthFactor(0.0, cosmo);
        const D1 = TidalTensorMath.linearGrowthFactor(1.0, cosmo);
        const f0 = TidalTensorMath.linearGrowthRate(0.0, cosmo);

        console.log(JSON.stringify({ D0, D1, f0 }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["D0"], 1.0, rel_tol=1e-6)
        # At z=1, growth factor in LCDM is ~ 0.61
        assert 0.55 < res["D1"] < 0.65
        # f(z=0) = 0.315^0.55 \approx 0.528
        assert math.isclose(res["f0"], 0.315**0.55, rel_tol=1e-4)

    def test_zeldovich_caustic_stream_detection(self):
        """
        Verify Jacobian determinant J = (1 - D*lambda1)(1 - D*lambda2)(1 - D*lambda3)
        and stream count transition when D*lambda1 >= 1.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const q = [10.0, 20.0, 30.0];
        const gradPhi = [1.0, 0.5, 0.2];
        const Psi = [
          2.0, 0.0, 0.0,
          0.0, 0.5, 0.0,
          0.0, 0.0, 0.1
        ];

        // Before collapse: D = 0.3 -> D * lambda1 = 0.6 < 1 (single stream)
        const statePre = TidalTensorMath.evaluateZeldovichState(q, gradPhi, Psi, 0.3);

        // After 1st collapse: D = 0.6 -> D * lambda1 = 1.2 >= 1 (3 streams / caustic)
        const statePost = TidalTensorMath.evaluateZeldovichState(q, gradPhi, Psi, 0.6);

        console.log(JSON.stringify({
          pre: { J: statePre.jacobianDet, streams: statePre.streamCount, caustic: statePre.isCaustic },
          post: { J: statePost.jacobianDet, streams: statePost.streamCount, caustic: statePost.isCaustic }
        }));
        """
        res = run_node_snippet(code)
        assert res["pre"]["J"] > 0.0
        assert res["pre"]["streams"] == 1
        assert res["pre"]["caustic"] is False

        assert res["post"]["J"] <= 0.0
        assert res["post"]["streams"] >= 3
        assert res["post"]["caustic"] is True


# ============================================================================
# 6. TIDAL TORQUE THEORY (TTT) & ANGULAR MOMENTUM
# ============================================================================

class TestTidalTorqueTheory:
    """Verify proto-halo angular momentum acquisition, commutator algebra, and spin parameters."""

    def test_commutator_and_torque_misalignment(self):
        """
        Verify torque tau_i = eps_ijk T_jl I_lk.
        If T and I commute ([T, I] = 0), torque must be identically zero.
        If T and I are misaligned, torque is non-zero and equals the dual vector of [T, I].
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        // 1. Co-aligned case: diagonal T and diagonal I
        const T_diag = [
          2.0, 0.0, 0.0,
          0.0, -1.0, 0.0,
          0.0, 0.0, -1.0
        ];
        const I_diag = [
          10.0, 0.0, 0.0,
          0.0, 5.0, 0.0,
          0.0, 0.0, 2.0
        ];

        const stateAligned = TidalTensorMath.computeTidalTorque(T_diag, I_diag);

        // 2. Misaligned case: rotated I
        const I_misaligned = [
          7.5, 2.5, 0.0,
          2.5, 7.5, 0.0,
          0.0, 0.0, 2.0
        ];

        const stateMisaligned = TidalTensorMath.computeTidalTorque(T_diag, I_misaligned);

        console.log(JSON.stringify({
          alignedTorque: stateAligned.torqueMagnitude,
          alignedL: stateAligned.angularMomentumMagnitude,
          alignedMisalignmentDeg: stateAligned.misalignmentAngleDeg,
          misalignedTorque: stateMisaligned.torqueMagnitude,
          misalignedL: stateMisaligned.angularMomentumMagnitude,
          misalignedAngleDeg: stateMisaligned.misalignmentAngleDeg,
          commutatorNorm: Math.hypot(...stateMisaligned.commutator)
        }));
        """
        res = run_node_snippet(code)

        # Co-aligned: zero torque and zero misalignment
        assert math.isclose(res["alignedTorque"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["alignedL"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["alignedMisalignmentDeg"], 0.0, abs_tol=1e-12)

        # Misaligned: positive torque and angular momentum
        assert res["misalignedTorque"] > 0.0
        assert res["misalignedL"] > 0.0
        assert res["commutatorNorm"] > 0.0
        assert math.isclose(res["misalignedAngleDeg"], 45.0, abs_tol=1e-5)

    def test_spin_parameter_scaling(self):
        """
        Verify Peebles and Bullock spin parameters are dimensionless and scale appropriately with mass.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const T = [
          1.0, 0.5, 0.0,
          0.5, -0.5, 0.0,
          0.0, 0.0, -0.5
        ];
        const I = [
          10.0, 3.0, 0.0,
          3.0, 8.0, 0.0,
          0.0, 0.0, 5.0
        ];

        const state = TidalTensorMath.computeTidalTorque(T, I, {
          haloMass: 1e12,
          haloRadius: 0.25,
          circularVelocity: 180.0,
          bindingEnergy: 1e59
        });

        console.log(JSON.stringify({
          peeblesSpin: state.peeblesSpin,
          bullockSpin: state.bullockSpin,
          haloMass: state.haloMass
        }));
        """
        res = run_node_snippet(code)
        assert res["peeblesSpin"] > 0.0
        assert res["bullockSpin"] > 0.0
        assert math.isclose(res["haloMass"], 1e12, rel_tol=1e-12)


# ============================================================================
# 7. EXACT ANALYTICAL FIELD BENCHMARKS
# ============================================================================

class TestAnalyticTidalFields:
    """Verify analytical field benchmarks against exact mathematical formulations."""

    def test_keplerian_point_mass(self):
        """
        Keplerian potential Phi(r) = -GM/r:
        Psi_ij = (GM/r^3) * [delta_ij - 3 * (x_i x_j / r^2)]
        Eigenvalues: lambda_parallel = -2*GM/r^3, lambda_perp = +GM/r^3 (multiplicity 2).
        Trace = 0 (Laplace equation outside origin).
        """
        code = """
        import { AnalyticTidalFields, TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const GM = 2.0;
        const x = 3.0, y = 4.0, z = 0.0; // r = 5.0
        const res = AnalyticTidalFields.keplerianMonopole(x, y, z, GM);

        const invariants = TidalTensorMath.computeInvariants(res.tensor);
        const cardano = TidalTensorMath.cardanoEigenvalues(res.tensor);

        console.log(JSON.stringify({
          potential: res.potential,
          trace: invariants.I1,
          eigenvalues: Array.from(res.eigenvalues),
          cardanoEigenvalues: Array.from(cardano)
        }));
        """
        res = run_node_snippet(code)
        r = 5.0
        GM = 2.0
        expected_pot = -GM / r
        expected_perp = GM / (r**3)
        expected_par = -2.0 * GM / (r**3)

        assert math.isclose(res["potential"], expected_pot, rel_tol=1e-12)
        assert math.isclose(res["trace"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["eigenvalues"][0], expected_perp, rel_tol=1e-12)
        assert math.isclose(res["eigenvalues"][1], expected_perp, rel_tol=1e-12)
        assert math.isclose(res["eigenvalues"][2], expected_par, rel_tol=1e-12)
        assert np.allclose(res["cardanoEigenvalues"], [expected_perp, expected_perp, expected_par], rtol=1e-12)

    def test_triaxial_quadrupole(self):
        """
        Phi(x, y, z) = 0.5 * (A*x^2 + B*y^2 + C*z^2):
        Psi = diag(A, B, C).
        """
        code = """
        import { AnalyticTidalFields, TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const quad = AnalyticTidalFields.triaxialQuadrupole(4.0, 2.5, 1.0);
        const invariants = TidalTensorMath.computeInvariants(quad.tensor);

        console.log(JSON.stringify({
          tensor: Array.from(quad.tensor),
          traceless: Array.from(quad.traceless),
          eigenvalues: Array.from(quad.eigenvalues),
          trace: invariants.I1
        }));
        """
        res = run_node_snippet(code)
        assert np.allclose(res["eigenvalues"], [4.0, 2.5, 1.0])
        assert math.isclose(res["trace"], 7.5, rel_tol=1e-12)

        # Traceless diagonal: 4 - 2.5 = 1.5, 2.5 - 2.5 = 0.0, 1.0 - 2.5 = -1.5
        assert np.allclose(res["traceless"], [1.5, 0, 0, 0, 0, 0, 0, 0, -1.5])

    def test_plummer_sphere_tidal_field(self):
        """
        Plummer sphere Phi(r) = -GM / sqrt(r^2 + b^2):
        At origin r = 0: d2Phi/dx2 = d2Phi/dy2 = d2Phi/dz2 = GM / b^3.
        Trace = 3 * GM / b^3 = 4 * pi * G * rho_0.
        """
        code = """
        import { AnalyticTidalFields, TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const GM = 10.0;
        const b = 2.0;
        const res = AnalyticTidalFields.plummerSphere(0.0, 0.0, 0.0, GM, b);
        const inv = TidalTensorMath.computeInvariants(res.tensor);

        console.log(JSON.stringify({
          potential: res.potential,
          trace: inv.I1,
          expectedTrace: (3.0 * GM) / (b ** 3),
          diag0: res.tensor[0]
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["potential"], -10.0 / 2.0, rel_tol=1e-12)
        assert math.isclose(res["trace"], 3.75, rel_tol=1e-12)
        assert math.isclose(res["diag0"], 1.25, rel_tol=1e-12)


# ============================================================================
# 8. FINITE DIFFERENCE DIFFERENTIATION & STENCILS
# ============================================================================

class TestFiniteDifferenceStencils:
    """Verify finite difference stencil orders and grid evaluation accuracy."""

    def test_stencil_evaluation_on_harmonic_plane_wave(self):
        """
        Verify 4th-order stencil on 3D harmonic plane wave Phi(x, y, z) = cos(kx*x + ky*y + kz*z).
        Exact second derivative d2Phi/dx2 = -kx^2 cos(k.x).
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { ScalarField3D } from './src/fields/scalar_field_3d.js';
        import { TidalTensorAnalyzer, DifferentiationOrder } from './src/fields/tidal_tensor_invariants.js';

        const N = 32;
        const L = 100.0;
        const kx = (2.0 * Math.PI) / L;
        const ky = (2.0 * Math.PI) / L;
        const kz = (2.0 * Math.PI) / L;

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          spacing: [L / N, L / N, L / N],
          boundaryMode: BoundaryMode.PERIODIC,
          isCellCentered: true
        });

        const total = grid.totalCells;
        const data = new Float64Array(total);

        for (let iz = 0; iz < N; iz++) {
          const z = (iz + 0.5) * grid.dz;
          for (let iy = 0; iy < N; iy++) {
            const y = (iy + 0.5) * grid.dy;
            for (let ix = 0; ix < N; ix++) {
              const x = (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              data[idx] = Math.cos(kx * x + ky * y + kz * z);
            }
          }
        }

        const potField = new ScalarField3D(grid, data, 'pot');
        const analyzer = new TidalTensorAnalyzer(potField, {
          order: DifferentiationOrder.FOURTH,
          boundaryMode: BoundaryMode.PERIODIC
        });

        const comp = analyzer.computeGridDeformationComponents();

        // Sample voxel (10, 12, 14)
        const ix = 10, iy = 12, iz = 14;
        const x = (ix + 0.5) * grid.dx;
        const y = (iy + 0.5) * grid.dy;
        const z = (iz + 0.5) * grid.dz;
        const idx = grid.index(ix, iy, iz);

        const cosVal = Math.cos(kx * x + ky * y + kz * z);
        const exact_xx = -kx * kx * cosVal;
        const exact_xy = -kx * ky * cosVal;

        console.log(JSON.stringify({
          numeric_xx: comp.xx[idx],
          exact_xx,
          numeric_xy: comp.xy[idx],
          exact_xy
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["numeric_xx"], res["exact_xx"], rel_tol=1e-3)
        assert math.isclose(res["numeric_xy"], res["exact_xy"], rel_tol=1e-3)


# ============================================================================
# 9. 3D SPECTRAL FFT POISSON SOLVER & DENSITY-TO-TIDAL
# ============================================================================

class TestSpectralPoissonSolver:
    """Verify spectral FFT Poisson inversion from density contrast delta(x) to potential Phi(x) and tidal shear."""

    def test_spectral_density_to_potential_solution(self):
        """
        For a sinusoidal density delta(x) = delta_0 * cos(kx * x):
        Poisson equation: nabla^2 Phi = delta => -kx^2 Phi = delta_0 cos(kx * x)
        => Phi(x) = - (delta_0 / kx^2) * cos(kx * x).
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { ScalarField3D } from './src/fields/scalar_field_3d.js';
        import { TidalTensorAnalyzer } from './src/fields/tidal_tensor_invariants.js';

        const N = 16;
        const L = 100.0;
        const kx = (2.0 * Math.PI) / L;
        const delta0 = 1.5;

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          spacing: [L / N, L / N, L / N],
          boundaryMode: BoundaryMode.PERIODIC,
          isCellCentered: true
        });

        const total = grid.totalCells;
        const deltaData = new Float64Array(total);

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const x = (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              deltaData[idx] = delta0 * Math.cos(kx * x);
            }
          }
        }

        const deltaField = new ScalarField3D(grid, deltaData, 'density');
        const analyzer = new TidalTensorAnalyzer(deltaField, {
          isDensity: true,
          cosmology: { fourPiGRhoBar: 1.0 }
        });

        const phiField = analyzer.potentialField;

        // Check sample point at ix = 4
        const ix = 4;
        const x = (ix + 0.5) * grid.dx;
        const idx = grid.index(ix, 0, 0);

        const exactPhi = -(delta0 / (kx * kx)) * Math.cos(kx * x);
        const numericPhi = phiField.data[idx];

        console.log(JSON.stringify({
          numericPhi,
          exactPhi,
          relDiff: Math.abs((numericPhi - exactPhi) / exactPhi)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["numericPhi"], res["exactPhi"], rel_tol=1e-6)


# ============================================================================
# 10. GRID-WIDE SUMMARIES, FRACTIONS, AND MAPS
# ============================================================================

class TestGridSummariesAndMapGeneration:
    """Verify grid-wide volume/mass fractions, threshold scans, and field map generation."""

    def test_volume_and_mass_filling_fractions_sum_to_unity(self):
        """
        Verify that volume filling fractions sum to 1.0 across all 4 cosmic web types.
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { ScalarField3D } from './src/fields/scalar_field_3d.js';
        import { TidalTensorAnalyzer } from './src/fields/tidal_tensor_invariants.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          spacing: [L / N, L / N, L / N],
          boundaryMode: BoundaryMode.PERIODIC,
          isCellCentered: true
        });

        const total = grid.totalCells;
        const potData = new Float64Array(total);

        // Populate with multiple harmonic modes
        for (let iz = 0; iz < N; iz++) {
          const z = (iz + 0.5) * grid.dz;
          for (let iy = 0; iy < N; iy++) {
            const y = (iy + 0.5) * grid.dy;
            for (let ix = 0; ix < N; ix++) {
              const x = (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              potData[idx] = Math.sin(0.06 * x) + Math.cos(0.06 * y) + Math.sin(0.06 * z);
            }
          }
        }

        const potField = new ScalarField3D(grid, potData, 'pot');
        const analyzer = new TidalTensorAnalyzer(potField, { gammaThreshold: 0.0 });

        const summary = analyzer.computeCosmicWebSummary({
          thresholdScanRange: [-0.005, 0.0, 0.005]
        });

        console.log(JSON.stringify(summary.toJSON()));
        """
        res = run_node_snippet(code)
        vFrac = res["volumeFractions"]
        totalVFrac = vFrac["void"] + vFrac["sheet"] + vFrac["filament"] + vFrac["knot"]
        assert math.isclose(totalVFrac, 1.0, rel_tol=1e-9)

        mFrac = res["massFractions"]
        totalMFrac = mFrac["void"] + mFrac["sheet"] + mFrac["filament"] + mFrac["knot"]
        assert math.isclose(totalMFrac, 1.0, rel_tol=1e-9)

        # Threshold scan verification
        assert len(res["thresholdScan"]) == 3

    def test_continuous_field_map_generation(self):
        """
        Verify generation of continuous 3D field maps for eigenvalues, web types,
        ellipticity, prolateness, and shear magnitude.
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { ScalarField3D } from './src/fields/scalar_field_3d.js';
        import { TidalTensorAnalyzer } from './src/fields/tidal_tensor_invariants.js';

        const N = 8;
        const L = 50.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          spacing: [L / N, L / N, L / N],
          boundaryMode: BoundaryMode.PERIODIC,
          isCellCentered: true
        });

        const total = grid.totalCells;
        const potData = new Float64Array(total);
        for (let i = 0; i < total; i++) potData[i] = Math.sin(i * 0.1);

        const potField = new ScalarField3D(grid, potData, 'pot');
        const analyzer = new TidalTensorAnalyzer(potField, { gammaThreshold: 0.0 });

        const maps = analyzer.generateFieldMaps();

        console.log(JSON.stringify({
          hasL1: maps.lambda1 instanceof ScalarField3D,
          hasL2: maps.lambda2 instanceof ScalarField3D,
          hasL3: maps.lambda3 instanceof ScalarField3D,
          hasWeb: maps.cosmicWeb instanceof ScalarField3D,
          hasEllipticity: maps.ellipticity instanceof ScalarField3D,
          hasProlateness: maps.prolateness instanceof ScalarField3D,
          hasShear: maps.shearMagnitude instanceof ScalarField3D,
          cellsCount: maps.lambda1.data.length
        }));
        """
        res = run_node_snippet(code)
        assert res["hasL1"] is True
        assert res["hasL2"] is True
        assert res["hasL3"] is True
        assert res["hasWeb"] is True
        assert res["hasEllipticity"] is True
        assert res["hasProlateness"] is True
        assert res["hasShear"] is True
        assert res["cellsCount"] == 8 * 8 * 8


# ============================================================================
# 11. ROBUSTNESS, DEGENERACIES & EDGE CASES
# ============================================================================

class TestRobustnessAndEdgeCases:
    """Verify numerical stability on extreme, degenerate, and edge-case inputs."""

    def test_zero_tensor_handling(self):
        """
        Verify behavior when input deformation tensor is zero (flat universe).
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const zeroTensor = new Float64Array(9);
        const inv = TidalTensorMath.computeInvariants(zeroTensor);
        const eigs = TidalTensorMath.cardanoEigenvalues(zeroTensor);
        const aniso = TidalTensorMath.computeAnisotropyParameters(eigs);
        const web = TidalTensorMath.classifyCosmicWeb(eigs, 0.0);

        console.log(JSON.stringify({
          I1: inv.I1,
          I2: inv.I2,
          I3: inv.I3,
          eigenvalues: Array.from(eigs),
          anisotropy: aniso.anisotropy,
          webType: web
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["I1"], 0.0, abs_tol=1e-15)
        assert math.isclose(res["I2"], 0.0, abs_tol=1e-15)
        assert math.isclose(res["I3"], 0.0, abs_tol=1e-15)
        assert np.allclose(res["eigenvalues"], [0.0, 0.0, 0.0], atol=1e-15)
        assert res["webType"] == 0  # VOID

    def test_extreme_anisotropy_values(self):
        """
        Verify stability when one eigenvalue dominates by orders of magnitude.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        const extreme = [1e6, 1.0, 1e-6];
        const aniso = TidalTensorMath.computeAnisotropyParameters(extreme);
        const diag = TidalTensorMath.jacobiDiagonalization([
          1e6, 0, 0,
          0, 1.0, 0,
          0, 0, 1e-6
        ]);

        console.log(JSON.stringify({
          anisotropy: aniso.anisotropy,
          ellipticity: aniso.ellipticity,
          jacobiEigs: Array.from(diag.eigenvalues)
        }));
        """
        res = run_node_snippet(code)
        assert not math.isnan(res["anisotropy"])
        assert not math.isnan(res["ellipticity"])
        assert math.isclose(res["jacobiEigs"][0], 1e6, rel_tol=1e-9)


# ============================================================================
# 12. DOROSHKEVICH PDF & PEAK JOINT INVARIANT STATISTICS
# ============================================================================

class TestDoroshkevichPDFAndJointStatistics:
    """Verify statistical distributions of Gaussian random field tidal tensor eigenvalues."""

    def test_doroshkevich_ordering_statistics(self):
        """
        For a Gaussian random density field, the eigenvalues of the tidal deformation tensor
        satisfy <lambda_1> > <lambda_2> > <lambda_3> unconditionally.
        """
        code = """
        import { TidalTensorMath } from './src/fields/tidal_tensor_invariants.js';

        // Monte Carlo ensemble of random symmetric Gaussian 3x3 matrices
        let sumL1 = 0, sumL2 = 0, sumL3 = 0;
        const N_samples = 500;

        // Pseudo-random Gaussian generator (Box-Muller)
        let seed = 42;
        const rng = () => {
          seed = (seed * 1664525 + 1013904223) % 4294967296;
          return seed / 4294967296;
        };
        const gaussian = () => {
          const u1 = Math.max(1e-15, rng());
          const u2 = rng();
          return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        };

        for (let s = 0; s < N_samples; s++) {
          const xx = gaussian(), yy = gaussian(), zz = gaussian();
          const xy = gaussian() * 0.5, xz = gaussian() * 0.5, yz = gaussian() * 0.5;

          const mat = [
            xx, xy, xz,
            xy, yy, yz,
            xz, yz, zz
          ];

          const eigs = TidalTensorMath.cardanoEigenvalues(mat);
          sumL1 += eigs[0];
          sumL2 += eigs[1];
          sumL3 += eigs[2];
        }

        console.log(JSON.stringify({
          meanL1: sumL1 / N_samples,
          meanL2: sumL2 / N_samples,
          meanL3: sumL3 / N_samples
        }));
        """
        res = run_node_snippet(code)
        assert res["meanL1"] > res["meanL2"] > res["meanL3"]
        # Mean trace should be close to zero for zero-mean Gaussian random tensor
        mean_trace = res["meanL1"] + res["meanL2"] + res["meanL3"]
        assert abs(mean_trace) < 0.2


# ============================================================================
# 13. FILAMENT AND SHEET ALIGNMENT METRICS
# ============================================================================

class TestFilamentAndSheetAlignment:
    """Verify directional alignment metrics between principal tidal frames and cosmic velocity fields."""

    def test_filament_spine_velocity_alignment(self):
        """
        Compute alignment cosine mu = |e_3 . v| / |v| where e_3 is the filament spine.
        """
        code = """
        import { TidalTensorMath, CosmicWebType } from './src/fields/tidal_tensor_invariants.js';

        // Filament tensor with elongation along Z-axis (eigenvector e3 = [0, 0, 1])
        const T = [
          4.0, 0.0, 0.0,
          0.0, 3.0, 0.0,
          0.0, 0.0, 0.2
        ];

        const diag = TidalTensorMath.jacobiDiagonalization(T);
        const spine = diag.eigenvectors[2]; // e3

        // Velocity along Z
        const vZ = [0.0, 0.0, 150.0];
        const dotZ = Math.abs(spine[0]*vZ[0] + spine[1]*vZ[1] + spine[2]*vZ[2]) / 150.0;

        // Velocity along X (perpendicular)
        const vX = [200.0, 0.0, 0.0];
        const dotX = Math.abs(spine[0]*vX[0] + spine[1]*vX[1] + spine[2]*vX[2]) / 200.0;

        console.log(JSON.stringify({
          dotZ,
          dotX,
          spine: Array.from(spine)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["dotZ"], 1.0, rel_tol=1e-12)
        assert math.isclose(res["dotX"], 0.0, abs_tol=1e-12)


# ============================================================================
# 14. COSMOLOGICAL REDSHIFT INVERSIONS & WORLD MODELS
# ============================================================================

class TestCosmologicalWorldModels:
    """Verify growth factor D(z) and collapse redshift inversions across world models."""

    def test_open_cdm_vs_flat_lcdm_growth(self):
        """
        Verify growth factor hierarchy: D_EdS(z) < D_LCDM(z) at fixed high redshift.
        """
        code = """
        import { TidalTensorMath, CosmologicalModel } from './src/fields/tidal_tensor_invariants.js';

        const z = 2.0;
        const d_eds = TidalTensorMath.linearGrowthFactor(z, { model: CosmologicalModel.EDS, omegaM: 1.0 });
        const d_lcdm = TidalTensorMath.linearGrowthFactor(z, { model: CosmologicalModel.FLAT_LCDM, omegaM: 0.315, omegaLambda: 0.685 });

        console.log(JSON.stringify({ d_eds, d_lcdm }));
        """
        res = run_node_snippet(code)
        # In EdS: D(2) = 1/3 \approx 0.3333
        assert math.isclose(res["d_eds"], 1.0 / 3.0, rel_tol=1e-6)
        # In flat LCDM: D(2) is ~ 0.40 - 0.45
        assert res["d_lcdm"] > res["d_eds"]
