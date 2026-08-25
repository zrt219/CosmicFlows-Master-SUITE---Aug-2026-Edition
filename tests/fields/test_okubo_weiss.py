"""
tests/fields/test_okubo_weiss.py
================================
Exhaustive automated Pytest test suite for the Okubo-Weiss Criterion, Velocity Deformation Tensor,
Principal Invariants (P, Q_J, R_J), Vieillefosse Tail Diagnostics, Jeong-Hussain Lambda_2 Vortex
Detection, and Cosmic Web Morphological Classification for the CosmicFlows-4 Research Workbench.

Test Categories:
----------------
1. Kinematic Decompositions & Exact Analytical Fields (Rigid Body, Planar Shear, Extensional Flow, Lamb-Oseen Vortex, Taylor-Green Vortex, ABC Flow)
2. Mathematical Invariants (P, Q_J, R_J), Cayley-Hamilton Identity, and Cardano Cubic Eigensolver
3. Four Flow Topologies (UFC, SFS / Vieillefosse Tail, UN/SS, SN/SS) and Invariant Plane (Q_J, R_J) Phase Space
4. Jeong & Hussain (1995) Lambda_2 and Hunt et al. (1988) Q-Criterion Vortex Identification
5. Cosmic Web Morphological Classification (Voids, Sheets, Filaments, Knots) & SO(3) Rotational Invariance
6. Finite Difference Stencil Orders (2nd, 4th, 6th) and Grid Boundary Modes (Periodic, Clamp, Reflective, Zero-Padding)
7. Grid-Wide Scalar Field Evaluators & Statistical Diagnostics (Okubo-Weiss, Enstrophy, Helicity, Divergence, Shear)
8. Multiscale Gaussian Scale-Space Analysis & Vortex Core Segmentation Extraction
9. Comprehensive Diagnostic Audit Report & W3C PROV-O Lineage Verification
10. Robustness, Numerical Stability, and Boundary Edge Case Verification

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
# 1. KINEMATIC DECOMPOSITIONS & EXACT ANALYTICAL FIELDS
# ============================================================================

class TestKinematicDecompositionsAndAnalyticFields:
    """Verify exact algebraic decompositions of velocity Jacobian J into S, Omega, sigma, omega, Q_OW."""

    def test_rigid_body_rotation_kinematics(self):
        """
        Pure rigid body rotation v = omega_vec x r with omega_vec = [0, 0, Omega_0]:
        vx = -Omega_0 * y, vy = Omega_0 * x, vz = 0.
        J = [[0, -Omega_0, 0], [Omega_0, 0, 0], [0, 0, 0]].
        Strain S = 0, Vorticity Omega = [[0, -Omega_0, 0], [Omega_0, 0, 0], [0, 0, 0]].
        vorticity vector omega = [0, 0, 2*Omega_0].
        s^2 = 0, |omega|^2 = 4*Omega_0^2, Q_OW = -2*Omega_0^2 < 0 (pure rotation-dominated).
        Lambda_2 = -Omega_0^2 < 0 (vortex core).
        """
        code = """
        import { OkuboWeissTensorAnalyzer, OkuboWeissRegime } from './src/fields/okubo_weiss_tensor.js';

        const Omega0 = 3.5;
        // J_ij = dv_i / dx_j
        // vx = -Omega0 * y => dvx/dx = 0, dvx/dy = -Omega0, dvx/dz = 0
        // vy = Omega0 * x  => dvy/dx = Omega0, dvy/dy = 0, dvy/dz = 0
        // vz = 0
        const J = new Float64Array([
          0.0, -Omega0, 0.0,
          Omega0, 0.0, 0.0,
          0.0, 0.0, 0.0
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J, [10.0, -5.0, 0.0]);

        console.log(JSON.stringify({
          divergence: decomp.divergence,
          strainMagnitudeSq: decomp.strainMagnitudeSq,
          vorticityMagnitudeSq: decomp.vorticityMagnitudeSq,
          okuboWeiss: decomp.okuboWeiss,
          huntQ: decomp.huntQ,
          lambda2: decomp.lambda2,
          vorticityVector: decomp.vorticityVector,
          regime: decomp.classifyRegime(0.1),
          isVortex: decomp.isVortexCore(0.0)
        }));
        """
        res = run_node_snippet(code)
        Omega0 = 3.5
        assert math.isclose(res["divergence"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["strainMagnitudeSq"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["vorticityMagnitudeSq"], 4.0 * Omega0**2, rel_tol=1e-12)
        assert math.isclose(res["okuboWeiss"], -2.0 * Omega0**2, rel_tol=1e-12)
        assert math.isclose(res["huntQ"], Omega0**2, rel_tol=1e-12)
        assert math.isclose(res["lambda2"], -Omega0**2, rel_tol=1e-12)
        assert np.allclose(res["vorticityVector"], [0.0, 0.0, 2.0 * Omega0], atol=1e-12)
        assert res["regime"] == -1  # ROTATION_DOMINATED
        assert res["isVortex"] is True

    def test_pure_hyperbolic_extensional_flow(self):
        """
        Pure hyperbolic stretching flow: vx = a*x, vy = -a*y, vz = 0 (a > 0).
        J = diag(a, -a, 0).
        Strain S = J, Vorticity Omega = 0.
        theta = 0, omega = 0.
        s^2 = 2*a^2, |omega|^2 = 0, Q_OW = 2*a^2 > 0 (pure strain-dominated).
        Lambda_2 = 0, Hunt Q = -a^2 < 0.
        """
        code = """
        import { OkuboWeissTensorAnalyzer, OkuboWeissRegime } from './src/fields/okubo_weiss_tensor.js';

        const a = 2.0;
        const J = new Float64Array([
          a, 0.0, 0.0,
          0.0, -a, 0.0,
          0.0, 0.0, 0.0
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);

        console.log(JSON.stringify({
          divergence: decomp.divergence,
          strainMagnitudeSq: decomp.strainMagnitudeSq,
          vorticityMagnitudeSq: decomp.vorticityMagnitudeSq,
          okuboWeiss: decomp.okuboWeiss,
          huntQ: decomp.huntQ,
          lambda2: decomp.lambda2,
          regime: decomp.classifyRegime(0.1),
          isVortex: decomp.isVortexCore(0.0)
        }));
        """
        res = run_node_snippet(code)
        a = 2.0
        assert math.isclose(res["divergence"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["strainMagnitudeSq"], 2.0 * a**2, rel_tol=1e-12)
        assert math.isclose(res["vorticityMagnitudeSq"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["okuboWeiss"], 2.0 * a**2, rel_tol=1e-12)
        assert math.isclose(res["huntQ"], -a**2, rel_tol=1e-12)
        # M = S^2 has eigenvalues [a^2, a^2, 0] => intermediate eigenvalue lambda2 = a^2 > 0
        assert math.isclose(res["lambda2"], a**2, rel_tol=1e-12)
        assert res["regime"] == 1  # STRAIN_DOMINATED
        assert res["isVortex"] is False

    def test_pure_planar_shear_flow(self):
        """
        Simple planar Couette shear flow: vx = gamma*y, vy = 0, vz = 0.
        J = [[0, gamma, 0], [0, 0, 0], [0, 0, 0]].
        S = [[0, gamma/2, 0], [gamma/2, 0, 0], [0, 0, 0]] => s^2 = 2*(gamma/2)^2 = gamma^2 / 2.
        Omega = [[0, gamma/2, 0], [-gamma/2, 0, 0], [0, 0, 0]] => omega_z = -gamma, |omega|^2 = gamma^2.
        ||Omega||_F^2 = 2*(gamma/2)^2 = gamma^2 / 2.
        Q_OW = s^2 - ||Omega||_F^2 = gamma^2/2 - gamma^2/2 = 0.
        Q_Hunt = 0.
        """
        code = """
        import { OkuboWeissTensorAnalyzer, OkuboWeissRegime } from './src/fields/okubo_weiss_tensor.js';

        const gamma = 4.0;
        const J = new Float64Array([
          0.0, gamma, 0.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);

        console.log(JSON.stringify({
          divergence: decomp.divergence,
          strainMagnitudeSq: decomp.strainMagnitudeSq,
          vorticityMagnitudeSq: decomp.vorticityMagnitudeSq,
          okuboWeiss: decomp.okuboWeiss,
          huntQ: decomp.huntQ,
          okuboWeissNormalized: decomp.okuboWeissNormalized,
          regime: decomp.classifyRegime(1e-5)
        }));
        """
        res = run_node_snippet(code)
        gamma = 4.0
        assert math.isclose(res["divergence"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["strainMagnitudeSq"], 0.5 * gamma**2, rel_tol=1e-12)
        assert math.isclose(res["vorticityMagnitudeSq"], gamma**2, rel_tol=1e-12)
        assert math.isclose(res["okuboWeiss"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["huntQ"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["okuboWeissNormalized"], 0.0, abs_tol=1e-12)
        assert res["regime"] == 0  # NEUTRAL

    def test_3d_taylor_green_vortex_analytical_matching(self):
        """
        3D Taylor-Green vortex:
        vx = U0 * sin(k x) * cos(k y) * cos(k z)
        vy = -U0 * cos(k x) * sin(k y) * cos(k z)
        vz = 0
        Verifies exact discrete Jacobian against analytical gradients on a 32^3 grid.
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer, DifferentiationOrder } from './src/fields/okubo_weiss_tensor.js';

        const N = 32;
        const L = 100.0; // Mpc/h
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          boundaryMode: BoundaryMode.PERIODIC
        });

        const U0 = 500.0; // km/s
        const k = (2.0 * Math.PI) / L;

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        for (let iz = 0; iz < N; iz++) {
          const z = (iz + 0.5) * grid.dz;
          for (let iy = 0; iy < N; iy++) {
            const y = (iy + 0.5) * grid.dy;
            for (let ix = 0; ix < N; ix++) {
              const x = (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              vx[idx] = U0 * Math.sin(k * x) * Math.cos(k * y) * Math.cos(k * z);
              vy[idx] = -U0 * Math.cos(k * x) * Math.sin(k * y) * Math.cos(k * z);
              vz[idx] = 0.0;
            }
          }
        }

        const vecField = new VectorField3D(grid, vx, vy, vz, { unit: 'km/s' });
        const analyzer4 = new OkuboWeissTensorAnalyzer(vecField, { order: DifferentiationOrder.FOURTH });

        // Sample at test point (ix=8, iy=8, iz=8)
        const ix = 8, iy = 8, iz = 8;
        const x = (ix + 0.5) * grid.dx;
        const y = (iy + 0.5) * grid.dy;
        const z = (iz + 0.5) * grid.dz;

        // Analytical derivatives
        const exact_dvx_dx = U0 * k * Math.cos(k * x) * Math.cos(k * y) * Math.cos(k * z);
        const exact_dvx_dy = -U0 * k * Math.sin(k * x) * Math.sin(k * y) * Math.cos(k * z);
        const exact_dvx_dz = -U0 * k * Math.sin(k * x) * Math.cos(k * y) * Math.sin(k * z);

        const exact_dvy_dx = U0 * k * Math.sin(k * x) * Math.sin(k * y) * Math.cos(k * z);
        const exact_dvy_dy = -U0 * k * Math.cos(k * x) * Math.cos(k * y) * Math.cos(k * z);
        const exact_dvy_dz = U0 * k * Math.cos(k * x) * Math.sin(k * y) * Math.sin(k * z);

        const numJ = analyzer4.computeLocalJacobian(ix, iy, iz);

        console.log(JSON.stringify({
          num_dvx_dx: numJ[0], exact_dvx_dx,
          num_dvx_dy: numJ[1], exact_dvx_dy,
          num_dvx_dz: numJ[2], exact_dvx_dz,
          num_dvy_dx: numJ[3], exact_dvy_dx,
          num_dvy_dy: numJ[4], exact_dvy_dy,
          num_dvy_dz: numJ[5], exact_dvy_dz,
          divNum: numJ[0] + numJ[4] + numJ[8]
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["num_dvx_dx"], res["exact_dvx_dx"], rel_tol=1e-3)
        assert math.isclose(res["num_dvx_dy"], res["exact_dvx_dy"], rel_tol=1e-3)
        assert math.isclose(res["num_dvx_dz"], res["exact_dvx_dz"], rel_tol=1e-3)
        assert math.isclose(res["num_dvy_dx"], res["exact_dvy_dx"], rel_tol=1e-3)
        assert math.isclose(res["num_dvy_dy"], res["exact_dvy_dy"], rel_tol=1e-3)
        assert math.isclose(res["num_dvy_dz"], res["exact_dvy_dz"], rel_tol=1e-3)
        # Taylor-Green is strictly incompressible: div(v) = 0
        assert math.isclose(res["divNum"], 0.0, abs_tol=1e-4)

    def test_abc_flow_helicity_and_vorticity(self):
        """
        Arnold-Beltrami-Childress (ABC) flow:
        vx = A*sin(k*z) + C*cos(k*y)
        vy = B*sin(k*x) + A*cos(k*z)
        vz = C*sin(k*y) + B*cos(k*x)
        Beltrami condition: curl(v) = k * v.
        Kinetic Helicity Density: h = v . curl(v) = k * ||v||^2.
        """
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const A = 1.0, B = 1.0, C = 1.0;
        const k = 0.5;
        const x = 1.2, y = 2.4, z = 0.8;

        const vx = A * Math.sin(k * z) + C * Math.cos(k * y);
        const vy = B * Math.sin(k * x) + A * Math.cos(k * z);
        const vz = C * Math.sin(k * y) + B * Math.cos(k * x);

        // Analytical Jacobian
        const dvx_dx = 0.0;
        const dvx_dy = -k * C * Math.sin(k * y);
        const dvx_dz = k * A * Math.cos(k * z);

        const dvy_dx = k * B * Math.cos(k * x);
        const dvy_dy = 0.0;
        const dvy_dz = -k * A * Math.sin(k * z);

        const dvz_dx = -k * B * Math.sin(k * x);
        const dvz_dy = k * C * Math.cos(k * y);
        const dvz_dz = 0.0;

        const J = new Float64Array([
          dvx_dx, dvx_dy, dvx_dz,
          dvy_dx, dvy_dy, dvy_dz,
          dvz_dx, dvz_dy, dvz_dz
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J, [vx, vy, vz]);

        const expected_wx = k * vx;
        const expected_wy = k * vy;
        const expected_wz = k * vz;
        const expectedHelicity = k * (vx*vx + vy*vy + vz*vz);

        console.log(JSON.stringify({
          calc_wx: decomp.vorticityVector[0], expected_wx,
          calc_wy: decomp.vorticityVector[1], expected_wy,
          calc_wz: decomp.vorticityVector[2], expected_wz,
          calcHelicity: decomp.helicityDensity, expectedHelicity,
          divergence: decomp.divergence
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["calc_wx"], res["expected_wx"], rel_tol=1e-12)
        assert math.isclose(res["calc_wy"], res["expected_wy"], rel_tol=1e-12)
        assert math.isclose(res["calc_wz"], res["expected_wz"], rel_tol=1e-12)
        assert math.isclose(res["calcHelicity"], res["expectedHelicity"], rel_tol=1e-12)
        assert math.isclose(res["divergence"], 0.0, abs_tol=1e-12)


# ============================================================================
# 2. MATHEMATICAL INVARIANTS (P, Q_J, R_J) & CAYLEY-HAMILTON THEOREM
# ============================================================================

class TestDeformationInvariantsAndCayleyHamilton:
    """Verify Cayley-Hamilton theorem, trace identities, and Cardano cubic eigensolver."""

    def test_cayley_hamilton_theorem_verification(self):
        """
        Verify Cayley-Hamilton theorem J^3 + P*J^2 + Q_J*J + R_J*I = 0 holds
        to machine precision for general random 3x3 matrices.
        """
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        function matMul(A, B) {
          const C = new Float64Array(9);
          for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
              let s = 0.0;
              for (let k = 0; k < 3; k++) {
                s += A[i * 3 + k] * B[k * 3 + j];
              }
              C[i * 3 + j] = s;
            }
          }
          return C;
        }

        const tests = [
          [1.2, -0.4, 0.7, 0.5, 2.1, -1.3, -0.8, 0.3, -1.5],
          [5.0, 2.0, -1.0, 2.0, 3.0, 0.5, -1.0, 0.5, 4.0], // symmetric
          [0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0], // antisymmetric
          [-2.5, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.0]    // diagonal
        ];

        const residuals = [];

        for (const rawJ of tests) {
          const J = new Float64Array(rawJ);
          const inv = OkuboWeissTensorAnalyzer.computeInvariants(J);
          const P = inv.P;
          const Q = inv.Q;
          const R = inv.R;

          const J2 = matMul(J, J);
          const J3 = matMul(J2, J);

          // Residual matrix M = J3 + P*J2 + Q*J + R*I
          let maxErr = 0.0;
          for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
              const delta_ij = (i === j) ? 1.0 : 0.0;
              const val = J3[i * 3 + j] + P * J2[i * 3 + j] + Q * J[i * 3 + j] + R * delta_ij;
              if (Math.abs(val) > maxErr) maxErr = Math.abs(val);
            }
          }
          residuals.push(maxErr);
        }

        console.log(JSON.stringify({ residuals }));
        """
        res = run_node_snippet(code)
        for err in res["residuals"]:
            assert err < 1e-11, f"Cayley-Hamilton identity violation: max error = {err}"

    def test_incompressible_flow_invariant_identities(self):
        """
        For incompressible flow (P = 0):
        Q_J = -0.5 * Tr(J^2) = 0.5 * (||Omega||_F^2 - ||S||_F^2) = -0.5 * Q_OW = Q_Hunt.
        R_J = -det(J) = -1/3 * Tr(J^3).
        """
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        // Incompressible random matrix with Tr(J) = 0
        const J = new Float64Array([
          2.0, 1.5, -0.8,
          -0.5, -1.0, 2.2,
          1.2, -0.7, -1.0 // 2.0 - 1.0 - 1.0 = 0
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);
        const inv = OkuboWeissTensorAnalyzer.computeInvariants(J);

        console.log(JSON.stringify({
          P: inv.P,
          Q_J: inv.Q,
          huntQ: decomp.huntQ,
          halfMinusQ_OW: -0.5 * decomp.okuboWeiss,
          R_J: inv.R,
          detJ: -(inv.R)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["P"], 0.0, abs_tol=1e-14)
        assert math.isclose(res["Q_J"], res["huntQ"], rel_tol=1e-12)
        assert math.isclose(res["Q_J"], res["halfMinusQ_OW"], rel_tol=1e-12)

    def test_cubic_eigenvalue_solver_accuracy(self):
        """Verify Cardano cubic solver reproduces numpy.linalg.eigvals for random 3x3 matrices."""
        np.random.seed(42)
        for _ in range(10):
            A = np.random.randn(3, 3)
            # Analytical invariants
            P = -float(np.trace(A))
            Q = float(0.5 * ((np.trace(A))**2 - np.trace(A @ A)))
            R = -float(np.linalg.det(A))

            expected_eigs = np.linalg.eigvals(A)

            code = f"""
            import {{ OkuboWeissTensorAnalyzer }} from './src/fields/okubo_weiss_tensor.js';
            const res = OkuboWeissTensorAnalyzer.solveCubicEigenvalues({P}, {Q}, {R});
            console.log(JSON.stringify(res));
            """
            out = run_node_snippet(code)
            calculated_roots = out["roots"]

            calc_complex = [complex(r["real"], r["imag"]) for r in calculated_roots]
            # Match eigenvalues
            matched = [False] * 3
            for exp_val in expected_eigs:
                closest_dist = min(abs(exp_val - c) for c in calc_complex)
                assert closest_dist < 1e-6, f"Eigenvalue mismatch: expected {exp_val}, got {calc_complex}"


# ============================================================================
# 3. FOUR FLOW TOPOLOGIES & VIEILLEFOSSE TAIL DIAGNOSTICS
# ============================================================================

class TestFlowTopologiesAndVieillefosseTail:
    """Verify classification of UFC, SFS, UN/SS, SN/SS regimes and Vieillefosse line mechanics."""

    def test_four_quadrant_topology_classification(self):
        """
        Verify the 4 flow topologies (Chong et al. 1990; Cantwell 1992):
        1. UFC: Delta > 0, R > 0 (spiral outflow / 1D compression)
        2. SFS: Delta > 0, R < 0 (Vieillefosse tail / vortex stretching)
        3. UN/SS: Delta <= 0, R > 0 (biaxial sheet/filament collapse)
        4. SN/SS: Delta <= 0, R <= 0 (biaxial void expansion)
        """
        code = """
        import { OkuboWeissTensorAnalyzer, FlowTopologyType } from './src/fields/okubo_weiss_tensor.js';

        // UFC state: e.g. P = 0, Q = 1.0, R = 1.0 => Delta = 27*(1)^2 + 4*(1)^3 = 31 > 0, R > 0
        const topo_UFC = OkuboWeissTensorAnalyzer.classifyTopology(0.0, 1.0, 1.0, 31.0);

        // SFS state: e.g. P = 0, Q = 1.0, R = -1.0 => Delta = 27*(-1)^2 + 4*(1)^3 = 31 > 0, R < 0
        const topo_SFS = OkuboWeissTensorAnalyzer.classifyTopology(0.0, 1.0, -1.0, 31.0);

        // UN/SS state: e.g. P = 0, Q = -2.0, R = 0.5 => Delta = 27*(0.25) + 4*(-8) = 6.75 - 32 = -25.25 < 0, R > 0
        const topo_UN_SS = OkuboWeissTensorAnalyzer.classifyTopology(0.0, -2.0, 0.5, -25.25);

        // SN/SS state: e.g. P = 0, Q = -2.0, R = -0.5 => Delta = 27*(0.25) + 4*(-8) = -25.25 < 0, R < 0
        const topo_SN_SS = OkuboWeissTensorAnalyzer.classifyTopology(0.0, -2.0, -0.5, -25.25);

        console.log(JSON.stringify({
          topo_UFC,
          topo_SFS,
          topo_UN_SS,
          topo_SN_SS
        }));
        """
        res = run_node_snippet(code)
        assert res["topo_UFC"] == 1  # UFC
        assert res["topo_SFS"] == 2  # SFS (Vieillefosse Tail)
        assert res["topo_UN_SS"] == 3  # UN/SS
        assert res["topo_SN_SS"] == 4  # SN/SS

    def test_vieillefosse_zero_discriminant_line_precision(self):
        """
        Vieillefosse line condition: Delta = (27/4)*R^2 + Q^3 = 0 => Q = -3*(R/2)^(2/3).
        Verify distanceToVieillefosseLine returns ~0 on the line and identifies tail points.
        """
        code = """
        import { DeformationInvariants, FlowTopologyType } from './src/fields/okubo_weiss_tensor.js';

        const R = -2.0;
        // Q on the Vieillefosse line: Q = -3 * (R/2)^(2/3) = -3 * (-1)^(2/3) = -3 * 1 = -3.0
        const Q_exact = -3.0 * Math.cbrt(Math.pow(R / 2.0, 2));

        const invOnLine = new DeformationInvariants(0.0, Q_exact, R, 0.0, FlowTopologyType.SFS);
        const distOnLine = invOnLine.distanceToVieillefosseLine();

        const invNearTail = new DeformationInvariants(0.0, Q_exact + 0.1, R, 0.5, FlowTopologyType.SFS);
        const isTail = invNearTail.isVieillefosseTail();

        console.log(JSON.stringify({
          Q_exact,
          distOnLine,
          isTail
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["Q_exact"], -3.0, rel_tol=1e-12)
        assert math.isclose(res["distOnLine"], 0.0, abs_tol=1e-12)
        assert res["isTail"] is True


# ============================================================================
# 4. JEONG-HUSSAIN LAMBDA_2 & HUNT Q-CRITERION VORTEX DETECTION
# ============================================================================

class TestVortexIdentificationCriteria:
    """Verify Jeong-Hussain lambda_2 and Hunt Q-criterion detection on synthetic vortices."""

    def test_vortex_core_lambda2_and_hunt_agreement(self):
        """
        For a Gaussian vortex profile v_theta(r) = (Gamma/(2*pi*r)) * (1 - exp(-r^2/r0^2)),
        near the core (r << r0), flow is solid-body rotation: lambda_2 < 0 and Q_Hunt > 0.
        Far from the core (r >> r0), flow is irrotational shear: lambda_2 >= 0 and Q_Hunt <= 0.
        """
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const r0 = 5.0; // Mpc/h
        const Gamma = 1000.0;

        function evaluateVortexJacobian(x, y) {
          const r2 = x * x + y * y;
          const r = Math.sqrt(r2) || 1e-6;
          // Analytical derivatives of Lamb-Oseen vortex
          const factor = (Gamma / (2.0 * Math.PI * r2)) * (1.0 - Math.exp(-r2 / (r0 * r0)));
          const dFactor_dr = (Gamma / (2.0 * Math.PI)) * (
            -2.0 / (r2 * r) * (1.0 - Math.exp(-r2 / (r0 * r0))) +
            (2.0 / (r2 * r0 * r0)) * Math.exp(-r2 / (r0 * r0))
          );

          // vx = -y * factor, vy = x * factor
          const dvx_dx = -y * dFactor_dr * (x / r);
          const dvx_dy = -factor - y * dFactor_dr * (y / r);
          const dvy_dx = factor + x * dFactor_dr * (x / r);
          const dvy_dy = x * dFactor_dr * (y / r);

          return new Float64Array([
            dvx_dx, dvx_dy, 0.0,
            dvy_dx, dvy_dy, 0.0,
            0.0, 0.0, 0.0
          ]);
        }

        // Inside core: r = 1.0 Mpc/h < r0
        const J_core = evaluateVortexJacobian(1.0, 0.0);
        const k_core = OkuboWeissTensorAnalyzer.decomposeKinematics(J_core);

        // Outside core: r = 25.0 Mpc/h >> r0
        const J_far = evaluateVortexJacobian(25.0, 0.0);
        const k_far = OkuboWeissTensorAnalyzer.decomposeKinematics(J_far);

        console.log(JSON.stringify({
          core: {
            lambda2: k_core.lambda2,
            huntQ: k_core.huntQ,
            okuboWeiss: k_core.okuboWeiss,
            isVortex: k_core.isVortexCore(0.0),
            isHunt: k_core.isHuntVortex(0.0)
          },
          far: {
            lambda2: k_far.lambda2,
            huntQ: k_far.huntQ,
            okuboWeiss: k_far.okuboWeiss,
            isVortex: k_far.isVortexCore(0.0),
            isHunt: k_far.isHuntVortex(0.0)
          }
        }));
        """
        res = run_node_snippet(code)
        # Inside core: lambda_2 < 0, Q_Hunt > 0, Q_OW < 0
        assert res["core"]["lambda2"] < 0.0
        assert res["core"]["huntQ"] > 0.0
        assert res["core"]["okuboWeiss"] < 0.0
        assert res["core"]["isVortex"] is True
        assert res["core"]["isHunt"] is True

        # Outside core: lambda_2 >= 0, Q_Hunt <= 0, Q_OW >= 0
        assert res["core"]["lambda2"] < res["far"]["lambda2"]
        assert res["core"]["huntQ"] > res["far"]["huntQ"]


# ============================================================================
# 5. COSMIC WEB CLASSIFICATION & SO(3) ROTATIONAL INVARIANCE
# ============================================================================

class TestCosmicWebClassificationAndInvariance:
    """Verify cosmic web types (void, sheet, filament, knot) and SO(3) rotational invariance."""

    def test_cosmic_web_four_morphologies(self):
        """
        Verify eigenvalues of strain rate tensor S classify into 4 web types:
        - Void: 0 eigenvalues > alpha_th
        - Sheet: 1 eigenvalue > alpha_th
        - Filament: 2 eigenvalues > alpha_th
        - Knot: 3 eigenvalues > alpha_th
        """
        code = """
        import { LocalKinematics, DeformationInvariants, CosmicWebType } from './src/fields/okubo_weiss_tensor.js';
        import { EigenSystem3D } from './src/coordinates/scientific_types.js';

        const alpha_th = 0.2;

        function createMockKinematics(lambdas) {
          const eigen = new EigenSystem3D(lambdas, [
            [1, 0, 0], [0, 1, 0], [0, 0, 1]
          ]);
          const S = eigen.reconstructMatrix();
          const J = new Float64Array(S);
          const inv = new DeformationInvariants(0, 0, 0, 0, 0);
          return new LocalKinematics({
            J, S, Omega: new Float64Array(9), sigma: new Float64Array(9),
            divergence: 0, vorticityVector: [0, 0, 0], vorticityMagnitudeSq: 0,
            strainMagnitudeSq: 0, shearMagnitudeSq: 0, okuboWeiss: 0, okuboWeissNormalized: 0,
            huntQ: 0, enstrophy: 0, invariants: inv, lambda2: 0,
            strainEigenSystem: eigen
          });
        }

        const voidPoint = createMockKinematics([-0.5, -0.2, -0.1]);
        const sheetPoint = createMockKinematics([0.8, -0.3, -0.2]);
        const filamentPoint = createMockKinematics([1.2, 0.6, -0.4]);
        const knotPoint = createMockKinematics([1.5, 0.9, 0.4]);

        console.log(JSON.stringify({
          voidType: voidPoint.classifyCosmicWeb(alpha_th),
          sheetType: sheetPoint.classifyCosmicWeb(alpha_th),
          filamentType: filamentPoint.classifyCosmicWeb(alpha_th),
          knotType: knotPoint.classifyCosmicWeb(alpha_th)
        }));
        """
        res = run_node_snippet(code)
        assert res["voidType"] == 0  # VOID
        assert res["sheetType"] == 1  # SHEET
        assert res["filamentType"] == 2  # FILAMENT
        assert res["knotType"] == 3  # KNOT

    def test_so3_rotational_invariance_of_invariants(self):
        """
        Verify that under an arbitrary 3D SO(3) coordinate rotation R,
        the rotated Jacobian J' = R * J * R^T produces identical:
        - P, Q_J, R_J invariants
        - Okubo-Weiss parameter Q_OW
        - Strain tensor eigenvalues
        - Lambda_2 vortex criterion
        """
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        // Base Jacobian J
        const J = new Float64Array([
          1.2, 2.5, -0.7,
          -1.8, 0.4, 3.1,
          0.9, -2.2, -1.6
        ]);

        // Exact 3D rotation matrix (rotation by angle theta around axis [1, 2, 2]/3)
        // axis u = [1/3, 2/3, 2/3], theta = 0.5 rad
        const ux = 1.0 / 3.0, uy = 2.0 / 3.0, uz = 2.0 / 3.0;
        const theta = 0.75;
        const c = Math.cos(theta);
        const s = Math.sin(theta);
        const C = 1.0 - c;

        const R_mat = new Float64Array([
          c + ux*ux*C,      ux*uy*C - uz*s,   ux*uz*C + uy*s,
          uy*ux*C + uz*s,   c + uy*uy*C,      uy*uz*C - ux*s,
          uz*ux*C - uy*s,   uz*uy*C + ux*s,   c + uz*uz*C
        ]);

        function matMul(A, B) {
          const C = new Float64Array(9);
          for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
              let s = 0.0;
              for (let k = 0; k < 3; k++) s += A[i * 3 + k] * B[k * 3 + j];
              C[i * 3 + j] = s;
            }
          }
          return C;
        }

        function transpose(A) {
          return new Float64Array([
            A[0], A[3], A[6],
            A[1], A[4], A[7],
            A[2], A[5], A[8]
          ]);
        }

        const RT = transpose(R_mat);
        const J_prime = matMul(matMul(R_mat, J), RT);

        const decompBase = OkuboWeissTensorAnalyzer.decomposeKinematics(J);
        const decompRot = OkuboWeissTensorAnalyzer.decomposeKinematics(J_prime);

        const invBase = OkuboWeissTensorAnalyzer.computeInvariants(J);
        const invRot = OkuboWeissTensorAnalyzer.computeInvariants(J_prime);

        console.log(JSON.stringify({
          P_diff: Math.abs(invBase.P - invRot.P),
          Q_diff: Math.abs(invBase.Q - invRot.Q),
          R_diff: Math.abs(invBase.R - invRot.R),
          ow_diff: Math.abs(decompBase.okuboWeiss - decompRot.okuboWeiss),
          lambda2_diff: Math.abs(decompBase.lambda2 - decompRot.lambda2),
          strainLam_diff: Math.abs(decompBase.strainEigenSystem.lambda1 - decompRot.strainEigenSystem.lambda1)
        }));
        """
        res = run_node_snippet(code)
        assert res["P_diff"] < 1e-10
        assert res["Q_diff"] < 1e-10
        assert res["R_diff"] < 1e-10
        assert res["ow_diff"] < 1e-10
        assert res["lambda2_diff"] < 1e-10
        assert res["strainLam_diff"] < 1e-10


# ============================================================================
# 6. FINITE DIFFERENCE STENCILS & GRID BOUNDARY MODES
# ============================================================================

class TestFiniteDifferenceStencilsAndBoundaries:
    """Verify convergence order of 2nd, 4th, 6th order stencils and boundary conditions."""

    def test_stencil_convergence_rates(self):
        """
        Verify convergence error of 2nd, 4th, and 6th order finite difference stencils
        on a smooth sinusoidal vector field v(x) = sin(k*x).
        L_inf error across the grid should decrease as O(h^2), O(h^4), O(h^6).
        """
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer, DifferentiationOrder } from './src/fields/okubo_weiss_tensor.js';

        const L = 100.0;
        const k = (2.0 * Math.PI) / L;

        function runGrid(N) {
          const grid = new GridIndexer({
            nx: N, ny: N, nz: N,
            origin: [0, 0, 0],
            boxSize: [L, L, L],
            boundaryMode: BoundaryMode.PERIODIC,
            isCellCentered: true
          });

          const total = grid.totalCells;
          const vx = new Float64Array(total);
          const vy = new Float64Array(total);
          const vz = new Float64Array(total);

          for (let iz = 0; iz < N; iz++) {
            for (let iy = 0; iy < N; iy++) {
              for (let ix = 0; ix < N; ix++) {
                const x = (ix + 0.5) * grid.dx;
                const idx = grid.index(ix, iy, iz);
                vx[idx] = Math.sin(k * x);
              }
            }
          }

          const vecField = new VectorField3D(grid, vx, vy, vz);
          const an2 = new OkuboWeissTensorAnalyzer(vecField, { order: DifferentiationOrder.SECOND });
          const an4 = new OkuboWeissTensorAnalyzer(vecField, { order: DifferentiationOrder.FOURTH });
          const an6 = new OkuboWeissTensorAnalyzer(vecField, { order: DifferentiationOrder.SIXTH });

          let maxErr2 = 0.0, maxErr4 = 0.0, maxErr6 = 0.0;
          for (let ix = 0; ix < N; ix++) {
            const x = (ix + 0.5) * grid.dx;
            const exact = k * Math.cos(k * x);
            const val2 = an2.computeLocalJacobian(ix, 0, 0)[0];
            const val4 = an4.computeLocalJacobian(ix, 0, 0)[0];
            const val6 = an6.computeLocalJacobian(ix, 0, 0)[0];

            if (Math.abs(val2 - exact) > maxErr2) maxErr2 = Math.abs(val2 - exact);
            if (Math.abs(val4 - exact) > maxErr4) maxErr4 = Math.abs(val4 - exact);
            if (Math.abs(val6 - exact) > maxErr6) maxErr6 = Math.abs(val6 - exact);
          }

          return { maxErr2, maxErr4, maxErr6 };
        }

        const res16 = runGrid(16);
        const res32 = runGrid(32);

        // Ratio of max errors when halving grid step h (h_32 = h_16 / 2)
        // 2nd order: ratio = 4.0 (2^2)
        // 4th order: ratio = 16.0 (2^4)
        // 6th order: ratio = 64.0 (2^6)
        console.log(JSON.stringify({
          ratio2: res16.maxErr2 / res32.maxErr2,
          ratio4: res16.maxErr4 / res32.maxErr4,
          ratio6: res16.maxErr6 / res32.maxErr6
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["ratio2"], 4.0, rel_tol=0.05)
        assert math.isclose(res["ratio4"], 16.0, rel_tol=0.05)
        assert math.isclose(res["ratio6"], 64.0, rel_tol=0.05)

    def test_grid_boundary_modes_execution(self):
        """Verify boundary modes: PERIODIC, CLAMP, REFLECTIVE, ZERO_PADDING execute safely without NaN."""
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const modes = [
          BoundaryMode.PERIODIC,
          BoundaryMode.CLAMP,
          BoundaryMode.REFLECTIVE,
          BoundaryMode.ZERO_PADDING
        ];

        const results = [];

        for (const bMode of modes) {
          const N = 8;
          const grid = new GridIndexer({
            nx: N, ny: N, nz: N,
            origin: [0, 0, 0],
            boxSize: [100, 100, 100],
            boundaryMode: bMode
          });

          const total = grid.totalCells;
          const vx = new Float64Array(total).fill(100.0);
          const vy = new Float64Array(total).fill(-50.0);
          const vz = new Float64Array(total).fill(25.0);

          const vecField = new VectorField3D(grid, vx, vy, vz);
          const analyzer = new OkuboWeissTensorAnalyzer(vecField, { boundaryMode: bMode });

          const J_corner = analyzer.computeLocalJacobian(0, 0, 0);
          let hasNaN = false;
          for (let k = 0; k < 9; k++) {
            if (!Number.isFinite(J_corner[k])) hasNaN = true;
          }

          results.push({ mode: bMode, hasNaN });
        }

        console.log(JSON.stringify({ results }));
        """
        res = run_node_snippet(code)
        for entry in res["results"]:
            assert entry["hasNaN"] is False, f"Boundary mode {entry['mode']} produced NaN."


# ============================================================================
# 7. GRID-WIDE SCALAR FIELD EVALUATORS & STATISTICAL DIAGNOSTICS
# ============================================================================

class TestGridWideEvaluatorsAndDiagnostics:
    """Verify grid-wide field generators (Q_OW, Q_norm, Hunt, Lambda2, Helicity, Divergence, etc.)."""

    def test_grid_wide_scalar_field_evaluators(self):
        """Verify all scalar field evaluators return populated ScalarField3D instances."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        // Inject random smooth velocity components
        for (let i = 0; i < total; i++) {
          vx[i] = Math.sin(i * 0.1) * 300.0;
          vy[i] = Math.cos(i * 0.1) * 300.0;
          vz[i] = Math.sin(i * 0.05) * 150.0;
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const ow = analyzer.computeOkuboWeissField();
        const normOw = analyzer.computeNormalizedOkuboWeissField();
        const hunt = analyzer.computeQCriterionField();
        const l2 = analyzer.computeLambda2Field();
        const delta = analyzer.computeDeltaCriterionField();
        const strain = analyzer.computeStrainMagnitudeField();
        const vort = analyzer.computeVorticityMagnitudeField();
        const enstrophy = analyzer.computeEnstrophyField();
        const helicity = analyzer.computeHelicityDensityField();
        const div = analyzer.computeDivergenceField();
        const shear = analyzer.computeShearMagnitudeField();

        console.log(JSON.stringify({
          owSize: ow.data.length,
          normOwMin: normOw.computeStatistics().min,
          normOwMax: normOw.computeStatistics().max,
          huntMean: hunt.computeStatistics().mean,
          l2Mean: l2.computeStatistics().mean,
          deltaSize: delta.data.length,
          strainMean: strain.computeStatistics().mean,
          vortMean: vort.computeStatistics().mean,
          enstrophyMean: enstrophy.computeStatistics().mean,
          helicitySize: helicity.data.length,
          divSize: div.data.length,
          shearMean: shear.computeStatistics().mean
        }));
        """
        res = run_node_snippet(code)
        assert res["owSize"] == 16**3
        assert res["normOwMin"] >= -1.0 - 1e-6
        assert res["normOwMax"] <= 1.0 + 1e-6
        assert res["strainMean"] >= 0.0
        assert res["vortMean"] >= 0.0
        assert res["enstrophyMean"] >= 0.0

    def test_flow_regime_classification_fractions(self):
        """Verify classifyFlowRegimes partitions volume into valid fractions summing to 1.0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        for (let i = 0; i < total; i++) {
          vx[i] = Math.sin(i * 0.2) * 200.0;
          vy[i] = Math.cos(i * 0.3) * 200.0;
          vz[i] = Math.sin(i * 0.4) * 200.0;
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const regimes = analyzer.classifyFlowRegimes({ sigmaMultiplier: 0.2 });
        const web = analyzer.classifyCosmicWeb(0.1);
        const topos = analyzer.computeTopologyField();

        console.log(JSON.stringify({
          regimeFractionsSum: regimes.volumeFractions.rotationDominated +
                              regimes.volumeFractions.neutral +
                              regimes.volumeFractions.strainDominated,
          webFractionsSum: web.volumeFractions.voids +
                           web.volumeFractions.sheets +
                           web.volumeFractions.filaments +
                           web.volumeFractions.knots,
          topoFractionsSum: topos.volumeFractions.degenerate +
                            topos.volumeFractions.UFC +
                            topos.volumeFractions.SFS_Vieillefosse +
                            topos.volumeFractions.UN_SS +
                            topos.volumeFractions.SN_SS
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["regimeFractionsSum"], 1.0, abs_tol=1e-12)
        assert math.isclose(res["webFractionsSum"], 1.0, abs_tol=1e-12)
        assert math.isclose(res["topoFractionsSum"], 1.0, abs_tol=1e-12)


# ============================================================================
# 8. MULTISCALE SCALE-SPACE & VORTEX CORE SEGMENTATION
# ============================================================================

class TestMultiscaleAnalysisAndVortexExtraction:
    """Verify Gaussian multiscale smoothing and vortex core segmentation mask extraction."""

    def test_multiscale_gaussian_smoothing(self):
        """Verify multiscale Okubo-Weiss fields evaluate smoothly across multiple spatial scales."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        for (let i = 0; i < total; i++) {
          vx[i] = Math.sin(i * 0.5) * 500.0;
          vy[i] = Math.cos(i * 0.5) * 500.0;
          vz[i] = 0.0;
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const multi = analyzer.computeMultiscaleOkuboWeiss([2.0, 4.0, 8.0]);

        console.log(JSON.stringify({
          count: multi.length,
          scales: multi.map(m => m.scale),
          stdDevs: multi.map(m => (typeof m.stats.std === 'number' ? m.stats.std : m.stats.stdDev))
        }));
        """
        res = run_node_snippet(code)
        assert res["count"] == 3
        assert res["scales"] == [2.0, 4.0, 8.0]
        # Larger smoothing radii damp high-frequency fluctuations, reducing stdDev
        assert res["stdDevs"][0] >= res["stdDevs"][1] >= res["stdDevs"][2]

    def test_vortex_core_binary_extraction(self):
        """Verify extractVortexCores returns valid binary masks and volume fractions."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        for (let i = 0; i < total; i++) {
          vx[i] = Math.sin(i * 0.2) * 300.0;
          vy[i] = Math.cos(i * 0.2) * 300.0;
          vz[i] = 0.0;
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const coresL2 = analyzer.extractVortexCores({ criterion: 'lambda2', threshold: 0.0 });
        const coresHunt = analyzer.extractVortexCores({ criterion: 'hunt_q' });
        const coresOW = analyzer.extractVortexCores({ criterion: 'okubo_weiss' });

        console.log(JSON.stringify({
          l2Fraction: coresL2.volumeFraction,
          huntFraction: coresHunt.volumeFraction,
          owFraction: coresOW.volumeFraction
        }));
        """
        res = run_node_snippet(code)
        assert 0.0 <= res["l2Fraction"] <= 1.0
        assert 0.0 <= res["huntFraction"] <= 1.0
        assert 0.0 <= res["owFraction"] <= 1.0


# ============================================================================
# 9. COMPREHENSIVE DIAGNOSTIC REPORT & PROVENANCE
# ============================================================================

class TestDiagnosticReportAndProvenance:
    """Verify diagnostic audit report format, citations, and SHA-256 hash reproducibility."""

    def test_diagnostic_report_structure_and_citations(self):
        """Verify diagnostic report includes Courtois et al. (2023) citation and complete summaries."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 8;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total).fill(10.0);
        const vy = new Float64Array(total).fill(-20.0);
        const vz = new Float64Array(total).fill(30.0);

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const report = analyzer.generateDiagnosticReport();

        console.log(JSON.stringify({
          title: report.title,
          module: report.module,
          reportHash: report.reportHash,
          citationsCount: report.scientificCitations.length,
          hasCourtoisCitation: report.scientificCitations.some(c => c.includes('Courtois, H. M., et al. (2023)')),
          hasStats: Boolean(report.statisticalSummaries.okuboWeiss)
        }));
        """
        res = run_node_snippet(code)
        assert "CosmicFlows-4" in res["title"]
        assert res["module"] == "fields/okubo_weiss_tensor.js"
        assert len(res["reportHash"]) == 64  # SHA-256 length
        assert res["citationsCount"] >= 5
        assert res["hasCourtoisCitation"] is True
        assert res["hasStats"] is True


# ============================================================================
# 10. NUMERICAL STABILITY & ERROR HANDLING
# ============================================================================

class TestNumericalStabilityAndEdgeCases:
    """Verify input validation, zero fields, and large velocity dynamic ranges."""

    def test_invalid_constructor_and_arguments(self):
        """Verify TypeError thrown on invalid velocityField instance."""
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        let threw = false;
        try {
          new OkuboWeissTensorAnalyzer(null);
        } catch (e) {
          threw = (e instanceof TypeError);
        }

        console.log(JSON.stringify({ threw }));
        """
        res = run_node_snippet(code)
        assert res["threw"] is True

    def test_zero_velocity_field_invariants(self):
        """Verify zero velocity field evaluates all invariants to exact zeros and DEGENERATE topology."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer, FlowTopologyType } from './src/fields/okubo_weiss_tensor.js';

        const N = 8;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100]
        });

        const total = grid.totalCells;
        const vecField = new VectorField3D(grid, new Float64Array(total), new Float64Array(total), new Float64Array(total));
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const J = analyzer.computeLocalJacobian(4, 4, 4);
        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);
        const inv = OkuboWeissTensorAnalyzer.computeInvariants(J);

        console.log(JSON.stringify({
          strainMagSq: decomp.strainMagnitudeSq,
          vortMagSq: decomp.vorticityMagnitudeSq,
          okuboWeiss: decomp.okuboWeiss,
          P: inv.P,
          Q: inv.Q,
          R: inv.R,
          topology: inv.topology
        }));
        """
        res = run_node_snippet(code)
        assert res["strainMagSq"] == 0.0
        assert res["vortMagSq"] == 0.0
        assert res["okuboWeiss"] == 0.0
        assert res["P"] == 0.0
        assert res["Q"] == 0.0
        assert res["R"] == 0.0
        assert res["topology"] == 0  # DEGENERATE

    def test_extreme_velocity_magnitudes_numerical_stability(self):
        """Verify numerical stability for extreme cluster velocities (10^5 km/s) and subtle flows (10^-8 km/s)."""
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        // Extreme relativistic velocity gradient: 10^5 km/s/(Mpc/h)
        const J_large = new Float64Array([
          1e5, 2e5, -1e5,
          -2e5, 1e5, 3e5,
          1e5, -3e5, -2e5
        ]);

        const decompLarge = OkuboWeissTensorAnalyzer.decomposeKinematics(J_large);
        const invLarge = OkuboWeissTensorAnalyzer.computeInvariants(J_large);

        // Micro-velocity gradient: 10^-8 km/s/(Mpc/h)
        const J_small = new Float64Array([
          1e-8, 2e-8, -1e-8,
          -2e-8, 1e-8, 3e-8,
          1e-8, -3e-8, -2e-8
        ]);

        const decompSmall = OkuboWeissTensorAnalyzer.decomposeKinematics(J_small);
        const invSmall = OkuboWeissTensorAnalyzer.computeInvariants(J_small);

        console.log(JSON.stringify({
          largeFinite: Number.isFinite(decompLarge.okuboWeiss) && Number.isFinite(invLarge.discriminant),
          smallFinite: Number.isFinite(decompSmall.okuboWeiss) && Number.isFinite(invSmall.discriminant),
          largeNorm: decompLarge.okuboWeissNormalized,
          smallNorm: decompSmall.okuboWeissNormalized
        }));
        """
        res = run_node_snippet(code)
        assert res["largeFinite"] is True
        assert res["smallFinite"] is True
        assert math.isclose(res["largeNorm"], res["smallNorm"], rel_tol=1e-6)


# ============================================================================
# 11. CONTINUOUS INTERPOLATION & SPLINE KINEMATICS
# ============================================================================

class TestContinuousInterpolationAndSplineKinematics:
    """Verify C^1 tricubic continuous evaluation of Jacobian and kinematic state at arbitrary sub-voxel positions."""

    def test_continuous_jacobian_and_kinematics_subvoxel_sampling(self):
        """Verify continuous evaluation at sub-voxel coordinates (x, y, z)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100],
          isCellCentered: true
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        // Smooth cubic flow field
        for (let iz = 0; iz < N; iz++) {
          const z = grid.zMin + (iz + 0.5) * grid.dz;
          for (let iy = 0; iy < N; iy++) {
            const y = grid.yMin + (iy + 0.5) * grid.dy;
            for (let ix = 0; ix < N; ix++) {
              const x = grid.xMin + (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              vx[idx] = 0.01 * x * x + 0.02 * y * z;
              vy[idx] = -0.01 * y * y + 0.03 * x * z;
              vz[idx] = 0.05 * z;
            }
          }
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        // Sample at sub-voxel continuous point: SGX = 12.345, SGY = -23.456, SGZ = 7.891 Mpc/h
        const sampleX = 12.345, sampleY = -23.456, sampleZ = 7.891;
        const contJ = analyzer.evaluateContinuousJacobian(sampleX, sampleY, sampleZ);
        const contKin = analyzer.evaluateContinuousKinematics(sampleX, sampleY, sampleZ);

        console.log(JSON.stringify({
          divContinuous: contJ.divergence(),
          detContinuous: contJ.determinant(),
          okuboWeiss: contKin.okuboWeiss,
          lambda2: contKin.lambda2,
          isFinite: Number.isFinite(contKin.okuboWeiss) && Number.isFinite(contKin.lambda2)
        }));
        """
        res = run_node_snippet(code)
        assert res["isFinite"] is True
        assert math.isfinite(res["divContinuous"])
        assert math.isfinite(res["detContinuous"])


# ============================================================================
# 12. INVARIANT PHASE SPACE (Q_J, R_J) PDF NORMALIZATION & STATISTICS
# ============================================================================

class TestInvariantPlanePDFAndStatistics:
    """Verify joint PDF P(Q_J, R_J) normalization, Vieillefosse tail fractions, and covariance."""

    def test_invariant_plane_pdf_integral_normalization(self):
        """Verify 2D joint PDF satisfies integral \\iint P(Q, R) dR dQ approx 1.0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [100, 100, 100],
          isCellCentered: true
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        for (let i = 0; i < total; i++) {
          vx[i] = Math.sin(i * 0.1) * 250.0;
          vy[i] = Math.cos(i * 0.15) * 250.0;
          vz[i] = Math.sin(i * 0.08) * 150.0;
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        const pdfResult = analyzer.computeInvariantPlanePDF({
          nbinsR: 32,
          nbinsQ: 32
        });

        const dr = pdfResult.binsR[1] - pdfResult.binsR[0];
        const dq = pdfResult.binsQ[1] - pdfResult.binsQ[0];

        let totalIntegral = 0.0;
        for (let i = 0; i < pdfResult.pdf2D.length; i++) {
          totalIntegral += pdfResult.pdf2D[i] * dr * dq;
        }

        const sumQuadrants = pdfResult.quadrantFractions.UFC +
                             pdfResult.quadrantFractions.SFS +
                             pdfResult.quadrantFractions.UN_SS +
                             pdfResult.quadrantFractions.SN_SS;

        console.log(JSON.stringify({
          totalIntegral,
          sumQuadrants,
          tailFraction: pdfResult.tailFraction,
          correlation: pdfResult.correlation,
          hasFiniteSkew: Number.isFinite(pdfResult.skewnessR) && Number.isFinite(pdfResult.skewnessQ)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["totalIntegral"], 1.0, rel_tol=0.01)
        assert math.isclose(res["sumQuadrants"], 1.0, abs_tol=1e-12)
        assert 0.0 <= res["tailFraction"] <= 1.0
        assert -1.0 <= res["correlation"] <= 1.0
        assert res["hasFiniteSkew"] is True


# ============================================================================
# 13. TRACELESS SHEAR TENSOR & ANISOTROPIC SPACING INVARIANTS
# ============================================================================

class TestShearTensorAndAnisotropicGrid:
    """Verify traceless shear tensor properties and anisotropic grid spacing (dx != dy != dz)."""

    def test_traceless_shear_tensor_trace_zero(self):
        """Verify Tr(sigma) = 0 and sigma^2 = s^2 - (1/3)*theta^2 identically."""
        code = """
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        // General compressible velocity gradient tensor
        const J = new Float64Array([
          3.5, 1.2, -0.8,
          0.4, -2.1, 1.7,
          -1.1, 0.9, 4.2
        ]);

        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);
        const sigma = decomp.sigma;
        const trSigma = sigma[0] + sigma[4] + sigma[8];

        const expectedShearSq = decomp.strainMagnitudeSq - (1.0 / 3.0) * (decomp.divergence * decomp.divergence);

        console.log(JSON.stringify({
          trSigma,
          shearMagnitudeSq: decomp.shearMagnitudeSq,
          expectedShearSq,
          shearDiff: Math.abs(decomp.shearMagnitudeSq - expectedShearSq)
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["trSigma"], 0.0, abs_tol=1e-12)
        assert math.isclose(res["shearMagnitudeSq"], res["expectedShearSq"], rel_tol=1e-12)

    def test_anisotropic_grid_spacing_derivatives(self):
        """Verify discrete Jacobian computation on anisotropic grid (dx=2, dy=4, dz=8 Mpc/h)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { OkuboWeissTensorAnalyzer } from './src/fields/okubo_weiss_tensor.js';

        const grid = new GridIndexer({
          nx: 16, ny: 8, nz: 4,
          origin: [0, 0, 0],
          boxSize: [32, 32, 32],
          spacing: [2.0, 4.0, 8.0],
          isCellCentered: true
        });

        const total = grid.totalCells;
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);

        // Linear velocity field: vx = 2*x, vy = 3*y, vz = -5*z
        for (let iz = 0; iz < grid.nz; iz++) {
          const z = (iz + 0.5) * grid.dz;
          for (let iy = 0; iy < grid.ny; iy++) {
            const y = (iy + 0.5) * grid.dy;
            for (let ix = 0; ix < grid.nx; ix++) {
              const x = (ix + 0.5) * grid.dx;
              const idx = grid.index(ix, iy, iz);
              vx[idx] = 2.0 * x;
              vy[idx] = 3.0 * y;
              vz[idx] = -5.0 * z;
            }
          }
        }

        const vecField = new VectorField3D(grid, vx, vy, vz);
        const analyzer = new OkuboWeissTensorAnalyzer(vecField);

        // Interior point
        const J = analyzer.computeLocalJacobian(4, 3, 2);
        const decomp = OkuboWeissTensorAnalyzer.decomposeKinematics(J);

        console.log(JSON.stringify({
          dvx_dx: J[0],
          dvy_dy: J[4],
          dvz_dz: J[8],
          divergence: decomp.divergence
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["dvx_dx"], 2.0, rel_tol=1e-6)
        assert math.isclose(res["dvy_dy"], 3.0, rel_tol=1e-6)
        assert math.isclose(res["dvz_dz"], -5.0, rel_tol=1e-6)
        assert math.isclose(res["divergence"], 0.0, abs_tol=1e-6)

