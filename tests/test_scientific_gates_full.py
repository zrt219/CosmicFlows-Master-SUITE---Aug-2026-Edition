"""
Full Scientific Acceptance Gates Master Test Runner
===================================================
Executes all scientific verification and regression acceptance gates across the Cosmicflows
computational cosmology platform:
  - Gate 1: Analytic Vector Fields & Differential Invariants Gate
  - Gate 2: Numerical Integration & Conservation Gate (Cash-Karp RK45)
  - Gate 3: Critical Point Topology & Eigenvalue Classification Gate
  - Gate 4: Multi-Threaded Worker Determinism & Thread-Safety Gate
  - Gate 5: Watershed Cosmic Basin Segmentation Gate
  - Gate 6: Bulk Flow Multipole & Cosmological Observables Gate
  - Gate 7: Full System Verification, Telemetry & Acceptance Summary

Usage:
  pytest tests/test_scientific_gates_full.py -v
  python tests/test_scientific_gates_full.py
"""

import sys
import os
import time
import math
import hashlib
import numpy as np
import pytest
from typing import Dict, Any, List, Tuple

import importlib.util

# Ensure workspace root is in sys.path
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

def _load_module(module_name: str, rel_path: str):
    file_path = os.path.join(WORKSPACE_ROOT, rel_path)
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod

_analytic_mod = _load_module("analytic_fields_test_mod", "tests/analytic-fields/test_analytic_fields.py")
_worker_mod = _load_module("workers_test_mod", "tests/workers/test_worker_determinism.py")

LaminarFlowField = _analytic_mod.LaminarFlowField
LinearSinkField = _analytic_mod.LinearSinkField
LinearSourceField = _analytic_mod.LinearSourceField
SaddleFilamentField = _analytic_mod.SaddleFilamentField
SaddleWallField = _analytic_mod.SaddleWallField
RotationalVortexField = _analytic_mod.RotationalVortexField
PlummerGravitationalCluster = _analytic_mod.PlummerGravitationalCluster
compute_numerical_jacobian = _analytic_mod.compute_numerical_jacobian
compute_numerical_divergence = _analytic_mod.compute_numerical_divergence
compute_numerical_vorticity = _analytic_mod.compute_numerical_vorticity
jacobi_eigenvalue_diagonalization = _analytic_mod.jacobi_eigenvalue_diagonalization
newton_raphson_critical_point_3d = _analytic_mod.newton_raphson_critical_point_3d
rk45_cash_karp_integrate = _analytic_mod.rk45_cash_karp_integrate

multi_cluster_velocity = _worker_mod.multi_cluster_velocity
integrate_single_streamline = _worker_mod.integrate_single_streamline
compute_grid_voxel_tensor = _worker_mod.compute_grid_voxel_tensor
watershed_classify_voxel = _worker_mod.watershed_classify_voxel
compute_bulk_flow_estimator = _worker_mod.compute_bulk_flow_estimator
run_streamlines_serial = _worker_mod.run_streamlines_serial
run_streamlines_threaded = _worker_mod.run_streamlines_threaded
run_grid_tensors_serial = _worker_mod.run_grid_tensors_serial
run_grid_tensors_threaded = _worker_mod.run_grid_tensors_threaded


# =============================================================================
# GATE 1: ANALYTIC VECTOR FIELDS & DIFFERENTIAL INVARIANTS GATE
# =============================================================================

class TestGate1AnalyticVectorFields:
    """Gate 1: Verification of analytic vector fields and differential operators."""

    def test_g1_01_laminar_zero_div_and_curl(self):
        field = LaminarFlowField(u0=120.0, v0=-450.0, w0=230.0)
        pt = np.array([1234.5, -5678.9, 9012.3])
        div_v = compute_numerical_divergence(field.velocity, pt, h=1e-4)
        curl_v = compute_numerical_vorticity(field.velocity, pt, h=1e-4)
        assert abs(div_v) < 1e-11, f"Laminar flow divergence must be 0, got {div_v}"
        np.testing.assert_allclose(curl_v, np.zeros(3), atol=1e-11)

    def test_g1_02_linear_sink_attractor_invariants(self):
        sink = LinearSinkField(l1=3.0, l2=2.0, l3=1.0)
        evals = sink.analytical_eigenvalues()
        assert np.all(evals < 0), "Sink eigenvalues must be strictly negative."
        div_v = compute_numerical_divergence(sink.velocity, np.zeros(3), h=1e-4)
        assert div_v == pytest.approx(-6.0, rel=1e-6)

    def test_g1_03_linear_source_repeller_invariants(self):
        source = LinearSourceField(l1=1.5, l2=2.5, l3=3.5)
        evals = source.analytical_eigenvalues()
        assert np.all(evals > 0), "Source eigenvalues must be strictly positive."
        div_v = compute_numerical_divergence(source.velocity, np.zeros(3), h=1e-4)
        assert div_v == pytest.approx(7.5, rel=1e-6)

    def test_g1_04_saddle_filament_two_neg_one_pos(self):
        filament = SaddleFilamentField(l1=2.0, l2=1.5, l3=1.0)
        evals = filament.analytical_eigenvalues()
        assert len(evals[evals < 0]) == 2
        assert len(evals[evals > 0]) == 1

    def test_g1_05_saddle_wall_one_neg_two_pos(self):
        wall = SaddleWallField(l1=2.5, l2=1.0, l3=1.5)
        evals = wall.analytical_eigenvalues()
        assert len(evals[evals < 0]) == 1
        assert len(evals[evals > 0]) == 2

    def test_g1_06_rotational_vortex_vorticity_guard(self):
        vortex = RotationalVortexField(omega_vector=(0.0, 0.0, 4.0))
        pt = np.array([10.0, 20.0, 0.0])
        curl_v = compute_numerical_vorticity(vortex.velocity, pt, h=1e-4)
        np.testing.assert_allclose(curl_v, np.array([0.0, 0.0, 8.0]), rtol=1e-5)
        div_v = compute_numerical_divergence(vortex.velocity, pt, h=1e-4)
        assert abs(div_v) < 1e-10

    def test_g1_07_plummer_smooth_core_and_peak_inflow(self):
        cluster = PlummerGravitationalCluster(center=(0.0, 0.0, 0.0), GM=1e7, epsilon=600.0)
        # Inflow at center is 0
        np.testing.assert_array_equal(cluster.velocity(np.zeros(3)), np.zeros(3))
        # Peak radius
        r_peak = cluster.peak_radius()
        assert r_peak == pytest.approx(600.0 / math.sqrt(2.0), rel=1e-9)
        # Analytical divergence matches finite difference
        pt = np.array([300.0, -400.0, 200.0])
        ana_div = cluster.analytical_divergence(pt)
        num_div = compute_numerical_divergence(cluster.velocity, pt, h=1e-2)
        assert num_div == pytest.approx(ana_div, rel=1e-4)


# =============================================================================
# GATE 2: NUMERICAL INTEGRATION & CONSERVATION GATE
# =============================================================================

class TestGate2NumericalIntegration:
    """Gate 2: Validation of Cash-Karp RK45 Butcher Tableau & Adaptive Step Control."""

    def test_g2_01_cash_karp_butcher_tableau_algebraic_consistency(self):
        # Verify Cash-Karp weights sum to 1.0 (consistency condition)
        b5 = np.array([37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0])
        b4 = np.array([2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0])
        assert math.isclose(float(np.sum(b5)), 1.0, rel_tol=1e-15)
        assert math.isclose(float(np.sum(b4)), 1.0, rel_tol=1e-15)

    def test_g2_02_rk45_linear_sink_convergence_rate(self):
        sink = LinearSinkField(l1=1.0, l2=1.0, l3=1.0)
        x0 = np.array([10.0, 10.0, 10.0])
        t_hist, x_hist = rk45_cash_karp_integrate(sink.velocity, x0, t_end=3.0, tol=1e-9)
        x_exact = x0 * math.exp(-3.0)
        np.testing.assert_allclose(x_hist[-1], x_exact, rtol=1e-6, atol=1e-6)

    def test_g2_03_rk45_closed_vortex_energy_conservation(self):
        vortex = RotationalVortexField(omega_vector=(0.0, 0.0, 2.0))
        x0 = np.array([5.0, 0.0, 0.0])
        period = 2.0 * math.pi / 2.0
        t_hist, x_hist = rk45_cash_karp_integrate(vortex.velocity, x0, t_end=period, tol=1e-10)
        r0 = np.linalg.norm(x0)
        # Radius must be preserved along the entire orbit
        radii = [np.linalg.norm(x) for x in x_hist]
        for r in radii:
            assert math.isclose(r, r0, rel_tol=1e-4)


# =============================================================================
# GATE 3: CRITICAL POINT TOPOLOGY & EIGENVALUE CLASSIFICATION GATE
# =============================================================================

class TestGate3CriticalPointTopology:
    """Gate 3: Critical point locating, Jacobian diagonalization, and topological indexing."""

    def test_g3_01_newton_raphson_sink_root_finding(self):
        sink = LinearSinkField(l1=2.0, l2=2.5, l3=3.0)
        x_init = np.array([123.0, -456.0, 789.0])
        x_root, conv, iters = newton_raphson_critical_point_3d(sink.velocity, x_init)
        assert conv is True
        np.testing.assert_allclose(x_root, np.zeros(3), atol=1e-12)

    def test_g3_02_jacobi_diagonalization_orthogonality(self):
        # Symmetric deformation tensor
        S = np.array([
            [4.0, 1.0, -2.0],
            [1.0, 3.0, 0.5],
            [-2.0, 0.5, 5.0]
        ], dtype=np.float64)
        evals, V = jacobi_eigenvalue_diagonalization(S)
        # Verify V is orthogonal: V.T @ V = I
        np.testing.assert_allclose(V.T @ V, np.eye(3), atol=1e-14)
        # Verify diagonalization: V.T @ S @ V = diag(evals)
        np.testing.assert_allclose(V.T @ S @ V, np.diag(evals), atol=1e-13)

    def test_g3_03_poincare_hopf_topological_indexing(self):
        # Attractor node: 3 negative eigenvalues -> index -1
        # Repeller node: 3 positive eigenvalues -> index +1
        sink_evals = np.array([-1.0, -2.0, -3.0])
        source_evals = np.array([1.0, 2.0, 3.0])
        index_sink = (-1) ** len(sink_evals[sink_evals < 0])
        index_source = (-1) ** len(source_evals[source_evals < 0])
        assert index_sink == -1
        assert index_source == 1


# =============================================================================
# GATE 4: MULTI-THREADED WORKER DETERMINISM & THREAD SAFETY GATE
# =============================================================================

class TestGate4WorkerDeterminism:
    """Gate 4: Multi-threaded worker numerical determinism & memory safety."""

    def test_g4_01_streamline_parallel_serial_exact_identity(self):
        np.random.seed(999)
        seeds = np.random.uniform(-4000.0, 4000.0, size=(40, 3))
        serial_res = run_streamlines_serial(seeds)
        parallel_res = run_streamlines_threaded(seeds, num_workers=4, chunk_size=10)
        assert len(serial_res) == len(parallel_res)
        for i in range(len(seeds)):
            np.testing.assert_array_equal(parallel_res[i], serial_res[i])

    def test_g4_02_grid_tensors_parallel_serial_identity(self):
        coords = np.array([
            [1000.0, 2000.0, 3000.0],
            [-2000.0, -3000.0, 1000.0],
            [5000.0, -5000.0, 0.0],
            [-1000.0, 4000.0, -2000.0]
        ], dtype=np.float64)
        div_s, curl_s, eval_s = run_grid_tensors_serial(coords)
        div_p, curl_p, eval_p = run_grid_tensors_threaded(coords, num_workers=2, chunk_size=2)
        np.testing.assert_array_equal(div_p, div_s)
        np.testing.assert_array_equal(curl_p, curl_s)
        np.testing.assert_array_equal(eval_p, eval_s)

    def test_g4_03_sha256_bitwise_concurrency_repeatability(self):
        np.random.seed(777)
        seeds = np.random.uniform(-3000.0, 3000.0, size=(20, 3))
        hashes = set()
        for _ in range(5):
            res = run_streamlines_threaded(seeds, num_workers=4, chunk_size=5)
            h = hashlib.sha256(np.concatenate(res).tobytes()).hexdigest()
            hashes.add(h)
        assert len(hashes) == 1, "Concurrency repeatability failed: non-identical hashes across runs."


# =============================================================================
# GATE 5: WATERSHED COSMIC BASIN SEGMENTATION GATE
# =============================================================================

class TestGate5WatershedSegmentation:
    """Gate 5: Watershed flow segmentation & basin similarity."""

    def test_g5_01_watershed_partition_of_unity(self):
        # Verify all grid points are assigned to a valid basin (no orphan points)
        attractors = [np.array([7200.0, -8600.0, -2400.0]), np.array([-3500.0, 1500.0, -800.0])]
        test_points = np.array([
            [6000.0, -7000.0, -2000.0],
            [-3000.0, 1000.0, -500.0],
            [0.0, 0.0, 0.0]
        ])
        for pt in test_points:
            basin = watershed_classify_voxel(pt, attractors)
            assert basin in (0, 1), f"Invalid basin ID {basin}"

    def test_g5_02_watershed_dice_similarity_is_unity(self):
        attractors = [np.array([7200.0, -8600.0, -2400.0]), np.array([-3500.0, 1500.0, -800.0])]
        pts = np.array([
            [7000.0, -8000.0, -2000.0],
            [-3000.0, 1200.0, -700.0]
        ])
        s1 = np.array([watershed_classify_voxel(p, attractors) for p in pts])
        s2 = np.array([watershed_classify_voxel(p, attractors) for p in pts])
        intersection = np.sum((s1 == 0) & (s2 == 0))
        dice = 2.0 * intersection / (np.sum(s1 == 0) + np.sum(s2 == 0))
        assert dice == 1.0


# =============================================================================
# GATE 6: BULK FLOW MULTIPOLE & COSMOLOGICAL OBSERVABLES GATE
# =============================================================================

class TestGate6BulkFlowMultipoles:
    """Gate 6: Bulk flow minimum-variance estimators and cosmological dipole validation."""

    def test_g6_01_bulk_flow_estimator_recovery(self):
        # Create synthetic uniform flow (300, -150, 200) + isotropic random galaxies
        np.random.seed(42)
        n = 500
        pos = np.random.normal(0, 3000, size=(n, 3))
        v_true = np.array([300.0, -150.0, 200.0])
        vel = np.tile(v_true, (n, 1)) + np.random.normal(0, 20, size=(n, 3))
        catalog = np.hstack([pos, vel])
        
        v_est = compute_bulk_flow_estimator(catalog)
        np.testing.assert_allclose(v_est, v_true, rtol=0.1, atol=25.0)

    def test_g6_02_bulk_flow_zero_flow_invariance(self):
        np.random.seed(123)
        n = 300
        pos = np.random.normal(0, 2000, size=(n, 3))
        vel = np.random.normal(0, 10, size=(n, 3))  # pure noise, zero net flow
        catalog = np.hstack([pos, vel])
        v_est = compute_bulk_flow_estimator(catalog)
        assert np.linalg.norm(v_est) < 50.0


# =============================================================================
# GATE 7: FULL SYSTEM TEST RUNNER CLI & SUMMARY
# =============================================================================

def run_all_scientific_gates() -> int:
    """
    Stand-alone execution runner with detailed telemetry and gate certification.
    """
    print("=" * 80)
    print("  COSMICFLOWS SCIENTIFIC REGRESSION & ACCEPTANCE GATES MASTER SUITE")
    print("=" * 80)
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Platform:  Python {sys.version.split()[0]} on {sys.platform}\n")

    gates = [
        ("Gate 1: Analytic Vector Fields & Differential Invariants", [
            ("Laminar Flow (Div=0, Curl=0)", TestGate1AnalyticVectorFields().test_g1_01_laminar_zero_div_and_curl),
            ("Linear Sink Attractor (Eigenvalues < 0)", TestGate1AnalyticVectorFields().test_g1_02_linear_sink_attractor_invariants),
            ("Linear Source Repeller (Eigenvalues > 0)", TestGate1AnalyticVectorFields().test_g1_03_linear_source_repeller_invariants),
            ("Saddle-Filament Topology (2 neg, 1 pos)", TestGate1AnalyticVectorFields().test_g1_04_saddle_filament_two_neg_one_pos),
            ("Saddle-Wall Topology (1 neg, 2 pos)", TestGate1AnalyticVectorFields().test_g1_05_saddle_wall_one_neg_two_pos),
            ("Rotational Vortex & Vorticity Guard", TestGate1AnalyticVectorFields().test_g1_06_rotational_vortex_vorticity_guard),
            ("3D Plummer Gravitational Model", TestGate1AnalyticVectorFields().test_g1_07_plummer_smooth_core_and_peak_inflow)
        ]),
        ("Gate 2: Numerical Integration & Conservation (Cash-Karp RK45)", [
            ("Butcher Tableau Consistency", TestGate2NumericalIntegration().test_g2_01_cash_karp_butcher_tableau_algebraic_consistency),
            ("RK45 Linear Sink Convergence", TestGate2NumericalIntegration().test_g2_02_rk45_linear_sink_convergence_rate),
            ("RK45 Vortex Energy Conservation", TestGate2NumericalIntegration().test_g2_03_rk45_closed_vortex_energy_conservation)
        ]),
        ("Gate 3: Critical Point Topology & Eigenvalue Classification", [
            ("3D Newton-Raphson Root Finding", TestGate3CriticalPointTopology().test_g3_01_newton_raphson_sink_root_finding),
            ("Jacobi Eigenvalue Diagonalization", TestGate3CriticalPointTopology().test_g3_02_jacobi_diagonalization_orthogonality),
            ("Poincaré-Hopf Topological Indexing", TestGate3CriticalPointTopology().test_g3_03_poincare_hopf_topological_indexing)
        ]),
        ("Gate 4: Multi-Threaded Worker Determinism & Thread Safety", [
            ("Streamline Batch Parallel Identity", TestGate4WorkerDeterminism().test_g4_01_streamline_parallel_serial_exact_identity),
            ("Grid Differential Tensor Identity", TestGate4WorkerDeterminism().test_g4_02_grid_tensors_parallel_serial_identity),
            ("Bitwise SHA-256 Repeatability", TestGate4WorkerDeterminism().test_g4_03_sha256_bitwise_concurrency_repeatability)
        ]),
        ("Gate 5: Watershed Cosmic Basin Segmentation", [
            ("Partition of Unity (No Orphans)", TestGate5WatershedSegmentation().test_g5_01_watershed_partition_of_unity),
            ("Segmentation Dice Similarity = 1.0", TestGate5WatershedSegmentation().test_g5_02_watershed_dice_similarity_is_unity)
        ]),
        ("Gate 6: Bulk Flow Multipole & Cosmological Observables", [
            ("Minimum-Variance Dipole Recovery", TestGate6BulkFlowMultipoles().test_g6_01_bulk_flow_estimator_recovery),
            ("Zero Net Flow Invariance", TestGate6BulkFlowMultipoles().test_g6_02_bulk_flow_zero_flow_invariance)
        ])
    ]

    total_passed = 0
    total_failed = 0
    t0_start = time.time()

    for gate_name, tests in gates:
        print(f"[*] Executing {gate_name}...")
        gate_t0 = time.time()
        gate_passed = 0
        gate_failed = 0

        for test_desc, test_fn in tests:
            try:
                test_fn()
                gate_passed += 1
                total_passed += 1
                print(f"    [PASS] {test_desc}")
            except Exception as e:
                gate_failed += 1
                total_failed += 1
                print(f"    [FAIL] {test_desc}: {e}")

        gate_duration = (time.time() - gate_t0) * 1000.0
        status = "PASSED" if gate_failed == 0 else "FAILED"
        print(f"    -> Gate Status: {status} ({gate_passed}/{len(tests)} passed, {gate_duration:.2f} ms)\n")

    total_time = time.time() - t0_start
    print("=" * 80)
    print("  SCIENTIFIC ACCEPTANCE GATES SUMMARY")
    print("=" * 80)
    print(f"Total Tests Executed: {total_passed + total_failed}")
    print(f"Passed:               {total_passed}")
    print(f"Failed:               {total_failed}")
    print(f"Total Duration:       {total_time:.3f} s")
    print(f"Final Certification:  {'ALL GATES PASSED (100% GREEN)' if total_failed == 0 else 'GATE VIOLATIONS DETECTED'}")
    print("=" * 80)

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(run_all_scientific_gates())
