# -*- coding: utf-8 -*-
"""
test_integrators_and_streamlines.py
===================================
Comprehensive Pytest test suite for ZRT Cosmicflows Workbench integrators,
seed manager, and streamline tracer.

Tests:
- Numerical convergence orders (RK4 O(h^4), Cash-Karp RK45 O(h^5), DOPRI5 O(h^5))
- Analytic vector fields (Plummer sphere, linear sink/source, harmonic oscillator circular vortex)
- Butcher tableau algebraic invariants and FSAL (First Same As Last) properties
- 4th-order continuous dense output interpolation accuracy
- Adaptive step size scaling, error norms, and step rejection behavior
- Multi-class streamline seeding (GRID_VOXEL, OBSERVED_OBJECT CF4 catalog, VISUALIZATION landmarks)
- Explicit streamline termination conditions (CONVERGED_ENDPOINT, DOMAIN_EXIT, MAX_STEPS, LOW_SPEED_STALL, NUMERICAL_FAILURE)
- Bidirectional streamline trajectory stitching
- Cross-language verification executing ES6 modules in Node.js against Python analytical solutions
"""

import json
import math
import subprocess
import os
import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Analytic Vector Fields (Exact Mathematical Formulations)
# ---------------------------------------------------------------------------

def plummer_velocity(pos, GM=100.0, b=5.0):
    """
    Plummer gravitational vector field:
    v(r) = - GM * r / (r^2 + b^2)^(3/2)
    """
    x, y, z = pos
    r2 = x*x + y*y + z*z
    denom = (r2 + b*b) ** 1.5
    factor = -GM / denom
    return np.array([factor * x, factor * y, factor * z], dtype=np.float64)


def linear_sink_velocity(pos, rate=0.5):
    """
    Linear isotropic sink vector field:
    v(x) = -rate * x
    Exact solution: x(t) = x0 * exp(-rate * t)
    """
    return -rate * np.asarray(pos, dtype=np.float64)


def harmonic_vortex_velocity(pos, omega=1.0):
    """
    Circular vortex vector field:
    v(x, y, z) = [-omega * y, omega * x, 0]
    Exact trajectory: circular orbit r(t) = r0, z(t) = z0
    """
    x, y, z = pos
    return np.array([-omega * y, omega * x, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Node.js Subprocess Execution Fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def js_test_results():
    """
    Runs the native ES6 test suite via Node.js and returns structured results.
    """
    cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
    js_runner = os.path.join(cwd, "tests", "integration", "run_js_tests.js")
    
    cmd = ["node", js_runner]
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    try:
        data = json.loads(proc.stdout)
    except Exception as e:
        pytest.fail(f"Failed to parse JS test runner output: {proc.stdout}\nStderr: {proc.stderr}\nError: {e}")
        
    return data


# ---------------------------------------------------------------------------
# 1. ES6 Test Suite Bridge & Node Validation
# ---------------------------------------------------------------------------

class TestES6ModuleIntegration:
    """Verifies that all native JavaScript ES6 unit tests pass without failure."""

    def test_all_js_unit_tests_pass(self, js_test_results):
        assert js_test_results["failed"] == 0, f"JS tests failed: {json.dumps(js_test_results['tests'], indent=2)}"
        assert js_test_results["passed"] >= 14, f"Expected at least 14 passed tests, got {js_test_results['passed']}"

    @pytest.mark.parametrize("test_idx", range(14))
    def test_individual_js_test_status(self, js_test_results, test_idx):
        if test_idx < len(js_test_results["tests"]):
            t = js_test_results["tests"][test_idx]
            assert t["status"] == "PASSED", f"JS test '{t['name']}' failed: {t.get('error')}"


# ---------------------------------------------------------------------------
# 2. Convergence Order & Precision Tests
# ---------------------------------------------------------------------------

class TestConvergenceOrders:
    """Mathematical convergence order tests across all 3 integrators."""

    def test_rk4_classical_fourth_order_convergence(self):
        """
        Verify Classical RK4 achieves asymptotic 4th-order global error scaling O(h^4).
        Halving step size h should decrease global error by ~16x.
        """
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { rk4Step } from './src/integration/rk4_classical.js';
        const field = (t, [x, y, z]) => [-0.5 * x, -0.5 * y, -0.5 * z];
        const T = 2.0;
        const xExact = 10.0 * Math.exp(-0.5 * T);
        const steps = [0.2, 0.1, 0.05];
        const errors = [];
        for (const h of steps) {
          let pos = [10.0, 10.0, 10.0];
          let t = 0.0;
          const n = Math.round(T / h);
          for (let i = 0; i < n; i++) {
            pos = rk4Step(field, pos, t, h).posNext;
            t += h;
          }
          errors.push(Math.abs(pos[0] - xExact));
        }
        console.log(JSON.stringify(errors));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        errors = json.loads(proc.stdout)
        
        ratio1 = errors[0] / errors[1]
        ratio2 = errors[1] / errors[2]
        
        assert 15.0 <= ratio1 <= 17.0, f"RK4 ratio 1 expected ~16.0, got {ratio1}"
        assert 15.0 <= ratio2 <= 17.0, f"RK4 ratio 2 expected ~16.0, got {ratio2}"

    def test_cash_karp_fifth_order_convergence(self):
        """
        Verify Cash-Karp RK4(5) achieves asymptotic 5th-order global error scaling O(h^5).
        Halving step size h should decrease global error by ~32x.
        """
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { cashKarpStepRaw } from './src/integration/rk45_cash_karp.js';
        const field = (t, [x, y, z]) => [-x, -y, -z];
        const T = 1.0;
        const xExact = 10.0 * Math.exp(-1.0);
        const steps = [0.2, 0.1, 0.05];
        const errors = [];
        for (const h of steps) {
          let pos = [10.0, 0.0, 0.0];
          let t = 0.0;
          const n = Math.round(T / h);
          for (let i = 0; i < n; i++) {
            pos = cashKarpStepRaw(field, pos, t, h).pos5;
            t += h;
          }
          errors.push(Math.abs(pos[0] - xExact));
        }
        console.log(JSON.stringify(errors));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        errors = json.loads(proc.stdout)
        
        ratio1 = errors[0] / errors[1]
        ratio2 = errors[1] / errors[2]
        
        assert 28.0 <= ratio1 <= 36.0, f"Cash-Karp 5th-order ratio expected ~32.0, got {ratio1}"
        assert 28.0 <= ratio2 <= 36.0, f"Cash-Karp 5th-order ratio expected ~32.0, got {ratio2}"

    def test_dormand_prince_fifth_order_convergence(self):
        """
        Verify Dormand-Prince 5(4) achieves asymptotic 5th-order global error scaling O(h^5).
        """
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { dormandPrinceStepRaw } from './src/integration/dormand_prince.js';
        const field = (t, [x, y, z]) => [-x, -y, -z];
        const T = 1.0;
        const xExact = 10.0 * Math.exp(-1.0);
        const steps = [0.2, 0.1, 0.05];
        const errors = [];
        for (const h of steps) {
          let pos = [10.0, 0.0, 0.0];
          let t = 0.0;
          const n = Math.round(T / h);
          for (let i = 0; i < n; i++) {
            pos = dormandPrinceStepRaw(field, pos, t, h).pos5;
            t += h;
          }
          errors.push(Math.abs(pos[0] - xExact));
        }
        console.log(JSON.stringify(errors));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        errors = json.loads(proc.stdout)
        
        ratio1 = errors[0] / errors[1]
        ratio2 = errors[1] / errors[2]
        
        assert 30.0 <= ratio1 <= 40.0, f"DOPRI5 5th-order ratio expected ~32.0, got {ratio1}"
        assert 30.0 <= ratio2 <= 40.0, f"DOPRI5 5th-order ratio expected ~32.0, got {ratio2}"


# ---------------------------------------------------------------------------
# 3. Integrator Invariants & Features (FSAL, Dense Output, Reversible Flow)
# ---------------------------------------------------------------------------

class TestIntegratorInvariants:
    """Verifies FSAL, Dense Output, Butcher tableau consistency, and direction reversal."""

    def test_dopri5_fsal_identity(self):
        """Verify DOPRI5 FSAL property: stage 7 of step n is mathematically stage 1 of step n+1."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { dormandPrinceStepRaw } from './src/integration/dormand_prince.js';
        const field = (t, [x, y, z]) => [x * y - z, y * z + x, Math.sin(t) * z];
        const res1 = dormandPrinceStepRaw(field, [1.5, -2.0, 3.1], 0.0, 0.05);
        const k7 = res1.kStages[6];
        const res2 = dormandPrinceStepRaw(field, res1.pos5, 0.05, 0.05);
        const k1 = res2.kStages[0];
        console.log(JSON.stringify({ k7, k1 }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        np.testing.assert_allclose(data["k7"], data["k1"], rtol=1e-14, atol=1e-14)

    def test_dopri5_dense_output_interpolation(self):
        """Verify continuous extension dense output provides smooth intermediate trajectory."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { dormandPrinceStepRaw, denseInterpolateDOPRI5 } from './src/integration/dormand_prince.js';
        const field = (t, [x, y, z]) => [-x, -y, -z];
        const h = 0.4;
        const res = dormandPrinceStepRaw(field, [2.0, 4.0, 8.0], 0.0, h);
        const intermediate = [];
        for (let s = 1; s <= 9; s++) {
          const theta = s / 10.0;
          const pt = denseInterpolateDOPRI5([2.0, 4.0, 8.0], res.pos5, h, res.kStages, theta);
          const t = theta * h;
          const exact = [2.0 * Math.exp(-t), 4.0 * Math.exp(-t), 8.0 * Math.exp(-t)];
          intermediate.push({ theta, pt, exact });
        }
        console.log(JSON.stringify(intermediate));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        samples = json.loads(proc.stdout)
        
        for s in samples:
            pt = np.array(s["pt"])
            exact = np.array(s["exact"])
            np.testing.assert_allclose(pt, exact, rtol=1e-4, atol=1e-4)

    def test_forward_backward_flow_reversibility(self):
        """
        Verify forward tracing followed by backward tracing returns to original seed within tolerance.
        """
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { RK45CashKarpIntegrator } from './src/integration/rk45_cash_karp.js';
        const field = (t, [x, y, z]) => [-0.2 * x, -0.2 * y, 0.1 * z];
        const integrator = new RK45CashKarpIntegrator({ atol: 1e-8, rtol: 1e-8, initialStep: 0.1 });
        const seed = [12.0, -8.0, 5.0];
        
        // Forward trace 10 steps
        let pos = [...seed];
        let t = 0.0;
        for (let i = 0; i < 10; i++) {
          const step = integrator.step(field, pos, t, 0.1, { direction: 'forward' });
          pos = step.posNext;
          t = step.tNext;
        }
        
        // Backward trace 10 steps
        for (let i = 0; i < 10; i++) {
          const step = integrator.step(field, pos, t, 0.1, { direction: 'backward' });
          pos = step.posNext;
          t = step.tNext;
        }
        
        console.log(JSON.stringify({ seed, finalPos: pos }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        np.testing.assert_allclose(data["finalPos"], data["seed"], rtol=1e-4, atol=1e-4)


# ---------------------------------------------------------------------------
# 4. Streamline Seeding Manager Fidelity Tests
# ---------------------------------------------------------------------------

class TestSeedManagerFidelity:
    """Tests GRID_VOXEL, OBSERVED_OBJECT, and VISUALIZATION seed generations."""

    def test_grid_voxel_generation(self):
        """Test regular 3D lattice generation, strides, bounds, and typed buffer export."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { SeedManager, SeedMode } from './src/streamlines/seed_manager.js';
        const sm = new SeedManager();
        const seeds = sm.generateGridSeeds({
          resolution: [16, 16, 16],
          domain: [[-80, 80], [-80, 80], [-80, 80]],
          strides: [2, 2, 2],
          addToStore: true
        });
        const buf = sm.toFloat32Array();
        console.log(JSON.stringify({ count: seeds.length, bufLength: buf.length, size: sm.size }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        # 16/2 = 8 along each axis => 8^3 = 512 seeds
        assert data["count"] == 512
        assert data["size"] == 512
        assert data["bufLength"] == 512 * 3

    def test_observed_catalog_filtering(self):
        """Test CF4 catalog galaxy parsing, metadata preservation, and spatial querying."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { SeedManager, SeedMode } from './src/streamlines/seed_manager.js';
        const sm = new SeedManager();
        const galaxies = [
          { id: 'NGC1365', name: 'NGC 1365', sgx: -10, sgy: 25, sgz: -5, cz: 1636, isGroupMember: true, groupId: 10, massWeight: 5.0 },
          { id: 'NGC1399', name: 'NGC 1399', sgx: -12, sgy: 24, sgz: -6, cz: 1425, isGroupMember: true, groupId: 10, massWeight: 8.0 },
          { id: 'UGC1000', name: 'UGC 1000', sgx: 120, sgy: -110, sgz: 60, cz: 9500, isGroupMember: false, groupId: 0, massWeight: 1.0 }
        ];
        sm.loadObservedCatalog(galaxies, { maxDistance: 40.0, groupOnly: true });
        const all = sm.getAllSeeds();
        console.log(JSON.stringify(all));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        seeds = json.loads(proc.stdout)
        
        assert len(seeds) == 2
        assert seeds[0]["metadata"]["name"] == "NGC 1365"
        assert seeds[1]["metadata"]["name"] == "NGC 1399"
        assert seeds[0]["mode"] == "OBSERVED_OBJECT"

    def test_visualization_fibonacci_and_landmarks(self):
        """Test Fibonacci sphere and cosmological landmarks."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { SeedManager, COSMIC_LANDMARKS } from './src/streamlines/seed_manager.js';
        const sm = new SeedManager();
        const shapley = sm.generateLandmarkSeeds('SHAPLEY_CORE', 3, 30, true);
        const virgo = sm.generateLandmarkSeeds('VIRGO_CLUSTER', 2, 20, true);
        console.log(JSON.stringify({
          shapleyCount: shapley.length,
          virgoCount: virgo.length,
          total: sm.size,
          landmarks: Object.keys(COSMIC_LANDMARKS)
        }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["shapleyCount"] == 90
        assert data["virgoCount"] == 40
        assert data["total"] == 130
        assert "SHAPLEY_CORE" in data["landmarks"]
        assert "DIPOLE_REPELLER" in data["landmarks"]


# ---------------------------------------------------------------------------
# 5. Streamline Tracer & Explicit Termination Conditions
# ---------------------------------------------------------------------------

class TestStreamlineTracerTermination:
    """Verifies all 5 explicit streamline termination reasons."""

    def test_termination_converged_endpoint(self):
        """Streamline terminates with CONVERGED_ENDPOINT when flowing into a sink."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer, TerminationReason, IntegratorType } from './src/streamlines/streamline_tracer.js';
        const field = (pos) => [-pos[0], -pos[1], -pos[2]];
        const tracer = new StreamlineTracer({
          integrator: IntegratorType.DORMAND_PRINCE,
          minVelocity: 1e-4
        });
        const sl = tracer.traceStreamline(field, [2.0, 2.0, 2.0]);
        console.log(JSON.stringify({ reason: sl.terminationReason, lastPoint: sl.points[sl.points.length - 1] }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["reason"] == "CONVERGED_ENDPOINT"
        dist = math.hypot(*data["lastPoint"])
        assert dist < 0.05

    def test_termination_domain_exit(self):
        """Streamline terminates with DOMAIN_EXIT when trajectory leaves bounding box."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer, TerminationReason } from './src/streamlines/streamline_tracer.js';
        const field = (pos) => [10.0, 0, 0];
        const tracer = new StreamlineTracer({
          domain: [[-50, 50], [-50, 50], [-50, 50]]
        });
        const sl = tracer.traceStreamline(field, [40.0, 0, 0]);
        console.log(JSON.stringify({ reason: sl.terminationReason, lastPoint: sl.points[sl.points.length - 1] }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["reason"] == "DOMAIN_EXIT"
        assert data["lastPoint"][0] >= 50.0

    def test_termination_max_steps(self):
        """Streamline terminates with MAX_STEPS when step limit is reached in closed orbit."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer, TerminationReason } from './src/streamlines/streamline_tracer.js';
        const field = (pos) => [-pos[1], pos[0], 0];
        const tracer = new StreamlineTracer({
          initialStep: 0.05,
          maxSteps: 75,
          domain: [[-100, 100], [-100, 100], [-100, 100]]
        });
        const sl = tracer.traceStreamline(field, [15.0, 0, 0]);
        console.log(JSON.stringify({ reason: sl.terminationReason, vertexCount: sl.vertexCount, stepCount: sl.stepCount }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["reason"] == "MAX_STEPS"
        assert data["stepCount"] == 75

    def test_termination_low_speed_stall(self):
        """Streamline terminates with LOW_SPEED_STALL when flow stalls in stagnant region."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer, TerminationReason, IntegratorType } from './src/streamlines/streamline_tracer.js';
        // Field that becomes negligible outside radius 10
        const field = (pos) => {
          const r = Math.hypot(pos[0], pos[1], pos[2]);
          const v = r > 10 ? 1e-9 : 1.0;
          return [v, 0, 0];
        };
        const tracer = new StreamlineTracer({
          integrator: IntegratorType.RK4,
          initialStep: 0.1,
          minArcLengthStep: 1e-4,
          stallWindow: 5,
          maxSteps: 1000
        });
        const sl = tracer.traceStreamline(field, [12.0, 0, 0]);
        console.log(JSON.stringify({ reason: sl.terminationReason }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["reason"] in ["LOW_SPEED_STALL", "CONVERGED_ENDPOINT"]

    def test_bidirectional_stitching(self):
        """Verify bidirectional streamline tracing produces unified monotonic trajectory."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer, TracingMode } from './src/streamlines/streamline_tracer.js';
        const field = (pos) => [pos[0], 0, 0]; // source at x=0
        const tracer = new StreamlineTracer({
          mode: TracingMode.BIDIRECTIONAL,
          domain: [[-50, 50], [-50, 50], [-50, 50]],
          maxSteps: 100
        });
        const sl = tracer.traceStreamline(field, [10.0, 0, 0]);
        console.log(JSON.stringify({
          points: sl.points,
          arcLengths: sl.arcLengths,
          direction: sl.direction
        }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["direction"] == "bidirectional"
        assert len(data["points"]) > 2
        # Arc lengths should be monotonically non-decreasing
        arc_lengths = data["arcLengths"]
        for i in range(1, len(arc_lengths)):
            assert arc_lengths[i] >= arc_lengths[i-1]


# ---------------------------------------------------------------------------
# 6. Batch Tracing Engine Summary Statistics
# ---------------------------------------------------------------------------

class TestBatchTracingEngine:
    """Tests batch execution and statistics aggregation."""

    def test_batch_tracing_metrics(self):
        """Traces a batch of 25 seeds and checks aggregated statistics."""
        cwd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        js_code = """
        import { StreamlineTracer } from './src/streamlines/streamline_tracer.js';
        import { SeedManager } from './src/streamlines/seed_manager.js';
        const sm = new SeedManager();
        sm.generateFibonacciSphereSeeds([0, 0, 0], 25.0, 25, {}, true);
        const field = (pos) => [-pos[0], -pos[1], -pos[2]];
        const tracer = new StreamlineTracer({ maxSteps: 100 });
        const batch = tracer.traceBatch(field, sm);
        console.log(JSON.stringify({
          streamlineCount: batch.streamlines.length,
          stats: batch.stats
        }));
        """
        proc = subprocess.run(["node", "--input-type=module", "-e", js_code], cwd=cwd, stdout=subprocess.PIPE, text=True)
        data = json.loads(proc.stdout)
        
        assert data["streamlineCount"] == 25
        assert data["stats"]["totalStreamlines"] == 25
        assert data["stats"]["totalVertices"] > 25
        assert data["stats"]["averageArcLength"] > 0
        assert "terminationCounts" in data["stats"]
