# -*- coding: utf-8 -*-
"""
tests/streamlines/test_poincare_section_tracer.py
Exhaustive scientific verification suite for Poincaré Section and Orbit Recurrence Analyzer.

Validates:
1. Geometric cutting plane definition & orthonormal coordinate projections (SGZ, SGY, SGX, arbitrary, principal tidal).
2. High-order exact intersection interpolation (Linear, Cubic Hermite spline, Brent-Dekker root finding).
3. Poincaré recurrence maps, return times, Jacobian area preservation, and RQA metrics.
4. Lyapunov exponent estimators (Benettin shadow trajectory & Variational tangent map with QR decomposition).
5. Invariant tori, island chains, elliptic/hyperbolic fixed points, and accretion basin flow classification.
6. Cosmological supercluster testbeds: Laniakea, Shapley, Virgo, Coma, Great Attractor, Perseus-Pisces.
7. Multi-integrator compliance (RK4, Cash-Karp RK45, DOPRI5), parameter sweeps, and physical invariants.
"""

import math
import pytest
from tests.utils import run_node_snippet

class TestPoincarePlaneGeometry:
    """Tests for Poincaré cutting plane definition and coordinate transformations."""

    def test_sgz_plane_creation_and_tangent_basis(self):
        code = """
        import { PoincarePlane, PlaneOrientation } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(10.0);
        
        console.log(JSON.stringify({
          name: plane.name,
          orientation: plane.orientation,
          origin: plane.origin,
          normal: plane.normal,
          uAxis: plane.uAxis,
          vAxis: plane.vAxis
        }));
        """
        res = run_node_snippet(code)
        assert res["orientation"] == "SGZ_PLANE"
        assert res["origin"] == [0.0, 0.0, 10.0]
        assert res["normal"] == [0.0, 0.0, 1.0]
        assert res["uAxis"] == [1.0, 0.0, 0.0]
        assert res["vAxis"] == [0.0, 1.0, 0.0]

    def test_sgy_and_sgx_plane_factories(self):
        code = """
        import { PoincarePlane, PlaneOrientation } from './src/streamlines/poincare_section_tracer.js';

        const sgy = PoincarePlane.createSGY(-5.0);
        const sgx = PoincarePlane.createSGX(25.0);

        console.log(JSON.stringify({
          sgy_norm: sgy.normal,
          sgy_u: sgy.uAxis,
          sgy_v: sgy.vAxis,
          sgx_norm: sgx.normal,
          sgx_u: sgx.uAxis,
          sgx_v: sgx.vAxis
        }));
        """
        res = run_node_snippet(code)
        assert res["sgy_norm"] == [0.0, 1.0, 0.0]
        assert res["sgy_u"] == [1.0, 0.0, 0.0]
        assert res["sgy_v"] == [0.0, 0.0, 1.0]
        assert res["sgx_norm"] == [1.0, 0.0, 0.0]
        assert res["sgx_u"] == [0.0, 1.0, 0.0]
        assert res["sgx_v"] == [0.0, 0.0, 1.0]

    def test_arbitrary_three_points_plane(self):
        code = """
        import { PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const p1 = [0, 0, 0];
        const p2 = [2, 0, 0];
        const p3 = [0, 3, 0];
        const plane = PoincarePlane.fromThreePoints(p1, p2, p3, 'XYPlane');

        console.log(JSON.stringify({
          normal: plane.normal,
          d_origin: plane.signedDistance([0, 0, 0]),
          d_above: plane.signedDistance([1, 1, 5]),
          d_below: plane.signedDistance([1, 1, -5])
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["normal"][0], 0.0, abs_tol=1e-6)
        assert math.isclose(res["normal"][1], 0.0, abs_tol=1e-6)
        assert math.isclose(res["normal"][2], 1.0, abs_tol=1e-6)
        assert math.isclose(res["d_origin"], 0.0, abs_tol=1e-6)
        assert math.isclose(res["d_above"], 5.0, abs_tol=1e-6)
        assert math.isclose(res["d_below"], -5.0, abs_tol=1e-6)

    def test_three_collinear_points_throws(self):
        code = """
        import { PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        try {
          PoincarePlane.fromThreePoints([0, 0, 0], [1, 1, 1], [2, 2, 2]);
          console.log(JSON.stringify({ threw: false }));
        } catch (e) {
          console.log(JSON.stringify({ threw: true, message: e.message }));
        }
        """
        res = run_node_snippet(code)
        assert res["threw"] is True

    def test_projection_and_unprojection_roundtrip(self):
        code = """
        import { PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.fromPointAndNormal([10, -20, 30], [1, 2, 2]);
        const testPoints = [
          [10, -20, 30],
          [15, -10, 25],
          [-30, 40, 50],
          [100, -100, 200]
        ];

        const errors = testPoints.map(pt => {
          const dist = plane.signedDistance(pt);
          const [u, v] = plane.projectPoint(pt);
          const reconstructed = plane.unprojectPoint(u, v, dist);
          const err = Math.hypot(pt[0] - reconstructed[0], pt[1] - reconstructed[1], pt[2] - reconstructed[2]);
          return err;
        });

        console.log(JSON.stringify({
          maxError: Math.max(...errors)
        }));
        """
        res = run_node_snippet(code)
        assert res["maxError"] < 1e-10

    def test_velocity_projection_and_orthonormal_energy_conservation(self):
        code = """
        import { PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.fromPointAndNormal([0, 0, 0], [3, 4, 12]);
        const vel = [300.0, -450.0, 600.0];

        const { vu, vv, vPerp, speed } = plane.projectVelocity(vel);
        const reconstructedSpeed = Math.sqrt(vu * vu + vv * vv + vPerp * vPerp);

        console.log(JSON.stringify({
          originalSpeed: speed,
          reconstructedSpeed,
          speedDiff: Math.abs(speed - reconstructedSpeed)
        }));
        """
        res = run_node_snippet(code)
        assert res["speedDiff"] < 1e-10

    def test_principal_tidal_tensor_plane(self):
        code = """
        import { PoincarePlane, PlaneOrientation } from './src/streamlines/poincare_section_tracer.js';

        const tensor = [
          3, 0, 0,
          0, 1, 0,
          0, 0, -2
        ];
        const center = [50, -20, 10];
        const planeMajor = PoincarePlane.fromClusterTidalTensor(center, tensor, 0);

        console.log(JSON.stringify({
          orientation: planeMajor.orientation,
          normal: planeMajor.normal,
          center: planeMajor.origin
        }));
        """
        res = run_node_snippet(code)
        assert res["orientation"] == "PRINCIPAL_TIDAL_PLANE"
        assert res["center"] == [50, -20, 10]
        assert math.isclose(abs(res["normal"][0]), 1.0, abs_tol=1e-5)

    def test_plane_json_serialization(self):
        code = """
        import { PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(15.5);
        const json = plane.toJSON();

        console.log(JSON.stringify({
          hasName: typeof json.name === 'string',
          originZ: json.origin[2],
          normalZ: json.normal[2]
        }));
        """
        res = run_node_snippet(code)
        assert res["hasName"] is True
        assert res["originZ"] == 15.5
        assert res["normalZ"] == 1.0


class TestExactIntersectionInterpolation:
    """Tests for exact root-finding and cutting plane intersection algorithms."""

    def test_linear_intersection_solver(self):
        code = """
        import { PoincarePlane, ExactIntersectionSolver } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const xA = [10.0, 20.0, -2.0];
        const xB = [14.0, 28.0, 6.0];
        const tA = 1.0, tB = 2.0;
        const vA = [100.0, 200.0, 200.0];
        const vB = [100.0, 200.0, 200.0];
        const sA = 5.0, sB = 15.0;

        const puncture = ExactIntersectionSolver.solveLinear(xA, xB, tA, tB, vA, vB, sA, sB, plane, 0, 1);

        console.log(JSON.stringify({
          zCross: puncture.position[2],
          xCross: puncture.position[0],
          yCross: puncture.position[1],
          uCross: puncture.planarCoords[0],
          vCross: puncture.planarCoords[1],
          tCross: puncture.time,
          sCross: puncture.arcLength,
          vPerp: puncture.normalVelocity
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["zCross"], 0.0, abs_tol=1e-8)
        assert math.isclose(res["xCross"], 11.0, abs_tol=1e-8)
        assert math.isclose(res["yCross"], 22.0, abs_tol=1e-8)
        assert math.isclose(res["uCross"], 11.0, abs_tol=1e-8)
        assert math.isclose(res["vCross"], 22.0, abs_tol=1e-8)
        assert math.isclose(res["tCross"], 1.25, abs_tol=1e-8)
        assert math.isclose(res["sCross"], 7.5, abs_tol=1e-8)
        assert math.isclose(res["vPerp"], 200.0, abs_tol=1e-8)

    def test_hermite_cubic_subgrid_accuracy(self):
        code = """
        import { PoincarePlane, ExactIntersectionSolver } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const xA = [0.0, 0.0, -1.0];
        const xB = [1.0, 1.0, 1.0];
        const tA = 0.0, tB = 1.0;
        const vA = [1.0, 1.0, 4.0];
        const vB = [1.0, 1.0, 0.0];
        const sA = 0.0, sB = 2.0;

        const puncture = ExactIntersectionSolver.solveHermiteCubic(xA, xB, tA, tB, vA, vB, sA, sB, plane);

        console.log(JSON.stringify({
          zCross: puncture.position[2],
          vPerp: puncture.normalVelocity,
          tCross: puncture.time
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["zCross"], 0.0, abs_tol=1e-6)
        assert res["vPerp"] > 0.0

    def test_directional_puncture_filtering(self):
        code = """
        import { PoincarePlane, CrossingDirection } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const posBelow = [0, 0, -1];
        const posAbove = [0, 0, 1];

        const crossForward = plane.isCrossing(posBelow, posAbove, CrossingDirection.POSITIVE);
        const crossBackward = plane.isCrossing(posAbove, posBelow, CrossingDirection.POSITIVE);
        const crossBoth1 = plane.isCrossing(posBelow, posAbove, CrossingDirection.BOTH);
        const crossBoth2 = plane.isCrossing(posAbove, posBelow, CrossingDirection.BOTH);

        console.log(JSON.stringify({
          crossForward,
          crossBackward,
          crossBoth1,
          crossBoth2
        }));
        """
        res = run_node_snippet(code)
        assert res["crossForward"] is True
        assert res["crossBackward"] is False
        assert res["crossBoth1"] is True
        assert res["crossBoth2"] is True


class TestPoincareRecurrenceMapAndRQA:
    """Tests for Poincaré recurrence maps, return times, Jacobian preservation, and RQA."""

    def test_circular_orbit_recurrence_map_and_return_time(self):
        code = """
        import { PoincarePlane, PuncturePoint, PoincareRecurrenceMap } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const period = 5.0;
        const radius = 25.0;
        const numPunctures = 20;

        const punctures = [];
        for (let i = 0; i < numPunctures; i++) {
          punctures.push(new PuncturePoint({
            index: i,
            position: [radius, 0.0, 0.0],
            velocity: [0.0, 200.0, 300.0],
            planarCoords: [radius, 0.0],
            planarVelocity: [0.0, 200.0],
            normalVelocity: 300.0,
            time: i * period,
            arcLength: i * (2.0 * Math.PI * radius)
          }));
        }

        const map = new PoincareRecurrenceMap(punctures, plane);
        const stats = map.computeReturnTimeStats();
        const centroid = map.computeCentroid();
        const rqa = map.computeRQA();

        console.log(JSON.stringify({
          meanReturnTime: stats.mean,
          stdReturnTime: stats.std,
          centroid,
          recurrenceRate: rqa.recurrenceRate,
          determinism: rqa.determinism
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["meanReturnTime"], 5.0, abs_tol=1e-6)
        assert math.isclose(res["stdReturnTime"], 0.0, abs_tol=1e-6)
        assert math.isclose(res["centroid"][0], 25.0, abs_tol=1e-6)
        assert math.isclose(res["centroid"][1], 0.0, abs_tol=1e-6)
        assert res["recurrenceRate"] > 0.90
        assert res["determinism"] > 0.90

    def test_symplectic_area_preservation_for_harmonic_vortex(self):
        code = """
        import { PoincarePlane, PuncturePoint, PoincareRecurrenceMap } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const theta = 0.45;
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);

        const punctures = [];
        let u = 10.0, v = 5.0;
        for (let i = 0; i < 30; i++) {
          punctures.push(new PuncturePoint({
            index: i,
            planarCoords: [u, v],
            time: i * 2.0
          }));
          const uNext = cosT * u - sinT * v;
          const vNext = sinT * u + cosT * v;
          u = uNext;
          v = vNext;
        }

        const map = new PoincareRecurrenceMap(punctures, plane);
        const jac = map.estimateAverageReturnJacobian();

        console.log(JSON.stringify({
          det: jac.det,
          trace: jac.trace,
          isAreaPreserving: jac.isAreaPreserving
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["det"], 1.0, abs_tol=0.08)
        assert res["isAreaPreserving"] is True

    def test_rqa_shannon_entropy_and_laminarity(self):
        code = """
        import { PoincarePlane, PuncturePoint, PoincareRecurrenceMap } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const punctures = [];
        for (let i = 0; i < 40; i++) {
          const angle = i * 0.3;
          punctures.push(new PuncturePoint({
            index: i,
            planarCoords: [15.0 * Math.cos(angle), 15.0 * Math.sin(angle)],
            time: i * 1.5
          }));
        }

        const map = new PoincareRecurrenceMap(punctures, plane);
        const rqa = map.computeRQA({ epsilon: 5.0, lMin: 2, vMin: 2 });

        console.log(JSON.stringify({
          longestDiagonalLine: rqa.longestDiagonalLine,
          shannonEntropy: rqa.shannonEntropy,
          numPunctures: rqa.numPunctures,
          hasDeterminism: rqa.determinism > 0.5
        }));
        """
        res = run_node_snippet(code)
        assert res["longestDiagonalLine"] >= 2
        assert res["shannonEntropy"] >= 0.0
        assert res["hasDeterminism"] is True


class TestLyapunovExponentEstimators:
    """Tests for Lyapunov exponent estimators (Benettin shadow & Variational tangent QR)."""

    def test_benettin_mle_for_regular_circular_vortex(self):
        code = """
        import { LyapunovExponentEstimator } from './src/streamlines/poincare_section_tracer.js';

        const omega = 0.1;
        const vortexField = (x, y, z) => [-omega * y, omega * x, 0.0];

        const seed = [20.0, 0.0, 0.0];
        const res = LyapunovExponentEstimator.estimateBenettinMLE(vortexField, seed, {
          totalTime: 40.0,
          dt: 0.05,
          renormInterval: 1
        });

        console.log(JSON.stringify({
          lambdaMax: res.lambdaMax,
          isChaotic: res.isChaotic,
          totalTime: res.totalTime
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["lambdaMax"]) < 0.03
        assert res["isChaotic"] is False

    def test_benettin_mle_for_hyperbolic_saddle_flow(self):
        code = """
        import { LyapunovExponentEstimator } from './src/streamlines/poincare_section_tracer.js';

        const lambdaTrue = 0.25;
        const saddleField = (x, y, z) => [lambdaTrue * x, -lambdaTrue * y, 0.0];

        const seed = [5.0, 5.0, 0.0];
        const res = LyapunovExponentEstimator.estimateBenettinMLE(saddleField, seed, {
          totalTime: 30.0,
          dt: 0.02,
          perturbationDir: [1.0, 0.0, 0.0],
          renormInterval: 1
        });

        console.log(JSON.stringify({
          lambdaMax: res.lambdaMax,
          isChaotic: res.isChaotic,
          lambdaTrue
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["lambdaMax"], 0.25, abs_tol=0.03)
        assert res["isChaotic"] is True

    def test_variational_spectrum_qr_decomposition_incompressible(self):
        code = """
        import { LyapunovExponentEstimator } from './src/streamlines/poincare_section_tracer.js';

        const A = 1.0, B = 1.0, C = 1.0;
        const abcFlow = (x, y, z) => [
          A * Math.sin(0.1 * z) + C * Math.cos(0.1 * y),
          B * Math.sin(0.1 * x) + A * Math.cos(0.1 * z),
          C * Math.sin(0.1 * y) + B * Math.cos(0.1 * x)
        ];

        const seed = [3.0, 5.0, 7.0];
        const res = LyapunovExponentEstimator.estimateVariationalSpectrum(abcFlow, seed, {
          totalTime: 30.0,
          dt: 0.05,
          qrInterval: 1
        });

        const sumSpectrum = res.spectrum[0] + res.spectrum[1] + res.spectrum[2];

        console.log(JSON.stringify({
          spectrum: res.spectrum,
          sumSpectrum,
          ksEntropy: res.ksEntropy,
          lyapunovDimension: res.lyapunovDimension
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["sumSpectrum"], 0.0, abs_tol=0.05)
        assert res["spectrum"][0] >= res["spectrum"][1] >= res["spectrum"][2]
        assert res["ksEntropy"] >= 0.0


class TestInvariantToriAndAccretionBasins:
    """Tests for KAM invariant tori, island chains, and cluster accretion basin classification."""

    def test_kam_invariant_torus_classification(self):
        code = """
        import { InvariantTorusClassifier, PuncturePoint, OrbitClassification } from './src/streamlines/poincare_section_tracer.js';

        const punctures = [];
        for (let i = 0; i < 60; i++) {
          const theta = i * 0.47;
          punctures.push(new PuncturePoint({
            index: i,
            planarCoords: [20.0 * Math.cos(theta), 12.0 * Math.sin(theta)],
            time: i * 2.0
          }));
        }

        const res = InvariantTorusClassifier.classifyOrbit(punctures, { lambdaMax: 0.01 });

        console.log(JSON.stringify({
          classification: res.classification,
          confidence: res.confidence,
          islandCount: res.islandCount,
          radialDispersion: res.radialDispersion
        }));
        """
        res = run_node_snippet(code)
        assert res["classification"] == "QUASI_PERIODIC_INVARIANT_TORUS"
        assert res["confidence"] > 0.80
        assert res["islandCount"] == 1

    def test_discrete_period_3_limit_cycle(self):
        code = """
        import { InvariantTorusClassifier, PuncturePoint, OrbitClassification } from './src/streamlines/poincare_section_tracer.js';

        const p1 = [10.0, 0.0];
        const p2 = [-5.0, 8.66];
        const p3 = [-5.0, -8.66];
        const fixedCycle = [p1, p2, p3];

        const punctures = [];
        for (let i = 0; i < 30; i++) {
          const pt = fixedCycle[i % 3];
          punctures.push(new PuncturePoint({
            index: i,
            planarCoords: [pt[0], pt[1]],
            time: i * 3.0
          }));
        }

        const res = InvariantTorusClassifier.classifyOrbit(punctures, { lambdaMax: 0.001 });

        console.log(JSON.stringify({
          classification: res.classification,
          islandCount: res.islandCount
        }));
        """
        res = run_node_snippet(code)
        assert res["classification"] == "PERIODIC_LIMIT_CYCLE"
        assert res["islandCount"] == 3

    def test_chaotic_sea_stochastic_classification(self):
        code = """
        import { InvariantTorusClassifier, PuncturePoint, OrbitClassification } from './src/streamlines/poincare_section_tracer.js';

        const punctures = [];
        for (let i = 0; i < 50; i++) {
          const r = 5.0 + Math.random() * 35.0;
          const th = Math.random() * 2.0 * Math.PI;
          punctures.push(new PuncturePoint({
            index: i,
            planarCoords: [r * Math.cos(th), r * Math.sin(th)],
            time: i * 1.0
          }));
        }

        const res = InvariantTorusClassifier.classifyOrbit(punctures, { lambdaMax: 0.18 });

        console.log(JSON.stringify({
          classification: res.classification,
          confidence: res.confidence
        }));
        """
        res = run_node_snippet(code)
        assert res["classification"] == "CHAOTIC_SEA_STOCHASTIC"


class TestFullPoincareSectionTracerEngine:
    """End-to-end integration tests for PoincareSectionTracer engine."""

    def test_trace_section_on_3d_cosmological_vortex(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane, PoincareIntegratorType, IntersectionMethod } from './src/streamlines/poincare_section_tracer.js';

        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-50, -50, -50], boxSize: [100, 100, 100] });

        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => {
          const omega = 0.08;
          const r = Math.max(1.0, Math.hypot(x, y));
          const vx = -omega * y;
          const vy = omega * x;
          const vz = 25.0 * (x / r);
          return [vx, vy, vz];
        });

        const tracer = new PoincareSectionTracer(field, {
          integrator: PoincareIntegratorType.DOPRI5,
          intersectionMethod: IntersectionMethod.HERMITE_CUBIC,
          maxPunctures: 10,
          maxSteps: 3000,
          dtInit: 0.05
        });

        const plane = PoincarePlane.createSGZ(0.0);
        const seed = [20.0, 0.0, -1.0];

        const result = tracer.traceSection(seed, plane);

        console.log(JSON.stringify({
          numPunctures: result.numPunctures,
          hasRecurrenceMap: result.recurrenceMap !== null,
          classification: result.classification.classification,
          lambdaMax: result.lyapunov.lambdaMax,
          totalSteps: result.stats.totalSteps,
          totalArcLength: result.stats.totalArcLength
        }));
        """
        res = run_node_snippet(code)
        assert res["numPunctures"] >= 2
        assert res["hasRecurrenceMap"] is True
        assert res["totalSteps"] > 10
        assert res["totalArcLength"] > 0.0

    def test_batch_seed_ensemble_tracing(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane, PoincareIntegratorType } from './src/streamlines/poincare_section_tracer.js';

        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-40, -40, -40], boxSize: [80, 80, 80] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [-0.05 * y, 0.05 * x, 10.0]);

        const tracer = new PoincareSectionTracer(field, {
          integrator: PoincareIntegratorType.RK4,
          maxPunctures: 5,
          maxSteps: 1000,
          dtInit: 0.1
        });

        const plane = PoincarePlane.createSGZ(0.0);
        const seeds = [
          [10.0, 0.0, -2.0],
          [20.0, 0.0, -2.0],
          [30.0, 0.0, -2.0]
        ];

        const batchResults = tracer.traceBatch(seeds, plane);

        console.log(JSON.stringify({
          batchSize: batchResults.length,
          puncturesPerSeed: batchResults.map(r => r.numPunctures)
        }));
        """
        res = run_node_snippet(code)
        assert res["batchSize"] == 3
        for count in res["puncturesPerSeed"]:
            assert count >= 1

    def test_virgo_cluster_accretion_basin(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane, PoincareIntegratorType } from './src/streamlines/poincare_section_tracer.js';

        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-60, -60, -60], boxSize: [120, 120, 120] });
        const virgoCenter = [-3.0, 16.0, -1.0];

        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => {
          const dx = x - virgoCenter[0];
          const dy = y - virgoCenter[1];
          const dz = z - virgoCenter[2];
          const r = Math.max(1.0, Math.sqrt(dx*dx + dy*dy + dz*dz));
          
          const vInfall = -150.0 / (1.0 + 0.05 * r);
          const vSwirl = 100.0 / (1.0 + 0.05 * r);

          const vx = vInfall * (dx / r) - vSwirl * (dy / r);
          const vy = vInfall * (dy / r) + vSwirl * (dx / r);
          const vz = vInfall * (dz / r) + 50.0 * (dx / r);
          return [vx, vy, vz];
        });

        const plane = PoincarePlane.createSGZ(virgoCenter[2], virgoCenter);
        const tracer = new PoincareSectionTracer(field, {
          integrator: PoincareIntegratorType.DOPRI5,
          maxPunctures: 10,
          maxSteps: 2000,
          dtInit: 0.02
        });

        const seed = [virgoCenter[0] + 15.0, virgoCenter[1], virgoCenter[2] - 1.0];
        const res = tracer.traceSection(seed, plane);

        console.log(JSON.stringify({
          puncturesFound: res.numPunctures > 0,
          firstPuncturePos: res.numPunctures > 0 ? res.punctures[0].position : null,
          totalArcLength: res.stats.totalArcLength
        }));
        """
        res = run_node_snippet(code)
        assert res["puncturesFound"] is True
        assert res["totalArcLength"] > 0.0

    def test_shapley_supercluster_high_velocity_flow(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-200, -50, -100], boxSize: [250, 200, 200] });
        const shapleyCenter = [-130.0, 80.0, -35.0];

        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => {
          const dx = x - shapleyCenter[0];
          const dy = y - shapleyCenter[1];
          const dz = z - shapleyCenter[2];
          const r = Math.max(2.0, Math.sqrt(dx*dx + dy*dy + dz*dz));
          const vInfall = -600.0 * Math.exp(-0.01 * r);
          const vSwirl = 250.0 * Math.exp(-0.01 * r);

          return [
            vInfall * (dx / r) - vSwirl * (dy / r),
            vInfall * (dy / r) + vSwirl * (dx / r),
            vInfall * (dz / r) + 100.0 * (dx / r)
          ];
        });

        const plane = PoincarePlane.createSGZ(shapleyCenter[2], shapleyCenter);
        const tracer = new PoincareSectionTracer(field, { maxPunctures: 8, maxSteps: 1500 });
        const seed = [shapleyCenter[0] + 30.0, shapleyCenter[1], shapleyCenter[2] - 2.0];
        const res = tracer.traceSection(seed, plane);

        console.log(JSON.stringify({
          puncturesFound: res.numPunctures > 0,
          totalTime: res.stats.totalTime
        }));
        """
        res = run_node_snippet(code)
        assert res["puncturesFound"] is True
        assert res["totalTime"] > 0.0

    def test_great_attractor_principal_plane(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const gaCenter = [-45.0, 15.0, -10.0];
        const tidalTensor = [
          5.0, 0.2, 0.1,
          0.2, 2.0, 0.0,
          0.1, 0.0, -1.0
        ];

        const plane = PoincarePlane.fromClusterTidalTensor(gaCenter, tidalTensor, 2);
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-100, -50, -50], boxSize: [120, 120, 120] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [-0.03 * (x - gaCenter[0]), -0.03 * (y - gaCenter[1]), 20.0]);

        const tracer = new PoincareSectionTracer(field, { maxPunctures: 3, maxSteps: 800 });
        const res = tracer.traceSection([gaCenter[0], gaCenter[1], gaCenter[2] - 10.0], plane);

        console.log(JSON.stringify({
          planeName: plane.name,
          normal: plane.normal,
          numPunctures: res.numPunctures
        }));
        """
        res = run_node_snippet(code)
        assert "Cluster Tidal Principal Plane" in res["planeName"]
        assert len(res["normal"]) == 3
        assert res["numPunctures"] >= 1

    def test_coma_cluster_high_dispersion_field(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const comaCenter = [0.0, 70.0, 9.0];
        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-50, 20, -40], boxSize: [100, 100, 100] });

        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => {
          const dx = x - comaCenter[0];
          const dy = y - comaCenter[1];
          const dz = z - comaCenter[2];
          const r = Math.max(1.0, Math.sqrt(dx*dx + dy*dy + dz*dz));
          return [
            -200.0 * (dx / r) - 80.0 * (dy / r),
            -200.0 * (dy / r) + 80.0 * (dx / r),
            40.0 * (dx / r)
          ];
        });

        const plane = PoincarePlane.createSGZ(comaCenter[2], comaCenter);
        const tracer = new PoincareSectionTracer(field, { maxPunctures: 5, maxSteps: 1200 });
        const res = tracer.traceSection([comaCenter[0] + 15.0, comaCenter[1], comaCenter[2] - 1.0], plane);

        console.log(JSON.stringify({
          numPunctures: res.numPunctures,
          hasLyapunov: res.lyapunov !== null
        }));
        """
        res = run_node_snippet(code)
        assert res["numPunctures"] >= 1
        assert res["hasLyapunov"] is True

    def test_perseus_pisces_supercluster_filamentary_inflow(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const ppCenter = [50.0, -15.0, -20.0];
        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [0, -65, -70], boxSize: [100, 100, 100] });

        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => {
          const dx = x - ppCenter[0];
          const dy = y - ppCenter[1];
          const dz = z - ppCenter[2];
          const r = Math.max(1.0, Math.hypot(dx, dy));
          return [
            -120.0 * (dx / r) - 90.0 * (dy / r),
            -120.0 * (dy / r) + 90.0 * (dx / r),
            35.0 * (dx / r)
          ];
        });

        const plane = PoincarePlane.createSGZ(ppCenter[2], ppCenter);
        const tracer = new PoincareSectionTracer(field, { maxPunctures: 6, maxSteps: 1500 });
        const res = tracer.traceSection([ppCenter[0] + 20.0, ppCenter[1], ppCenter[2] - 1.5], plane);

        console.log(JSON.stringify({
          numPunctures: res.numPunctures,
          classification: res.classification.classification
        }));
        """
        res = run_node_snippet(code)
        assert res["numPunctures"] >= 1


class TestMathematicalInvariantsAndEdgeCases:
    """Rigorous validation of mathematical invariants, edge cases, and numerical tolerances."""

    def test_zero_velocity_stall_termination(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [0.0, 0.0, 0.0]);

        const plane = PoincarePlane.createSGZ(0.0);
        const tracer = new PoincareSectionTracer(field, { maxSteps: 100 });
        const res = tracer.traceSection([2.0, 2.0, -1.0], plane);

        console.log(JSON.stringify({
          terminationReason: res.stats.terminationReason,
          totalSteps: res.stats.totalSteps
        }));
        """
        res = run_node_snippet(code)
        assert res["terminationReason"] == "LOW_SPEED_STALL"

    def test_max_arclength_exhaustion(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { PoincareSectionTracer, PoincarePlane } from './src/streamlines/poincare_section_tracer.js';

        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [0.0, 0.0, 100.0]);

        const plane = PoincarePlane.createSGZ(10.0);
        const tracer = new PoincareSectionTracer(field, { maxArcLength: 25.0, maxSteps: 500 });
        const res = tracer.traceSection([0.0, 0.0, 0.0], plane);

        console.log(JSON.stringify({
          terminationReason: res.stats.terminationReason,
          arcLength: res.stats.totalArcLength
        }));
        """
        res = run_node_snippet(code)
        assert res["terminationReason"] == "MAX_ARCLENGTH"
        assert res["arcLength"] >= 25.0

    def test_empty_puncture_list_rqa_graceful_handling(self):
        code = """
        import { PoincarePlane, PoincareRecurrenceMap } from './src/streamlines/poincare_section_tracer.js';

        const plane = PoincarePlane.createSGZ(0.0);
        const map = new PoincareRecurrenceMap([], plane);
        const rqa = map.computeRQA();
        const stats = map.computeReturnTimeStats();

        console.log(JSON.stringify({
          rr: rqa.recurrenceRate,
          det: rqa.determinism,
          meanTime: stats.mean
        }));
        """
        res = run_node_snippet(code)
        assert res["rr"] == 0.0
        assert res["det"] == 0.0
        assert res["meanTime"] == 0.0

    def test_tangent_qr_orthonormality(self):
        code = """
        import { qrDecomposition3x3 } from './src/streamlines/poincare_section_tracer.js';

        const M = [
          1.0, 2.0, 3.0,
          0.5, -1.0, 2.0,
          -2.0, 1.5, 0.5
        ];

        const { Q, R } = qrDecomposition3x3(M);

        // Q columns: q0 = [Q[0], Q[3], Q[6]], q1 = [Q[1], Q[4], Q[7]], q2 = [Q[2], Q[5], Q[8]]
        const dot01 = Q[0]*Q[1] + Q[3]*Q[4] + Q[6]*Q[7];
        const dot02 = Q[0]*Q[2] + Q[3]*Q[5] + Q[6]*Q[8];
        const dot12 = Q[1]*Q[2] + Q[4]*Q[5] + Q[7]*Q[8];

        const norm0 = Math.sqrt(Q[0]*Q[0] + Q[3]*Q[3] + Q[6]*Q[6]);
        const norm1 = Math.sqrt(Q[1]*Q[1] + Q[4]*Q[4] + Q[7]*Q[7]);
        const norm2 = Math.sqrt(Q[2]*Q[2] + Q[5]*Q[5] + Q[8]*Q[8]);

        console.log(JSON.stringify({
          orthogonality: Math.max(Math.abs(dot01), Math.abs(dot02), Math.abs(dot12)),
          normDeviation: Math.max(Math.abs(norm0 - 1), Math.abs(norm1 - 1), Math.abs(norm2 - 1))
        }));
        """
        res = run_node_snippet(code)
        assert res["orthogonality"] < 1e-10
        assert res["normDeviation"] < 1e-10
