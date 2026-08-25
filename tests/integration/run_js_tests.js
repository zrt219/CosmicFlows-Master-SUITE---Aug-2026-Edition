/**
 * @file run_js_tests.js
 * @description Node.js test harness for validating ES6 integrators, seed manager, and streamline tracer.
 * Returns structured JSON test results for Pytest integration.
 */

import { rk4Step, traceRK4Streamline, RK4Integrator, TerminationReason, IntegrationDirection } from '../../src/integration/rk4_classical.js';
import { rk45CashKarpStep, cashKarpStepRaw, traceCashKarpStreamline, RK45CashKarpIntegrator, CASH_KARP_COEFFS } from '../../src/integration/rk45_cash_karp.js';
import { dormandPrinceStep, dormandPrinceStepRaw, traceDormandPrinceStreamline, denseInterpolateDOPRI5, DormandPrinceIntegrator, DOPRI5_COEFFS } from '../../src/integration/dormand_prince.js';
import { SeedManager, StreamlineSeed, SeedMode, COSMIC_LANDMARKS } from '../../src/streamlines/seed_manager.js';
import { StreamlineTracer, Streamline, IntegratorType, TracingMode } from '../../src/streamlines/streamline_tracer.js';

function assert(condition, message) {
  if (!condition) {
    throw new Error(`Assertion failed: ${message}`);
  }
}

function assertClose(a, b, tol = 1e-5, message = '') {
  const diff = Math.abs(a - b);
  if (diff > tol) {
    throw new Error(`Assertion failed: ${a} not close to ${b} (diff=${diff}, tol=${tol}). ${message}`);
  }
}

export async function runAllTests() {
  const results = {
    passed: 0,
    failed: 0,
    tests: []
  };

  function test(name, fn) {
    try {
      fn();
      results.passed++;
      results.tests.push({ name, status: 'PASSED' });
    } catch (err) {
      results.failed++;
      results.tests.push({ name, status: 'FAILED', error: err.message, stack: err.stack });
    }
  }

  // =========================================================================
  // 1. Classical RK4 Tests
  // =========================================================================
  test('RK4: Exact Linear Field Convergence Order O(h^4)', () => {
    // dx/dt = -0.5 * x, x(0) = 10, exact x(t) = 10 * exp(-0.5 * t)
    const field = (t, [x, y, z]) => [-0.5 * x, -0.5 * y, -0.5 * z];
    const T = 2.0;
    const xExact = 10.0 * Math.exp(-0.5 * T);

    const errors = [];
    const steps = [0.2, 0.1, 0.05];

    for (const h of steps) {
      let x = 10.0;
      let y = 10.0;
      let z = 10.0;
      let t = 0.0;
      const n = Math.round(T / h);
      for (let i = 0; i < n; i++) {
        const res = rk4Step(field, [x, y, z], t, h);
        x = res.posNext[0];
        y = res.posNext[1];
        z = res.posNext[2];
        t = res.tNext;
      }
      errors.push(Math.abs(x - xExact));
    }

    const ratio1 = errors[0] / errors[1]; // should be ~ 2^4 = 16
    const ratio2 = errors[1] / errors[2];
    assert(ratio1 >= 15.0 && ratio1 <= 17.0, `RK4 convergence ratio 1 expected ~16, got ${ratio1}`);
    assert(ratio2 >= 15.0 && ratio2 <= 17.0, `RK4 convergence ratio 2 expected ~16, got ${ratio2}`);
  });

  test('RK4: Harmonic Oscillator Circular Orbit Invariance', () => {
    // v = [-y, x, 0], r = sqrt(x^2 + y^2) = const, exact orbit radius = 5
    const field = (t, [x, y, z]) => [-y, x, 0];
    const res = traceRK4Streamline(field, [5.0, 0.0, 0.0], {
      dt: 0.01,
      maxSteps: 628 // ~ 2*pi / 0.01
    });

    assert(res.points.length > 500, `Expected >500 points, got ${res.points.length}`);
    const lastP = res.points[res.points.length - 1];
    const rLast = Math.hypot(lastP[0], lastP[1]);
    assertClose(rLast, 5.0, 1e-4, 'Radius must remain conserved in circular vortex');
  });

  // =========================================================================
  // 2. Cash-Karp RK45 Tests
  // =========================================================================
  test('Cash-Karp: Error Estimation and Adaptation on Plummer Field', () => {
    // Plummer potential: v = -GM*r / (r^2 + b^2)^(3/2)
    const GM = 100.0;
    const b = 5.0;
    const plummer = (pos) => {
      const r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
      const factor = -GM / Math.pow(r2 + b*b, 1.5);
      return [factor * pos[0], factor * pos[1], factor * pos[2]];
    };

    const integrator = new RK45CashKarpIntegrator({ atol: 1e-7, rtol: 1e-7, initialStep: 0.1 });
    const trace = integrator.trace(plummer, [50.0, 0.0, 0.0], { maxSteps: 500 });

    assert(trace.points.length > 10, 'Streamline should have multiple points');
    assert(trace.stats.acceptedSteps > 0, 'Should have accepted steps');
    // Radial straight line toward center: y and z must stay 0
    for (const p of trace.points) {
      assertClose(p[1], 0.0, 1e-6, 'Plummer y-coordinate must be 0');
      assertClose(p[2], 0.0, 1e-6, 'Plummer z-coordinate must be 0');
    }
  });

  test('Cash-Karp: Convergence Order O(h^5) for 5th-order output', () => {
    const field = (t, [x, y, z]) => [-x, -y, -z];
    const T = 1.0;
    const xExact = 10.0 * Math.exp(-1.0);

    const steps = [0.2, 0.1, 0.05];
    const errors = [];

    for (const h of steps) {
      let pos = [10.0, 0, 0];
      let t = 0.0;
      const n = Math.round(T / h);
      for (let i = 0; i < n; i++) {
        const raw = cashKarpStepRaw(field, pos, t, h);
        pos = raw.pos5;
        t += h;
      }
      errors.push(Math.abs(pos[0] - xExact));
    }

    const ratio = errors[0] / errors[1]; // ~ 2^5 = 32
    assert(ratio >= 28.0 && ratio <= 36.0, `Cash-Karp 5th order ratio expected ~32, got ${ratio}`);
  });

  // =========================================================================
  // 3. Dormand-Prince DOPRI5 Tests
  // =========================================================================
  test('DOPRI5: FSAL Property Verification', () => {
    const field = (t, [x, y, z]) => [Math.cos(t) - x, Math.sin(t) - y, -z];
    const step1 = dormandPrinceStep(field, [1.0, 2.0, 3.0], 0.0, 0.1);
    assert(step1.accepted, 'Step 1 must be accepted');

    // Run step 2 with cached FSAL stage
    const step2 = dormandPrinceStep(field, step1.posNext, step1.tNext, 0.1, {}, step1.fsalNext);
    assert(step2.accepted, 'Step 2 must be accepted');
    // First stage of step 2 should match fsalNext of step 1
    assertClose(step2.velocity[0], step1.fsalNext[0], 1e-14, 'FSAL k1 must match previous k7');
    assertClose(step2.velocity[1], step1.fsalNext[1], 1e-14, 'FSAL k1 must match previous k7');
    assertClose(step2.velocity[2], step1.fsalNext[2], 1e-14, 'FSAL k1 must match previous k7');
  });

  test('DOPRI5: Dense Output Continuity and Accuracy', () => {
    // Exact ODE dx/dt = -x, x(t) = exp(-t)
    const field = (t, [x, y, z]) => [-x, 0, 0];
    const h = 0.5;
    const raw = dormandPrinceStepRaw(field, [1.0, 0, 0], 0.0, h);

    // Test intermediate points theta = 0.25, 0.5, 0.75
    for (const theta of [0.25, 0.5, 0.75]) {
      const denseP = denseInterpolateDOPRI5([1.0, 0, 0], raw.pos5, h, raw.kStages, theta);
      const exactX = Math.exp(-theta * h);
      assertClose(denseP[0], exactX, 1e-4, `Dense output at theta=${theta} matches exact solution`);
    }
  });

  // =========================================================================
  // 4. Seed Manager Tests
  // =========================================================================
  test('SeedManager: GRID_VOXEL lattice, strides, and bounding box', () => {
    const sm = new SeedManager();
    const seeds = sm.generateGridSeeds({
      resolution: [8, 8, 8],
      domain: [[-50, 50], [-50, 50], [-50, 50]],
      strides: [2, 2, 2],
      addToStore: true
    });

    // 8 / 2 = 4 along each dimension => 4^3 = 64 seeds
    assert(seeds.length === 64, `Expected 64 seeds, got ${seeds.length}`);
    assert(sm.size === 64, `Expected sm.size=64, got ${sm.size}`);

    const arr = sm.toFloat32Array();
    assert(arr.length === 64 * 3, `Expected Float32Array of length 192, got ${arr.length}`);
  });

  test('SeedManager: OBSERVED_OBJECT CF4 galaxy catalog ingestion and filtering', () => {
    const sm = new SeedManager();
    const mockCatalog = [
      { id: 'NGC1365', name: 'NGC 1365', sgx: -10.0, sgy: 25.0, sgz: -5.0, cz: 1636, isGroupMember: true, groupId: 101, massWeight: 10.5 },
      { id: 'ESO358-G006', name: 'ESO 358-G006', sgx: -12.0, sgy: 24.0, sgz: -6.0, cz: 1550, isGroupMember: true, groupId: 101, massWeight: 8.2 },
      { id: 'UGC12345', name: 'UGC 12345', sgx: 80.0, sgy: -90.0, sgz: 40.0, cz: 7500, isGroupMember: false, groupId: 0, massWeight: 4.1 }
    ];

    sm.loadObservedCatalog(mockCatalog, { maxDistance: 50.0 });
    assert(sm.size === 2, `Expected 2 galaxies within 50 Mpc/h, got ${sm.size}`);

    const seeds = sm.getAllSeeds();
    assert(seeds[0].metadata.name === 'NGC 1365', 'Galaxy metadata preserved');
    assert(seeds[0].mode === SeedMode.OBSERVED_OBJECT, 'Mode should be OBSERVED_OBJECT');
  });

  test('SeedManager: VISUALIZATION Fibonacci sphere and landmark presets', () => {
    const sm = new SeedManager();
    const fibSeeds = sm.generateFibonacciSphereSeeds([0, 0, 0], 20.0, 50, {}, false);
    assert(fibSeeds.length === 50, `Expected 50 seeds, got ${fibSeeds.length}`);

    // Verify all points are on sphere of radius 20
    for (const s of fibSeeds) {
      const r = Math.hypot(s.x, s.y, s.z);
      assertClose(r, 20.0, 1e-4, 'Point must lie on sphere surface');
    }

    // Landmark generation
    const shapleySeeds = sm.generateLandmarkSeeds('SHAPLEY_CORE', 2, 20, true);
    assert(shapleySeeds.length === 40, `Expected 40 landmark seeds, got ${shapleySeeds.length}`);
  });

  // =========================================================================
  // 5. Streamline Tracer & Termination Conditions Tests
  // =========================================================================
  test('StreamlineTracer: Termination on CONVERGED_ENDPOINT', () => {
    // Sinking field toward (0,0,0)
    const field = (pos) => [-2.0 * pos[0], -2.0 * pos[1], -2.0 * pos[2]];
    const tracer = new StreamlineTracer({
      integrator: IntegratorType.DORMAND_PRINCE,
      minVelocity: 1e-4,
      maxSteps: 2000
    });

    const sl = tracer.traceStreamline(field, [5.0, 5.0, 5.0]);
    assert(sl.terminationReason === TerminationReason.CONVERGED_ENDPOINT, `Expected CONVERGED_ENDPOINT, got ${sl.terminationReason}`);
    const lastP = sl.points[sl.points.length - 1];
    assert(Math.hypot(lastP[0], lastP[1], lastP[2]) < 0.1, 'Streamline should terminate near origin sink');
  });

  test('StreamlineTracer: Termination on DOMAIN_EXIT', () => {
    // Outward repelling field
    const field = (pos) => [pos[0] + 1.0, pos[1] + 1.0, pos[2] + 1.0];
    const tracer = new StreamlineTracer({
      domain: [[-20, 20], [-20, 20], [-20, 20]],
      maxSteps: 1000
    });

    const sl = tracer.traceStreamline(field, [15.0, 15.0, 15.0]);
    assert(sl.terminationReason === TerminationReason.DOMAIN_EXIT, `Expected DOMAIN_EXIT, got ${sl.terminationReason}`);
  });

  test('StreamlineTracer: Termination on MAX_STEPS', () => {
    // Constant circular vortex that never leaves domain and never stalls
    const field = (pos) => [-pos[1], pos[0], 0];
    const tracer = new StreamlineTracer({
      initialStep: 0.1,
      maxSteps: 50,
      domain: [[-100, 100], [-100, 100], [-100, 100]]
    });

    const sl = tracer.traceStreamline(field, [10.0, 0.0, 0.0]);
    assert(sl.terminationReason === TerminationReason.MAX_STEPS, `Expected MAX_STEPS, got ${sl.terminationReason}`);
    assert(sl.points.length >= 50, 'Should reach max steps');
  });

  test('StreamlineTracer: BIDIRECTIONAL Tracing Seamless Stitching', () => {
    const field = (pos) => [-pos[0], -pos[1], 0];
    const tracer = new StreamlineTracer({
      mode: TracingMode.BIDIRECTIONAL,
      domain: [[-50, 50], [-50, 50], [-50, 50]],
      maxSteps: 200
    });

    const sl = tracer.traceStreamline(field, [10.0, 5.0, 0.0]);
    assert(sl.direction === TracingMode.BIDIRECTIONAL, 'Direction should be bidirectional');
    assert(sl.points.length > 2, 'Should contain backward and forward points');

    // Arc length should be monotonically increasing
    for (let i = 1; i < sl.arcLengths.length; i++) {
      assert(sl.arcLengths[i] >= sl.arcLengths[i - 1], `Arc length must be monotonically increasing at ${i}`);
    }
  });

  test('StreamlineTracer: Batch Tracing & Progress Reporting', () => {
    const field = (pos) => [-pos[0], -pos[1], -pos[2]];
    const sm = new SeedManager();
    sm.generateFibonacciSphereSeeds([0, 0, 0], 30.0, 10, {}, true);

    const tracer = new StreamlineTracer({ maxSteps: 100 });
    let progressCalls = 0;

    const batch = tracer.traceBatch(field, sm, {}, (p) => {
      progressCalls++;
      assert(p.percentage >= 0 && p.percentage <= 100, 'Progress percentage within [0, 100]');
    });

    assert(batch.streamlines.length === 10, 'All 10 streamlines traced');
    assert(progressCalls === 10, `Expected 10 progress calls, got ${progressCalls}`);
    assert(batch.stats.totalVertices > 0, 'Vertices counted');
  });

  return results;
}

// Execute and output JSON
runAllTests().then(res => {
  console.log(JSON.stringify(res, null, 2));
  if (res.failed > 0) {
    process.exit(1);
  }
}).catch(err => {
  console.error(err);
  process.exit(1);
});
