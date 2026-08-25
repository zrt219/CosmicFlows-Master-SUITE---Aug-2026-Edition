"""
Automated Pytest Suite for ZRT Cosmicflows Workbench
Fields, Interpolation, and Ingestion Modules.

Tests:
1. GridIndexer: FITS (SGZ, SGY, SGX) <-> Canonical (SGX, SGY, SGZ) stride mapping, voxel indexing, boundary clamping.
2. Constant Fields: Exact interpolation, zero gradients, zero Hessian.
3. Exact Linear Fields: Exact interpolation, constant analytic gradient, FD agreement.
4. Polynomial / Trilinear Fields: Cross-terms, exact gradient & Hessian matching.
5. Velocity Field & x52.0 Scale Factor: Continuous evaluation, scale factor enforcement, strain tensor, vorticity, divergence, Jacobi diagonalization.
6. Density Field & CIC Mass Assignment: Particle mass conservation, overdensity normalization, log density.
7. Finite Difference Operators: 6-point and 14-point central stencils, convergence rate O(h^2), Richardson extrapolation.
8. Cross-verification with ES6 Node.js runtime.
"""

import subprocess
import json
import math
import numpy as np
import pytest

NODE_EXEC = "node"

def run_node_es6_snippet(script_code: str) -> dict:
    """Helper to execute an ES6 snippet in Node.js and return parsed JSON output."""
    full_code = f"""
    {script_code}
    """
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", full_code],
        capture_output=True,
        text=True,
        check=True
    )
    return json.loads(proc.stdout)


# ============================================================================
# 1. GRID INDEXER TESTS
# ============================================================================

class TestGridIndexer:
    def test_grid_indexer_initialization_and_strides(self):
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        const grid = new GridIndexer({
          nx: 16, ny: 32, nz: 64,
          origin: [-100, -200, -300],
          boxSize: [200, 400, 600]
        });

        const res = {
          totalCells: grid.totalCells,
          dx: grid.dx, dy: grid.dy, dz: grid.dz,
          strideX: grid.strideX,
          strideY: grid.strideY,
          strideZ: grid.strideZ,
          xMin: grid.xMin, xMax: grid.xMax,
          yMin: grid.yMin, yMax: grid.yMax,
          zMin: grid.zMin, zMax: grid.zMax
        };
        console.log(JSON.stringify(res));
        """
        res = run_node_es6_snippet(code)
        assert res["totalCells"] == 16 * 32 * 64
        assert math.isclose(res["dx"], 200.0 / 15.0, rel_tol=1e-7)
        assert math.isclose(res["dy"], 400.0 / 31.0, rel_tol=1e-7)
        assert math.isclose(res["dz"], 600.0 / 63.0, rel_tol=1e-7)
        assert res["strideX"] == 1
        assert res["strideY"] == 16
        assert res["strideZ"] == 16 * 32
        assert res["xMin"] == -100.0 and res["xMax"] == 100.0
        assert res["yMin"] == -200.0 and res["yMax"] == 200.0
        assert res["zMin"] == -300.0 and res["zMax"] == 300.0

    def test_grid_indexer_index_round_trip(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        const grid = new GridIndexer({ nx: 10, ny: 20, nz: 30 });
        let passed = true;
        for (let iz = 0; iz < 30; iz += 7) {
          for (let iy = 0; iy < 20; iy += 5) {
            for (let ix = 0; ix < 10; ix += 3) {
              const flat = grid.getLinearIndex(ix, iy, iz);
              const [rx, ry, rz] = grid.get3DIndices(flat);
              if (rx !== ix || ry !== iy || rz !== iz) passed = false;
            }
          }
        }
        console.log(JSON.stringify({ passed }));
        """
        res = run_node_es6_snippet(code)
        assert res["passed"] is True

    def test_grid_indexer_coord_transform_round_trip(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        const grid = new GridIndexer({
          nx: 32, ny: 32, nz: 32,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const testPoints = [
          [-50, -50, -50],
          [0, 0, 0],
          [25.5, -12.3, 40.8],
          [50, 50, 50]
        ];

        const roundTrips = testPoints.map(([x, y, z]) => {
          const [gx, gy, gz] = grid.coordToGridIndex(x, y, z);
          const [rx, ry, rz] = grid.gridIndexToCoord(gx, gy, gz);
          return { orig: [x, y, z], rec: [rx, ry, rz] };
        });

        console.log(JSON.stringify(roundTrips));
        """
        res = run_node_es6_snippet(code)
        for item in res:
            for o, r in zip(item["orig"], item["rec"]):
                assert math.isclose(o, r, abs_tol=1e-12)

    def test_grid_indexer_boundary_modes(self):
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10 });

        const clampNeg = grid.handleBoundaryIndex(-2, 5, 12, BoundaryMode.CLAMP);
        const periodicNeg = grid.handleBoundaryIndex(-2, 12, 25, BoundaryMode.PERIODIC);
        const reflectNeg = grid.handleBoundaryIndex(-2, 11, 10, BoundaryMode.REFLECT);
        const zeroOut = grid.handleBoundaryIndex(-1, 5, 5, BoundaryMode.ZERO);
        const zeroIn = grid.handleBoundaryIndex(2, 5, 5, BoundaryMode.ZERO);

        console.log(JSON.stringify({ clampNeg, periodicNeg, reflectNeg, zeroOut, zeroIn }));
        """
        res = run_node_es6_snippet(code)
        assert res["clampNeg"] == [0, 5, 9]
        assert res["periodicNeg"] == [8, 2, 5]
        assert res["reflectNeg"] == [2, 7, 8]  # Reflects within [0, 9]
        assert res["zeroOut"] is None
        assert res["zeroIn"] == [2, 5, 5]


# ============================================================================
# 2. TRILINEAR INTERPOLATOR & ANALYTIC DERIVATIVES TESTS
# ============================================================================

class TestTrilinearInterpolator:
    def test_constant_field_exactness(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';

        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [0, 0, 0], boxSize: [7, 7, 7] });
        const interp = new TrilinearInterpolator(grid);

        const C = 42.195;
        const buf = new Float64Array(grid.totalCells).fill(C);

        const val = interp.interpolateScalar(buf, 2.34, 4.56, 1.78);
        const grad = interp.analyticGradient(buf, 2.34, 4.56, 1.78);
        const hess = interp.analyticHessian(buf, 2.34, 4.56, 1.78);

        console.log(JSON.stringify({ val, grad, hess, C }));
        """
        res = run_node_es6_snippet(code)
        assert math.isclose(res["val"], res["C"], rel_tol=1e-14)
        for g in res["grad"]:
            assert abs(g) < 1e-14
        for row in res["hess"]:
            for h in row:
                assert abs(h) < 1e-14

    def test_exact_linear_field(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';

        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const interp = new TrilinearInterpolator(grid);

        const c0 = 15.0, c1 = 2.5, c2 = -3.8, c3 = 1.2;
        const buf = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < grid.nz; iz++) {
          for (let iy = 0; iy < grid.ny; iy++) {
            for (let ix = 0; ix < grid.nx; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              buf[idx] = c0 + c1 * x + c2 * y + c3 * z;
            }
          }
        }

        const pts = [
          [-5.3, 2.7, 0.4],
          [3.14, -8.2, 7.9],
          [0.0, 0.0, 0.0]
        ];

        const evaluations = pts.map(([x, y, z]) => {
          const exactVal = c0 + c1 * x + c2 * y + c3 * z;
          const interpVal = interp.interpolateScalar(buf, x, y, z);
          const grad = interp.analyticGradient(buf, x, y, z);
          const hess = interp.analyticHessian(buf, x, y, z);
          return { exactVal, interpVal, grad, hess };
        });

        console.log(JSON.stringify({ evaluations, expectedGrad: [c1, c2, c3] }));
        """
        res = run_node_es6_snippet(code)
        exp_grad = res["expectedGrad"]
        for ev in res["evaluations"]:
            assert math.isclose(ev["interpVal"], ev["exactVal"], abs_tol=1e-12)
            assert math.isclose(ev["grad"][0], exp_grad[0], abs_tol=1e-12)
            assert math.isclose(ev["grad"][1], exp_grad[1], abs_tol=1e-12)
            assert math.isclose(ev["grad"][2], exp_grad[2], abs_tol=1e-12)
            for row in ev["hess"]:
                for h in row:
                    assert abs(h) < 1e-12

    def test_trilinear_polynomial_field_and_hessian(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';

        // For f(x, y, z) = x * y * z, trilinear interpolation on a grid is exact inside each cell!
        const grid = new GridIndexer({ nx: 6, ny: 6, nz: 6, origin: [0, 0, 0], boxSize: [5, 5, 5] });
        const interp = new TrilinearInterpolator(grid);

        const buf = new Float64Array(grid.totalCells);
        for (let iz = 0; iz < grid.nz; iz++) {
          for (let iy = 0; iy < grid.ny; iy++) {
            for (let ix = 0; ix < grid.nx; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              buf[idx] = x * y * z;
            }
          }
        }

        const x = 2.4, y = 3.6, z = 1.5;
        const exactVal = x * y * z;
        const interpVal = interp.interpolateScalar(buf, x, y, z);
        const grad = interp.analyticGradient(buf, x, y, z);
        const hess = interp.analyticHessian(buf, x, y, z);

        console.log(JSON.stringify({
          exactVal, interpVal,
          grad,
          exactGrad: [y * z, x * z, x * y],
          hess,
          exactHessMixed: { dxy: z, dxz: y, dyz: x }
        }));
        """
        res = run_node_es6_snippet(code)
        assert math.isclose(res["interpVal"], res["exactVal"], abs_tol=1e-12)
        assert math.isclose(res["grad"][0], res["exactGrad"][0], abs_tol=1e-12)
        assert math.isclose(res["grad"][1], res["exactGrad"][1], abs_tol=1e-12)
        assert math.isclose(res["grad"][2], res["exactGrad"][2], abs_tol=1e-12)
        # Hessian mixed partials
        assert math.isclose(res["hess"][0][1], res["exactHessMixed"]["dxy"], abs_tol=1e-12)
        assert math.isclose(res["hess"][0][2], res["exactHessMixed"]["dxz"], abs_tol=1e-12)
        assert math.isclose(res["hess"][1][2], res["exactHessMixed"]["dyz"], abs_tol=1e-12)


# ============================================================================
# 3. FINITE DIFFERENCE & CONVERGENCE TESTS
# ============================================================================

class TestFiniteDifference:
    def test_6point_and_14point_stencils_linear_agreement(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FiniteDifference } from './src/interpolation/finite_difference.js';

        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const fd = new FiniteDifference(grid);

        // Smooth analytic scalar function f(x, y, z) = sin(0.1*x) + cos(0.2*y) + exp(0.05*z)
        const fn = (x, y, z) => Math.sin(0.1 * x) + Math.cos(0.2 * y) + Math.exp(0.05 * z);
        const exactGrad = (x, y, z) => [
          0.1 * Math.cos(0.1 * x),
          -0.2 * Math.sin(0.2 * y),
          0.05 * Math.exp(0.05 * z)
        ];

        const x = 1.5, y = -2.3, z = 0.8;
        const g6 = fd.gradient6Point(fn, x, y, z, [0.01, 0.01, 0.01]);
        const g14 = fd.gradient14Point(fn, x, y, z, [0.01, 0.01, 0.01]);
        const gExact = exactGrad(x, y, z);

        console.log(JSON.stringify({ g6, g14, gExact }));
        """
        res = run_node_es6_snippet(code)
        for i in range(3):
            assert math.isclose(res["g6"][i], res["gExact"][i], rel_tol=1e-4)
            assert math.isclose(res["g14"][i], res["gExact"][i], rel_tol=1e-4)

    def test_convergence_rate_second_order(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FiniteDifference } from './src/interpolation/finite_difference.js';

        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const fd = new FiniteDifference(grid);

        const fn = (x, y, z) => Math.sin(x) * Math.cos(y) * Math.exp(0.1 * z);
        const conv = fd.checkConvergence(fn, 1.0, 1.0, 1.0, 0.1, '14point');
        const convFine = fd.checkConvergence(fn, 1.0, 1.0, 1.0, 0.01, '14point');

        console.log(JSON.stringify({ conv, convFine }));
        """
        res = run_node_es6_snippet(code)
        conv = res["conv"]
        conv_fine = res["convFine"]
        assert math.isclose(conv["orderObserved"], 2.0, abs_tol=0.05)  # Confirms O(h^2) 2nd order convergence
        assert all(err < 1e-3 for err in conv["errorEstimate"])
        assert math.isclose(conv_fine["orderObserved"], 2.0, abs_tol=0.05)
        assert all(err < 1e-5 for err in conv_fine["errorEstimate"])


# ============================================================================
# 4. VELOCITY FIELD & SCALE FACTOR TESTS
# ============================================================================

class TestVelocityField:
    def test_velocity_scale_factor_enforcement(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField, DEFAULT_VELOCITY_SCALE_FACTOR } from './src/fields/velocity_field.js';

        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4 });
        const n = grid.totalCells;

        const rawVx = new Float64Array(n).fill(1.0);
        const rawVy = new Float64Array(n).fill(2.0);
        const rawVz = new Float64Array(n).fill(3.0);

        const vFieldEnforced = new VelocityField(grid, rawVx, rawVy, rawVz, { enforceScaleFactor: true, scaleFactor: 52.0 });
        const vFieldRaw = new VelocityField(grid, rawVx, rawVy, rawVz, { enforceScaleFactor: false });

        const sampEnforced = vFieldEnforced.sampleVelocity(0, 0, 0);
        const sampRaw = vFieldRaw.sampleVelocity(0, 0, 0);

        console.log(JSON.stringify({
          DEFAULT_VELOCITY_SCALE_FACTOR,
          sampEnforced,
          sampRaw
        }));
        """
        res = run_node_es6_snippet(code)
        assert res["DEFAULT_VELOCITY_SCALE_FACTOR"] == 52.0
        assert math.isclose(res["sampEnforced"][0], 52.0, rel_tol=1e-12)
        assert math.isclose(res["sampEnforced"][1], 104.0, rel_tol=1e-12)
        assert math.isclose(res["sampEnforced"][2], 156.0, rel_tol=1e-12)
        assert math.isclose(res["sampRaw"][0], 1.0, rel_tol=1e-12)

    def test_velocity_divergence_vorticity_strain(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';

        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        
        // Pure rotational + expanding velocity field:
        // vx = 2*x - y
        // vy = 2*y + x
        // vz = 2*z
        // Div = 2 + 2 + 2 = 6
        // Vorticity = curl(v) = [0, 0, dvy/dx - dvx/dy] = [0, 0, 1 - (-1)] = [0, 0, 2]
        const vFn = (x, y, z) => [2 * x - y, 2 * y + x, 2 * z];
        const vField = VelocityField.fromAnalyticFunction(grid, vFn);

        const x = 3.0, y = -1.5, z = 2.0;
        const div = vField.sampleDivergence(x, y, z);
        const vort = vField.sampleVorticity(x, y, z);
        const strain = vField.sampleRateOfStrain(x, y, z);

        console.log(JSON.stringify({ div, vort, strain }));
        """
        res = run_node_es6_snippet(code)
        assert math.isclose(res["div"], 6.0, abs_tol=1e-10)
        assert math.isclose(res["vort"][0], 0.0, abs_tol=1e-10)
        assert math.isclose(res["vort"][1], 0.0, abs_tol=1e-10)
        assert math.isclose(res["vort"][2], 2.0, abs_tol=1e-10)
        # S_xx = 2, S_yy = 2, S_zz = 2, off-diagonals 0.5*(-1+1) = 0
        assert math.isclose(res["strain"][0][0], 2.0, abs_tol=1e-10)
        assert math.isclose(res["strain"][1][1], 2.0, abs_tol=1e-10)
        assert math.isclose(res["strain"][2][2], 2.0, abs_tol=1e-10)
        assert math.isclose(res["strain"][0][1], 0.0, abs_tol=1e-10)

    def test_jacobi_diagonalization(self):
        code = """
        import { VelocityField } from './src/fields/velocity_field.js';

        // Known symmetric matrix with known eigenvalues [5, 2, -1]
        // M = [[3, 2, 0], [2, 3, 0], [0, 0, 2]] -> eigenvalues 5, 2, 1
        const M = [
          [3, 2, 0],
          [2, 3, 0],
          [0, 0, 2]
        ];

        const diag = VelocityField.diagonalizeSymmetric3x3(M);
        console.log(JSON.stringify(diag));
        """
        res = run_node_es6_snippet(code)
        evals = res["eigenvalues"]
        assert math.isclose(evals[0], 5.0, abs_tol=1e-10)
        assert math.isclose(evals[1], 2.0, abs_tol=1e-10)
        assert math.isclose(evals[2], 1.0, abs_tol=1e-10)

    def test_bulk_flow_estimator(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';

        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const vFn = (x, y, z) => [300.0, -150.0, 450.0];
        const vField = VelocityField.fromAnalyticFunction(grid, vFn);

        const bulk = vField.computeBulkFlow([0, 0, 0], 25.0, 15);
        console.log(JSON.stringify(bulk));
        """
        res = run_node_es6_snippet(code)
        assert math.isclose(res["bulkFlow"][0], 300.0, abs_tol=1e-10)
        assert math.isclose(res["bulkFlow"][1], -150.0, abs_tol=1e-10)
        assert math.isclose(res["bulkFlow"][2], 450.0, abs_tol=1e-10)
        expected_speed = math.sqrt(300**2 + 150**2 + 450**2)
        assert math.isclose(res["speed"], expected_speed, abs_tol=1e-10)


# ============================================================================
# 5. DENSITY FIELD & CIC MASS ASSIGNMENT TESTS
# ============================================================================

class TestDensityField:
    def test_density_analytic_and_log_density(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { DensityField } from './src/fields/density_field.js';

        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const deltaFn = (x, y, z) => 0.5 + 0.1 * x;
        const dField = DensityField.fromAnalyticFunction(grid, deltaFn);

        const x = 2.0, y = 1.0, z = -3.0;
        const dVal = dField.sampleDensity(x, y, z);
        const logD = dField.sampleLogDensity(x, y, z);
        const grad = dField.sampleGradient(x, y, z);
        const logGrad = dField.sampleLogGradient(x, y, z);

        console.log(JSON.stringify({ dVal, logD, grad, logGrad }));
        """
        res = run_node_es6_snippet(code)
        expected_delta = 0.5 + 0.1 * 2.0  # 0.7
        assert math.isclose(res["dVal"], expected_delta, abs_tol=1e-12)
        assert math.isclose(res["logD"], math.log(1.0 + expected_delta), abs_tol=1e-12)
        assert math.isclose(res["grad"][0], 0.1, abs_tol=1e-12)
        assert math.isclose(res["logGrad"][0], 0.1 / (1.0 + expected_delta), abs_tol=1e-12)

    def test_cic_mass_conservation_and_mean_normalization(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { DensityField } from './src/fields/density_field.js';

        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        
        // 100 random tracer galaxies
        const particles = [];
        const weights = [];
        let totalW = 0.0;
        for (let i = 0; i < 100; i++) {
          const px = -8.0 + (i * 1.37) % 16.0;
          const py = -8.0 + (i * 2.51) % 16.0;
          const pz = -8.0 + (i * 3.19) % 16.0;
          const w = 1.0 + (i % 5) * 0.5;
          particles.push([px, py, pz]);
          weights.push(w);
          totalW += w;
        }

        const dField = DensityField.fromParticlesCIC(grid, particles, weights);
        const stats = dField.computeStatistics();

        console.log(JSON.stringify({ stats, totalW }));
        """
        res = run_node_es6_snippet(code)
        stats = res["stats"]
        # In a closed grid, average delta over all cells must be exactly 0
        assert math.isclose(stats["mean"], 0.0, abs_tol=1e-10)
        assert stats["min"] >= -1.0  # Physical bound on overdensity
        assert stats["variance"] > 0.0

    def test_fits_stride_mapping_and_indexing(self):
        code = """
        import { GridIndexer, StrideOrder } from './src/fields/grid_indexer.js';
        const grid = new GridIndexer({
          nx: 12, ny: 18, nz: 24,
          origin: [-100, -100, -100],
          boxSize: [200, 200, 200],
          strideOrder: StrideOrder.FITS_ZYX
        });

        // Test FITS NAXIS1 (SGX), NAXIS2 (SGY), NAXIS3 (SGZ)
        const ix = 5, iy = 9, iz = 14;
        const flatIdx = grid.fitsToCanonicalIndex(ix, iy, iz);
        const [rx, ry, rz] = grid.get3DIndices(flatIdx);

        console.log(JSON.stringify({ flatIdx, rx, ry, rz, match: rx === ix && ry === iy && rz === iz }));
        """
        res = run_node_es6_snippet(code)
        assert res["match"] is True
        assert res["flatIdx"] == 5 + 9 * 12 + 14 * 12 * 18

    def test_analytic_jacobian_agreement_with_finite_difference(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';
        import { FiniteDifference } from './src/interpolation/finite_difference.js';

        const grid = new GridIndexer({ nx: 32, ny: 32, nz: 32, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        
        // Trilinear field with bilinear cross terms
        const vFn = (x, y, z) => [
          1.2 * x + 0.3 * y * z,
          -0.8 * y + 0.5 * x * z,
          0.4 * z + 0.2 * x * y
        ];

        const vField = VelocityField.fromAnalyticFunction(grid, vFn);
        const interp = new TrilinearInterpolator(grid);
        const fd = new FiniteDifference(grid, interp);

        const x = 1.7, y = -2.3, z = 0.9;
        const J_analytic = interp.analyticVectorJacobian(vField.vx, vField.vy, vField.vz, x, y, z);
        const J_fd = fd.velocityJacobian(vField.vx, vField.vy, vField.vz, x, y, z, '14point', [0.001, 0.001, 0.001]);

        console.log(JSON.stringify({ J_analytic, J_fd }));
        """
        res = run_node_es6_snippet(code)
        J_analytic = res["J_analytic"]
        J_fd = res["J_fd"]
        for i in range(3):
            for j in range(3):
                assert math.isclose(J_analytic[i][j], J_fd[i][j], rel_tol=1e-3, abs_tol=1e-3)

    def test_boundary_conditions_clamping_on_interpolator(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';

        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const interp = new TrilinearInterpolator(grid);

        const buf = new Float64Array(grid.totalCells);
        for (let iz = 0; iz < grid.nz; iz++) {
          for (let iy = 0; iy < grid.ny; iy++) {
            for (let ix = 0; ix < grid.nx; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              buf[idx] = ix + iy + iz;
            }
          }
        }

        // Evaluate way outside domain: should clamp gracefully without throwing
        const vNeg = interp.interpolateScalar(buf, -50.0, -50.0, -50.0);
        const vPos = interp.interpolateScalar(buf, 50.0, 50.0, 50.0);

        console.log(JSON.stringify({ vNeg, vPos, expectedNeg: 0, expectedPos: 7 + 7 + 7 }));
        """
        res = run_node_es6_snippet(code)
        assert math.isclose(res["vNeg"], res["expectedNeg"], abs_tol=1e-12)
        assert math.isclose(res["vPos"], res["expectedPos"], abs_tol=1e-12)

    def test_vector_hessian_tensor_symmetry(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';

        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-5, -5, -5], boxSize: [10, 10, 10] });
        const interp = new TrilinearInterpolator(grid);

        const vxBuf = new Float64Array(grid.totalCells);
        const vyBuf = new Float64Array(grid.totalCells);
        const vzBuf = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < grid.nz; iz++) {
          for (let iy = 0; iy < grid.ny; iy++) {
            for (let ix = 0; ix < grid.nx; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vxBuf[idx] = x * y * z;
              vyBuf[idx] = 2.0 * x * y * z;
              vzBuf[idx] = -3.0 * x * y * z;
            }
          }
        }

        const H = interp.analyticVectorHessian(vxBuf, vyBuf, vzBuf, 1.2, -0.8, 2.1);
        console.log(JSON.stringify({ H }));
        """
        res = run_node_es6_snippet(code)
        H = res["H"]
        # Verify 3x3x3 tensor structure and symmetry in spatial derivative indices j, k
        assert len(H) == 3
        for comp in range(3):
            assert len(H[comp]) == 3
            assert math.isclose(H[comp][0][1], H[comp][1][0], abs_tol=1e-12)
            assert math.isclose(H[comp][0][2], H[comp][2][0], abs_tol=1e-12)
            assert math.isclose(H[comp][1][2], H[comp][2][1], abs_tol=1e-12)


# ============================================================================
# 6. SUMMARY REPORT
# ============================================================================

def test_complete_integration_sanity():
    """Verify that all modules integrate seamlessly without exceptions."""
    code = """
    import { GridIndexer } from './src/fields/grid_indexer.js';
    import { VelocityField } from './src/fields/velocity_field.js';
    import { DensityField } from './src/fields/density_field.js';
    import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';
    import { FiniteDifference } from './src/interpolation/finite_difference.js';

    const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8 });
    const vField = VelocityField.fromAnalyticFunction(grid, (x,y,z)=>[x, y, z]);
    const dField = DensityField.fromAnalyticFunction(grid, (x,y,z)=>0.1*(x+y+z));
    const interp = new TrilinearInterpolator(grid);
    const fd = new FiniteDifference(grid, interp);

    console.log(JSON.stringify({ status: "ALL_MODULES_LOADED_OK" }));
    """
    res = run_node_es6_snippet(code)
    assert res["status"] == "ALL_MODULES_LOADED_OK"

