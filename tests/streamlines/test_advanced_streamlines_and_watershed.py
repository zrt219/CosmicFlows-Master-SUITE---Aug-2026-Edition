import os
import pytest
import math
from tests.utils import run_node_snippet

class TestMassiveStreamlinesAndWatershed:
    def test_massive_streamline_streamer_chunking(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { MassiveStreamlineStreamer } from './src/streamlines/massive_streamline_streamer.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [-0.05 * x, -0.05 * y, -0.05 * z]);
        
        const streamer = new MassiveStreamlineStreamer(field, { chunkSize: 64, maxSteps: 50 });
        const seeds = streamer.generateGridSeeds(2); // 4x4x4 = 64 seeds
        
        let chunkCalls = 0;
        const result = streamer.streamBatch(seeds, () => { chunkCalls++; });
        
        console.log(JSON.stringify({
          totalStreamlines: result.totalStreamlines,
          chunkCalls,
          hasSignature: result.signature.length === 64
        }));
        """
        res = run_node_snippet(code)
        assert res["totalStreamlines"] == 64
        assert res["chunkCalls"] == 1
        assert res["hasSignature"] is True

    def test_published_watershed_reader_and_geometry(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PublishedWatershedReader } from './src/watershed/published_watershed_reader.js';
        import { BasinGeometryEngine } from './src/watershed/basin_geometry_engine.js';
        
        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
        const total = grid.totalCells;
        
        // Synthetic segmentation: half Laniakea (1), half Shapley (6)
        const labelBuf = new Uint8Array(total);
        for (let i = 0; i < total; i++) {
          labelBuf[i] = (i < total / 2) ? 1 : 6;
        }
        
        const reader = new PublishedWatershedReader({ grid });
        const parsed = reader.parseLabelVolume(labelBuf);
        
        const metrics = BasinGeometryEngine.computeBasinMetrics(grid, labelBuf, 1);
        
        console.log(JSON.stringify({
          uniqueBasins: parsed.uniqueBasins,
          laniakeaPct: parsed.basins.find(b => b.id === 1).volumePercentage,
          shapleyPct: parsed.basins.find(b => b.id === 6).volumePercentage,
          metricsVoxelCount: metrics.voxelCount,
          equivalentRadius: metrics.equivalentRadius
        }));
        """
        res = run_node_snippet(code)
        assert res["uniqueBasins"] == 2
        assert res["laniakeaPct"] == 50.0
        assert res["shapleyPct"] == 50.0
        assert res["metricsVoxelCount"] == 2048
        assert res["equivalentRadius"] > 0.0

class TestAdvancedRootsAndReferenceSolvers:
    def test_advanced_root_solver_linear_sink(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { AdvancedRootSolver } from './src/topology/advanced_root_solver.js';
        
        const N = 16;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const field = VelocityField.fromAnalyticFunction(grid, (x, y, z) => [-0.1 * (x - 5.0), -0.1 * (y + 3.0), -0.1 * (z - 2.0)]);
        
        const solver = new AdvancedRootSolver(field);
        const res = solver.solveFromSeed(1.0, 1.0, 1.0);
        
        console.log(JSON.stringify({
          converged: res.converged,
          pos: res.position,
          residual: res.residual
        }));
        """
        res = run_node_snippet(code)
        assert res["converged"] is True
        assert abs(res["pos"][0] - 5.0) < 1e-4
        assert abs(res["pos"][1] - (-3.0)) < 1e-4
        assert abs(res["pos"][2] - 2.0) < 1e-4
        assert res["residual"] < 1e-6

    def test_reference_trilinear_interpolator_identity(self):
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TrilinearInterpolator } from './src/interpolation/trilinear_interpolator.js';
        import { ReferenceTrilinearInterpolator } from './src/reference/independent_reference_solvers.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const buf = new Float32Array(N * N * N);
        for (let i = 0; i < buf.length; i++) buf[i] = i * 1.5;
        
        const interp = new TrilinearInterpolator(grid);
        const valOpt = interp.interpolateScalar(buf, 2.5, -3.5, 1.5);
        const valRef = ReferenceTrilinearInterpolator.eval(grid, buf, 2.5, -3.5, 1.5);
        
        console.log(JSON.stringify({
          valOpt,
          valRef,
          diff: Math.abs(valOpt - valRef)
        }));
        """
        res = run_node_snippet(code)
        assert res["diff"] < 1e-5
