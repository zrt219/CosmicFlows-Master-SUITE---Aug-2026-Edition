# -*- coding: utf-8 -*-
"""
tests/uncertainty/test_cf4pp_and_multipolar_flows.py
Automated test suite verifying multipolar bulk flows, Monte Carlo bootstrap,
and CF4++ posterior sampling.
"""

import pytest
from tests.utils import run_node_snippet

class TestCF4ppAndMultipolarFlows:
    def test_multipolar_bulk_flow_dipole_recovery(self):
        """Verifies dipole estimator recovers known constant bulk flow vector."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { MultipolarBulkFlow } from './src/bulk-flow/multipolar_bulk_flow.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-20, -20, -20], boxSize: [40, 40, 40] });
        const total = grid.totalCells;
        
        // Constant flow: vx = 300, vy = 400, vz = 0 (|v| = 500 km/s)
        const vx = new Float32Array(total).fill(300 / 52.0);
        const vy = new Float32Array(total).fill(400 / 52.0);
        const vz = new Float32Array(total).fill(0.0);
        
        const vField = new VelocityField(grid, vx, vy, vz, { scaleApplied: false });
        const mbf = new MultipolarBulkFlow(vField);
        const res = mbf.computeDipole(20.0);
        
        console.log(JSON.stringify({
          magnitude: res.magnitude,
          SGL: res.SGL,
          SGB: res.SGB
        }));
        """
        res = run_node_snippet(code)
        assert res["magnitude"] == pytest.approx(500.0, rel=1e-3)
        assert res["SGL"] == pytest.approx(53.13, abs=0.1)
        assert res["SGB"] == pytest.approx(0.0, abs=0.1)

    def test_monte_carlo_bootstrap_mean_and_ci(self):
        """Verifies statistical bootstrap calculates empirical 95% confidence intervals."""
        code = """
        import { StatisticalResampling } from './src/statistics/statistical_resampling.js';
        
        // Gaussian distributed sample centered at 100 with std ~ 10
        const data = [90, 95, 100, 105, 110, 98, 102, 104, 96, 101];
        const meanFn = (arr) => arr.reduce((a, b) => a + b, 0) / arr.length;
        
        const bs = StatisticalResampling.bootstrap(data, meanFn, 500, 12345);
        
        console.log(JSON.stringify({
          mean: bs.mean,
          std: bs.std,
          ciLow: bs.ci95[0],
          ciHigh: bs.ci95[1]
        }));
        """
        res = run_node_snippet(code)
        assert res["mean"] == pytest.approx(100.1, abs=0.5)
        assert res["ciLow"] < res["mean"] < res["ciHigh"]

    def test_cf4pp_posterior_sampler_dispersion(self):
        """Verifies CF4++ posterior sampler generates Gaussian perturbations around mean."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { CF4PosteriorSampler } from './src/uncertainty/cf4pp_posterior_sampler.js';
        import { MT19937 } from './src/statistics/statistical_resampling.js';
        
        const N = 4;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const total = grid.totalCells;
        
        const vMean = new VelocityField(grid, new Float32Array(total).fill(100 / 52.0), new Float32Array(total), new Float32Array(total), { scaleApplied: false });
        const vRms = new VelocityField(grid, new Float32Array(total).fill(20 / 52.0), new Float32Array(total).fill(20 / 52.0), new Float32Array(total).fill(20 / 52.0), { scaleApplied: false });
        
        const sampler = new CF4PosteriorSampler(vMean, vRms);
        const prng = new MT19937(42);
        
        const samples = [];
        for (let i = 0; i < 1000; i++) {
          const s = sampler.sampleVelocity(0, 0, 0, prng);
          samples.push(s[0]);
        }
        
        const sampleMean = samples.reduce((a, b) => a + b, 0) / 1000;
        let sumSq = 0;
        for (let s of samples) sumSq += (s - sampleMean) ** 2;
        const sampleStd = Math.sqrt(sumSq / 999);
        
        console.log(JSON.stringify({
          mean: sampleMean,
          std: sampleStd
        }));
        """
        res = run_node_snippet(code)
        assert res["mean"] == pytest.approx(100.0, abs=2.0)
        assert res["std"] == pytest.approx(20.0, abs=2.0)
