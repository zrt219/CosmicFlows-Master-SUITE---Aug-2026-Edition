# -*- coding: utf-8 -*-
"""
tests/topology/test_advanced_roots_and_separatrices.py
Automated test suite verifying 3D Jacobian eigensystems, multi-scale watersheds,
event detectors, parameter sweeps, and adversarial guards.
"""

import pytest
from tests.utils import run_node_snippet

class TestJacobianEigensystemAndSweeps:
    def test_jacobian_eigensystem_attractor_classification(self):
        """Verifies true 3D sink/attractor has 3 negative eigenvalues and classifies as NODE_SINK."""
        code = """
        import { JacobianEigensystemSolver, CriticalPointDynamicalType } from './src/topology/jacobian_eigensystem.js';
        
        // J = diag(-0.5, -0.4, -0.6)
        const J = new Float64Array([
          -0.5, 0.0, 0.0,
          0.0, -0.4, 0.0,
          0.0, 0.0, -0.6
        ]);
        
        const analysis = JacobianEigensystemSolver.analyzeJacobian(J);
        
        console.log(JSON.stringify({
          dynType: analysis.dynType,
          divergence: analysis.divergence,
          isAttractor: analysis.dynamicalType === CriticalPointDynamicalType.NODE_SINK
        }));
        """
        res = run_node_snippet(code)
        assert res["isAttractor"] is True
        assert res["divergence"] == pytest.approx(-1.5)

    def test_streamline_event_detector_boundary_and_stall(self):
        """Verifies event detector flags domain exits and low speed stalls."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { StreamlineEventDetector, StreamlineEventType } from './src/streamlines/streamline_event_detector.js';
        
        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const detector = new StreamlineEventDetector(grid, { vMin: 5.0, maxArcLength: 200.0 });
        
        // 1. Boundary exit test
        const evtExit = detector.checkStep([60.0, 0.0, 0.0], [10.0, 0.0, 0.0], 10.0, []);
        
        // 2. Low speed stall test
        const evtStall = detector.checkStep([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 10.0, []);
        
        console.log(JSON.stringify({
          exitType: evtExit.type,
          stallType: evtStall.type
        }));
        """
        res = run_node_snippet(code)
        assert res["exitType"] == "DOMAIN_EXIT"
        assert res["stallType"] == "LOW_SPEED_STALL"

    def test_multi_scale_watershed_peak_detection(self):
        """Verifies 26-neighbor peak detection finds genuine endpoint concentration maxima."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField } from './src/fields/velocity_field.js';
        import { MultiScaleWatershed } from './src/watershed/multi_scale_watershed.js';
        
        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, origin: [-20, -20, -20], boxSize: [40, 40, 40] });
        const total = grid.totalCells;
        const density = new Float32Array(total);
        
        // Place isolated peak at center (ix=4, iy=4, iz=4)
        const centerIdx = grid.getLinearIndex(4, 4, 4);
        density[centerIdx] = 50.0;
        
        const ws = new MultiScaleWatershed({ grid }, { minPeakDensity: 10 });
        const peaks = ws.findLocalMaxima(density);
        
        console.log(JSON.stringify({
          peakCount: peaks.length,
          peakIdx: peaks[0] ? peaks[0].index : null,
          peakVal: peaks[0] ? peaks[0].peakValue : null
        }));
        """
        res = run_node_snippet(code)
        assert res["peakCount"] == 1
        assert res["peakVal"] == 50.0

    def test_adversarial_science_guard_boundary_exclusion(self):
        """Verifies adversarial guard rejects critical points within boundary layer."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { AdversarialScienceGuard } from './src/validation/adversarial_science_guard.js';
        
        const grid = new GridIndexer({ nx: 10, ny: 10, nz: 10, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const guard = new AdversarialScienceGuard(grid);
        
        const testPoints = [
          { position: [0.0, 0.0, 0.0], residual: 1e-8, type: 'Attractor' }, // Valid interior root
          { position: [49.0, 0.0, 0.0], residual: 1e-8, type: 'Attractor' }, // Boundary artifact
          { position: [0.0, 0.0, 0.0], residual: 0.1, type: 'Attractor' }    // High residual pseudo-root
        ];
        
        const audit = guard.auditCriticalPoints(testPoints);
        
        console.log(JSON.stringify({
          validCount: audit.validRoots.length,
          rejectedCount: audit.rejectedRoots.length
        }));
        """
        res = run_node_snippet(code)
        assert res["validCount"] == 1
        assert res["rejectedCount"] == 2

    def test_publication_tooling_latex_and_attribution(self):
        """Verifies LaTeX table generation and mandatory publication attributions."""
        code = """
        import { PublicationTooling } from './src/export/publication_tooling.js';
        
        const basins = [
          { id: 1, name: 'Laniakea', volumePercentage: 12.5, equivalentRadius: 82.3, centroid: [-12.0, 4.0, -8.0] },
          { id: 6, name: 'Shapley', volumePercentage: 24.1, equivalentRadius: 110.5, centroid: [-140.0, 80.0, -20.0] }
        ];
        
        const latex = PublicationTooling.generateBasinTableLaTeX(basins);
        const method = PublicationTooling.getStandardMethodologyText();
        
        console.log(JSON.stringify({
          hasLaniakea: latex.includes('Laniakea'),
          hasCourtois2023: method.includes('Courtois et al. 2023') && method.includes('A&A 670, L15'),
          hasDupuy2023: method.includes('Dupuy & Courtois (2023') && method.includes('A&A 678, A176')
        }));
        """
        res = run_node_snippet(code)
        assert res["hasLaniakea"] is True
        assert res["hasCourtois2023"] is True
        assert res["hasDupuy2023"] is True
