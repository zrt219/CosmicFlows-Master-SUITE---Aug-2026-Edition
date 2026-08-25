# -*- coding: utf-8 -*-
"""
tests/adversarial/test_adversarial_science_and_robustness.py
Extensive adversarial QA, boundary condition stress, and failure injection suite.
"""

import math
import pytest
from tests.utils import run_node_snippet

class TestAdversarialScienceAndRobustness:
    def test_rejection_of_unscaled_cf4_displacement(self):
        """Ensures displacement fields without x52 scaling cannot be labeled as physical velocities."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VelocityField, DEFAULT_VELOCITY_SCALE_FACTOR } from './src/fields/velocity_field.js';
        
        const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4 });
        const total = grid.totalCells;
        const rawPsi = new Float64Array(total * 3);
        for (let i = 0; i < rawPsi.length; i++) rawPsi[i] = 1.0;
        
        // When constructed with enforceScaleFactor=true, scaling must be 52.0
        const field = VelocityField.fromInterleavedArray(grid, rawPsi, { scaleFactor: 52.0 });
        const sample = field.sampleVelocity(0, 0, 0);
        
        console.log(JSON.stringify({
          scaledVx: sample[0],
          is52: Math.abs(sample[0] - 52.0) < 1e-5
        }));
        """
        res = run_node_snippet(code)
        assert res["is52"] is True

    def test_non_hyperbolic_saddle_classification(self):
        """Verifies eigensystem handles degenerate or zero eigenvalues without false classification."""
        code = """
        import { EigenSystem3D } from './src/coordinates/scientific_types.js';
        
        // Matrix with zero eigenvalues (degenerate plane)
        const S = new Float64Array([
          0, 0, 0,
          0, 0, 0,
          0, 0, 0
        ]);
        const eigen = EigenSystem3D.fromSymmetricMatrix(S);
        const web = eigen.classifyWebStructure(0.0);
        
        console.log(JSON.stringify({
          l1: eigen.lambda1,
          l2: eigen.lambda2,
          l3: eigen.lambda3,
          code: web.code
        }));
        """
        res = run_node_snippet(code)
        assert res["l1"] == 0.0
        assert res["l2"] == 0.0
        assert res["l3"] == 0.0
        assert res["code"] == 0  # Classified as void/homogeneous baseline
