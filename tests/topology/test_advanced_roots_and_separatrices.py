# -*- coding: utf-8 -*-
"""
Automated Pytest Suite for Advanced Vector Roots, Jacobian Eigensystems, and Manifold Separatrix Tracking.
"""

import pytest
import math
import numpy as np
from tests.utils import run_node_snippet

class TestAdvancedRootsAndSeparatrices:
    def test_separatrix_manifold_tracer_instantiation(self):
        res = run_node_snippet("""
        import { SeparatrixManifoldTracer } from './src/topology/manifold_separatrix_tracer.js';
        const mockField = {
            evaluate: (pos, outV) => {
                outV[0] = -pos[0];
                outV[1] = -pos[1];
                outV[2] = 2.0 * pos[2];
            }
        };
        const tracer = new SeparatrixManifoldTracer(mockField);
        console.log(JSON.stringify({
            hasTracer: tracer !== null && typeof tracer.traceSaddleManifolds === 'function',
            stepSize: tracer.stepSize,
            tolerance: tracer.tolerance
        }));
        """)
        assert res["hasTracer"] is True
        assert res["stepSize"] > 0
        assert res["tolerance"] > 0

    def test_saddle_filament_manifold_tracing(self):
        res = run_node_snippet("""
        import { SeparatrixManifoldTracer } from './src/topology/manifold_separatrix_tracer.js';
        const mockField = {
            evaluate: (pos, outV) => {
                outV[0] = -pos[0];
                outV[1] = -pos[1];
                outV[2] = 2.0 * pos[2];
            }
        };
        const tracer = new SeparatrixManifoldTracer(mockField, { stepSize: 0.2, maxSteps: 100 });
        const saddle = {
            id: 'saddle_origin',
            type: 'SADDLE_FILAMENT',
            position: [0.0, 0.0, 0.0],
            eigenvalues: [-1.0, -1.0, 2.0],
            eigenvectors: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]
            ]
        };
        const attractor = {
            id: 'attractor_north',
            type: 'NODE_SINK',
            position: [0.0, 0.0, 10.0]
        };
        const manifolds = tracer.traceSaddleManifolds(saddle, [attractor]);
        const posPts = manifolds.branchPositive.waypoints;
        const negPts = manifolds.branchNegative.waypoints;
        console.log(JSON.stringify({
            isFilament: manifolds.isFilamentManifold,
            posBranchLen: posPts.length,
            negBranchLen: negPts.length,
            posZ: posPts[posPts.length - 1][2],
            negZ: negPts[negPts.length - 1][2]
        }));
        """)
        assert res["isFilament"] is True
        assert res["posBranchLen"] > 5
        assert res["negBranchLen"] > 5
        assert res["posZ"] > 0.0
        assert res["negZ"] < 0.0


    def test_topology_graph_construction(self):
        res = run_node_snippet("""
        import { SeparatrixManifoldTracer } from './src/topology/manifold_separatrix_tracer.js';
        const mockField = {
            evaluate: (pos, outV) => {
                outV[0] = -pos[0];
                outV[1] = -pos[1];
                outV[2] = pos[2] > 0 ? (10.0 - pos[2]) : (-10.0 - pos[2]);
            }
        };
        const tracer = new SeparatrixManifoldTracer(mockField, { stepSize: 0.5, captureRadiusMpc: 2.0 });
        const criticalPoints = [
            {
                id: 'sink_1',
                type: 'NODE_SINK',
                position: [0.0, 0.0, 10.0],
                eigenvalues: [-2.0, -2.0, -1.0]
            },
            {
                id: 'sink_2',
                type: 'NODE_SINK',
                position: [0.0, 0.0, -10.0],
                eigenvalues: [-2.0, -2.0, -1.0]
            },
            {
                id: 'saddle_mid',
                type: 'SADDLE_FILAMENT',
                position: [0.0, 0.0, 0.0],
                eigenvalues: [-1.0, -1.0, 2.0],
                eigenvectors: [[1,0,0], [0,1,0], [0,0,1]]
            }
        ];
        const graph = tracer.buildTopologyGraph(criticalPoints);
        console.log(JSON.stringify({
            nodeCount: graph.nodeCount,
            edgeCount: graph.edgeCount,
            filamentCount: graph.filamentCount,
            hasMetadata: !!graph.metadata.citation
        }));
        """)
        assert res["nodeCount"] == 3
        assert res["edgeCount"] >= 2
        assert res["filamentCount"] >= 2
        assert res["hasMetadata"] is True

    def test_jacobian_dynamical_classification(self):
        res = run_node_snippet("""
        import { JacobianEigensystemSolver, CriticalPointDynamicalType } from './src/topology/jacobian_eigensystem.js';
        
        const J_sink = [-2, 0, 0, 0, -3, 0, 0, 0, -1];
        const resSink = JacobianEigensystemSolver.analyzeJacobian(J_sink);
        
        const J_source = [2, 0, 0, 0, 3, 0, 0, 0, 1];
        const resSource = JacobianEigensystemSolver.analyzeJacobian(J_source);
        
        const J_filament = [2, 0, 0, 0, -3, 0, 0, 0, -1];
        const resFilament = JacobianEigensystemSolver.analyzeJacobian(J_filament);
        
        const J_wall = [2, 0, 0, 0, 3, 0, 0, 0, -1];
        const resWall = JacobianEigensystemSolver.analyzeJacobian(J_wall);

        console.log(JSON.stringify({
            sinkType: resSink.dynamicalType,
            sourceType: resSource.dynamicalType,
            filamentType: resFilament.dynamicalType,
            wallType: resWall.dynamicalType,
            sinkDiv: resSink.divergence,
            sourceDiv: resSource.divergence
        }));
        """)
        assert res["sinkType"] == "NODE_SINK"
        assert res["sourceType"] == "NODE_SOURCE"
        assert res["filamentType"] == "SADDLE_FILAMENT"
        assert res["wallType"] == "SADDLE_WALL"
        assert res["sinkDiv"] < 0
        assert res["sourceDiv"] > 0

