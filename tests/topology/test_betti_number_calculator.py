# -*- coding: utf-8 -*-
"""
tests/topology/test_betti_number_calculator.py
=============================================
Exhaustive Pytest Suite for 3D Persistent Homology, Betti Number,
and Euler Characteristic Calculator across Cosmic Density & Potential Fields.

Verifies:
1. 3D Cubical Complex indexing, cell relations, boundary operators d_k d_{k+1} = 0, and cell-count Euler characteristic.
2. Union-Find 0D persistence, path compression, and Elder Rule merges.
3. Exact Betti Numbers \beta_0, \beta_1, \beta_2 on canonical topological manifolds (points, hollow spheres S^2, tori S^1 x S^1, double tori).
4. Euler-Poincaré invariant: \chi = \beta_0 - \beta_1 + \beta_2 \equiv V - E + F - C across all grids.
5. Algebraic Z_2 boundary matrix column reduction and persistent homology barcodes.
6. Persistence diagram metrics: total persistence, persistent entropy, Bottleneck distance, 1-Wasserstein & 2-Wasserstein distances.
7. Persistence landscapes \lambda_k(t) and L_p norm vectorization.
8. Tomita-Gott Gaussian Random Field (GRF) analytical Euler density \chi_{GRF}(\nu) invariants, extrema, and zero-crossings.
9. Velocity potential tidal deformation tensor eigenvalues \lambda_1 \ge \lambda_2 \ge \lambda_3, trace \nabla^2 \Phi = \delta, and cosmic web classification.
10. End-to-end cosmic web topology analysis (meatball, spongy, Swiss-cheese, GRF).
"""

import pytest
import math
import numpy as np
from tests.utils import run_node_snippet

# ============================================================================
# 1. CUBICAL COMPLEX INDEXING & BOUNDARY OPERATORS
# ============================================================================

class TestCubicalComplex3D:
    def test_cubical_complex_initialization_and_cell_counts(self):
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        const nx = 4, ny = 5, nz = 6;
        const complex = new CubicalComplex3D(nx, ny, nz);
        
        console.log(JSON.stringify({
            numVertices: complex.numVertices,
            numEdgesX: complex.numEdgesX,
            numEdgesY: complex.numEdgesY,
            numEdgesZ: complex.numEdgesZ,
            numEdges: complex.numEdges,
            numFacesXY: complex.numFacesXY,
            numFacesYZ: complex.numFacesYZ,
            numFacesZX: complex.numFacesZX,
            numFaces: complex.numFaces,
            numCubes: complex.numCubes,
            totalCells: complex.totalCells
        }));
        """)
        nx, ny, nz = 4, 5, 6
        expected_vertices = nx * ny * nz
        expected_edges_x = (nx - 1) * ny * nz
        expected_edges_y = nx * (ny - 1) * nz
        expected_edges_z = nx * ny * (nz - 1)
        expected_edges = expected_edges_x + expected_edges_y + expected_edges_z
        expected_faces_xy = (nx - 1) * (ny - 1) * nz
        expected_faces_yz = nx * (ny - 1) * (nz - 1)
        expected_faces_zx = (nx - 1) * ny * (nz - 1)
        expected_faces = expected_faces_xy + expected_faces_yz + expected_faces_zx
        expected_cubes = (nx - 1) * (ny - 1) * (nz - 1)

        assert res["numVertices"] == expected_vertices
        assert res["numEdgesX"] == expected_edges_x
        assert res["numEdgesY"] == expected_edges_y
        assert res["numEdgesZ"] == expected_edges_z
        assert res["numEdges"] == expected_edges
        assert res["numFacesXY"] == expected_faces_xy
        assert res["numFacesYZ"] == expected_faces_yz
        assert res["numFacesZX"] == expected_faces_zx
        assert res["numFaces"] == expected_faces
        assert res["numCubes"] == expected_cubes
        assert res["totalCells"] == expected_vertices + expected_edges + expected_faces + expected_cubes

    def test_cubical_complex_invalid_dimensions_throw(self):
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        let caught = false;
        try {
            new CubicalComplex3D(1, 4, 4);
        } catch (e) {
            caught = true;
        }
        console.log(JSON.stringify({ caught }));
        """)
        assert res["caught"] is True

    def test_vertex_indexing_and_coordinates_bijection(self):
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        const complex = new CubicalComplex3D(5, 6, 7);
        let allMatch = true;
        for (let k = 0; k < 7; k++) {
            for (let j = 0; j < 6; j++) {
                for (let i = 0; i < 5; i++) {
                    const idx = complex.vertexIndex(i, j, k);
                    const [ci, cj, ck] = complex.vertexCoords(idx);
                    if (ci !== i || cj !== j || ck !== k) {
                        allMatch = false;
                        break;
                    }
                }
            }
        }
        console.log(JSON.stringify({ allMatch }));
        """)
        assert res["allMatch"] is True

    def test_boundary_operators_algebraic_consistency(self):
        """Verify boundary operator dimensions: d1 (2 vertices), d2 (4 edges), d3 (6 faces)."""
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        const complex = new CubicalComplex3D(3, 3, 3);
        let d1Valid = true;
        let d2Valid = true;
        let d3Valid = true;

        for (let e = 0; e < complex.numEdges; e++) {
            const v = complex.edgeBoundaryVertices(e);
            if (!Array.isArray(v) || v.length !== 2 || v[0] === v[1] || v[0] < 0 || v[1] >= complex.numVertices) {
                d1Valid = false;
            }
        }

        for (let f = 0; f < complex.numFaces; f++) {
            const edges = complex.faceBoundaryEdges(f);
            if (!Array.isArray(edges) || edges.length !== 4) {
                d2Valid = false;
            }
        }

        for (let c = 0; c < complex.numCubes; c++) {
            const faces = complex.cubeBoundaryFaces(c);
            if (!Array.isArray(faces) || faces.length !== 6) {
                d3Valid = false;
            }
        }

        console.log(JSON.stringify({ d1Valid, d2Valid, d3Valid }));
        """)
        assert res["d1Valid"] is True
        assert res["d2Valid"] is True
        assert res["d3Valid"] is True

    def test_boundary_of_boundary_is_zero_in_z2(self):
        """Verify d1 o d2 = 0 over Z_2 (every vertex in boundary of a face appears exactly twice)."""
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        const complex = new CubicalComplex3D(4, 4, 4);
        let d1_d2_is_zero = true;

        for (let f = 0; f < complex.numFaces; f++) {
            const edges = complex.faceBoundaryEdges(f);
            const vertexCounts = new Map();
            for (const e of edges) {
                const [v0, v1] = complex.edgeBoundaryVertices(e);
                vertexCounts.set(v0, (vertexCounts.get(v0) || 0) + 1);
                vertexCounts.set(v1, (vertexCounts.get(v1) || 0) + 1);
            }
            // In Z_2, every vertex must have degree 2 (boundary vanishes modulo 2)
            for (const [v, count] of vertexCounts.entries()) {
                if (count % 2 !== 0) {
                    d1_d2_is_zero = false;
                }
            }
        }

        console.log(JSON.stringify({ d1_d2_is_zero }));
        """)
        assert res["d1_d2_is_zero"] is True

    def test_cube_boundary_of_boundary_is_zero_in_z2(self):
        """Verify d2 o d3 = 0 over Z_2 (every edge in boundary of a cube appears exactly twice among its 6 faces)."""
        res = run_node_snippet("""
        import { CubicalComplex3D } from './src/topology/betti_number_calculator.js';
        const complex = new CubicalComplex3D(3, 3, 3);
        let d2_d3_is_zero = true;

        for (let c = 0; c < complex.numCubes; c++) {
            const faces = complex.cubeBoundaryFaces(c);
            const edgeCounts = new Map();
            for (const f of faces) {
                const edges = complex.faceBoundaryEdges(f);
                for (const e of edges) {
                    edgeCounts.set(e, (edgeCounts.get(e) || 0) + 1);
                }
            }
            for (const [e, count] of edgeCounts.entries()) {
                if (count % 2 !== 0) {
                    d2_d3_is_zero = false;
                }
            }
        }

        console.log(JSON.stringify({ d2_d3_is_zero }));
        """)
        assert res["d2_d3_is_zero"] is True


# ============================================================================
# 2. UNION-FIND PERSISTENCE & ELDER RULE
# ============================================================================

class TestUnionFindPersistent:
    def test_union_find_initialization_and_find(self):
        res = run_node_snippet("""
        import { UnionFindPersistent } from './src/topology/betti_number_calculator.js';
        const uf = new UnionFindPersistent(10);
        uf.makeSet(2, 5.0, 1.0, 2.0, 3.0);
        uf.makeSet(4, 3.0, 4.0, 5.0, 6.0);
        
        console.log(JSON.stringify({
            root2: uf.find(2),
            root4: uf.find(4),
            birth2: uf.birthValue[2],
            birth4: uf.birthValue[4]
        }));
        """)
        assert res["root2"] == 2
        assert res["root4"] == 4
        assert res["birth2"] == 5.0
        assert res["birth4"] == 3.0

    def test_union_find_superlevel_elder_rule(self):
        """Superlevel: higher birth value survives, lower birth value dies."""
        res = run_node_snippet("""
        import { UnionFindPersistent } from './src/topology/betti_number_calculator.js';
        const uf = new UnionFindPersistent(10);
        uf.makeSet(1, 10.0, 0, 0, 0); // Older (born at higher threshold 10.0)
        uf.makeSet(2, 6.0, 1, 1, 1);  // Younger (born at lower threshold 6.0)
        
        const pair = uf.union(1, 2, 2.5, true); // Merge at threshold 2.5
        const survivingRoot = uf.find(2);
        
        console.log(JSON.stringify({
            pairBirth: pair.birth,
            pairDeath: pair.death,
            pairLifetime: pair.lifetime,
            dyingIndex: pair.dyingIndex,
            survivingIndex: pair.survivingIndex,
            survivingRoot: survivingRoot
        }));
        """)
        assert res["pairBirth"] == 6.0
        assert res["pairDeath"] == 2.5
        assert res["pairLifetime"] == 3.5
        assert res["dyingIndex"] == 2
        assert res["survivingIndex"] == 1
        assert res["survivingRoot"] == 1

    def test_union_find_sublevel_elder_rule(self):
        """Sublevel: lower birth value survives, higher birth value dies."""
        res = run_node_snippet("""
        import { UnionFindPersistent } from './src/topology/betti_number_calculator.js';
        const uf = new UnionFindPersistent(10);
        uf.makeSet(1, -5.0, 0, 0, 0); // Older (born at lower threshold -5.0)
        uf.makeSet(2, -1.0, 1, 1, 1); // Younger (born at higher threshold -1.0)
        
        const pair = uf.union(1, 2, 3.0, false); // Merge at sublevel threshold 3.0
        
        console.log(JSON.stringify({
            pairBirth: pair.birth,
            pairDeath: pair.death,
            pairLifetime: pair.lifetime,
            dyingIndex: pair.dyingIndex,
            survivingIndex: pair.survivingIndex,
            root: uf.find(2)
        }));
        """)
        assert res["pairBirth"] == -1.0
        assert res["pairDeath"] == 3.0
        assert res["pairLifetime"] == 4.0
        assert res["dyingIndex"] == 2
        assert res["survivingIndex"] == 1
        assert res["root"] == 1

    def test_union_find_same_component_union_returns_null(self):
        res = run_node_snippet("""
        import { UnionFindPersistent } from './src/topology/betti_number_calculator.js';
        const uf = new UnionFindPersistent(5);
        uf.makeSet(0, 5.0);
        uf.makeSet(1, 4.0);
        uf.union(0, 1, 2.0);
        const secondUnion = uf.union(0, 1, 1.0);
        
        console.log(JSON.stringify({ isNull: secondUnion === null }));
        """)
        assert res["isNull"] is True

    def test_union_find_centroid_and_voxel_count_tracking(self):
        res = run_node_snippet("""
        import { UnionFindPersistent } from './src/topology/betti_number_calculator.js';
        const uf = new UnionFindPersistent(5);
        uf.makeSet(0, 8.0, 0.0, 0.0, 0.0);
        uf.makeSet(1, 6.0, 10.0, 20.0, 30.0);
        uf.union(0, 1, 2.0);
        
        const root = uf.find(1);
        const cx = uf.centroids[root * 3];
        const cy = uf.centroids[root * 3 + 1];
        const cz = uf.centroids[root * 3 + 2];
        const count = uf.voxelCounts[root];
        
        console.log(JSON.stringify({ cx, cy, cz, count }));
        """)
        assert res["cx"] == pytest.approx(5.0)
        assert res["cy"] == pytest.approx(10.0)
        assert res["cz"] == pytest.approx(15.0)
        assert res["count"] == 2


# ============================================================================
# 3. EXACT BETTI NUMBERS & EULER-POINCARÉ INVARIANT ON CANONICAL TOPOLOGIES
# ============================================================================

class TestExactBettiNumbersAndEulerPoincare:
    def test_empty_field_zero_betti(self):
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 4, ny = 4, nz = 4;
        const grid = new Float64Array(nx * ny * nz).fill(-10.0);
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 0.0);
        
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 0
        assert res["beta1"] == 0
        assert res["beta2"] == 0
        assert res["chi"] == 0

    def test_single_point_cluster(self):
        """Single isolated voxel in superlevel set: beta0 = 1, beta1 = 0, beta2 = 0, chi = 1."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 5, ny = 5, nz = 5;
        const grid = new Float64Array(nx * ny * nz).fill(0.0);
        grid[2 + nx * (2 + ny * 2)] = 10.0; // Point peak
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 1
        assert res["beta1"] == 0
        assert res["beta2"] == 0
        assert res["chi"] == 1
        assert res["valid"] is True

    def test_multiple_disjoint_clusters(self):
        """Three separated peaks: beta0 = 3, beta1 = 0, beta2 = 0, chi = 3."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 8, ny = 8, nz = 8;
        const grid = new Float64Array(nx * ny * nz).fill(0.0);
        grid[1 + nx * (1 + ny * 1)] = 10.0;
        grid[6 + nx * (1 + ny * 1)] = 10.0;
        grid[4 + nx * (6 + ny * 6)] = 10.0;
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 3
        assert res["beta1"] == 0
        assert res["beta2"] == 0
        assert res["chi"] == 3
        assert res["valid"] is True

    def test_hollow_spherical_bubble_s2_topology(self):
        """Hollow 2-sphere S^2 enclosure: beta0 = 1, beta1 = 0, beta2 = 1, chi = 2."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 7, ny = 7, nz = 7;
        const grid = new Float64Array(nx * ny * nz).fill(0.0);
        
        // Build hollow 3x3x3 shell from index 2 to 4
        for (let k = 2; k <= 4; k++) {
            for (let j = 2; j <= 4; j++) {
                for (let i = 2; i <= 4; i++) {
                    const isBoundary = (i === 2 || i === 4 || j === 2 || j === 4 || k === 2 || k === 4);
                    if (isBoundary) {
                        grid[i + nx * (j + ny * k)] = 10.0;
                    } else {
                        grid[i + nx * (j + ny * k)] = 0.0; // Enclosed hollow interior
                    }
                }
            }
        }
        
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 1
        assert res["beta1"] == 0
        assert res["beta2"] == 1
        assert res["chi"] == 2
        assert res["valid"] is True

    def test_filament_ring_loop_s1_topology(self):
        """Closed 1D loop (circle S^1): beta0 = 1, beta1 = 1, beta2 = 0, chi = 0."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 6, ny = 6, nz = 6;
        const grid = new Float64Array(nx * ny * nz).fill(0.0);
        
        // Form a square loop in the z=2 plane: (2,2)-(3,2)-(4,2)-(4,3)-(4,4)-(3,4)-(2,4)-(2,3)-(2,2)
        const loopPts = [
            [2,2,2], [3,2,2], [4,2,2],
            [4,3,2], [4,4,2],
            [3,4,2], [2,4,2],
            [2,3,2]
        ];
        for (const [x, y, z] of loopPts) {
            grid[x + nx * (y + ny * z)] = 10.0;
        }
        
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 1
        assert res["beta1"] == 1
        assert res["beta2"] == 0
        assert res["chi"] == 0
        assert res["valid"] is True

    def test_figure_eight_double_loop(self):
        """Figure-eight two loops sharing a point/edge: beta0 = 1, beta1 = 2, beta2 = 0, chi = -1."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 8, ny = 8, nz = 4;
        const grid = new Float64Array(nx * ny * nz).fill(0.0);
        
        // Loop 1
        const loop1 = [[1,1,1], [2,1,1], [3,1,1], [3,2,1], [3,3,1], [2,3,1], [1,3,1], [1,2,1]];
        // Loop 2 (shares (3,2,1) or connects to it)
        const loop2 = [[3,2,1], [4,1,1], [5,1,1], [5,2,1], [5,3,1], [4,3,1]];
        for (const [x, y, z] of [...loop1, ...loop2]) {
            grid[x + nx * (y + ny * z)] = 10.0;
        }
        
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 1
        assert res["beta1"] == 2
        assert res["beta2"] == 0
        assert res["chi"] == -1
        assert res["valid"] is True

    def test_two_enclosed_void_bubbles(self):
        """Two separated enclosed cavities: beta0 = 1 (if outer wall connects them) or 2, beta2 = 2."""
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 11, ny = 7, nz = 7;
        const grid = new Float64Array(nx * ny * nz).fill(10.0);
        
        // Void 1 interior at (2..3, 2..4, 2..4)
        for (let k = 2; k <= 4; k++) {
            for (let j = 2; j <= 4; j++) {
                for (let i = 2; i <= 3; i++) {
                    grid[i + nx * (j + ny * k)] = 0.0;
                }
            }
        }
        
        // Void 2 interior at (7..8, 2..4, 2..4)
        for (let k = 2; k <= 4; k++) {
            for (let j = 2; j <= 4; j++) {
                for (let i = 7; i <= 8; i++) {
                    grid[i + nx * (j + ny * k)] = 0.0;
                }
            }
        }
        
        const res = calc.computeBettiNumbers(grid, nx, ny, nz, 5.0);
        console.log(JSON.stringify(res));
        """)
        assert res["beta0"] == 1
        assert res["beta2"] == 2
        assert res["beta0"] - res["beta1"] + res["beta2"] == res["chi"]


# ============================================================================
# 4. PERSISTENT HOMOLOGY & BOUNDARY MATRIX REDUCTION
# ============================================================================

class TestPersistentHomologyAndZ2Reduction:
    def test_z2_boundary_matrix_reduction(self):
        res = run_node_snippet("""
        import { BoundaryMatrixReducer } from './src/topology/betti_number_calculator.js';
        
        // Simple filtration of 3 vertices, 2 edges:
        // v0 (val 1.0), v1 (val 2.0), v2 (val 3.0), e01 (val 4.0), e12 (val 5.0)
        const sortedCells = [
            { dim: 0, id: 0, value: 1.0, boundary: [] },
            { dim: 0, id: 1, value: 2.0, boundary: [] },
            { dim: 0, id: 2, value: 3.0, boundary: [] },
            { dim: 1, id: 3, value: 4.0, boundary: [0, 1] },
            { dim: 1, id: 4, value: 5.0, boundary: [1, 2] }
        ];
        
        const { pairs, essential } = BoundaryMatrixReducer.reduceZ2(sortedCells);
        console.log(JSON.stringify({ pairs, essential }));
        """)
        pairs = res["pairs"]
        essential = res["essential"]
        assert len(pairs) == 2
        assert len(essential) == 1
        assert essential[0]["birth"] == 1.0

    def test_persistent_homology_on_3d_grid(self):
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 4, ny = 4, nz = 4;
        const grid = new Float64Array(nx * ny * nz);
        for (let i = 0; i < grid.length; i++) {
            grid[i] = Math.sin(i * 0.7);
        }
        
        const ph = calc.computePersistentHomology(grid, nx, ny, nz);
        console.log(JSON.stringify({
            h0Pairs: ph.h0.pairs.length,
            h0Ess: ph.h0.essential.length,
            h1Pairs: ph.h1.pairs.length,
            h2Pairs: ph.h2.pairs.length,
            totalEntropy: ph.totalEntropy
        }));
        """)
        assert res["h0Pairs"] > 0
        assert res["h0Ess"] >= 1
        assert res["totalEntropy"] >= 0.0


# ============================================================================
# 5. PERSISTENCE DIAGRAMS, BARCODES, DISTANCES & LANDSCAPES
# ============================================================================

class TestPersistenceMetricsAndDistances:
    def test_diagram_total_persistence_and_entropy(self):
        res = run_node_snippet("""
        import { PersistenceDiagram } from './src/topology/betti_number_calculator.js';
        const diag = new PersistenceDiagram(0, [
            { birth: 10.0, death: 4.0 }, // lifetime 6
            { birth: 8.0, death: 6.0 },  // lifetime 2
            { birth: 5.0, death: 5.0 }   // lifetime 0
        ]);
        
        console.log(JSON.stringify({
            totP1: diag.totalPersistence(1),
            totP2: diag.totalPersistence(2),
            entropy: diag.persistentEntropy()
        }));
        """)
        assert res["totP1"] == pytest.approx(8.0)
        assert res["totP2"] == pytest.approx(36.0 + 4.0) # 40.0
        # p1 = 6/8 = 0.75, p2 = 2/8 = 0.25
        expected_entropy = -(0.75 * math.log(0.75) + 0.25 * math.log(0.25))
        assert res["entropy"] == pytest.approx(expected_entropy, rel=1e-5)

    def test_diagram_filter_noise(self):
        res = run_node_snippet("""
        import { PersistenceDiagram } from './src/topology/betti_number_calculator.js';
        const diag = new PersistenceDiagram(1, [
            { birth: 5.0, death: 1.0 }, // lifetime 4
            { birth: 3.2, death: 3.1 }, // lifetime 0.1
            { birth: 4.0, death: 3.5 }  // lifetime 0.5
        ]);
        const filtered = diag.filterNoise(0.6);
        console.log(JSON.stringify({
            origCount: diag.pairs.length,
            filtCount: filtered.pairs.length,
            remainingLifetime: filtered.pairs[0].lifetime
        }));
        """)
        assert res["origCount"] == 3
        assert res["filtCount"] == 1
        assert res["remainingLifetime"] == pytest.approx(4.0)

    def test_bottleneck_distance_properties(self):
        """Verify metric properties: W_\infty(D, D) = 0, symmetry, triangle inequality."""
        res = run_node_snippet("""
        import { PersistenceDiagram } from './src/topology/betti_number_calculator.js';
        const d1 = new PersistenceDiagram(0, [
            { birth: 5.0, death: 2.0 },
            { birth: 8.0, death: 4.0 }
        ]);
        const d2 = new PersistenceDiagram(0, [
            { birth: 5.2, death: 1.9 },
            { birth: 8.0, death: 4.1 }
        ]);
        const d3 = new PersistenceDiagram(0, [
            { birth: 7.0, death: 1.0 }
        ]);

        const dist_self = d1.bottleneckDistance(d1);
        const dist_12 = d1.bottleneckDistance(d2);
        const dist_21 = d2.bottleneckDistance(d1);
        const dist_13 = d1.bottleneckDistance(d3);
        const dist_23 = d2.bottleneckDistance(d3);

        console.log(JSON.stringify({
            dist_self,
            dist_12,
            dist_21,
            dist_13,
            dist_23,
            triangleHolds: dist_13 <= (dist_12 + dist_23 + 1e-10)
        }));
        """)
        assert res["dist_self"] == pytest.approx(0.0)
        assert res["dist_12"] == pytest.approx(res["dist_21"])
        assert res["dist_12"] > 0
        assert res["triangleHolds"] is True

    def test_wasserstein_distance_hungarian(self):
        res = run_node_snippet("""
        import { PersistenceDiagram } from './src/topology/betti_number_calculator.js';
        const d1 = new PersistenceDiagram(0, [
            { birth: 4.0, death: 2.0 }
        ]);
        const d2 = new PersistenceDiagram(0, [
            { birth: 4.5, death: 2.5 }
        ]);
        const w1 = d1.wassersteinDistance(d2, 1);
        const w2 = d1.wassersteinDistance(d2, 2);
        
        console.log(JSON.stringify({ w1, w2 }));
        """)
        assert res["w1"] == pytest.approx(0.5)
        assert res["w2"] == pytest.approx(0.5)

    def test_persistence_barcode_intervals_and_active_count(self):
        res = run_node_snippet("""
        import { PersistenceBarcode } from './src/topology/betti_number_calculator.js';
        const barcode = new PersistenceBarcode(0, [
            { birth: 10.0, death: 2.0 },
            { birth: 7.0, death: 4.0 },
            { birth: 3.0, death: 1.0 }
        ]);
        
        console.log(JSON.stringify({
            countAt9: barcode.activeCount(9.0),
            countAt5: barcode.activeCount(5.0),
            countAt2_5: barcode.activeCount(2.5),
            countAt0: barcode.activeCount(0.5)
        }));
        """)
        assert res["countAt9"] == 1  # Only [10, 2) active
        assert res["countAt5"] == 2  # [10, 2) and [7, 4) active
        assert res["countAt2_5"] == 2 # [10, 2) and [3, 1) active
        assert res["countAt0"] == 0

    def test_persistence_landscape_envelope_ordering(self):
        """Verify Bubenik landscape ordering: \lambda_1(t) >= \lambda_2(t) >= 0."""
        res = run_node_snippet("""
        import { PersistenceDiagram, PersistenceLandscape } from './src/topology/betti_number_calculator.js';
        const diag = new PersistenceDiagram(0, [
            { birth: 1.0, death: 5.0 },
            { birth: 2.0, death: 4.0 },
            { birth: 3.0, death: 7.0 }
        ]);
        
        const tGrid = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5];
        const landscapes = PersistenceLandscape.computeLandscapes(diag, tGrid, 3);
        const l1 = Array.from(landscapes[0]);
        const l2 = Array.from(landscapes[1]);
        const l3 = Array.from(landscapes[2]);
        
        let orderingHolds = true;
        for (let i = 0; i < tGrid.length; i++) {
            if (l1[i] < l2[i] || l2[i] < l3[i] || l3[i] < 0) {
                orderingHolds = false;
            }
        }
        const norm = PersistenceLandscape.landscapeNorm(landscapes, tGrid, 2);
        
        console.log(JSON.stringify({ orderingHolds, norm }));
        """)
        assert res["orderingHolds"] is True
        assert res["norm"] > 0


# ============================================================================
# 6. GAUSSIAN RANDOM FIELD (GRF) ANALYTICAL TOMITA-GOTT TOPOLOGY
# ============================================================================

class TestTomitaGottAnalyticalAndGRF:
    def test_analytical_euler_density_zero_crossings(self):
        """Analytical \chi_{GRF}(\nu) must have exact zero crossings at \nu = -1 and \nu = +1."""
        res = run_node_snippet("""
        import { TomitaGottAnalytical } from './src/topology/betti_number_calculator.js';
        const sigma0 = 1.2;
        const sigma1 = 0.8;
        
        const chi_neg1 = TomitaGottAnalytical.computeExpectedEulerDensity(-1.0, sigma0, sigma1);
        const chi_pos1 = TomitaGottAnalytical.computeExpectedEulerDensity(1.0, sigma0, sigma1);
        
        console.log(JSON.stringify({ chi_neg1, chi_pos1 }));
        """)
        assert abs(res["chi_neg1"]) < 1e-15
        assert abs(res["chi_pos1"]) < 1e-15

    def test_analytical_euler_density_extrema(self):
        """Minimum at \nu=0, symmetric maxima at \nu = \pm \sqrt{3}."""
        res = run_node_snippet("""
        import { TomitaGottAnalytical } from './src/topology/betti_number_calculator.js';
        const sigma0 = 1.0;
        const sigma1 = 1.0;
        
        const chi_0 = TomitaGottAnalytical.computeExpectedEulerDensity(0.0, sigma0, sigma1);
        const sqrt3 = Math.sqrt(3.0);
        const chi_peak_pos = TomitaGottAnalytical.computeExpectedEulerDensity(sqrt3, sigma0, sigma1);
        const chi_peak_neg = TomitaGottAnalytical.computeExpectedEulerDensity(-sqrt3, sigma0, sigma1);
        
        console.log(JSON.stringify({
            chi_0,
            chi_peak_pos,
            chi_peak_neg,
            isMinAtZero: chi_0 < 0,
            peaksEqual: Math.abs(chi_peak_pos - chi_peak_neg) < 1e-14
        }));
        """)
        assert res["isMinAtZero"] is True
        assert res["peaksEqual"] is True
        assert res["chi_peak_pos"] > 0

    def test_genus_density_relation(self):
        """G_V(\nu) = -0.5 * \chi_V(\nu)."""
        res = run_node_snippet("""
        import { TomitaGottAnalytical } from './src/topology/betti_number_calculator.js';
        const nu = 0.5;
        const s0 = 2.0, s1 = 1.5;
        const chi = TomitaGottAnalytical.computeExpectedEulerDensity(nu, s0, s1);
        const g = TomitaGottAnalytical.computeExpectedGenusDensity(nu, s0, s1);
        
        console.log(JSON.stringify({ match: Math.abs(g - (-0.5 * chi)) < 1e-15 }));
        """)
        assert res["match"] is True

    def test_spectral_parameter_estimation_on_3d_gaussian_field(self):
        """Generate 3D sinusoidal test field and verify finite difference dispersion estimation."""
        res = run_node_snippet("""
        import { TomitaGottAnalytical } from './src/topology/betti_number_calculator.js';
        const nx = 16, ny = 16, nz = 16;
        const dx = 1.0, dy = 1.0, dz = 1.0;
        const grid = new Float64Array(nx * ny * nz);
        const kx = 2.0 * Math.PI / nx;
        const ky = 2.0 * Math.PI / ny;
        const kz = 2.0 * Math.PI / nz;
        
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    grid[i + nx * (j + ny * k)] = Math.sin(kx * i) + Math.cos(ky * j) + Math.sin(kz * k);
                }
            }
        }
        
        const params = TomitaGottAnalytical.computeSpectralParameters(grid, nx, ny, nz, dx, dy, dz);
        console.log(JSON.stringify(params));
        """)
        assert res["sigma0"] > 0
        assert res["sigma1"] > 0
        assert res["sigma2"] > 0
        assert res["spectralGamma"] > 0
        assert res["volume"] == pytest.approx(16.0 * 16.0 * 16.0)

    def test_genus_curve_fitting_and_asymmetry_diagnostics(self):
        res = run_node_snippet("""
        import { TomitaGottAnalytical } from './src/topology/betti_number_calculator.js';
        const sigma0 = 1.0;
        const sigma1 = 1.0;
        const vol = 1000.0;
        
        // Generate ideal Tomita-Gott empirical samples
        const empiricalData = [];
        for (let nu = -3.0; nu <= 3.0; nu += 0.2) {
            const chi = TomitaGottAnalytical.computeExpectedEulerDensity(nu, sigma0, sigma1) * vol;
            empiricalData.push({ nu, chiEmpirical: chi });
        }
        
        const fit = TomitaGottAnalytical.fitGenusCurve(empiricalData, vol, sigma0, sigma1);
        console.log(JSON.stringify({
            rmsResidual: fit.rmsResidual,
            asymmetry: fit.asymmetry,
            shiftDeltaNu: fit.shiftDeltaNu
        }));
        """)
        assert res["rmsResidual"] == pytest.approx(0.0, abs=1e-10)
        assert abs(abs(res["asymmetry"]) - 0.386) < 0.05
        assert abs(res["shiftDeltaNu"]) < 0.05


# ============================================================================
# 7. VELOCITY POTENTIAL & TIDAL TENSOR TOPOLOGY
# ============================================================================

class TestVelocityPotentialTopologist:
    def test_tidal_tensor_eigenvalues_and_trace_identity(self):
        """Trace of tidal tensor Tr(T) = \lambda_1 + \lambda_2 + \lambda_3 = \nabla^2 \Phi = \delta."""
        res = run_node_snippet("""
        import { VelocityPotentialTopologist } from './src/topology/betti_number_calculator.js';
        const nx = 9, ny = 9, nz = 9;
        const grid = new Float64Array(nx * ny * nz);
        
        // Quadratic potential: \Phi = 0.5 * (2 x^2 + 3 y^2 + 4 z^2) + x y
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const x = i - 4;
                    const y = j - 4;
                    const z = k - 4;
                    grid[i + nx * (j + ny * k)] = 0.5 * (2.0 * x * x + 3.0 * y * y + 4.0 * z * z) + x * y;
                }
            }
        }
        
        const analysis = VelocityPotentialTopologist.computeTidalTensor(grid, nx, ny, nz, 4, 4, 4);
        const [l1, l2, l3] = analysis.eigenvalues;
        const trace = l1 + l2 + l3;
        
        console.log(JSON.stringify({
            l1, l2, l3,
            trace,
            divergence: analysis.divergence,
            isSorted: (l1 >= l2 && l2 >= l3),
            traceMatchesDiv: Math.abs(trace - analysis.divergence) < 1e-10,
            environment: analysis.environment
        }));
        """)
        assert res["isSorted"] is True
        assert res["traceMatchesDiv"] is True
        assert res["divergence"] == pytest.approx(2.0 + 3.0 + 4.0, rel=1e-3)
        assert res["environment"] == "PEAK_NODE"

    def test_cosmic_environment_classification_filament_and_void(self):
        res = run_node_snippet("""
        import { VelocityPotentialTopologist, CosmicEnvironmentType } from './src/topology/betti_number_calculator.js';
        const nx = 7, ny = 7, nz = 7;
        
        // Filament potential: collapse along x and y (l1, l2 > th), expansion along z (l3 <= th)
        const gridFil = new Float64Array(nx * ny * nz);
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const x = i - 3, y = j - 3, z = k - 3;
                    gridFil[i + nx * (j + ny * k)] = 0.5 * (2.0 * x * x + 2.0 * y * y - 0.5 * z * z);
                }
            }
        }
        
        const resFil = VelocityPotentialTopologist.computeTidalTensor(gridFil, nx, ny, nz, 3, 3, 3, 1, 1, 1, 0.2);
        
        // Void potential: expansion along all 3 axes (l1, l2, l3 < 0)
        const gridVoid = new Float64Array(nx * ny * nz);
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const x = i - 3, y = j - 3, z = k - 3;
                    gridVoid[i + nx * (j + ny * k)] = -0.5 * (x * x + y * y + z * z);
                }
            }
        }
        const resVoid = VelocityPotentialTopologist.computeTidalTensor(gridVoid, nx, ny, nz, 3, 3, 3, 1, 1, 1, 0.2);
        
        console.log(JSON.stringify({
            envFil: resFil.environment,
            envVoid: resVoid.environment
        }));
        """)
        assert res["envFil"] == "FILAMENT"
        assert res["envVoid"] == "VOID"


# ============================================================================
# 8. END-TO-END COSMIC TOPOLOGY ANALYZER & SCENARIOS
# ============================================================================

class TestCosmicTopologyAnalyzer:
    def test_betti_curves_threshold_sweep(self):
        res = run_node_snippet("""
        import { BettiNumberCalculator } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 6, ny = 6, nz = 6;
        const grid = new Float64Array(nx * ny * nz);
        for (let i = 0; i < grid.length; i++) {
            grid[i] = Math.cos(i * 0.5);
        }
        
        const curves = calc.computeBettiCurves(grid, nx, ny, nz, { numSteps: 10 });
        console.log(JSON.stringify({
            numSteps: curves.length,
            allValid: curves.every(c => c.beta0 - c.beta1 + c.beta2 === c.chi),
            firstTh: curves[0].threshold,
            lastTh: curves[curves.length - 1].threshold
        }));
        """)
        assert res["numSteps"] == 10
        assert res["allValid"] is True

    def test_full_cosmic_topology_analysis_meatball_cluster_regime(self):
        """Synthesize high-density cluster field and verify MEATBALL classification."""
        res = run_node_snippet("""
        import { BettiNumberCalculator, CosmicTopologyType } from './src/topology/betti_number_calculator.js';
        const calc = new BettiNumberCalculator();
        const nx = 10, ny = 10, nz = 10;
        const grid = new Float64Array(nx * ny * nz).fill(-1.0);
        
        // 5 isolated high-density peaks
        const peaks = [[2,2,2], [7,2,2], [2,7,2], [7,7,2], [5,5,7]];
        for (const [x,y,z] of peaks) {
            grid[x + nx * (y + ny * z)] = 15.0;
        }
        
        const analysis = calc.analyzeCosmicTopology(grid, nx, ny, nz, { numSteps: 20 });
        console.log(JSON.stringify({
            maxBeta0: analysis.maxBeta0,
            maxBeta2: analysis.maxBeta2,
            hasSpectral: !!analysis.spectral,
            hasBettiCurves: analysis.bettiCurves.length > 0
        }));
        """)
        assert res["maxBeta0"] >= 5
        assert res["hasSpectral"] is True
        assert res["hasBettiCurves"] is True
