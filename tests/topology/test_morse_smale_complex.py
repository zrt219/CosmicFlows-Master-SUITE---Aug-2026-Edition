"""
tests/topology/test_morse_smale_complex.py
Exhaustive Pytest verification suite for:
- 3D Morse-Smale Complex (MSC) decomposition & cell hierarchy (0-cells, 1-cells, 2-cells, 3-cells)
- Boundary operators \\partial_1, \\partial_2, \\partial_3 and topological homology chain groups C_0, C_1, C_2, C_3
- Exact boundary operator nilpotency verification (\\partial_1 \\circ \\partial_2 = 0 and \\partial_2 \\circ \\partial_3 = 0)
- Betti numbers b_0 (connected components), b_1 (tunnels / loops), b_2 (enclosed cosmic voids), b_3
- Euler-Poincare characteristic \\chi = b_0 - b_1 + b_2 - b_3 = |C_0| - |C_1| + |C_2| - |C_3|
- Topological simplification via Morse cancellation of persistence pairs with threshold \\epsilon_p
- Cosmic Web Graph network model (filaments, wall sheets, void crystals, Dijkstra routing, GeoJSON)
- 3D Discrete Morse-Smale Complex extraction on scalar cosmological grids
"""

import json
import subprocess
import os
import math
import pytest
import numpy as np

NODE_EXEC = "node"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def run_node_snippet(code: str) -> dict:
    """Executes an inline ES6 Node.js snippet and returns parsed JSON output."""
    wrapped_code = f"""
    {code}
    """
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", wrapped_code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed with code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")
    return json.loads(proc.stdout)


# ============================================================================
# TIER 1: Vector3D & Matrix Linear Algebra Engine Tests
# ============================================================================

def test_vector3d_algebraic_operations():
    """Verify Vector3D arithmetic, dot product, cross product, norm, and distance."""
    code = """
    import { Vector3D } from './src/topology/morse_smale_complex.js';

    const v1 = new Vector3D(1.0, 2.0, 3.0);
    const v2 = new Vector3D(4.0, -5.0, 6.0);

    const sum = v1.add(v2);
    const diff = v1.sub(v2);
    const scaled = v1.scale(2.5);
    const dot = v1.dot(v2);
    const cross = v1.cross(v2);
    const norm1 = v1.norm();
    const unit1 = v1.normalize();
    const dist = v1.distanceTo(v2);

    console.log(JSON.stringify({
      sum: sum.toArray(),
      diff: diff.toArray(),
      scaled: scaled.toArray(),
      dot: dot,
      cross: cross.toArray(),
      norm1: norm1,
      unit1: unit1.toArray(),
      unitNorm: unit1.norm(),
      dist: dist
    }));
    """
    res = run_node_snippet(code)
    assert res["sum"] == [5.0, -3.0, 9.0]
    assert res["diff"] == [-3.0, 7.0, -3.0]
    assert res["scaled"] == [2.5, 5.0, 7.5]
    assert math.isclose(res["dot"], 1.0 * 4.0 + 2.0 * -5.0 + 3.0 * 6.0) # 4 - 10 + 18 = 12
    assert res["cross"] == [27.0, 6.0, -13.0]
    assert math.isclose(res["norm1"], math.sqrt(1 + 4 + 9))
    assert math.isclose(res["unitNorm"], 1.0, rel_tol=1e-9)
    assert math.isclose(res["dist"], math.sqrt(9 + 49 + 9))


def test_linear_algebra_jacobi_eigenvalues():
    """Verify Jacobi 3x3 eigensolver accurately diagonalizes symmetric tensor and preserves invariants."""
    code = """
    import { LinearAlgebraEngine } from './src/topology/morse_smale_complex.js';

    // Symmetric matrix with known eigenvalues [1, 2, 5]
    const mat = [
      [3.0, 2.0, 0.0],
      [2.0, 3.0, 0.0],
      [0.0, 0.0, 2.0]
    ];

    const eigen = LinearAlgebraEngine.jacobiEigen3x3(mat);
    const [l1, l2, l3] = eigen.eigenvalues;
    const [v1, v2, v3] = eigen.eigenvectors;

    // Check orthogonality: v1 . v2 == 0
    const dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    const dot13 = v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2];
    const dot23 = v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2];

    console.log(JSON.stringify({
      eigenvalues: [l1, l2, l3],
      trace: l1 + l2 + l3,
      det: l1 * l2 * l3,
      orthogonality: [dot12, dot13, dot23]
    }));
    """
    res = run_node_snippet(code)
    evals = sorted(res["eigenvalues"])
    assert math.isclose(evals[0], 1.0, abs_tol=1e-7)
    assert math.isclose(evals[1], 2.0, abs_tol=1e-7)
    assert math.isclose(evals[2], 5.0, abs_tol=1e-7)
    assert math.isclose(res["trace"], 8.0, abs_tol=1e-7) # 3 + 3 + 2 = 8
    assert math.isclose(res["det"], 10.0, abs_tol=1e-7)   # (9-4)*2 = 10
    for dot in res["orthogonality"]:
        assert math.isclose(dot, 0.0, abs_tol=1e-7)


def test_linear_algebra_matrix_rank_gaussian():
    """Verify Gaussian elimination rank solver on full rank, rank-deficient, and rectangular matrices."""
    code = """
    import { LinearAlgebraEngine } from './src/topology/morse_smale_complex.js';

    // Full rank 3x3
    const fullRank3x3 = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];

    // Rank 2 3x3 (row3 = row1 + row2)
    const rank2 = [
      [1, 2, 3],
      [4, 5, 6],
      [5, 7, 9]
    ];

    // Rectangular 2x4 with rank 2
    const rect2x4 = [
      [1, 0, 2, -1],
      [0, 1, 3, 4]
    ];

    // Zero matrix 3x3
    const zeroMat = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0]
    ];

    console.log(JSON.stringify({
      rankFull: LinearAlgebraEngine.computeRank(fullRank3x3),
      rank2: LinearAlgebraEngine.computeRank(rank2),
      rankRect: LinearAlgebraEngine.computeRank(rect2x4),
      rankZero: LinearAlgebraEngine.computeRank(zeroMat)
    }));
    """
    res = run_node_snippet(code)
    assert res["rankFull"] == 3
    assert res["rank2"] == 2
    assert res["rankRect"] == 2
    assert res["rankZero"] == 0


def test_linear_algebra_rank_z2():
    """Verify GF(2) / Z_2 rank computation on binary incidence matrices."""
    code = """
    import { LinearAlgebraEngine } from './src/topology/morse_smale_complex.js';

    // In Z_2: row1 + row2 = [1, 1, 0] + [0, 1, 1] = [1, 0, 1]
    const depZ2 = [
      [1, 1, 0],
      [0, 1, 1],
      [1, 0, 1]
    ];

    const indepZ2 = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ];

    console.log(JSON.stringify({
      rankDepZ2: LinearAlgebraEngine.computeRankZ2(depZ2),
      rankIndepZ2: LinearAlgebraEngine.computeRankZ2(indepZ2)
    }));
    """
    res = run_node_snippet(code)
    assert res["rankDepZ2"] == 2
    assert res["rankIndepZ2"] == 3


# ============================================================================
# TIER 2: Morse Critical Points & Cosmological Scientific Invariants
# ============================================================================

def test_morse_critical_point_indices_and_types():
    """Verify critical points for all 4 Morse indices (0: Void, 1: Wall, 2: Filament, 3: Cluster)."""
    code = """
    import { MorseCriticalPoint, MorseIndex, CosmicWebElementType } from './src/topology/morse_smale_complex.js';

    const minVoid = new MorseCriticalPoint({
      id: 'void_01',
      index: MorseIndex.MINIMUM,
      position: [-20.0, 10.0, 5.0],
      value: -0.95,
      eigenvalues: [0.02, 0.05, 0.08],
      velocity: [-50.0, 30.0, 10.0],
      name: 'Sculptor Void Core'
    });

    const wallSad = new MorseCriticalPoint({
      id: 'wall_01',
      index: MorseIndex.WALL_SADDLE,
      position: [-10.0, 20.0, 15.0],
      value: -0.15,
      eigenvalues: [-0.04, 0.03, 0.06]
    });

    const filSad = new MorseCriticalPoint({
      id: 'fil_01',
      index: MorseIndex.FILAMENT_SADDLE,
      position: [15.0, -30.0, 25.0],
      value: 0.65,
      eigenvalues: [-0.08, -0.05, 0.07]
    });

    const maxCluster = new MorseCriticalPoint({
      id: 'cluster_01',
      index: MorseIndex.MAXIMUM,
      position: [-35.0, 15.0, -10.0],
      value: 2.85,
      eigenvalues: [-0.15, -0.11, -0.08],
      name: 'Laniakea Core'
    });

    console.log(JSON.stringify({
      void: minVoid.toJSON(),
      wall: wallSad.toJSON(),
      fil: filSad.toJSON(),
      cluster: maxCluster.toJSON()
    }));
    """
    res = run_node_snippet(code)
    assert res["void"]["index"] == 0
    assert res["void"]["type"] == "MINIMUM_VOID"
    assert res["void"]["cosmicType"] == "VOID_CORE"
    assert res["void"]["positionMpc"] == [-20.0, 10.0, 5.0]

    assert res["wall"]["index"] == 1
    assert res["wall"]["type"] == "SADDLE_WALL"
    assert res["wall"]["cosmicType"] == "WALL_HUB"

    assert res["fil"]["index"] == 2
    assert res["fil"]["type"] == "SADDLE_FILAMENT"
    assert res["fil"]["cosmicType"] == "FILAMENT_HUB"

    assert res["cluster"]["index"] == 3
    assert res["cluster"]["type"] == "MAXIMUM_CLUSTER"
    assert res["cluster"]["cosmicType"] == "CLUSTER_NODE"


def test_morse_critical_point_bounds_and_error_handling():
    """Verify invalid Morse index or missing mandatory parameters throw TypeError / RangeError."""
    code = """
    import { MorseCriticalPoint } from './src/topology/morse_smale_complex.js';

    const errors = [];
    try {
      new MorseCriticalPoint({ id: 'bad1', index: 4, value: 1.0 }); // index > 3
    } catch (e) {
      errors.push(e.name);
    }

    try {
      new MorseCriticalPoint({ id: 'bad2', index: -1, value: 1.0 }); // index < 0
    } catch (e) {
      errors.push(e.name);
    }

    try {
      new MorseCriticalPoint({ id: 'bad3' }); // missing index and value
    } catch (e) {
      errors.push(e.name);
    }

    console.log(JSON.stringify({ errorCount: errors.length, errors: errors }));
    """
    res = run_node_snippet(code)
    assert res["errorCount"] == 3
    assert "RangeError" in res["errors"]
    assert "TypeError" in res["errors"]


# ============================================================================
# TIER 3: Cell Hierarchy (0-Cells, 1-Cells, 2-Cells, 3-Cells)
# ============================================================================

def test_morse_1cell_filament_spine_geometry():
    """Verify Morse1Cell filament spine path, arc length, tortuosity, and linear mass density."""
    code = """
    import { MorseCriticalPoint, Morse1Cell, MorseIndex } from './src/topology/morse_smale_complex.js';

    const saddle = new MorseCriticalPoint({
      id: 'sad1',
      index: MorseIndex.FILAMENT_SADDLE,
      position: [0.0, 0.0, 0.0],
      value: 0.4
    });

    const maxNode = new MorseCriticalPoint({
      id: 'max1',
      index: MorseIndex.MAXIMUM,
      position: [30.0, 40.0, 0.0], // Chord length = 50.0 Mpc/h
      value: 2.5
    });

    // S-curved polyline path: (0,0,0) -> (10, 25, 5) -> (20, 15, -5) -> (30, 40, 0)
    const path = [
      [0.0, 0.0, 0.0],
      [10.0, 25.0, 5.0],
      [20.0, 15.0, -5.0],
      [30.0, 40.0, 0.0]
    ];

    const fil = new Morse1Cell({
      id: 'fil_01',
      source: saddle,
      target: maxNode,
      path: path,
      linearMassDensity: 3.2e13
    });

    console.log(JSON.stringify({
      lengthMpc: fil.lengthMpc,
      tortuosity: fil.getTortuosity(),
      boundary: fil.getBoundary(),
      meanDensity: fil.meanDensityContrast,
      linearMass: fil.linearMassDensity
    }));
    """
    res = run_node_snippet(code)
    assert res["lengthMpc"] > 50.0 # Curving path exceeds Euclidean chord 50.0
    assert res["tortuosity"] > 1.0 # Tortuosity > 1
    assert res["boundary"] == [
        {"cellId": "max1", "orientation": 1},
        {"cellId": "sad1", "orientation": -1}
    ]
    assert math.isclose(res["meanDensity"], 1.45, rel_tol=1e-5) # 0.5 * (0.4 + 2.5)
    assert math.isclose(res["linearMass"], 3.2e13)


def test_morse_2cell_wall_sheet_surface():
    """Verify Morse2Cell cosmic wall sheet surface area and boundary 1-cells."""
    code = """
    import { MorseCriticalPoint, Morse1Cell, Morse2Cell, MorseIndex } from './src/topology/morse_smale_complex.js';

    const p0 = new MorseCriticalPoint({ id: 'p0', index: 1, position: [0, 0, 0], value: 0 });
    const p1 = new MorseCriticalPoint({ id: 'p1', index: 0, position: [100, 0, 0], value: -1 });
    const p2 = new MorseCriticalPoint({ id: 'p2', index: 0, position: [100, 100, 0], value: -1 });

    const e01 = new Morse1Cell({ id: 'e01', source: p0, target: p1 });
    const e12 = new Morse1Cell({ id: 'e12', source: p1, target: p2 });
    const e20 = new Morse1Cell({ id: 'e20', source: p2, target: p0 });

    const wall = new Morse2Cell({
      id: 'wall_sheet_01',
      boundary1Cells: [e01, e12, e20],
      boundaryOrientations: [1, 1, 1],
      surfaceVertices: [p0.position, p1.position, p2.position],
      normalVector: [0, 0, 1]
    });

    console.log(JSON.stringify({
      area: wall.areaMpc2,
      normal: wall.normalVector.toArray(),
      boundary: wall.getBoundary()
    }));
    """
    res = run_node_snippet(code)
    assert math.isclose(res["area"], 5000.0, rel_tol=1e-3) # Triangle (0,0) (100,0) (100,100) area = 0.5 * 100 * 100 = 5000
    assert res["normal"] == [0.0, 0.0, 1.0]
    assert len(res["boundary"]) == 3


def test_morse_3cell_void_crystal_volume_and_sphericity():
    """Verify Morse3Cell void crystal volume, effective radius, and sphericity."""
    code = """
    import { MorseCriticalPoint, Morse2Cell, Morse3Cell, MorseIndex } from './src/topology/morse_smale_complex.js';

    const minP = new MorseCriticalPoint({ id: 'void_min', index: MorseIndex.MINIMUM, position: [0, 0, 0], value: -0.9 });
    const maxP = new MorseCriticalPoint({ id: 'cluster_max', index: MorseIndex.MAXIMUM, position: [50, 50, 50], value: 2.1 });

    // Sphere equivalent volume V = 33510.32 Mpc^3 => R_eff = 20.0 Mpc/h
    const V = (4.0 / 3.0) * Math.PI * Math.pow(20.0, 3.0);
    const surfaceArea = 4.0 * Math.PI * Math.pow(20.0, 2.0); // 5026.55

    // Dummy face with known area
    const face = new Morse2Cell({ id: 'face1', areaMpc2: surfaceArea });

    const crystal = new Morse3Cell({
      id: 'void_crystal_01',
      minimum: minP,
      maximum: maxP,
      boundary2Cells: [face],
      boundaryOrientations: [1],
      volumeMpc3: V
    });

    console.log(JSON.stringify({
      volume: crystal.volumeMpc3,
      rEff: crystal.getEffectiveRadiusMpc(),
      sphericity: crystal.getSphericity(),
      centroid: crystal.centroid.toArray()
    }));
    """
    res = run_node_snippet(code)
    assert math.isclose(res["rEff"], 20.0, rel_tol=1e-3)
    assert math.isclose(res["sphericity"], 1.0, rel_tol=1e-3) # Perfect sphere sphericity = 1.0
    assert res["centroid"] == [25.0, 25.0, 25.0]


# ============================================================================
# TIER 4: Boundary Operator Nilpotency (\partial \circ \partial = 0)
# ============================================================================

def test_tetrahedron_complex_nilpotency_exactness():
    """Verify \partial_1 \circ \partial_2 = 0 and \partial_2 \circ \partial_3 = 0 on solid tetrahedron."""
    code = """
    import { buildTetrahedralComplex } from './src/topology/morse_smale_complex.js';

    const tet = buildTetrahedralComplex();
    const nilpotency = tet.chainComplex.verifyNilpotency();
    const betti = tet.chainComplex.computeBettiNumbers();

    console.log(JSON.stringify({
      nilpotency: nilpotency,
      betti: betti
    }));
    """
    res = run_node_snippet(code)
    assert res["nilpotency"]["d1_d2_isZero"] is True
    assert res["nilpotency"]["d2_d3_isZero"] is True
    assert res["nilpotency"]["maxNorm_d1_d2"] < 1e-12
    assert res["nilpotency"]["maxNorm_d2_d3"] < 1e-12
    assert res["nilpotency"]["isValid"] is True

    # Homology of solid tetrahedron (contractible 3-ball B^3):
    # b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0, Euler characteristic \chi = 1
    assert res["betti"]["b0"] == 1
    assert res["betti"]["b1"] == 0
    assert res["betti"]["b2"] == 0
    assert res["betti"]["b3"] == 0
    assert res["betti"]["eulerPoincare"] == 1
    assert res["betti"]["isEulerConsistent"] is True


def test_enclosed_void_bubble_betti_numbers():
    """Verify hollow 2-sphere void bubble has b_0 = 1, b_1 = 0, b_2 = 1, b_3 = 0, and \chi = 2."""
    code = """
    import { buildEnclosedVoidBubbleComplex } from './src/topology/morse_smale_complex.js';

    const bubble = buildEnclosedVoidBubbleComplex();
    const summary = bubble.getTopologySummary();

    console.log(JSON.stringify(summary));
    """
    res = run_node_snippet(code)
    assert res["nilpotencyVerified"] is True
    assert res["bettiNumbers"]["b0"] == 1
    assert res["bettiNumbers"]["b1"] == 0
    assert res["bettiNumbers"]["b2"] == 1 # Enclosed 2D cavity!
    assert res["bettiNumbers"]["b3"] == 0
    assert res["eulerCharacteristic"] == 2 # S^2 Euler characteristic = 2


# ============================================================================
# TIER 5: Multi-Loop & Tunnel Betti Numbers (b_1 > 0)
# ============================================================================

def test_filament_loop_tunnel_topology():
    """Verify a network with 2 independent cosmic filament loops has b_1 = 2."""
    code = """
    import {
      MorseCriticalPoint,
      Morse1Cell,
      Morse2Cell,
      Morse3Cell,
      MorseSmaleComplex3D,
      MorseIndex
    } from './src/topology/morse_smale_complex.js';

    // Graph with 2 independent cycles (e.g. Figure-8 filament loop)
    // Vertices: A, B, C, D, E
    const vA = new MorseCriticalPoint({ id: 'vA', index: 3, position: [0, 0, 0], value: 2.0 });
    const vB = new MorseCriticalPoint({ id: 'vB', index: 2, position: [10, 0, 0], value: 0.5 });
    const vC = new MorseCriticalPoint({ id: 'vC', index: 3, position: [10, 10, 0], value: 2.0 });
    const vD = new MorseCriticalPoint({ id: 'vD', index: 2, position: [-10, 0, 0], value: 0.5 });
    const vE = new MorseCriticalPoint({ id: 'vE', index: 3, position: [-10, 10, 0], value: 2.0 });

    // Loop 1: A-B-C-A (3 edges)
    const eAB = new Morse1Cell({ id: 'eAB', source: vA, target: vB });
    const eBC = new Morse1Cell({ id: 'eBC', source: vB, target: vC });
    const eCA = new Morse1Cell({ id: 'eCA', source: vC, target: vA });

    // Loop 2: A-D-E-A (3 edges)
    const eAD = new Morse1Cell({ id: 'eAD', source: vA, target: vD });
    const eDE = new Morse1Cell({ id: 'eDE', source: vD, target: vE });
    const eEA = new Morse1Cell({ id: 'eEA', source: vE, target: vA });

    const msc = new MorseSmaleComplex3D({
      criticalPoints: [vA, vB, vC, vD, vE],
      cells1: [eAB, eBC, eCA, eAD, eDE, eEA],
      cells2: [], // No 2-cells filling the loops => 2 independent 1D tunnels
      cells3: []
    });

    const summary = msc.getTopologySummary();
    console.log(JSON.stringify(summary));
    """
    res = run_node_snippet(code)
    assert res["bettiNumbers"]["b0"] == 1 # 1 connected component
    assert res["bettiNumbers"]["b1"] == 2 # 2 independent filament loops (tunnels)
    assert res["bettiNumbers"]["b2"] == 0
    assert res["eulerCharacteristic"] == 1 - 2 # \chi = b_0 - b_1 = -1
    assert res["cellCounts"]["cells0"] == 5
    assert res["cellCounts"]["cells1"] == 6 # V - E = 5 - 6 = -1


# ============================================================================
# TIER 6: Topological Simplification via Morse Persistence Cancellation
# ============================================================================

def test_morse_persistence_pair_extraction_and_simplification():
    """Verify persistence pair extraction and cancellation of low-persistence noise below epsilon_p."""
    code = """
    import {
      MorseCriticalPoint,
      Morse1Cell,
      MorseIndex,
      MorsePersistenceEngine,
      MorseSmaleComplex3D
    } from './src/topology/morse_smale_complex.js';

    // True dominant void and noise void
    const vDeep = new MorseCriticalPoint({ id: 'v_deep', index: MorseIndex.MINIMUM, position: [0, 0, 0], value: -0.90 });
    const vShallow = new MorseCriticalPoint({ id: 'v_noise', index: MorseIndex.MINIMUM, position: [10, 0, 0], value: -0.15 });

    // Saddles
    const sSpurious = new MorseCriticalPoint({ id: 's_noise', index: MorseIndex.WALL_SADDLE, position: [5, 0, 0], value: -0.10 });
    const sDominant = new MorseCriticalPoint({ id: 's_wall', index: MorseIndex.WALL_SADDLE, position: [25, 0, 0], value: 0.10 });

    // Maxima
    const mCore = new MorseCriticalPoint({ id: 'm_core', index: MorseIndex.MAXIMUM, position: [50, 0, 0], value: 2.50 });

    const eNoise = new Morse1Cell({ id: 'e_noise', source: sSpurious, target: vShallow });
    const eDom1 = new Morse1Cell({ id: 'e_dom1', source: sDominant, target: vDeep });

    const msc = new MorseSmaleComplex3D({
      criticalPoints: [vDeep, vShallow, sSpurious, sDominant, mCore],
      cells1: [eNoise, eDom1]
    });

    const pairs = MorsePersistenceEngine.extractPersistencePairs(msc.criticalPoints, msc.cells1);

    // Simplify with threshold epsilon_p = 0.10 (should cancel the spurious pair with persistence 0.05)
    const simplified = msc.simplify(0.10);

    console.log(JSON.stringify({
      pairs: pairs.map(p => p.toJSON()),
      survivingPoints: simplified.criticalPoints.map(p => p.id),
      survivingPointCount: simplified.criticalPoints.length
    }));
    """
    res = run_node_snippet(code)
    # Check persistence pairs
    pairs = res["pairs"]
    assert len(pairs) >= 1
    # Noise pair |(-0.10) - (-0.15)| = 0.05
    noisePair = next(p for p in pairs if p["creatorId"] == "v_noise" or p["destroyerId"] == "s_noise")
    assert math.isclose(noisePair["persistence"], 0.05, rel_tol=1e-3)

    # Post-simplification: v_noise and s_noise should be canceled
    assert "v_noise" not in res["survivingPoints"]
    assert "s_noise" not in res["survivingPoints"]
    assert "v_deep" in res["survivingPoints"]
    assert "m_core" in res["survivingPoints"]


def test_synthetic_cosmic_web_simplification_sweep():
    """Verify persistence threshold sweep simplifies synthetic Cosmic Web while preserving major attractors."""
    code = """
    import { buildSyntheticCosmicWeb } from './src/topology/morse_smale_complex.js';

    const web = buildSyntheticCosmicWeb();
    const origCount = web.criticalPoints.length;

    // Threshold 0: no simplification
    const s0 = web.simplify(0.0);
    // Threshold 0.5: cancels subtle saddles
    const s1 = web.simplify(0.5);
    // Threshold 1.5: cancels medium structures
    const s2 = web.simplify(1.5);

    console.log(JSON.stringify({
      origCount: origCount,
      countS0: s0.criticalPoints.length,
      countS1: s1.criticalPoints.length,
      countS2: s2.criticalPoints.length,
      shapleySurvivesS1: s1.criticalPoints.some(p => p.id === 'max_shapley'),
      shapleySurvivesS2: s2.criticalPoints.some(p => p.id === 'max_shapley')
    }));
    """
    res = run_node_snippet(code)
    assert res["origCount"] == 10
    assert res["countS0"] == 10
    assert res["countS1"] <= res["origCount"]
    assert res["countS2"] <= res["countS1"]
    assert res["shapleySurvivesS1"] is True
    assert res["shapleySurvivesS2"] is True # Shapley has delta = 3.40 (highest persistence)


# ============================================================================
# TIER 7: Cosmic Web Graph & Network Analysis
# ============================================================================

def test_cosmic_web_graph_dijkstra_routing():
    """Verify Dijkstra shortest path routing along connected cosmic filaments in Supergalactic space."""
    code = """
    import { buildSyntheticCosmicWeb } from './src/topology/morse_smale_complex.js';

    const web = buildSyntheticCosmicWeb();
    const graph = web.cosmicWebGraph;

    // Route between Laniakea and Shapley through Centaurus-Shapley filament saddle
    const routeLanShapley = graph.findShortestPath('max_laniakea', 'max_shapley');
    // Route to Perseus
    const routeLanPerseus = graph.findShortestPath('max_laniakea', 'max_perseus');
    // Route to disconnected/non-adjacent component
    const routeVoid = graph.findShortestPath('max_laniakea', 'void_bootes');

    console.log(JSON.stringify({
      routeLanShapley: routeLanShapley,
      routeLanPerseus: routeLanPerseus,
      totalLength: graph.getTotalFilamentLengthMpc(),
      totalVoidVolume: graph.getTotalVoidVolumeMpc3(),
      meanTortuosity: graph.getMeanTortuosity()
    }));
    """
    res = run_node_snippet(code)
    routeLS = res["routeLanShapley"]
    assert routeLS["found"] is True
    assert "max_laniakea" in routeLS["pathNodeIds"]
    assert "saddle_fil_1" in routeLS["pathNodeIds"]
    assert "max_shapley" in routeLS["pathNodeIds"]
    assert math.isclose(routeLS["distanceMpc"], 42.5 + 48.2, rel_tol=1e-2)

    routeLP = res["routeLanPerseus"]
    assert routeLP["found"] is True
    assert math.isclose(routeLP["distanceMpc"], 55.4 + 41.0, rel_tol=1e-2)

    assert res["totalLength"] > 0
    assert res["totalVoidVolume"] > 0
    assert res["meanTortuosity"] >= 1.0


def test_cosmic_web_graph_geojson_export():
    """Verify GeoJSON FeatureCollection export containing Point nodes and LineString edges."""
    code = """
    import { buildSyntheticCosmicWeb } from './src/topology/morse_smale_complex.js';

    const web = buildSyntheticCosmicWeb();
    const geojson = web.cosmicWebGraph.toGeoJSON();

    const points = geojson.features.filter(f => f.geometry.type === 'Point');
    const lines = geojson.features.filter(f => f.geometry.type === 'LineString');

    console.log(JSON.stringify({
      type: geojson.type,
      totalFeatures: geojson.features.length,
      pointCount: points.length,
      lineCount: lines.length,
      samplePoint: points[0],
      sampleLine: lines[0]
    }));
    """
    res = run_node_snippet(code)
    assert res["type"] == "FeatureCollection"
    assert res["pointCount"] == 10
    assert res["lineCount"] == 8
    assert "morseIndex" in res["samplePoint"]["properties"]
    assert "lengthMpc" in res["sampleLine"]["properties"]


# ============================================================================
# TIER 8: 3D Discrete Morse-Smale Complex Extractor on Scalar Grids
# ============================================================================

def test_discrete_msc_extractor_synthetic_gaussian_grid():
    """Verify DiscreteMorseSmaleComplexExtractor detects local minima, maxima, and saddles on 3D grid."""
    code = """
    import { DiscreteMorseSmaleComplexExtractor } from './src/topology/morse_smale_complex.js';

    const nx = 16;
    const ny = 16;
    const nz = 16;
    const data = new Float64Array(nx * ny * nz);

    // Create 1 dominant maximum at (12, 12, 12) and 1 minimum at (4, 4, 4)
    for (let iz = 0; iz < nz; iz++) {
      for (let iy = 0; iy < ny; iy++) {
        for (let ix = 0; ix < nx; ix++) {
          const idx = (iz * ny + iy) * nx + ix;
          // Gaussian peak
          const dMax = Math.hypot(ix - 12, iy - 12, iz - 12);
          const dMin = Math.hypot(ix - 4, iy - 4, iz - 4);
          data[idx] = 2.0 * Math.exp(-dMax*dMax / 8.0) - 1.0 * Math.exp(-dMin*dMin / 8.0);
        }
      }
    }

    const extractor = new DiscreteMorseSmaleComplexExtractor({
      gridData: data,
      nx: nx,
      ny: ny,
      nz: nz,
      boundsMpc: [-100, 100, -100, 100, -100, 100]
    });

    const cps = extractor.extractCriticalPoints();
    const complex = extractor.extractComplex();

    const maxCount = cps.filter(p => p.index === 3).length;
    const minCount = cps.filter(p => p.index === 0).length;

    console.log(JSON.stringify({
      totalCPs: cps.length,
      maxCount: maxCount,
      minCount: minCount,
      cell1Count: complex.cells1.length,
      cell2Count: complex.cells2.length
    }));
    """
    res = run_node_snippet(code)
    assert res["maxCount"] >= 1
    assert res["minCount"] >= 1
    assert res["totalCPs"] >= 2
    assert res["cell1Count"] >= 0


def test_boundary_matrix_rank_and_subspaces():
    """Verify BoundaryMatrix algebraic operations, entry access, and dense array exports."""
    code = """
    import { BoundaryMatrix } from './src/topology/morse_smale_complex.js';

    const bm = new BoundaryMatrix(1, ['v0', 'v1', 'v2'], ['e0', 'e1', 'e2']);
    // Triangle boundary
    // e0: v1 - v0
    // e1: v2 - v1
    // e2: v0 - v2
    bm.setEntry('v0', 'e0', -1);
    bm.setEntry('v1', 'e0', +1);

    bm.setEntry('v1', 'e1', -1);
    bm.setEntry('v2', 'e1', +1);

    bm.setEntry('v2', 'e2', -1);
    bm.setEntry('v0', 'e2', +1);

    console.log(JSON.stringify({
      rank: bm.getRank(),
      rankZ2: bm.getRankZ2(),
      dense: bm.toDenseArray(),
      v1_e0: bm.getEntry('v1', 'e0'),
      v2_e0: bm.getEntry('v2', 'e0')
    }));
    """
    res = run_node_snippet(code)
    assert res["rank"] == 2 # 3 edges forming closed cycle has rank 2
    assert res["rankZ2"] == 2
    assert res["v1_e0"] == 1.0
    assert res["v2_e0"] == 0.0
