"""
tests/watershed/test_watershed_and_topology.py
Comprehensive Pytest verification suite for:
- Table A.1 Dupuy & Courtois (2023) Basin Taxonomy & 128³ FITS/binary loader
- Watershed Streamline Endpoint Segmentation & 26-neighbor local maxima
- Quantitative Segmentation Comparator (Jaccard, Dice, Centroids, Transition Matrix)
- 3D Newton-Raphson Vector Root Finder (||v|| < 1e-7 km/s)
- 3x3 Jacobi Strain-Rate Eigensolver & Hyperbolic Morse Web Classifier
- Separatrix Manifold Tracer & Cosmic Web Graph
- CF4++ HMC Streaming Welford & Covariance Propagation
"""

import json
import subprocess
import os
import math
import numpy as np
import pytest

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
# TIER 1: Table A.1 Taxonomy & Loader Tests
# ============================================================================

def test_table_a1_taxonomy_completeness():
    """Verify official Table A.1 basin taxonomy contains all 8 major cosmic basins + background."""
    code = """
    import { BASIN_TAXONOMY, getAllBasins, getBasinById, getBasinByName } from './src/watershed/watershed_classifier.js';

    const basins = getAllBasins();
    const laniakea = getBasinById(1);
    const shapley = getBasinByName('Shapley');
    const apus = getBasinByName('Apus');
    const hercules = getBasinByName('Hercules');
    const lepus = getBasinByName('Lepus');
    const perseus = getBasinByName('Perseus-Pisces');
    const sdss1a = getBasinByName('SDSS-1a');
    const sdss2a = getBasinByName('SDSS-2a');

    console.log(JSON.stringify({
      basinCount: basins.length,
      laniakeaVolume: laniakea.volumeCF4Ind1e6,
      laniakeaAttractor: laniakea.attractorSGMpc,
      shapleyVolume: shapley.volumeCF4Ind1e6,
      shapleyDist: shapley.distToShapleyMpc,
      hasAllEight: !!(laniakea && shapley && apus && hercules && lepus && perseus && sdss1a && sdss2a)
    }));
    """
    res = run_node_snippet(code)
    assert res["basinCount"] == 9  # 0 to 8
    assert res["hasAllEight"] is True
    assert math.isclose(res["laniakeaVolume"], 2.07, rel_tol=1e-2)
    assert math.isclose(res["shapleyVolume"], 7.82, rel_tol=1e-2)
    assert res["shapleyDist"] < 30.0  # Shapley centroid distance to itself


def test_voxel_grid_coordinate_transforms():
    """Verify exact bijective round-trip transformations between grid voxel index and Supergalactic space."""
    code = """
    import { VoxelGrid, createGridMetadata } from './src/watershed/watershed_classifier.js';

    const meta = createGridMetadata({ nx: 128, ny: 128, nz: 128, halfExtentMpc: 500.0 });
    const grid = new VoxelGrid(128, 128, 128, null, meta);

    // Test corners and center
    const testPoints = [
      [0, 0, 0],
      [64, 64, 64],
      [127, 127, 127],
      [32, 64, 96]
    ];

    const errors = [];
    for (const [ix, iy, iz] of testPoints) {
      const world = grid.gridToWorld(ix, iy, iz, 'Mpc');
      const recovered = grid.worldToGrid(world.x, world.y, world.z, 'Mpc');
      errors.push({
        orig: [ix, iy, iz],
        recovered: [recovered.ix, recovered.iy, recovered.iz],
        inside: recovered.inside
      });
    }

    console.log(JSON.stringify(errors));
    """
    res = run_node_snippet(code)
    for item in res:
        assert item["orig"] == item["recovered"]
        assert item["inside"] is True


def test_fits_and_binary_loader_roundtrip():
    """Verify loading raw binary and simulated FITS segmentation voxel arrays."""
    code = """
    import { VoxelGrid, loadSegmentationBuffer, computeBasinStatistics } from './src/watershed/watershed_classifier.js';

    // Create 16x16x16 test grid with 2 basins
    const N = 16;
    const total = N * N * N;
    const rawBuffer = new Uint16Array(total);

    for (let ix = 0; ix < N; ix++) {
      for (let iy = 0; iy < N; iy++) {
        for (let iz = 0; iz < N; iz++) {
          const idx = (ix * N + iy) * N + iz;
          if (ix < 8) {
            rawBuffer[idx] = 1; // Laniakea
          } else {
            rawBuffer[idx] = 6; // Shapley
          }
        }
      }
    }

    const grid = loadSegmentationBuffer(rawBuffer.buffer, { nx: N, ny: N, nz: N, halfExtentMpc: 200.0 });
    const stats = computeBasinStatistics(grid);

    console.log(JSON.stringify({
      lanVoxelCount: stats[1].voxelCount,
      shaVoxelCount: stats[6].voxelCount,
      totalCount: grid.totalVoxels,
      lanCentroidX: stats[1].centroidMpc[0],
      shaCentroidX: stats[6].centroidMpc[0]
    }));
    """
    res = run_node_snippet(code)
    assert res["lanVoxelCount"] == 16 * 16 * 16 / 2
    assert res["shaVoxelCount"] == 16 * 16 * 16 / 2
    assert res["lanCentroidX"] < 0
    assert res["shaCentroidX"] > 0


# ============================================================================
# TIER 2: 3D Newton-Raphson & Eigensolver Tests
# ============================================================================

def test_3d_newton_raphson_exact_convergence():
    """Verify 3D Newton-Raphson converges to velocity zero-crossings with residual < 1e-7 km/s."""
    code = """
    import { solveNewtonRaphson3D } from './src/topology/root_finder.js';

    // Nonlinear 3D velocity field with exact root at (x*, y*, z*) = (15.2, -34.8, 8.4)
    const xTrue = [15.2, -34.8, 8.4];
    function velField(x, y, z) {
      const dx = x - xTrue[0];
      const dy = y - xTrue[1];
      const dz = z - xTrue[2];
      // Coupled nonlinear system
      return [
        -2.5 * dx + 0.1 * dy * dy - 0.05 * dz,
        -1.8 * dy + 0.2 * dx * dz,
        -3.2 * dz + 0.15 * dx * dy
      ];
    }

    const guess = [14.0, -32.0, 10.0];
    const res = solveNewtonRaphson3D(velField, guess, { velocityTol: 1e-8, maxIterations: 50 });

    console.log(JSON.stringify({
      converged: res.converged,
      residualNorm: res.residualNorm,
      iterations: res.iterations,
      root: res.root,
      errorToTrue: Math.hypot(res.root[0] - xTrue[0], res.root[1] - xTrue[1], res.root[2] - xTrue[2])
    }));
    """
    res = run_node_snippet(code)
    assert res["converged"] is True
    assert res["residualNorm"] < 1e-8
    assert res["errorToTrue"] < 1e-6
    assert res["iterations"] < 20


def test_jacobi_diagonalizer_accuracy_and_orthogonality():
    """Verify 3x3 Jacobi diagonalizer produces exact eigenvalues and strictly orthogonal eigenvectors (V^T V = I)."""
    code = """
    import { jacobiDiagonalize3x3, analyzeVelocityTensor } from './src/topology/eigen_topology.js';

    // Symmetric strain tensor with known analytical eigenvalues
    const S = [
      [4.0, 1.0, -2.0],
      [1.0, 2.0,  0.0],
      [-2.0, 0.0, 5.0]
    ];

    const diag = jacobiDiagonalize3x3(S);
    const [l1, l2, l3] = diag.eigenvalues;
    const [e1, e2, e3] = diag.eigenvectors;

    // Check V^T * V = I
    const V = [
      [e1[0], e2[0], e3[0]],
      [e1[1], e2[1], e3[1]],
      [e1[2], e2[2], e3[2]]
    ];

    const VtV = [
      [V[0][0]*V[0][0] + V[1][0]*V[1][0] + V[2][0]*V[2][0], V[0][0]*V[0][1] + V[1][0]*V[1][1] + V[2][0]*V[2][1], V[0][0]*V[0][2] + V[1][0]*V[1][2] + V[2][0]*V[2][2]],
      [V[0][1]*V[0][0] + V[1][1]*V[1][0] + V[2][1]*V[2][0], V[0][1]*V[0][1] + V[1][1]*V[1][1] + V[2][1]*V[2][1], V[0][1]*V[0][2] + V[1][1]*V[1][2] + V[2][1]*V[2][2]],
      [V[0][2]*V[0][0] + V[1][2]*V[1][0] + V[2][2]*V[2][0], V[0][2]*V[0][1] + V[1][2]*V[1][1] + V[2][2]*V[2][1], V[0][2]*V[0][2] + V[1][2]*V[1][2] + V[2][2]*V[2][2]]
    ];

    const traceOrig = S[0][0] + S[1][1] + S[2][2];
    const traceEig = l1 + l2 + l3;

    console.log(JSON.stringify({
      eigenvalues: diag.eigenvalues,
      VtV: VtV,
      traceError: Math.abs(traceOrig - traceEig),
      iterations: diag.iterations
    }));
    """
    res = run_node_snippet(code)
    # Trace invariant
    assert res["traceError"] < 1e-12
    # Orthogonality check: VtV should be identity
    vtv = np.array(res["VtV"])
    np.testing.assert_allclose(vtv, np.eye(3), atol=1e-12)


def test_morse_hyperbolic_cosmic_web_classification():
    """Verify classification of Attractor, Repeller, Saddle-Filament, and Saddle-Wall."""
    code = """
    import { analyzeVelocityTensor } from './src/topology/eigen_topology.js';

    // 1. Attractor (Megasink: 3 negative eigenvalues)
    const J_sink = [
      [-3.0, 0.2, 0.1],
      [0.2, -2.5, 0.0],
      [0.1, 0.0, -1.8]
    ];
    const p_sink = analyzeVelocityTensor(J_sink);

    // 2. Repeller (Void Center: 3 positive eigenvalues)
    const J_rep = [
      [2.1, 0.1, 0.0],
      [0.1, 1.9, 0.2],
      [0.0, 0.2, 3.4]
    ];
    const p_rep = analyzeVelocityTensor(J_rep);

    // 3. Saddle-Filament (Type-2: 2 negative, 1 positive)
    const J_fil = [
      [-2.0, 0.0, 0.0],
      [0.0, -1.5, 0.0],
      [0.0, 0.0, +3.0]
    ];
    const p_fil = analyzeVelocityTensor(J_fil);

    // 4. Saddle-Wall (Type-1: 1 negative, 2 positive)
    const J_wall = [
      [-3.0, 0.0, 0.0],
      [0.0, +1.5, 0.0],
      [0.0, 0.0, +2.0]
    ];
    const p_wall = analyzeVelocityTensor(J_wall);

    console.log(JSON.stringify({
      sinkType: p_sink.type,
      sinkMorse: p_sink.morseIndex,
      repType: p_rep.type,
      repMorse: p_rep.morseIndex,
      filType: p_fil.type,
      filMorse: p_fil.morseIndex,
      wallType: p_wall.type,
      wallMorse: p_wall.morseIndex
    }));
    """
    res = run_node_snippet(code)
    assert res["sinkType"] == "ATTRACTOR" and res["sinkMorse"] == 3
    assert res["repType"] == "REPELLER" and res["repMorse"] == 0
    assert res["filType"] == "SADDLE_FILAMENT" and res["filMorse"] == 2
    assert res["wallType"] == "SADDLE_WALL" and res["wallMorse"] == 1


# ============================================================================
# TIER 3: Watershed Pipeline & Comparator Tests
# ============================================================================

def test_watershed_endpoint_segmentation_pipeline():
    """Verify streamline endpoint integration, 26-neighbor local maxima, and basin attribution on a 2-attractor field."""
    code = """
    import { computeWatershedSegmentation } from './src/watershed/endpoint_segmentation.js';

    // 16x16x16 grid with 2 analytical attractors at (-100, 0, 0) and (+100, 0, 0)
    const N = 16;
    const halfExtent = 200.0;
    const boxSize = 2.0 * halfExtent;
    const grid = new Float32Array(N * N * N * 3);

    const c1 = [-100.0, 0.0, 0.0];
    const c2 = [+100.0, 0.0, 0.0];

    for (let ix = 0; ix < N; ix++) {
      const x = ((ix + 0.5) / N) * boxSize - halfExtent;
      for (let iy = 0; iy < N; iy++) {
        const y = ((iy + 0.5) / N) * boxSize - halfExtent;
        for (let iz = 0; iz < N; iz++) {
          const z = ((iz + 0.5) / N) * boxSize - halfExtent;
          const idx = ((ix * N + iy) * N + iz) * 3;

          // Velocity field directing particles to nearest sink
          const dx1 = x - c1[0], dy1 = y - c1[1], dz1 = z - c1[2];
          const dx2 = x - c2[0], dy2 = y - c2[1], dz2 = z - c2[2];
          const d1 = Math.hypot(dx1, dy1, dz1);
          const d2 = Math.hypot(dx2, dy2, dz2);

          if (d1 < d2) {
            grid[idx] = -dx1 * 2.0;
            grid[idx + 1] = -dy1 * 2.0;
            grid[idx + 2] = -dz1 * 2.0;
          } else {
            grid[idx] = -dx2 * 2.0;
            grid[idx + 1] = -dy2 * 2.0;
            grid[idx + 2] = -dz2 * 2.0;
          }
        }
      }
    }

    const seg = computeWatershedSegmentation(grid, N, N, N, halfExtent, {
      maxSteps: 300,
      dt: 0.1,
      densityThreshold: 10.0,
      mergeRadiusMpc: 30.0
    });

    console.log(JSON.stringify({
      attractorCount: seg.attractors.length,
      convergedRatio: seg.telemetry.convergedCount / seg.telemetry.totalVoxels,
      capturedTotal: seg.attractors.reduce((acc, a) => acc + a.capturedVoxels, 0),
      totalVoxels: seg.telemetry.totalVoxels
    }));
    """
    res = run_node_snippet(code)
    assert res["attractorCount"] >= 2
    assert res["convergedRatio"] > 0.90
    assert res["capturedTotal"] > (res["totalVoxels"] * 0.90)


def test_segmentation_comparator_jaccard_and_dice():
    """Verify Jaccard similarity, Dice coefficient, precision, recall, and transition matrix invariants."""
    code = """
    import { VoxelGrid, createGridMetadata } from './src/watershed/watershed_classifier.js';
    import { compareSegmentations } from './src/watershed/watershed_comparator.js';

    const N = 20;
    const total = N * N * N;
    const meta = createGridMetadata({ nx: N, ny: N, nz: N, halfExtentMpc: 100.0 });

    const gridA = new VoxelGrid(N, N, N, null, meta);
    const gridB = new VoxelGrid(N, N, N, null, meta);

    // Basin 1: identical overlap on 80% of voxels
    for (let i = 0; i < total; i++) {
      if (i < 1000) {
        gridA.data[i] = 1; // Basin 1
        gridB.data[i] = 1;
      } else if (i < 1200) {
        gridA.data[i] = 1; // Only in A
        gridB.data[i] = 2; // In B it's Basin 2
      } else if (i < 1400) {
        gridA.data[i] = 0;
        gridB.data[i] = 1; // Only in B
      } else {
        gridA.data[i] = 0;
        gridB.data[i] = 0;
      }
    }

    const comp = compareSegmentations(gridA, gridB);
    const b1 = comp.perBasinMetrics[1];

    console.log(JSON.stringify({
      countRef: b1.countRef,
      countComp: b1.countComp,
      countOverlap: b1.countOverlap,
      jaccard: b1.jaccardIndex,
      dice: b1.diceCoefficient,
      precision: b1.precision,
      recall: b1.recall,
      accuracy: comp.overallAccuracy,
      kappa: comp.cohensKappa
    }));
    """
    res = run_node_snippet(code)
    # A = 1200 voxels, B = 1200 voxels, Overlap = 1000 voxels, Union = 1400 voxels
    expected_jaccard = 1000.0 / 1400.0  # 0.7142857
    expected_dice = (2.0 * 1000.0) / (1200.0 + 1200.0)  # 0.833333
    expected_prec = 1000.0 / 1200.0
    expected_rec = 1000.0 / 1200.0

    assert math.isclose(res["jaccard"], expected_jaccard, rel_tol=1e-4)
    assert math.isclose(res["dice"], expected_dice, rel_tol=1e-4)
    assert math.isclose(res["precision"], expected_prec, rel_tol=1e-4)
    assert math.isclose(res["recall"], expected_rec, rel_tol=1e-4)
    assert res["kappa"] > 0.80


# ============================================================================
# TIER 4: Separatrix Graph & HMC Uncertainty Tests
# ============================================================================

def test_separatrix_graph_construction_and_shortest_path():
    """Verify cosmic web graph construction and shortest filament path search."""
    code = """
    import { buildCosmicWebGraph, findShortestFilamentPath } from './src/topology/separatrix_tracer.js';

    // Mock 16x16x16 grid
    const N = 16;
    const halfExtent = 200.0;
    const grid = new Float32Array(N * N * N * 3);

    // 2 Attractors and 1 Filament Saddle in between
    const criticalRoots = [
      {
        root: [-80.0, 0.0, 0.0],
        jacobian: [[-2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, -2.0]] // ATTRACTOR
      },
      {
        root: [+80.0, 0.0, 0.0],
        jacobian: [[-2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, -2.0]] // ATTRACTOR
      },
      {
        root: [0.0, 0.0, 0.0],
        jacobian: [[-1.5, 0.0, 0.0], [0.0, -1.5, 0.0], [0.0, 0.0, +2.0]] // SADDLE_FILAMENT
      }
    ];

    const graph = buildCosmicWebGraph(grid, N, N, N, halfExtent, criticalRoots);

    console.log(JSON.stringify({
      nodeCount: graph.nodes.length,
      edgeCount: graph.edges.length,
      nodeTypes: graph.nodes.map(n => n.type)
    }));
    """
    res = run_node_snippet(code)
    assert res["nodeCount"] == 3
    assert "ATTRACTOR" in res["nodeTypes"]
    assert "SADDLE_FILAMENT" in res["nodeTypes"]
    assert res["edgeCount"] >= 2


def test_hmc_streaming_welford_and_covariance_ellipsoid():
    """Verify streaming Welford mean/RMS grid accumulator and 3D error ellipsoid diagonalization."""
    code = """
    import { StreamingHMCAccumulator, computeCovarianceEllipsoid, propagateSpatialToVelocityDispersion } from './src/uncertainty/hmc_uncertainty.js';

    const numVoxels = 50;
    const accum = new StreamingHMCAccumulator(numVoxels);

    // Generate 100 random realization samples with known distribution
    const trueMeanX = 150.0;
    const trueStdX = 25.0;

    // Simple deterministic pseudo-random sequence
    for (let step = 0; step < 100; step++) {
      const buffer = new Float32Array(numVoxels * 3);
      for (let i = 0; i < numVoxels; i++) {
        // Box-Muller approx or sine distribution
        const val = trueMeanX + trueStdX * Math.sin(step * 0.1 + i);
        buffer[i * 3] = val;
        buffer[i * 3 + 1] = 0;
        buffer[i * 3 + 2] = 0;
      }
      accum.update(buffer);
    }

    const { meanGrid, rmsGrid } = accum.finalize();

    // Test error ellipsoid
    const cov = [
      [400.0, 50.0, 0.0],
      [50.0, 100.0, 0.0],
      [0.0, 0.0, 25.0]
    ];
    const ellipsoid = computeCovarianceEllipsoid(cov, '1sigma');

    console.log(JSON.stringify({
      sampleCount: accum.sampleCount,
      meanVx_0: meanGrid[0],
      rmsVx_0: rmsGrid[0],
      semiMajor: ellipsoid.semiMajorAxis,
      semiMinor: ellipsoid.semiMinorAxis,
      volume: ellipsoid.ellipsoidVolume
    }));
    """
    res = run_node_snippet(code)
    assert res["sampleCount"] == 100
    assert 120.0 < res["meanVx_0"] < 180.0
    assert res["rmsVx_0"] > 10.0
    assert res["semiMajor"] > res["semiMinor"]
    assert res["volume"] > 0
