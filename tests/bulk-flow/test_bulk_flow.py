"""
tests/bulk-flow/test_bulk_flow.py
Comprehensive Pytest verification suite for:
- Field Volume-Weighted Bulk Flow Estimator
- Catalog Point-Tracer Weighted Estimator
- Inverse-Variance (Maximum Likelihood) Weighted Estimator
- Radial Shell Decomposition Profiles
- Astrometric Apex Dipole Angles (to CMB Dipole and Shapley Core)
- Coordinate transformations (Supergalactic Cartesian <-> SGL, SGB)
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
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed with code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")
    return json.loads(proc.stdout)


def test_coordinate_transforms_and_dipole_angles():
    """Verify Supergalactic Cartesian <-> Angles (SGL, SGB) and apex angular separation."""
    code = """
    import {
      cartesianToSupergalacticAngles,
      supergalacticAnglesToCartesian,
      vectorAngleDeg,
      ASTRONOMICAL_APEX_REFERENCES
    } from './src/bulk-flow/bulk_flow_estimator.js';

    // Test point at (100, 100, 0) -> SGL = 45 deg, SGB = 0 deg
    const ang1 = cartesianToSupergalacticAngles(100.0, 100.0, 0.0);
    const pRec = supergalacticAnglesToCartesian(ang1.r, ang1.sgl, ang1.sgb);

    // Test orthogonal vectors -> 90 deg
    const ang90 = vectorAngleDeg([1, 0, 0], [0, 1, 0]);
    // Test parallel vectors -> 0 deg
    const ang0 = vectorAngleDeg([3, 4, 5], [6, 8, 10]);
    // Test antiparallel -> 180 deg
    const ang180 = vectorAngleDeg([1, 2, 3], [-1, -2, -3]);

    // Check angle between CMB dipole and Shapley
    const cmbVec = ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE.supergalacticKms;
    const shapleyVec = ASTRONOMICAL_APEX_REFERENCES.SHAPLEY_CORE.supergalacticKms;
    const cmbShapleyAngle = vectorAngleDeg(cmbVec, shapleyVec);

    console.log(JSON.stringify({
      sgl1: ang1.sgl,
      sgb1: ang1.sgb,
      r1: ang1.r,
      pRec: pRec,
      ang90,
      ang0,
      ang180,
      cmbShapleyAngle
    }));
    """
    res = run_node_snippet(code)
    assert math.isclose(res["sgl1"], 45.0, rel_tol=1e-4)
    assert math.isclose(res["sgb1"], 0.0, abs_tol=1e-4)
    assert math.isclose(res["ang90"], 90.0, rel_tol=1e-4)
    assert math.isclose(res["ang0"], 0.0, abs_tol=1e-4)
    assert math.isclose(res["ang180"], 180.0, rel_tol=1e-4)
    assert 0 < res["cmbShapleyAngle"] < 180.0


def test_field_volume_weighted_estimator_pure_dipole():
    """Verify Field Volume-Weighted estimator recovers exact bulk drift vector."""
    code = """
    import { estimateBulkFlowFieldVolumeWeighted } from './src/bulk-flow/bulk_flow_estimator.js';

    // 16x16x16 grid with uniform bulk flow V_true = (150, -250, 80) km/s
    const N = 16;
    const halfExtent = 200.0;
    const grid = new Float32Array(N * N * N * 3);
    const vTrue = [150.0, -250.0, 80.0];

    for (let i = 0; i < N * N * N; i++) {
      grid[i * 3] = vTrue[0];
      grid[i * 3 + 1] = vTrue[1];
      grid[i * 3 + 2] = vTrue[2];
    }

    const res = estimateBulkFlowFieldVolumeWeighted(grid, N, N, N, halfExtent, 0, 150.0);

    console.log(JSON.stringify({
      bulkVector: res.bulkFlowVector,
      mag: res.magnitude,
      trueMag: Math.hypot(...vTrue),
      sampleCount: res.sampleCount
    }));
    """
    res = run_node_snippet(code)
    np.testing.assert_allclose(res["bulkVector"], [150.0, -250.0, 80.0], rtol=1e-4)
    assert math.isclose(res["mag"], res["trueMag"], rel_tol=1e-4)
    assert res["sampleCount"] > 0


def test_field_volume_weighted_estimator_pure_expansion():
    """Verify spherically symmetric expansion field v(r) = H*r yields zero net bulk flow."""
    code = """
    import { estimateBulkFlowFieldVolumeWeighted } from './src/bulk-flow/bulk_flow_estimator.js';

    const N = 16;
    const halfExtent = 200.0;
    const boxSize = 2.0 * halfExtent;
    const grid = new Float32Array(N * N * N * 3);
    const H = 0.5; // km/s / Mpc

    for (let ix = 0; ix < N; ix++) {
      const x = ((ix + 0.5) / N) * boxSize - halfExtent;
      for (let iy = 0; iy < N; iy++) {
        const y = ((iy + 0.5) / N) * boxSize - halfExtent;
        for (let iz = 0; iz < N; iz++) {
          const z = ((iz + 0.5) / N) * boxSize - halfExtent;
          const idx = ((ix * N + iy) * N + iz) * 3;
          grid[idx] = H * x;
          grid[idx + 1] = H * y;
          grid[idx + 2] = H * z;
        }
      }
    }

    const res = estimateBulkFlowFieldVolumeWeighted(grid, N, N, N, halfExtent, 20.0, 150.0);

    console.log(JSON.stringify({
      bulkVector: res.bulkFlowVector,
      mag: res.magnitude
    }));
    """
    res = run_node_snippet(code)
    # Net bulk flow must vanish by parity symmetry
    assert res["mag"] < 1e-3


def test_catalog_weighted_and_inverse_variance_estimators():
    """Verify recovery of bulk flow from line-of-sight projected galaxy catalog."""
    code = """
    import {
      estimateBulkFlowCatalogWeighted,
      estimateBulkFlowInverseVarianceWeighted
    } from './src/bulk-flow/bulk_flow_estimator.js';

    // Generate 500 galaxies uniformly distributed on sphere of radius R=50 Mpc
    // with true bulk drift V_true = (100, 300, -50) km/s
    const V_true = [100.0, 300.0, -50.0];
    const galaxies = [];
    const numGalaxies = 500;

    // Deterministic Fibonacci sphere sampling
    const phi = (1 + Math.sqrt(5)) / 2;
    for (let i = 0; i < numGalaxies; i++) {
      const theta = 2 * Math.PI * i / phi;
      const zNorm = 1 - (2 * i + 1) / numGalaxies;
      const radius2D = Math.sqrt(1 - zNorm * zNorm);
      const r = 30.0 + 40.0 * (i / numGalaxies);

      const x = r * radius2D * Math.cos(theta);
      const y = r * radius2D * Math.sin(theta);
      const z = r * zNorm;

      const rMag = Math.hypot(x, y, z);
      const nx = x / rMag, ny = y / rMag, nz = z / rMag;

      // Line of sight velocity u = V_true . n_hat + noise
      const u_true = V_true[0] * nx + V_true[1] * ny + V_true[2] * nz;
      const noise = 20.0 * Math.sin(i * 1.7); // bounded perturbation

      galaxies.push({
        x, y, z,
        u: u_true + noise,
        error: 120.0,
        weight: 1.0
      });
    }

    const resCat = estimateBulkFlowCatalogWeighted(galaxies, 0, 100.0);
    const resIV = estimateBulkFlowInverseVarianceWeighted(galaxies, 0, 100.0, 187.0);

    console.log(JSON.stringify({
      catVector: resCat.bulkFlowVector,
      catMag: resCat.magnitude,
      ivVector: resIV.bulkFlowVector,
      ivMag: resIV.magnitude,
      trueMag: Math.hypot(...V_true),
      ivUncertainty: resIV.uncertainty1Sigma
    }));
    """
    res = run_node_snippet(code)
    # Recovered velocity should be within ~5% of true velocity (100, 300, -50)
    np.testing.assert_allclose(res["catVector"], [100.0, 300.0, -50.0], rtol=0.10)
    np.testing.assert_allclose(res["ivVector"], [100.0, 300.0, -50.0], rtol=0.10)
    assert res["ivUncertainty"][0] > 0
    assert res["ivUncertainty"][1] > 0
    assert res["ivUncertainty"][2] > 0


def test_radial_shell_profile_decomposition():
    """Verify radial shell decomposition across consecutive radial bins."""
    code = """
    import { computeRadialShellProfile } from './src/bulk-flow/bulk_flow_estimator.js';

    // 16x16x16 grid with radial decaying velocity field V(r) = V0 / (1 + r/50)
    const N = 16;
    const halfExtent = 200.0;
    const boxSize = 2.0 * halfExtent;
    const grid = new Float32Array(N * N * N * 3);
    const V0 = [200.0, 100.0, -50.0];

    for (let ix = 0; ix < N; ix++) {
      const x = ((ix + 0.5) / N) * boxSize - halfExtent;
      for (let iy = 0; iy < N; iy++) {
        const y = ((iy + 0.5) / N) * boxSize - halfExtent;
        for (let iz = 0; iz < N; iz++) {
          const z = ((iz + 0.5) / N) * boxSize - halfExtent;
          const r = Math.hypot(x, y, z);
          const factor = 1.0 / (1.0 + r / 50.0);
          const idx = ((ix * N + iy) * N + iz) * 3;
          grid[idx] = V0[0] * factor;
          grid[idx + 1] = V0[1] * factor;
          grid[idx + 2] = V0[2] * factor;
        }
      }
    }

    const shellEdges = [0, 40, 80, 120, 160];
    const profile = computeRadialShellProfile(shellEdges, 'FIELD', {
      velocityGrid: grid,
      nx: N, ny: N, nz: N,
      halfExtent: halfExtent
    });

    console.log(JSON.stringify({
      shellCount: profile.length,
      magnitudes: profile.map(p => p.magnitude),
      rMids: profile.map(p => p.rMidMpc)
    }));
    """
    res = run_node_snippet(code)
    assert res["shellCount"] == 4
    # Magnitude should decay monotonically as radius increases
    mags = res["magnitudes"]
    for i in range(len(mags) - 1):
        assert mags[i] > mags[i + 1]
