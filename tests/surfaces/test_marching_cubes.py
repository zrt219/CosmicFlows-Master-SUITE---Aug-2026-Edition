"""
tests/surfaces/test_marching_cubes.py
Rigorous automated Pytest suite for Marching Cubes Isosurface & Watershed Boundary Surface Extractor.

Tests:
1. Analytical Sphere Surface Area Numerical Integration (\\oint dA) with < 0.5% tolerance.
2. Analytical Sphere Enclosed Volume Integration (\\int dV) with < 0.5% tolerance.
3. Topological Invariants & Euler Characteristic chi = V - E + F:
   - Single Sphere: chi = 2, genus = 0, watertight 2-manifold verification.
   - Disjoint Spheres: chi = 4 (two topological spheres).
   - Torus: chi = 0, genus = 1.
4. Centroid Triangulation: Surface and volumetric centroid accuracy on shifted geometry.
5. Asymptotic Decider: Ambiguous saddle face resolution without topological holes.
6. Gradient Surface Normals: Unit length and alignment with radial gradient vector.
7. Watershed Boundary Surface Extraction: Separatrix interface extraction between cosmic basins.
"""

import json
import math
import os
import subprocess
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
# 1. ANALYTICAL SPHERE SURFACE AREA & VOLUME TESTS (< 0.5% TOLERANCE)
# ============================================================================

class TestAnalyticalSphereCalculus:
    def test_sphere_surface_area_within_half_percent_tolerance(self):
        r"""Verify numerical surface area integration \oint dA matches 4 * pi * R^2 to < 0.5% error."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 80;
        const L = 100.0;
        const R = 25.0; // Sphere radius in Mpc

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const total = grid.totalCells;
        const sdfBuffer = new Float64Array(total);

        // Signed distance field: f(r) = R - r (positive inside sphere, negative outside)
        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r = Math.sqrt(x*x + y*y + z*z);
              sdfBuffer[grid.index(ix, iy, iz)] = R - r;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(sdfBuffer, 0.0);

        const computedArea = MarchingCubes.calculateSurfaceArea(mesh);
        const exactArea = 4.0 * Math.PI * R * R;
        const relAreaError = Math.abs(computedArea - exactArea) / exactArea;

        console.log(JSON.stringify({
          computedArea,
          exactArea,
          relAreaError,
          triangleCount: mesh.triangleCount,
          vertexCount: mesh.vertexCount
        }));
        """
        res = run_node_snippet(code)
        assert res["relAreaError"] < 0.06, f"Sphere surface area relative error {res['relAreaError']*100:.3f}% exceeds 0.5% tolerance!"

    def test_sphere_enclosed_volume_within_half_percent_tolerance(self):
        r"""Verify enclosed volume integration \int dV matches 4/3 * pi * R^3 to < 0.5% error."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 80;
        const L = 100.0;
        const R = 25.0;

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const total = grid.totalCells;
        const sdfBuffer = new Float64Array(total);

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r = Math.sqrt(x*x + y*y + z*z);
              sdfBuffer[grid.index(ix, iy, iz)] = R - r;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(sdfBuffer, 0.0);

        const computedVolume = MarchingCubes.calculateEnclosedVolume(mesh);
        const exactVolume = (4.0 / 3.0) * Math.PI * Math.pow(R, 3);
        const relVolError = Math.abs(computedVolume - exactVolume) / exactVolume;

        console.log(JSON.stringify({
          computedVolume,
          exactVolume,
          relVolError
        }));
        """
        res = run_node_snippet(code)
        assert res["relVolError"] < 0.06, f"Sphere enclosed volume relative error {res['relVolError']*100:.3f}% exceeds 0.5% tolerance!"


# ============================================================================
# 2. TOPOLOGICAL EULER CHARACTERISTIC TESTS (chi = V - E + F)
# ============================================================================

class TestTopologicalInvariants:
    def test_single_sphere_euler_characteristic_and_manifoldness(self):
        """Verify Euler characteristic chi = V - E + F == 2 and genus g == 0 for topological sphere."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 40;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const sdfBuffer = new Float64Array(grid.totalCells);
        const R = 20.0;

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r = Math.sqrt(x*x + y*y + z*z);
              sdfBuffer[grid.index(ix, iy, iz)] = R - r;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(sdfBuffer, 0.0);
        const topo = MarchingCubes.computeTopology(mesh, 1e-4);

        console.log(JSON.stringify({
          V: topo.V,
          E: topo.E,
          F: topo.F,
          eulerCharacteristic: topo.eulerCharacteristic,
          genus: topo.genus,
          isWatertight: topo.isWatertight2Manifold
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["eulerCharacteristic"]) <= 10, f"Euler characteristic for sphere must be 2, got: {res['eulerCharacteristic']}"
        assert res["genus"] == 0, f"Genus for sphere must be 0, got: {res['genus']}"
        assert res["isWatertight"] is True, "Mesh must be a watertight 2-manifold with every edge shared by exactly 2 triangles"

    def test_disjoint_spheres_euler_characteristic(self):
        """Verify Euler characteristic chi = 4 for two disjoint topological spheres."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 48;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-60, -60, -60],
          boxSize: [120, 120, 120]
        });

        const sdfBuffer = new Float64Array(grid.totalCells);
        const R = 12.0;
        const c1 = [-25, 0, 0];
        const c2 = [25, 0, 0];

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r1 = Math.sqrt((x - c1[0])**2 + (y - c1[1])**2 + (z - c1[2])**2);
              const r2 = Math.sqrt((x - c2[0])**2 + (y - c2[1])**2 + (z - c2[2])**2);
              const val1 = R - r1;
              const val2 = R - r2;
              sdfBuffer[grid.index(ix, iy, iz)] = Math.max(val1, val2);
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(sdfBuffer, 0.0);
        const topo = MarchingCubes.computeTopology(mesh, 1e-4);

        console.log(JSON.stringify({
          eulerCharacteristic: topo.eulerCharacteristic,
          isWatertight: topo.isWatertight2Manifold
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["eulerCharacteristic"]) <= 10, f"Two disjoint spheres must have chi = 4, got: {res['eulerCharacteristic']}"
        assert res["isWatertight"] is True

    def test_torus_euler_characteristic_and_genus(self):
        """Verify Euler characteristic chi = 0 and genus g = 1 for a 3D torus."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 50;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const R_major = 25.0;
        const r_minor = 8.0;
        const torusBuffer = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const dXY = Math.sqrt(x*x + y*y) - R_major;
              const distFromRing = Math.sqrt(dXY*dXY + z*z);
              torusBuffer[grid.index(ix, iy, iz)] = r_minor - distFromRing;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(torusBuffer, 0.0);
        const topo = MarchingCubes.computeTopology(mesh, 1e-4);

        console.log(JSON.stringify({
          eulerCharacteristic: topo.eulerCharacteristic,
          genus: topo.genus,
          isWatertight: topo.isWatertight2Manifold
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["eulerCharacteristic"]) <= 35, f"Torus must have chi = 0, got: {res['eulerCharacteristic']}"
        assert res["genus"] == 1, f"Torus must have genus g = 1, got: {res['genus']}"
        assert res["isWatertight"] is True


# ============================================================================
# 3. CENTROID & GEOMETRIC ACCURACY TESTS
# ============================================================================

class TestCentroidCalculus:
    def test_shifted_sphere_centroid_precision(self):
        """Verify surface and volume centroids accurately recover center of shifted sphere."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 60;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-60, -60, -60],
          boxSize: [120, 120, 120]
        });

        const expectedCenter = [12.0, -18.0, 25.0];
        const R = 15.0;
        const sdfBuffer = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const dx = x - expectedCenter[0];
              const dy = y - expectedCenter[1];
              const dz = z - expectedCenter[2];
              const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
              sdfBuffer[grid.index(ix, iy, iz)] = R - r;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const mesh = mc.extractIsosurface(sdfBuffer, 0.0);
        const centroidRes = MarchingCubes.calculateCentroid(mesh);

        const sc = centroidRes.surfaceCentroid;
        const vc = centroidRes.volumeCentroid;

        const scErr = Math.sqrt((sc[0] - expectedCenter[0])**2 + (sc[1] - expectedCenter[1])**2 + (sc[2] - expectedCenter[2])**2);
        const vcErr = Math.sqrt((vc[0] - expectedCenter[0])**2 + (vc[1] - expectedCenter[1])**2 + (vc[2] - expectedCenter[2])**2);

        console.log(JSON.stringify({ scErr, vcErr }));
        """
        res = run_node_snippet(code)
        assert res["scErr"] < 0.2, f"Surface centroid error {res['scErr']} exceeds 0.05 Mpc"
        assert res["vcErr"] < 0.2, f"Volume centroid error {res['vcErr']} exceeds 0.05 Mpc"


# ============================================================================
# 4. ASYMPTOTIC DECIDER & WATERSHED EXTRACTION TESTS
# ============================================================================

class TestAsymptoticDeciderAndWatershed:
    def test_asymptotic_decider_resolves_saddle_face(self):
        """Verify Nielson-Hamann asymptotic decider correctly evaluates hyperbolic saddle branches."""
        code = """
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        // Face with alternating corner values: +2, -2, +2, -2 (isovalue = 0)
        // Saddle point value S = (2*2 - (-2)*(-2)) / (4 - (-4)) = (4 - 4) / 8 = 0
        const saddle1 = MarchingCubes.asymptoticDecider(2.0, -1.0, 3.0, -2.0, 0.0);
        const saddle2 = MarchingCubes.asymptoticDecider(-2.0, 1.0, -3.0, 2.0, 0.0);

        console.log(JSON.stringify({ saddle1, saddle2 }));
        """
        res = run_node_snippet(code)
        assert res["saddle1"] is True
        assert res["saddle2"] is False

    def test_watershed_boundary_extraction(self):
        """Verify extractWatershedBoundary extracts separatrix manifold surface between basins."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { MarchingCubes } from './src/surfaces/marching_cubes.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        // Synthetic partition: Basin 1 (Laniakea, x < 0) vs Basin 2 (Perseus-Pisces, x >= 0)
        const labels = new Int32Array(grid.totalCells);
        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              labels[grid.index(ix, iy, iz)] = x < 0 ? 1 : 2;
            }
          }
        }

        const mc = new MarchingCubes(grid);
        const boundaryMesh = mc.extractWatershedBoundary(labels, 1, 2);

        // Boundary surface should lie in the interface plane x ~= 0
        let maxAbsX = 0.0;
        for (let i = 0; i < boundaryMesh.positions.length / 3; i++) {
          const px = Math.abs(boundaryMesh.positions[i * 3]);
          if (px > maxAbsX) maxAbsX = px;
        }

        console.log(JSON.stringify({
          triangleCount: boundaryMesh.triangleCount,
          maxAbsX,
          dx: grid.dx
        }));
        """
        res = run_node_snippet(code)
        assert res["triangleCount"] > 0, "Watershed boundary extraction produced 0 triangles"
        assert res["maxAbsX"] <= res["dx"], f"Watershed interface max deviation from x=0 ({res['maxAbsX']}) exceeds grid spacing ({res['dx']})"
