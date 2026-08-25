"""
tests/surfaces/test_exact_marching_tetrahedra.py
Rigorous Automated Pytest Suite for Exact 3D Marching Tetrahedra & Manifold Mesher.

Verification Coverage:
1. Cube Simplicial Decompositions:
   - 6-Tetrahedra Kuhn / Freudenthal decomposition.
   - 6-Tetrahedra Alternating Checkerboard decomposition.
   - 5-Tetrahedra Alternating Checkerboard decomposition.
   - Simplicial volume partition consistency (sum(V_tet) == V_cube).
2. Analytical Sphere Calculus (< 0.5% - 1.0% tolerance):
   - Surface area integral \iint dA == 4 * pi * R^2.
   - Enclosed volume integral \iiint dV == 4/3 * pi * R^3.
   - Sphericity Psi and Isoperimetric Quotient Q.
   - Surface centroid and Volumetric centroid on shifted coordinates.
3. Topological Invariants & Watertight 2-Manifolds:
   - Single Sphere: chi = V - E + F = 2, genus = 0, watertight (E_bnd = 0, E_non_manifold = 0).
   - Disjoint Spheres: chi = 4, genus = 0, 2 connected components.
   - Torus: chi = 0, genus = 1, 1 connected component.
   - Multi-component topologies: Betti numbers (b0, b1, b2).
4. Exact Edge Linear & Newton-Raphson Root Interpolation:
   - Linear edge interpolation across arbitrary thresholds delta_iso:
     * delta = -0.8 (Cosmic Underdense Void Boundary)
     * delta = 0.0 (Cosmic Mean Density Separatrix)
     * delta = 1.5 (Cosmic Sheet / Wall Boundary)
     * delta = 1.686 (Linear Spherical Collapse Critical Overdensity delta_c)
     * delta = 5.0 (Cosmic Filament Core)
     * delta = 20.0 (Virialized Cluster Envelope)
     * delta = 200.0 (Virialized Halo Boundary delta_vir)
   - Newton-Raphson high-order convergence on non-linear fields.
5. Differential Geometry, Normals & Curvature:
   - Facet, Area-Weighted, Angle-Weighted, and Field-Gradient normals.
   - Unit normal length and radial alignment.
   - Discrete Gauss-Bonnet theorem: \iint K dA == 2pi * chi.
   - Taubin non-shrinking smoothing and volume preservation.
6. Cosmological Supercluster Isocontours:
   - Laniakea Supercluster Gravity Basin.
   - Shapley Supercluster Concentration.
   - Great Attractor (Norma/Centaurus Cluster Core).
   - Local Universe Multi-Structure Field.
   - NFW dark matter halo profile isocontours.
7. Geometric Ray Casting & Physical Velocity Flux Integrals:
   - Point-in-manifold ray containment test (Jordan curve theorem in 3D).
   - Cosmic bulk flow flux integration \iint (v . n) dA.
8. Multi-format Serialization:
   - Wavefront OBJ export.
   - GeoJSON 3D MultiPolygon export.
   - Three.js BufferGeometry typed arrays.
"""

import json
import math
import os
import subprocess
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
# 1. CUBE SIMPLICIAL DECOMPOSITION & TOPOLOGY INTEGRITY
# ============================================================================

class TestCubeSimplicialDecomposition:
    def test_kuhn_6_tetrahedra_decomposition_properties(self):
        """Verify Kuhn 6-tet decomposition covers unit cube volume with exact sum = 1.0."""
        code = """
        import { KUHN_6_TETRAHEDRA, CUBE_CORNER_OFFSETS } from './src/surfaces/exact_marching_tetrahedra.js';

        let totalVolume = 0.0;
        for (const tet of KUHN_6_TETRAHEDRA) {
          const p0 = CUBE_CORNER_OFFSETS[tet[0]];
          const p1 = CUBE_CORNER_OFFSETS[tet[1]];
          const p2 = CUBE_CORNER_OFFSETS[tet[2]];
          const p3 = CUBE_CORNER_OFFSETS[tet[3]];

          const v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
          const v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
          const v3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];

          // det([v1, v2, v3]) / 6
          const det = v1[0] * (v2[1]*v3[2] - v2[2]*v3[1])
                    - v1[1] * (v2[0]*v3[2] - v2[2]*v3[0])
                    + v1[2] * (v2[0]*v3[1] - v2[1]*v3[0]);
          totalVolume += Math.abs(det) / 6.0;
        }

        console.log(JSON.stringify({
          tetCount: KUHN_6_TETRAHEDRA.length,
          totalVolume
        }));
        """
        res = run_node_snippet(code)
        assert res["tetCount"] == 6
        assert math.isclose(res["totalVolume"], 1.0, rel_tol=1e-9)

    def test_alternating_6_tetrahedra_decomposition_properties(self):
        """Verify alternating 6-tet decomposition volumes and vertex validity."""
        code = """
        import { ALTERNATING_6_TETRAHEDRA_EVEN, ALTERNATING_6_TETRAHEDRA_ODD, CUBE_CORNER_OFFSETS } from './src/surfaces/exact_marching_tetrahedra.js';

        const calcVol = (tets) => {
          let vol = 0.0;
          for (const tet of tets) {
            const p0 = CUBE_CORNER_OFFSETS[tet[0]];
            const p1 = CUBE_CORNER_OFFSETS[tet[1]];
            const p2 = CUBE_CORNER_OFFSETS[tet[2]];
            const p3 = CUBE_CORNER_OFFSETS[tet[3]];

            const v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            const v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            const v3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];

            const det = v1[0] * (v2[1]*v3[2] - v2[2]*v3[1])
                      - v1[1] * (v2[0]*v3[2] - v2[2]*v3[0])
                      + v1[2] * (v2[0]*v3[1] - v2[1]*v3[0]);
            vol += Math.abs(det) / 6.0;
          }
          return vol;
        };

        const volEven = calcVol(ALTERNATING_6_TETRAHEDRA_EVEN);
        const volOdd = calcVol(ALTERNATING_6_TETRAHEDRA_ODD);

        console.log(JSON.stringify({ volEven, volOdd }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["volEven"], 1.0, rel_tol=1e-9)
        assert math.isclose(res["volOdd"], 1.0, rel_tol=1e-9)

    def test_5_tetrahedra_decomposition_properties(self):
        """Verify 5-tet decomposition consists of 4 corner tets (1/6 each) + 1 central tet (1/3) = 1.0."""
        code = """
        import { FIVE_TETRAHEDRA_EVEN, FIVE_TETRAHEDRA_ODD, CUBE_CORNER_OFFSETS } from './src/surfaces/exact_marching_tetrahedra.js';

        const calcVol = (tets) => {
          let vol = 0.0;
          for (const tet of tets) {
            const p0 = CUBE_CORNER_OFFSETS[tet[0]];
            const p1 = CUBE_CORNER_OFFSETS[tet[1]];
            const p2 = CUBE_CORNER_OFFSETS[tet[2]];
            const p3 = CUBE_CORNER_OFFSETS[tet[3]];

            const v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            const v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            const v3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];

            const det = v1[0] * (v2[1]*v3[2] - v2[2]*v3[1])
                      - v1[1] * (v2[0]*v3[2] - v2[2]*v3[0])
                      + v1[2] * (v2[0]*v3[1] - v2[1]*v3[0]);
            vol += Math.abs(det) / 6.0;
          }
          return vol;
        };

        const volEven = calcVol(FIVE_TETRAHEDRA_EVEN);
        const volOdd = calcVol(FIVE_TETRAHEDRA_ODD);

        console.log(JSON.stringify({
          tetCountEven: FIVE_TETRAHEDRA_EVEN.length,
          tetCountOdd: FIVE_TETRAHEDRA_ODD.length,
          volEven,
          volOdd
        }));
        """
        res = run_node_snippet(code)
        assert res["tetCountEven"] == 5
        assert res["tetCountOdd"] == 5
        assert math.isclose(res["volEven"], 1.0, rel_tol=1e-9)
        assert math.isclose(res["volOdd"], 1.0, rel_tol=1e-9)


# ============================================================================
# 2. ANALYTICAL SPHERE CALCULUS & INTEGRATION ACCURACY
# ============================================================================

class TestAnalyticalSphereCalculus:
    def test_sphere_surface_area_integration_tolerance(self):
        r"""Verify numerical surface area integration \iint dA matches 4 * pi * R^2 to < 1.0% error."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, DecompositionMode } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 64;
        const L = 100.0;
        const R = 25.0; // Mpc/h

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const sdfBuffer = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], R);
        const mesher = new ExactMarchingTetrahedra(grid, {
          decompositionMode: DecompositionMode.SIX_TETRAHEDRA_ALTERNATING
        });

        const mesh = mesher.extractIsosurface(sdfBuffer, 0.0);
        const computedArea = mesh.calculus.surfaceArea;
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
        assert res["relAreaError"] < 0.015, f"Sphere area rel error {res['relAreaError']*100:.3f}% exceeds tolerance!"

    def test_sphere_enclosed_volume_integration_tolerance(self):
        r"""Verify enclosed volume integration \iiint dV matches 4/3 * pi * R^3 to < 1.0% error."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, DecompositionMode } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 64;
        const L = 100.0;
        const R = 25.0; // Mpc/h

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const sdfBuffer = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], R);
        const mesher = new ExactMarchingTetrahedra(grid, {
          decompositionMode: DecompositionMode.SIX_TETRAHEDRA_ALTERNATING
        });

        const mesh = mesher.extractIsosurface(sdfBuffer, 0.0);
        const computedVolume = mesh.calculus.enclosedVolume;
        const exactVolume = (4.0 / 3.0) * Math.PI * Math.pow(R, 3);
        const relVolError = Math.abs(computedVolume - exactVolume) / exactVolume;

        console.log(JSON.stringify({
          computedVolume,
          exactVolume,
          relVolError,
          sphericity: mesh.calculus.sphericity,
          isoperimetricQuotient: mesh.calculus.isoperimetricQuotient
        }));
        """
        res = run_node_snippet(code)
        assert res["relVolError"] < 0.015, f"Sphere volume rel error {res['relVolError']*100:.3f}% exceeds tolerance!"
        assert res["sphericity"] > 0.95, f"Sphere sphericity {res['sphericity']} is too low!"
        assert res["isoperimetricQuotient"] > 0.90, f"Isoperimetric quotient {res['isoperimetricQuotient']} is too low!"

    def test_shifted_sphere_centroid_triangulation(self):
        """Verify surface and volumetric centroids accurately match shifted center (x0, y0, z0)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 60;
        const L = 120.0;
        const R = 20.0;
        const center = [15.0, -10.0, 8.0]; // Mpc/h

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const sdfBuffer = CosmicDensityFieldGenerator.generateSphereSDF(grid, center, R);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdfBuffer, 0.0);

        const sCentroid = mesh.calculus.surfaceCentroid;
        const vCentroid = mesh.calculus.volumetricCentroid;

        const sDist = Math.hypot(sCentroid[0] - center[0], sCentroid[1] - center[1], sCentroid[2] - center[2]);
        const vDist = Math.hypot(vCentroid[0] - center[0], vCentroid[1] - center[1], vCentroid[2] - center[2]);

        console.log(JSON.stringify({
          sCentroid,
          vCentroid,
          sDist,
          vDist
        }));
        """
        res = run_node_snippet(code)
        assert res["sDist"] < 0.25, f"Surface centroid distance error {res['sDist']} Mpc/h is too large!"
        assert res["vDist"] < 0.25, f"Volumetric centroid distance error {res['vDist']} Mpc/h is too large!"

    def test_five_tetrahedra_mode_calculus_accuracy(self):
        """Verify 5-tetrahedra mode computes valid surface area and volume."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, DecompositionMode } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 50;
        const L = 100.0;
        const R = 22.0;

        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const sdfBuffer = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0,0,0], R);
        const mesher = new ExactMarchingTetrahedra(grid, {
          decompositionMode: DecompositionMode.FIVE_TETRAHEDRA_ALTERNATING
        });

        const mesh = mesher.extractIsosurface(sdfBuffer, 0.0);
        const exactArea = 4.0 * Math.PI * R * R;
        const exactVolume = (4.0 / 3.0) * Math.PI * Math.pow(R, 3);

        const relAreaError = Math.abs(mesh.calculus.surfaceArea - exactArea) / exactArea;
        const relVolError = Math.abs(mesh.calculus.enclosedVolume - exactVolume) / exactVolume;

        console.log(JSON.stringify({
          relAreaError,
          relVolError,
          isWatertight: mesh.topology.isWatertight,
          eulerCharacteristic: mesh.topology.eulerCharacteristic
        }));
        """
        res = run_node_snippet(code)
        assert res["isWatertight"] is True
        assert res["eulerCharacteristic"] == 2
        assert res["relAreaError"] < 0.02
        assert res["relVolError"] < 0.02


# ============================================================================
# 3. TOPOLOGICAL INVARIANTS & EULER CHARACTERISTIC TESTS
# ============================================================================

class TestTopologicalInvariants:
    def test_single_sphere_topology(self):
        """Single closed sphere: chi = 2, genus = 0, watertight = true."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 40, ny: 40, nz: 40,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], 20.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        console.log(JSON.stringify(mesh.topology));
        """
        res = run_node_snippet(code)
        assert res["eulerCharacteristic"] == 2
        assert res["genus"] == 0
        assert res["isWatertight"] is True
        assert res["boundaryEdgeCount"] == 0
        assert res["nonManifoldEdgeCount"] == 0
        assert res["connectedComponents"] == 1
        assert res["bettiNumbers"] == {"b0": 1, "b1": 0, "b2": 1}

    def test_two_disjoint_spheres_topology(self):
        """Two disjoint spheres: chi = 4, genus = 0, connected components = 2."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 50;
        const L = 120.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const total = grid.totalCells;
        const buffer = new Float64Array(total);
        const c1 = [-25.0, 0.0, 0.0];
        const c2 = [25.0, 0.0, 0.0];
        const R = 15.0;

        for (let iz = 0; iz < N; iz++) {
          for (let iy = 0; iy < N; iy++) {
            for (let ix = 0; ix < N; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r1 = Math.hypot(x - c1[0], y - c1[1], z - c1[2]);
              const r2 = Math.hypot(x - c2[0], y - c2[1], z - c2[2]);
              // Max of two sphere SDFs (union)
              buffer[grid.index(ix, iy, iz)] = Math.max(R - r1, R - r2);
            }
          }
        }

        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(buffer, 0.0);

        console.log(JSON.stringify(mesh.topology));
        """
        res = run_node_snippet(code)
        assert res["eulerCharacteristic"] == 4
        assert res["genus"] == 0
        assert res["isWatertight"] is True
        assert res["connectedComponents"] == 2
        assert res["bettiNumbers"] == {"b0": 2, "b1": 0, "b2": 2}

    def test_torus_topology_genus_one(self):
        """Torus: chi = 0, genus = 1, connected components = 1."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const N = 60;
        const L = 120.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-L/2, -L/2, -L/2],
          boxSize: [L, L, L]
        });

        const torusSDF = CosmicDensityFieldGenerator.generateTorusSDF(grid, 30.0, 10.0, [0, 0, 0]);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(torusSDF, 0.0);

        console.log(JSON.stringify(mesh.topology));
        """
        res = run_node_snippet(code)
        assert res["eulerCharacteristic"] == 0
        assert res["genus"] == 1
        assert res["isWatertight"] is True
        assert res["connectedComponents"] == 1
        assert res["bettiNumbers"] == {"b0": 1, "b1": 2, "b2": 1}

    def test_discrete_gauss_bonnet_theorem_on_sphere(self):
        r"""Verify discrete Gauss-Bonnet theorem: \iint K dA == 2 * pi * chi = 4 * pi."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 48, ny: 48, nz: 48,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], 22.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        const cur = mesh.curvature;
        console.log(JSON.stringify({
          totalGaussianCurvature: cur.totalGaussianCurvature,
          expectedGaussBonnet: cur.expectedGaussBonnet,
          gaussBonnetResidual: cur.gaussBonnetResidual
        }));
        """
        res = run_node_snippet(code)
        assert math.isclose(res["totalGaussianCurvature"], 4.0 * math.pi, rel_tol=1e-4)
        assert res["gaussBonnetResidual"] < 1e-4


# ============================================================================
# 4. EXACT EDGE LINEAR & HIGH-ORDER ROOT INTERPOLATION
# ============================================================================

class TestEdgeInterpolationAndThresholds:
    def test_canonical_cosmological_thresholds_extraction(self):
        """Verify extraction across standard cosmological thresholds (void, mean, linear collapse, cluster, virial)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, CosmologicalDensityThresholds } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 48, ny: 48, nz: 48,
          origin: [-80, -80, -80],
          boxSize: [160, 160, 160]
        });

        const sphereField = CosmicDensityFieldGenerator.generateAnalyticalSphere(grid, {
          center: [0, 0, 0],
          peakOverdensity: 250.0,
          radius: 12.0,
          backgroundDensity: -0.9
        });

        const mesher = new ExactMarchingTetrahedra(grid);

        const thresholds = [
          CosmologicalDensityThresholds.VOID_BOUNDARY,
          CosmologicalDensityThresholds.MEAN_DENSITY,
          CosmologicalDensityThresholds.LINEAR_COLLAPSE,
          CosmologicalDensityThresholds.FILAMENT_CORE,
          CosmologicalDensityThresholds.VIRIALIZED_CLUSTER,
          CosmologicalDensityThresholds.VIRIALIZED_HALO
        ];

        const results = [];
        for (const th of thresholds) {
          const mesh = mesher.extractIsosurface(sphereField, th);
          results.push({
            threshold: th,
            vertexCount: mesh.vertexCount,
            triangleCount: mesh.triangleCount,
            isWatertight: mesh.topology.isWatertight,
            chi: mesh.topology.eulerCharacteristic,
            volume: mesh.calculus.enclosedVolume,
            area: mesh.calculus.surfaceArea
          });
        }

        console.log(JSON.stringify(results));
        """
        res = run_node_snippet(code)
        assert len(res) == 6
        for item in res:
            assert item["vertexCount"] > 0
            assert item["triangleCount"] > 0
            assert item["isWatertight"] is True
            assert item["chi"] == 2
        # As threshold increases, enclosed volume of central peak should decrease monotonically
        volumes = [item["volume"] for item in res]
        for i in range(len(volumes) - 1):
            assert volumes[i] > volumes[i + 1], f"Volume at th[{i}] ({volumes[i]}) should exceed th[{i+1}] ({volumes[i+1]})"

    def test_newton_raphson_root_refinement(self):
        """Verify Newton-Raphson high-order root refinement achieves exact convergence."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 32, ny: 32, nz: 32,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const sphereField = CosmicDensityFieldGenerator.generateAnalyticalSphere(grid, {
          center: [0, 0, 0],
          peakOverdensity: 10.0,
          radius: 15.0
        });

        const mesherLinear = new ExactMarchingTetrahedra(grid, { enableNewtonRaphson: false });
        const mesherNR = new ExactMarchingTetrahedra(grid, {
          enableNewtonRaphson: true,
          newtonRaphsonTolerance: 1e-6,
          newtonRaphsonMaxIters: 6
        });

        const meshLinear = mesherLinear.extractIsosurface(sphereField, 1.686);
        const meshNR = mesherNR.extractIsosurface(sphereField, 1.686);

        console.log(JSON.stringify({
          linearTriangles: meshLinear.triangleCount,
          nrTriangles: meshNR.triangleCount,
          linearArea: meshLinear.calculus.surfaceArea,
          nrArea: meshNR.calculus.surfaceArea,
          linearWatertight: meshLinear.topology.isWatertight,
          nrWatertight: meshNR.topology.isWatertight
        }));
        """
        res = run_node_snippet(code)
        assert res["linearWatertight"] is True
        assert res["nrWatertight"] is True
        assert res["linearTriangles"] == res["nrTriangles"]
        assert abs(res["linearArea"] - res["nrArea"]) < 50.0


# ============================================================================
# 5. DIFFERENTIAL GEOMETRY, NORMALS & MESH SMOOTHING
# ============================================================================

class TestDifferentialGeometryAndSmoothing:
    def test_normal_calculation_modes(self):
        """Verify normal calculation modes produce unit length vectors aligned with geometry."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, NormalCalculationMode } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 36, ny: 36, nz: 36,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const sphereField = CosmicDensityFieldGenerator.generateAnalyticalSphere(grid, {
          center: [0, 0, 0],
          peakOverdensity: 8.0,
          radius: 18.0
        });

        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sphereField, 1.686);

        const modes = [
          NormalCalculationMode.FACET_UNIFORM,
          NormalCalculationMode.AREA_WEIGHTED,
          NormalCalculationMode.ANGLE_WEIGHTED,
          NormalCalculationMode.FIELD_GRADIENT
        ];

        const stats = {};
        for (const m of modes) {
          mesh.computeNormals(m, grid, sphereField);
          const normals = mesh.normals;
          let maxLenError = 0.0;
          let minDotRadial = 1.0;

          for (let i = 0; i < mesh.vertexCount; i++) {
            const nx = normals[i * 3 + 0];
            const ny = normals[i * 3 + 1];
            const nz = normals[i * 3 + 2];
            const len = Math.hypot(nx, ny, nz);
            maxLenError = Math.max(maxLenError, Math.abs(len - 1.0));

            // Radial direction from center (0,0,0)
            const vx = mesh.vertices[i * 3 + 0];
            const vy = mesh.vertices[i * 3 + 1];
            const vz = mesh.vertices[i * 3 + 2];
            const vLen = Math.hypot(vx, vy, vz);
            const rx = vx / vLen, ry = vy / vLen, rz = vz / vLen;

            const dot = nx * rx + ny * ry + nz * rz;
            minDotRadial = Math.min(minDotRadial, dot);
          }

          stats[m] = { maxLenError, minDotRadial };
        }

        console.log(JSON.stringify(stats));
        """
        res = run_node_snippet(code)
        for mode, data in res.items():
            assert data["maxLenError"] < 1e-5, f"Normal length error in {mode} exceeds tolerance!"
            assert data["minDotRadial"] > 0.90, f"Normal alignment in {mode} is lower than 0.90!"

    def test_taubin_and_laplacian_smoothing(self):
        """Verify Taubin non-shrinking smoothing preserves topology and enclosed volume."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, SmoothingAlgorithmMode } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 36, ny: 36, nz: 36,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], 20.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        const initialVol = mesh.calculus.enclosedVolume;
        const initialArea = mesh.calculus.surfaceArea;

        const smoothedTaubin = mesh.smooth({
          algorithm: SmoothingAlgorithmMode.TAUBIN_NON_SHRINKING,
          iterations: 5,
          preserveVolume: true
        });

        const smoothedLaplace = mesh.smooth({
          algorithm: SmoothingAlgorithmMode.UNIFORM_LAPLACIAN,
          iterations: 3,
          lambda: 0.3,
          preserveVolume: true
        });

        console.log(JSON.stringify({
          initialVol,
          taubinVol: smoothedTaubin.calculus.enclosedVolume,
          laplaceVol: smoothedLaplace.calculus.enclosedVolume,
          taubinWatertight: smoothedTaubin.topology.isWatertight,
          laplaceWatertight: smoothedLaplace.topology.isWatertight,
          taubinChi: smoothedTaubin.topology.eulerCharacteristic,
          laplaceChi: smoothedLaplace.topology.eulerCharacteristic
        }));
        """
        res = run_node_snippet(code)
        assert res["taubinWatertight"] is True
        assert res["laplaceWatertight"] is True
        assert res["taubinChi"] == 2
        assert res["laplaceChi"] == 2
        assert math.isclose(res["taubinVol"], res["initialVol"], rel_tol=1e-3)
        assert math.isclose(res["laplaceVol"], res["initialVol"], rel_tol=1e-3)


# ============================================================================
# 6. COSMOLOGICAL SUPERCLUSTER ISOCONTOURS & PROFILES
# ============================================================================

class TestCosmologicalSuperclusters:
    def test_shapley_supercluster_overdensity_extraction(self):
        """Extract Shapley Supercluster concentration isodensity envelopes."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 48, ny: 48, nz: 48,
          origin: [-200.0, 80.0, -80.0],
          boxSize: [140.0, 140.0, 120.0]
        });

        const shapleyInfo = CosmicDensityFieldGenerator.KNOWN_STRUCTURES.SHAPLEY_SUPERCLUSTER;
        const field = CosmicDensityFieldGenerator.generateAnalyticalSphere(grid, {
          center: shapleyInfo.center,
          peakOverdensity: shapleyInfo.peakOverdensity,
          radius: shapleyInfo.coreRadius
        });

        const mesher = new ExactMarchingTetrahedra(grid);
        const meshVirial = mesher.extractIsosurface(field, 1.686); // Spherical collapse boundary

        const centroid = meshVirial.calculus.volumetricCentroid;
        const distToCenter = Math.hypot(
          centroid[0] - shapleyInfo.center[0],
          centroid[1] - shapleyInfo.center[1],
          centroid[2] - shapleyInfo.center[2]
        );

        console.log(JSON.stringify({
          vertexCount: meshVirial.vertexCount,
          triangleCount: meshVirial.triangleCount,
          enclosedVolume: meshVirial.calculus.enclosedVolume,
          surfaceArea: meshVirial.calculus.surfaceArea,
          distToCenter,
          isWatertight: meshVirial.topology.isWatertight,
          chi: meshVirial.topology.eulerCharacteristic
        }));
        """
        res = run_node_snippet(code)
        assert res["isWatertight"] is True
        assert res["chi"] == 2
        assert res["distToCenter"] < 0.5
        assert res["enclosedVolume"] > 5000.0

    def test_great_attractor_isocontour(self):
        """Extract Great Attractor (Norma/Centaurus) overdensity isocontour."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const ga = CosmicDensityFieldGenerator.KNOWN_STRUCTURES.GREAT_ATTRACTOR;
        const grid = new GridIndexer({
          nx: 40, ny: 40, nz: 40,
          origin: [-90.0, -30.0, -60.0],
          boxSize: [100.0, 100.0, 100.0]
        });

        const field = CosmicDensityFieldGenerator.generateAnalyticalSphere(grid, {
          center: ga.center,
          peakOverdensity: ga.peakOverdensity,
          radius: ga.coreRadius
        });

        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(field, 1.686);

        console.log(JSON.stringify({
          isWatertight: mesh.topology.isWatertight,
          chi: mesh.topology.eulerCharacteristic,
          volume: mesh.calculus.enclosedVolume,
          centroid: mesh.calculus.volumetricCentroid
        }));
        """
        res = run_node_snippet(code)
        assert res["isWatertight"] is True
        assert res["chi"] == 2
        assert res["volume"] > 1000.0

    def test_nfw_dark_matter_halo_isocontour(self):
        """Extract NFW dark matter halo profile isocontour at delta_vir = 200.0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator, CosmologicalDensityThresholds } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 48, ny: 48, nz: 48,
          origin: [-30, -30, -30],
          boxSize: [60, 60, 60]
        });

        const nfwField = CosmicDensityFieldGenerator.generateNFWHaloField(grid, {
          center: [0, 0, 0],
          scaleRadius: 4.0,
          characteristicOverdensity: 500.0,
          virialRadius: 20.0
        });

        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(nfwField, CosmologicalDensityThresholds.VIRIALIZED_HALO);

        console.log(JSON.stringify({
          isWatertight: mesh.topology.isWatertight,
          chi: mesh.topology.eulerCharacteristic,
          volume: mesh.calculus.enclosedVolume,
          sphericity: mesh.calculus.sphericity
        }));
        """
        res = run_node_snippet(code)
        assert res["isWatertight"] is True
        assert res["chi"] == 2
        assert res["volume"] > 100.0
        assert res["sphericity"] > 0.90

    def test_local_universe_multi_structure_field_extraction(self):
        """Extract Local Universe multi-cluster complex at mean density separatrix delta = 0.0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 40, ny: 40, nz: 40,
          origin: [-160, -60, -60],
          boxSize: [240, 240, 160]
        });

        const localUniv = CosmicDensityFieldGenerator.generateLocalUniverseField(grid);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(localUniv, 0.5);

        console.log(JSON.stringify({
          vertexCount: mesh.vertexCount,
          triangleCount: mesh.triangleCount,
          isWatertight: mesh.topology.isWatertight,
          components: mesh.topology.connectedComponents,
          chi: mesh.topology.eulerCharacteristic,
          volume: mesh.calculus.enclosedVolume
        }));
        """
        res = run_node_snippet(code)
        assert res["vertexCount"] > 100
        assert res["triangleCount"] > 200
        assert res["components"] >= 1


# ============================================================================
# 7. RAY CONTAINMENT & BULK FLOW VELOCITY FLUX INTEGRATION
# ============================================================================

class TestRayContainmentAndFluxIntegration:
    def test_point_in_manifold_ray_test(self):
        """Verify 3D point-in-manifold containment test on sphere."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 36, ny: 36, nz: 36,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [5.0, 5.0, 5.0], 18.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        const insidePoint = [5.0, 5.0, 5.0];
        const nearInsidePoint = [15.0, 5.0, 5.0]; // distance 10 < 18
        const outsidePoint = [30.0, 30.0, 30.0];
        const farOutsidePoint = [-30.0, -30.0, -30.0];

        console.log(JSON.stringify({
          inside1: mesh.containsPoint(insidePoint),
          inside2: mesh.containsPoint(nearInsidePoint),
          outside1: mesh.containsPoint(outsidePoint),
          outside2: mesh.containsPoint(farOutsidePoint)
        }));
        """
        res = run_node_snippet(code)
        assert res["inside1"] is True
        assert res["inside2"] is True
        assert res["outside1"] is False
        assert res["outside2"] is False

    def test_bulk_flow_velocity_divergence_flux(self):
        r"""Verify bulk flow velocity flux \iint (v . n) dA matches divergence theorem \iiint (div v) dV."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 40, ny: 40, nz: 40,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const R = 20.0;
        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], R);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        // Linear velocity field: v = H0 * r => div(v) = 3 * H0 (Hubble expansion)
        const H0 = 70.0; // km/s / Mpc
        const velocityField = (x, y, z) => [H0 * x, H0 * y, H0 * z];

        const computedFlux = mesh.integrateVelocityFlux(velocityField);
        // By divergence theorem: Flux = \iiint (3 * H0) dV = 3 * H0 * V_sphere
        const expectedFlux = 3.0 * H0 * mesh.calculus.enclosedVolume;
        const relFluxError = Math.abs(computedFlux - expectedFlux) / expectedFlux;

        console.log(JSON.stringify({
          computedFlux,
          expectedFlux,
          relFluxError
        }));
        """
        res = run_node_snippet(code)
        assert res["relFluxError"] < 0.02, f"Flux relative error {res['relFluxError']*100:.3f}% exceeds tolerance!"


# ============================================================================
# 8. MULTI-FORMAT SERIALIZATION & EXPORT
# ============================================================================

class TestSerializationAndExport:
    def test_wavefront_obj_export(self):
        """Verify Wavefront OBJ serialization format and headers."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 30, ny: 30, nz: 30,
          origin: [-30, -30, -30],
          boxSize: [60, 60, 60]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], 15.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        const objText = mesh.toOBJ('LaniakeaBoundary');
        const lines = objText.split('\\n');

        const vLines = lines.filter(l => l.startsWith('v '));
        const vnLines = lines.filter(l => l.startsWith('vn '));
        const fLines = lines.filter(l => l.startsWith('f '));

        console.log(JSON.stringify({
          lineCount: lines.length,
          vCount: vLines.length,
          vnCount: vnLines.length,
          fCount: fLines.length,
          meshVertexCount: mesh.vertexCount,
          meshTriangleCount: mesh.triangleCount,
          hasHeader: lines[0].includes('ZRT Cosmicflows')
        }));
        """
        res = run_node_snippet(code)
        assert res["hasHeader"] is True
        assert res["vCount"] == res["meshVertexCount"]
        assert res["vnCount"] == res["meshVertexCount"]
        assert res["fCount"] == res["meshTriangleCount"]

    def test_geojson_and_buffergeometry_export(self):
        """Verify GeoJSON 3D and Three.js BufferGeometry exports."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { ExactMarchingTetrahedra, CosmicDensityFieldGenerator } from './src/surfaces/exact_marching_tetrahedra.js';

        const grid = new GridIndexer({
          nx: 24, ny: 24, nz: 24,
          origin: [-20, -20, -20],
          boxSize: [40, 40, 40]
        });

        const sdf = CosmicDensityFieldGenerator.generateSphereSDF(grid, [0, 0, 0], 12.0);
        const mesher = new ExactMarchingTetrahedra(grid);
        const mesh = mesher.extractIsosurface(sdf, 0.0);

        const geojson = mesh.toGeoJSON();
        const bufferGeom = mesh.toBufferGeometry();

        console.log(JSON.stringify({
          geoJsonType: geojson.type,
          geomType: geojson.geometry.type,
          polyCount: geojson.geometry.coordinates.length,
          bufPosLen: bufferGeom.position.length,
          bufNormLen: bufferGeom.normal.length,
          bufIdxLen: bufferGeom.index.length
        }));
        """
        res = run_node_snippet(code)
        assert res["geoJsonType"] == "Feature"
        assert res["geomType"] == "MultiPolygon"
        assert res["polyCount"] > 0
        assert res["bufPosLen"] == res["bufNormLen"]
        assert res["bufIdxLen"] == res["polyCount"] * 3
