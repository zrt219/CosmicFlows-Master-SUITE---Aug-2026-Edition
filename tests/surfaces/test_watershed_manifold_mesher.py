# -*- coding: utf-8 -*-
"""
Automated Pytest Suite for Watershed Manifold Boundary Meshing and Inter-Basin Hydrodynamic Flux Integrals.
"""

import pytest
import numpy as np
from tests.utils import run_node_snippet

class TestWatershedManifoldMesher:
    def test_mesher_instantiation(self):
        res = run_node_snippet("""
        import { WatershedManifoldMesher } from './src/surfaces/watershed_manifold_mesher.js';
        const mesher = new WatershedManifoldMesher({ boxSizeMpc: 500.0 });
        console.log(JSON.stringify({
            hasMesher: !!mesher,
            boxSize: mesher.boxSizeMpc
        }));
        """)
        assert res["hasMesher"] is True
        assert res["boxSize"] == 500.0

    def test_boundary_surface_area_and_flux_planar_slab(self):
        res = run_node_snippet("""
        import { WatershedManifoldMesher } from './src/surfaces/watershed_manifold_mesher.js';
        const N = 16;
        const grid = new Uint16Array(N * N * N);
        // Basin 1 on left half (x < 8), Basin 2 on right half (x >= 8)
        for (let iz = 0; iz < N; iz++) {
            for (let iy = 0; iy < N; iy++) {
                for (let ix = 0; ix < N; ix++) {
                    const idx = ix + iy * N + iz * N * N;
                    grid[idx] = (ix < N / 2) ? 1 : 2;
                }
            }
        }

        const mockVelocity = {
            evaluate: (pos, outV) => {
                outV[0] = 300.0; // uniform flow along +x from Basin 1 to Basin 2
                outV[1] = 0.0;
                outV[2] = 0.0;
            }
        };

        const mesher = new WatershedManifoldMesher({ boxSizeMpc: 160.0 });
        const mesh = mesher.extractBoundaryInterface(grid, N, 1, 2, mockVelocity);
        
        console.log(JSON.stringify({
            triCount: mesh.triangleCount,
            area: mesh.totalAreaMpc2,
            netFlux: mesh.netFluxKmsMpc2,
            hasPositions: mesh.positions.length > 0
        }));
        """)
        # Domain is 160 x 160 x 160 Mpc/h, cross-section is (N-1)*(N-1)*dx^2 = 15*15*(10)^2 = 22500
        assert res["triCount"] > 0
        assert res["area"] > 0
        assert res["netFlux"] > 0, "Flow from Basin 1 to Basin 2 should have positive flux"
        assert res["hasPositions"] is True
