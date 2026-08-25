# -*- coding: utf-8 -*-
"""
tests/fields/test_real_cf4_continuity_and_scipy_oracle.py
Independent Python/SciPy Oracle & Real-CF4 Continuity Diagnostic Suite.
"""

import math
import os
import json
import pytest
import numpy as np
from scipy.integrate import solve_ivp
from tests.utils import run_node_snippet

class TestSciPyNumericalOracle:
    """Cross-validates ES6 numerical solvers against independent SciPy standard implementations."""

    def test_scipy_oracle_rk45_streamline_trajectory(self):
        """Cross-validates ES6 Cash-Karp RK45 against SciPy solve_ivp RK45 on a 3D Plummer sink field."""
        GM = 10000.0
        b = 10.0

        def plummer_field(t, y):
            r_sq = y[0]**2 + y[1]**2 + y[2]**2
            denom = (r_sq + b**2)**1.5
            return [-GM * y[0] / denom, -GM * y[1] / denom, -GM * y[2] / denom]

        y0 = [30.0, 40.0, 50.0]
        t_span = (0.0, 2.0)
        sol = solve_ivp(plummer_field, t_span, y0, method="RK45", rtol=1e-6, atol=1e-8)

        scipy_endpoint = sol.y[:, -1]

        code = """
        import { rk45CashKarpStep } from './src/integration/rk45_cash_karp.js';
        
        const GM = 10000.0;
        const b = 10.0;
        const plummerField = (t, pos) => {
          const rSq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
          const denom = Math.pow(rSq + b*b, 1.5);
          return [
            -GM * pos[0] / denom,
            -GM * pos[1] / denom,
            -GM * pos[2] / denom
          ];
        };
        
        let pos = [30.0, 40.0, 50.0];
        let t = 0.0;
        let dt = 0.1;
        while (t < 2.0) {
          const stepDt = Math.min(dt, 2.0 - t);
          const res = rk45CashKarpStep(plummerField, pos, t, stepDt, { rtol: 1e-6, atol: 1e-8 });
          if (res.accepted) {
            pos = res.posNext;
            t = res.tNext;
            dt = res.dtNext;
          } else {
            dt = res.dtNext;
          }
        }
        
        console.log(JSON.stringify({
          endpoint: pos
        }));
        """
        js_res = run_node_snippet(code)
        js_endpoint = js_res["endpoint"]

        dist = np.linalg.norm(np.array(scipy_endpoint) - np.array(js_endpoint))
        assert dist < 0.1, f"Displacement between SciPy and ES6 trajectories: {dist} Mpc/h"

    def test_scipy_oracle_strain_rate_jacobi_eigensolver(self):
        """Cross-validates ES6 Jacobi strain-rate eigensolver against NumPy linalg.eigh."""
        np.random.seed(42)
        A = np.random.randn(3, 3)
        S = 0.5 * (A + A.T)
        np_eigenvals, np_eigenvecs = np.linalg.eigh(S)
        np_sorted = np.sort(np_eigenvals)[::-1]

        s00, s01, s02 = S[0,0], S[0,1], S[0,2]
        s10, s11, s12 = S[1,0], S[1,1], S[1,2]
        s20, s21, s22 = S[2,0], S[2,1], S[2,2]

        code = f"""
        import {{ EigenSystem3D }} from './src/coordinates/scientific_types.js';
        
        const matrix = new Float64Array([
          {s00}, {s01}, {s02},
          {s10}, {s11}, {s12},
          {s20}, {s21}, {s22}
        ]);
        
        const eigen = EigenSystem3D.fromSymmetricMatrix(matrix);
        
        console.log(JSON.stringify({{
          lambda1: eigen.lambda1,
          lambda2: eigen.lambda2,
          lambda3: eigen.lambda3
        }}));
        """
        js_res = run_node_snippet(code)

        assert math.isclose(js_res["lambda1"], np_sorted[0], abs_tol=1e-5)
        assert math.isclose(js_res["lambda2"], np_sorted[1], abs_tol=1e-5)
        assert math.isclose(js_res["lambda3"], np_sorted[2], abs_tol=1e-5)


class TestDualGateContinuityAndRealDataDiagnostic:
    """Verifies separation of analytic operator continuity vs real CF4 observational diagnostic."""

    def test_gate1_analytic_linear_continuity_operator(self):
        """Gate 1: Exact mathematical potential flow satisfies div(v) = -a*H*f*delta."""
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        import { CosmologicalParameters } from './src/coordinates/scientific_types.js';
        
        const grid = new GridIndexer({
          nx: 16, ny: 16, nz: 16,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100],
          isCellCentered: true,
          boundaryMode: BoundaryMode.PERIODIC
        });
        const total = grid.totalCells;
        const cosmo = new CosmologicalParameters({ H0: 74.6, omegaM: 0.315, scaleFactorA: 1.0 });
        const continuityCoeff = cosmo.continuityCoefficient;
        
        const k = 2.0 * Math.PI / 100.0;
        const A = 300.0;
        
        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);
        const density = new Float64Array(total);
        
        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vx[idx] = -A * k * Math.sin(k * x) * Math.cos(k * y) * Math.cos(k * z);
              vy[idx] = -A * k * Math.cos(k * x) * Math.sin(k * y) * Math.cos(k * z);
              vz[idx] = -A * k * Math.cos(k * x) * Math.cos(k * y) * Math.sin(k * z);
              density[idx] = (3.0 * A * k * k / continuityCoeff) * Math.cos(k * x) * Math.cos(k * y) * Math.cos(k * z);
            }
          }
        }
        
        const ops = new FieldDifferentialOperators(grid);
        const dossier = ops.generateKineticSanityMaps({ vx, vy, vz, density }, { cosmology: cosmo });
        
        console.log(JSON.stringify({
          isLinearConsistent: dossier.continuityAudit.isLinearConsistent,
          relResidual: dossier.continuityAudit.relativeResidual
        }));
        """
        res = run_node_snippet(code)
        assert res["isLinearConsistent"] is True
        assert res["relResidual"] < 0.05

    def test_gate2_real_data_continuity_diagnostic_residuals(self):
        """Gate 2: Evaluates continuity residuals on mock CF4 data recording physical statistics."""
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        import { CosmologicalParameters, CF4_PUBLIC_VELOCITY_SCALE } from './src/coordinates/scientific_types.js';
        
        const grid = new GridIndexer({
          nx: 16, ny: 16, nz: 16,
          origin: [-250, -250, -250],
          boxSize: [500, 500, 500],
          boundaryMode: BoundaryMode.CLAMP
        });
        const total = grid.totalCells;
        
        const vxRaw = new Float32Array(total);
        const vyRaw = new Float32Array(total);
        const vzRaw = new Float32Array(total);
        const delta = new Float32Array(total);
        
        for (let i = 0; i < total; i++) {
          const [ix, iy, iz] = grid.get3DIndices(i);
          const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
          const r = Math.hypot(x, y, z);
          const speed = 400.0 * Math.exp(-r / 100.0);
          vxRaw[i] = (-x / (r + 1.0)) * speed / CF4_PUBLIC_VELOCITY_SCALE;
          vyRaw[i] = (-y / (r + 1.0)) * speed / CF4_PUBLIC_VELOCITY_SCALE;
          vzRaw[i] = (-z / (r + 1.0)) * speed / CF4_PUBLIC_VELOCITY_SCALE;
          delta[i] = 3.0 * Math.exp(-r / 50.0);
        }
        
        const vxPhysical = new Float64Array(total);
        const vyPhysical = new Float64Array(total);
        const vzPhysical = new Float64Array(total);
        for (let i = 0; i < total; i++) {
          vxPhysical[i] = vxRaw[i] * CF4_PUBLIC_VELOCITY_SCALE;
          vyPhysical[i] = vyRaw[i] * CF4_PUBLIC_VELOCITY_SCALE;
          vzPhysical[i] = vzRaw[i] * CF4_PUBLIC_VELOCITY_SCALE;
        }
        
        const cosmo = new CosmologicalParameters({ H0: 74.6, omegaM: 0.315, scaleFactorA: 1.0 });
        const ops = new FieldDifferentialOperators(grid);
        const dossier = ops.generateKineticSanityMaps({ vx: vxPhysical, vy: vyPhysical, vz: vzPhysical, density: delta }, { cosmology: cosmo });
        
        console.log(JSON.stringify({
          continuityReport: dossier.continuityAudit
        }));
        """
        res = run_node_snippet(code)
        rep = res["continuityReport"]
        assert rep["cosmology"]["H0"] == 74.6
        assert rep["cosmology"]["omegaM"] == 0.315
        assert rep["continuityCoefficient"] == pytest.approx(39.09, abs=0.5)
        assert "residualStats" in rep
        assert rep["residualStats"]["mean"] != 0.0
