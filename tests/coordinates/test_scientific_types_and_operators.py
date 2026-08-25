# -*- coding: utf-8 -*-
"""
Automated Pytest Suite for Scientific Type System, Units, Coordinates,
Tricubic Interpolator, Differential Operators, and Kinetic Sanity Maps.

Phases 2, 3, 4, 5, 6, 10, 11, 12, 13, 14.

Scientific Truth Rules Enforced:
1. SGZ, SGY, SGX -> SGX, SGY, SGZ stride handling.
2. Coordinates in Mpc/h, Velocities in km/s. Incompatible operations strictly throw TypeError.
3. Exact x52 velocity scale factor on CF4 velocity/error grids.
4. Table A.1 official watershed taxonomy (Basins 1-8: Laniakea, Apus, Hercules, Lepus, Perseus-Pisces, Shapley, SDSS-1a, SDSS-2a).
5. Exact 3D Newton-Raphson vector root finding and Jacobi strain rate tensor diagonalization.
6. W3C PROV-O cryptographic lineage and NIST SHA-256 manifests.
"""

import json
import math
import os
import subprocess
import pytest
import numpy as np

WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def run_node_eval(script_code: str):
    """Executes an ES6 JavaScript code snippet inside Node.js and returns JSON parsed output."""
    proc = subprocess.run(
        ["node", "--input-type=module", "-e", script_code],
        cwd=WORKSPACE_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8"
    )
    if proc.returncode != 0:
        pytest.fail(f"Node execution failed with exit code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")
    
    stdout = proc.stdout.strip()
    if not stdout:
        return None
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return stdout


# ============================================================================
# TIER 1: SCIENTIFIC TYPE SYSTEM & RUNTIME SAFETY
# ============================================================================

class TestSupergalacticPosition:
    """1. SupergalacticPosition validation, operations, dimensional rules, and conversions."""

    def test_position_instantiation_and_properties(self):
        res = run_node_eval("""
        import { SupergalacticPosition, ScientificUnits } from './src/coordinates/scientific_types.js';
        const pos = new SupergalacticPosition(-35.5, 12.0, -8.25, { h: 0.74 });
        console.log(JSON.stringify({
          sgx: pos.sgx,
          sgy: pos.sgy,
          sgz: pos.sgz,
          x: pos.x,
          y: pos.y,
          z: pos.z,
          unit: pos.unit,
          h: pos.h,
          type: pos.scientificType,
          arr: pos.toArray()
        }));
        """)
        assert res["sgx"] == -35.5
        assert res["sgy"] == 12.0
        assert res["sgz"] == -8.25
        assert res["x"] == -35.5
        assert res["unit"] == "Mpc/h"
        assert res["h"] == 0.74
        assert res["type"] == "SupergalacticPosition"
        assert res["arr"] == [-35.5, 12.0, -8.25]

    def test_position_pairwise_and_radial_distance(self):
        res = run_node_eval("""
        import { SupergalacticPosition } from './src/coordinates/scientific_types.js';
        const p1 = new SupergalacticPosition(10.0, 20.0, 30.0);
        const p2 = new SupergalacticPosition(13.0, 24.0, 30.0); // dx=3, dy=4, dz=0 -> dist=5
        const dist = p1.distanceTo(p2);
        const r1 = p1.radialDistance();
        console.log(JSON.stringify({
          distValue: dist.value,
          distUnit: dist.unit,
          r1Value: r1.value,
          expectedR1: Math.sqrt(10*10 + 20*20 + 30*30)
        }));
        """)
        assert math.isclose(res["distValue"], 5.0, abs_tol=1e-12)
        assert res["distUnit"] == "Mpc/h"
        assert math.isclose(res["r1Value"], res["expectedR1"], abs_tol=1e-12)

    def test_position_displacement_addition(self):
        res = run_node_eval("""
        import { SupergalacticPosition } from './src/coordinates/scientific_types.js';
        const p = new SupergalacticPosition(10.0, -20.0, 5.0);
        const translated = p.addDisplacement([5.0, 10.0, -15.0]);
        console.log(JSON.stringify({
          sgx: translated.sgx,
          sgy: translated.sgy,
          sgz: translated.sgz
        }));
        """)
        assert res["sgx"] == 15.0
        assert res["sgy"] == -10.0
        assert res["sgz"] == -10.0

    def test_position_dimensional_error_on_velocity_addition(self):
        """Rule: Adding VelocityVector to SupergalacticPosition must throw TypeError."""
        res = run_node_eval("""
        import { SupergalacticPosition, VelocityVector } from './src/coordinates/scientific_types.js';
        const pos = new SupergalacticPosition(10.0, 20.0, 30.0);
        const vel = new VelocityVector(500.0, -200.0, 150.0);
        let errorCaught = false;
        let errorMessage = '';
        try {
          pos.add(vel);
        } catch (err) {
          errorCaught = true;
          errorMessage = err.message;
        }
        console.log(JSON.stringify({ errorCaught, errorMessage }));
        """)
        assert res["errorCaught"] is True
        assert "Dimensional Error" in res["errorMessage"]

    def test_position_unit_conversions_with_h_scaling(self):
        res = run_node_eval("""
        import { SupergalacticPosition, ScientificUnits } from './src/coordinates/scientific_types.js';
        // 100 Mpc/h with h=0.7 -> physical 100 / 0.7 Mpc = 142.85714 Mpc
        const posMpcH = new SupergalacticPosition(100.0, 0.0, 0.0, { h: 0.7 });
        const posMpc = posMpcH.toUnit(ScientificUnits.MPC);
        const posKpcH = posMpcH.toUnit(ScientificUnits.KPC_OVER_H);
        const posRoundTrip = posMpc.toUnit(ScientificUnits.MPC_OVER_H);
        console.log(JSON.stringify({
          mpcVal: posMpc.sgx,
          kpcHVal: posKpcH.sgx,
          roundTripVal: posRoundTrip.sgx
        }));
        """)
        assert math.isclose(res["mpcVal"], 100.0 / 0.7, rel_tol=1e-6)
        assert math.isclose(res["kpcHVal"], 100000.0, rel_tol=1e-6)
        assert math.isclose(res["roundTripVal"], 100.0, rel_tol=1e-6)

    def test_position_spherical_round_trip(self):
        res = run_node_eval("""
        import { SupergalacticPosition } from './src/coordinates/scientific_types.js';
        const original = new SupergalacticPosition(30.0, 40.0, 50.0);
        const sph = original.toSpherical();
        const reconstructed = SupergalacticPosition.fromSpherical(sph.sgl, sph.sgb, sph.distance);
        console.log(JSON.stringify({
          errX: Math.abs(reconstructed.sgx - original.sgx),
          errY: Math.abs(reconstructed.sgy - original.sgy),
          errZ: Math.abs(reconstructed.sgz - original.sgz),
          distance: sph.distance
        }));
        """)
        assert res["errX"] < 1e-10
        assert res["errY"] < 1e-10
        assert res["errZ"] < 1e-10
        assert math.isclose(res["distance"], math.sqrt(30*30 + 40*40 + 50*50), rel_tol=1e-10)

    def test_position_box_containment(self):
        res = run_node_eval("""
        import { SupergalacticPosition } from './src/coordinates/scientific_types.js';
        const pIn = new SupergalacticPosition(10.0, -20.0, 35.0);
        const pOut = new SupergalacticPosition(110.0, -20.0, 35.0);
        const minB = [-100.0, -100.0, -100.0];
        const maxB = [100.0, 100.0, 100.0];
        console.log(JSON.stringify({
          inBox: pIn.isWithinBox(minB, maxB),
          outBox: pOut.isWithinBox(minB, maxB)
        }));
        """)
        assert res["inBox"] is True
        assert res["outBox"] is False


class TestVelocityVector:
    """2. VelocityVector kinematics, vector algebra, matrix rotation, and CF4 displacement scaling."""

    def test_velocity_instantiation_and_magnitude(self):
        res = run_node_eval("""
        import { VelocityVector } from './src/coordinates/scientific_types.js';
        const v = new VelocityVector(300.0, -400.0, 0.0);
        console.log(JSON.stringify({
          vx: v.vx,
          vy: v.vy,
          vz: v.vz,
          mag: v.magnitude(),
          magSq: v.magnitudeSquared(),
          dir: v.direction()
        }));
        """)
        assert res["vx"] == 300.0
        assert res["vy"] == -400.0
        assert res["vz"] == 0.0
        assert math.isclose(res["mag"], 500.0, abs_tol=1e-12)
        assert math.isclose(res["magSq"], 250000.0, abs_tol=1e-12)
        assert math.isclose(res["dir"][0], 0.6, abs_tol=1e-12)
        assert math.isclose(res["dir"][1], -0.8, abs_tol=1e-12)
        assert math.isclose(res["dir"][2], 0.0, abs_tol=1e-12)

    def test_velocity_vector_algebra(self):
        res = run_node_eval("""
        import { VelocityVector } from './src/coordinates/scientific_types.js';
        const v1 = new VelocityVector(100.0, 200.0, 300.0);
        const v2 = new VelocityVector(50.0, -50.0, 100.0);
        const sum = v1.add(v2);
        const diff = v1.subtract(v2);
        const scaled = v1.scale(2.5);
        const dot = v1.dot(v2);
        const cross = v1.cross(v2);
        console.log(JSON.stringify({
          sum: sum.toArray(),
          diff: diff.toArray(),
          scaled: scaled.toArray(),
          dot: dot,
          cross: cross.toArray()
        }));
        """)
        assert res["sum"] == [150.0, 150.0, 400.0]
        assert res["diff"] == [50.0, 250.0, 200.0]
        assert res["scaled"] == [250.0, 500.0, 750.0]
        assert res["dot"] == 100*50 + 200*(-50) + 300*100  # 5000 - 10000 + 30000 = 25000
        # Cross product: [200*100 - 300*(-50), 300*50 - 100*100, 100*(-50) - 200*50]
        # = [20000 + 15000, 15000 - 10000, -5000 - 10000] = [35000, 5000, -15000]
        assert res["cross"] == [35000.0, 5000.0, -15000.0]

    def test_velocity_dimensional_rejection_on_position_subtraction(self):
        res = run_node_eval("""
        import { VelocityVector, SupergalacticPosition } from './src/coordinates/scientific_types.js';
        const vel = new VelocityVector(300.0, 400.0, 0.0);
        const pos = new SupergalacticPosition(10.0, 20.0, 30.0);
        let errorCaught = false;
        try {
          vel.subtract(pos);
        } catch (err) {
          errorCaught = true;
        }
        console.log(JSON.stringify({ errorCaught }));
        """)
        assert res["errorCaught"] is True

    def test_velocity_to_spatial_displacement(self):
        """Rule: 1 km/s integrated over 1 Gyr ~= 1.022712 Mpc/h (for h=1.0)."""
        res = run_node_eval("""
        import { VelocityVector } from './src/coordinates/scientific_types.js';
        const v = new VelocityVector(1000.0, 0.0, 0.0); // 1000 km/s along X
        const disp = v.toDisplacement(1.0, 1.0); // 1 Gyr, h=1.0
        console.log(JSON.stringify({
          dx: disp[0],
          dy: disp[1],
          dz: disp[2]
        }));
        """)
        assert math.isclose(res["dx"], 1022.712165, rel_tol=1e-4)
        assert res["dy"] == 0.0
        assert res["dz"] == 0.0

    def test_velocity_wiener_filter_scale_factor_x52(self):
        """Rule: CF4 displacement field is converted via exact 52.0 scale factor."""
        res = run_node_eval("""
        import { VelocityVector, CF4_VELOCITY_SCALE_FACTOR } from './src/coordinates/scientific_types.js';
        const psi = [10.0, -5.0, 2.5];
        const v = VelocityVector.fromDisplacementField(psi[0], psi[1], psi[2]);
        console.log(JSON.stringify({
          scaleFactor: CF4_VELOCITY_SCALE_FACTOR,
          vx: v.vx,
          vy: v.vy,
          vz: v.vz
        }));
        """)
        assert res["scaleFactor"] == 52.0
        assert res["vx"] == 520.0
        assert res["vy"] == -260.0
        assert res["vz"] == 130.0

    def test_velocity_matrix_transformation(self):
        res = run_node_eval("""
        import { VelocityVector } from './src/coordinates/scientific_types.js';
        // 90 deg rotation around Z: [0, -1, 0, 1, 0, 0, 0, 0, 1]
        const rotZ90 = [0, -1, 0, 1, 0, 0, 0, 0, 1];
        const v = new VelocityVector(100.0, 0.0, 50.0);
        const vRot = v.transformWithMatrix(rotZ90);
        console.log(JSON.stringify({
          vx: vRot.vx,
          vy: vRot.vy,
          vz: vRot.vz,
          origMag: v.magnitude(),
          rotMag: vRot.magnitude()
        }));
        """)
        assert math.isclose(v_x := res["vx"], 0.0, abs_tol=1e-10)
        assert math.isclose(v_y := res["vy"], 100.0, abs_tol=1e-10)
        assert math.isclose(v_z := res["vz"], 50.0, abs_tol=1e-10)
        assert math.isclose(res["origMag"], res["rotMag"], rel_tol=1e-10)


class TestDensityContrastAndVelocityError:
    """3. DensityContrast physical bounds and VelocityError propagation."""

    def test_density_contrast_physical_bounds(self):
        """Rule: delta >= -1.0 strictly. delta < -1 throws RangeError unless clamped."""
        res = run_node_eval("""
        import { DensityContrast } from './src/coordinates/scientific_types.js';
        const dValid = new DensityContrast(2.5);
        let errorCaught = false;
        try {
          new DensityContrast(-1.5);
        } catch (err) {
          errorCaught = true;
        }
        const dClamped = new DensityContrast(-1.5, { clampToPhysical: true });
        console.log(JSON.stringify({
          validDelta: dValid.delta,
          errorCaught,
          clampedDelta: dClamped.delta,
          ratioValid: dValid.densityRatio(),
          ratioClamped: dClamped.densityRatio()
        }));
        """)
        assert res["validDelta"] == 2.5
        assert res["errorCaught"] is True
        assert res["clampedDelta"] == -1.0
        assert res["ratioValid"] == 3.5
        assert res["ratioClamped"] == 0.0

    def test_density_contrast_regime_classification(self):
        res = run_node_eval("""
        import { DensityContrast } from './src/coordinates/scientific_types.js';
        const dLinearOver = new DensityContrast(0.5);
        const dNonlinearOver = new DensityContrast(3.5);
        const dVoid = new DensityContrast(-0.8);
        console.log(JSON.stringify({
          d1Over: dLinearOver.isOverdense(),
          d1NonLin: dLinearOver.isNonLinear(),
          d2NonLin: dNonlinearOver.isNonLinear(),
          d3Under: dVoid.isUnderdense()
        }));
        """)
        assert res["d1Over"] is True
        assert res["d1NonLin"] is False
        assert res["d2NonLin"] is True
        assert res["d3Under"] is True

    def test_velocity_error_isotropic_and_scale_factor(self):
        res = run_node_eval("""
        import { VelocityError, CF4_VELOCITY_SCALE_FACTOR } from './src/coordinates/scientific_types.js';
        const err = VelocityError.isotropic(10.0);
        const scaled = err.scaleWithFactor(CF4_VELOCITY_SCALE_FACTOR);
        console.log(JSON.stringify({
          evx: err.evx,
          totalErr: err.totalError(),
          meanErr: err.meanError(),
          scaledEvx: scaled.evx,
          scaledTotal: scaled.totalError()
        }));
        """)
        assert res["evx"] == 10.0
        assert math.isclose(res["totalErr"], math.sqrt(300.0), rel_tol=1e-10)
        assert math.isclose(res["meanErr"], 10.0, rel_tol=1e-10)
        assert res["scaledEvx"] == 520.0
        assert math.isclose(res["scaledTotal"], math.sqrt(300.0) * 52.0, rel_tol=1e-10)


class TestGridIndexAndPhysicalDistance:
    """4. GridIndex reversible strides and PhysicalDistance conversions."""

    def test_grid_index_linear_stride_mapping(self):
        res = run_node_eval("""
        import { GridIndex } from './src/coordinates/scientific_types.js';
        const nx = 32, ny = 64, nz = 128;
        const g = new GridIndex(5, 12, 45);
        const linear = g.toLinearIndex(nx, ny);
        const reconstructed = GridIndex.fromLinearIndex(linear, nx, ny);
        console.log(JSON.stringify({
          linear,
          reconstructed: reconstructed.toArray(),
          valid: g.isValid(nx, ny, nz)
        }));
        """)
        assert res["linear"] == 5 + 32 * (12 + 64 * 45)
        assert res["reconstructed"] == [5, 12, 45]
        assert res["valid"] is True

    def test_grid_index_neighbors(self):
        res = run_node_eval("""
        import { GridIndex } from './src/coordinates/scientific_types.js';
        const g = new GridIndex(10, 10, 10);
        const nbrs = g.get6Neighbors().map(n => n.toArray());
        console.log(JSON.stringify({ nbrs }));
        """)
        assert len(res["nbrs"]) == 6
        assert [9, 10, 10] in res["nbrs"]
        assert [11, 10, 10] in res["nbrs"]
        assert [10, 9, 10] in res["nbrs"]
        assert [10, 11, 10] in res["nbrs"]
        assert [10, 10, 9] in res["nbrs"]
        assert [10, 10, 11] in res["nbrs"]

    def test_physical_distance_hubble_inversion(self):
        res = run_node_eval("""
        import { PhysicalDistance, ScientificUnits } from './src/coordinates/scientific_types.js';
        const d = new PhysicalDistance(100.0, ScientificUnits.MPC_OVER_H, { h: 0.75 });
        const cz = d.toHubbleVelocity(100.0);
        const dMpc = d.toMpc();
        console.log(JSON.stringify({
          cz,
          dMpc
        }));
        """)
        assert math.isclose(res["cz"], 10000.0, rel_tol=1e-10)  # 100 Mpc/h * 100 km/s/(Mpc/h) = 10000 km/s
        assert math.isclose(res["dMpc"], 100.0 / 0.75, rel_tol=1e-10)


class TestBasinTaxonomy:
    """5. Table A.1 official Cosmicflows watershed taxonomy (Basins 1 through 8)."""

    def test_all_official_basins_table_a1(self):
        """Verify Basins 1-8: Laniakea, Apus, Hercules, Lepus, Perseus-Pisces, Shapley, SDSS-1a, SDSS-2a."""
        res = run_node_eval("""
        import { BasinID, BASIN_TAXONOMY } from './src/coordinates/scientific_types.js';
        const basins = BasinID.getAllBasins();
        const info = basins.map(b => ({
          id: b.id,
          name: b.name,
          slug: b.slug,
          attractor: b.attractorName,
          center: b.nominalCenterSG,
          posType: b.getCenterPosition().scientificType
        }));
        console.log(JSON.stringify({ count: basins.length, info }));
        """)
        assert res["count"] == 8
        slugs = [item["slug"] for item in res["info"]]
        expected_slugs = [
            'laniakea', 'apus', 'hercules', 'lepus',
            'perseus-pisces', 'shapley', 'sdss-1a', 'sdss-2a'
        ]
        assert slugs == expected_slugs
        for item in res["info"]:
            assert item["posType"] == "SupergalacticPosition"
            assert len(item["center"]) == 3

    def test_basin_rejection_of_invalid_ids(self):
        res = run_node_eval("""
        import { BasinID } from './src/coordinates/scientific_types.js';
        const invalidIds = [0, 9, -1, 100, 'andromeda'];
        const results = invalidIds.map(id => {
          try {
            new BasinID(id);
            return { id, rejected: false };
          } catch (e) {
            return { id, rejected: true };
          }
        });
        console.log(JSON.stringify({ results }));
        """)
        for item in res["results"]:
            assert item["rejected"] is True


class TestJacobianAndEigenSystem:
    """6. Jacobian3x3 strain/vorticity decomposition, Jacobi diagonalization, and cosmic web classification."""

    def test_jacobian_kinematic_decomposition(self):
        res = run_node_eval("""
        import { Jacobian3x3 } from './src/coordinates/scientific_types.js';
        // J = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        const J = new Jacobian3x3([1, 2, 3, 4, 5, 6, 7, 8, 9]);
        const div = J.divergence();
        const S = J.symmetricStrainRate();
        const W = J.antisymmetricVorticity();
        const vort = J.vorticityVector();
        const sigma = J.shearTensor();
        console.log(JSON.stringify({
          div,
          S: Array.from(S),
          W: Array.from(W),
          vort,
          sigma: Array.from(sigma)
        }));
        """)
        assert res["div"] == 1 + 5 + 9  # 15
        # Symmetric strain: S[1] = S[3] = (2+4)/2 = 3; S[2] = S[6] = (3+7)/2 = 5; S[5] = S[7] = (6+8)/2 = 7
        assert res["S"] == [1, 3, 5, 3, 5, 7, 5, 7, 9]
        # Antisymmetric vorticity: W[1] = (2-4)/2 = -1; W[2] = (3-7)/2 = -2; W[5] = (6-8)/2 = -1
        assert res["W"] == [0, -1, -2, 1, 0, -1, 2, 1, 0]
        # Vorticity vector: [8 - 6, 3 - 7, 4 - 2] = [2, -4, 2]
        assert res["vort"] == [2.0, -4.0, 2.0]
        # Traceless shear: diagonal minus div/3 (15/3 = 5): [1-5, 5-5, 9-5] = [-4, 0, 4]
        assert res["sigma"][0] == -4.0
        assert res["sigma"][4] == 0.0
        assert res["sigma"][8] == 4.0

    def test_jacobian_linear_solve_cramer(self):
        res = run_node_eval("""
        import { Jacobian3x3 } from './src/coordinates/scientific_types.js';
        // Invertible matrix: J = diag(2, 3, 4)
        const J = new Jacobian3x3([2, 0, 0, 0, 3, 0, 0, 0, 4]);
        const b = [6.0, 12.0, 20.0];
        const x = J.solveLinear(b);
        console.log(JSON.stringify({ x }));
        """)
        assert res["x"] == [3.0, 4.0, 5.0]

    def test_jacobi_eigensystem_diagonalization(self):
        """Exact 3D Jacobi diagonalization of real symmetric matrix."""
        res = run_node_eval("""
        import { EigenSystem3D } from './src/coordinates/scientific_types.js';
        // Symmetric matrix A = [[2, 1, 0], [1, 3, 1], [0, 1, 2]]
        const A = [2, 1, 0, 1, 3, 1, 0, 1, 2];
        const eigen = EigenSystem3D.fromSymmetricMatrix(A);
        const isOrtho = eigen.isOrthonormal(1e-12);
        const reconstructed = eigen.reconstructMatrix();
        let maxDiff = 0.0;
        for (let i = 0; i < 9; i++) {
          maxDiff = Math.max(maxDiff, Math.abs(reconstructed[i] - A[i]));
        }
        console.log(JSON.stringify({
          eigenvalues: eigen.eigenvalues,
          isOrtho,
          maxDiff
        }));
        """)
        assert res["isOrtho"] is True
        assert res["maxDiff"] < 1e-12
        # Analytical eigenvalues for this tridiagonal symmetric matrix: lambda = 4, 2, 1
        assert math.isclose(res["eigenvalues"][0], 4.0, abs_tol=1e-10)
        assert math.isclose(res["eigenvalues"][1], 2.0, abs_tol=1e-10)
        assert math.isclose(res["eigenvalues"][2], 1.0, abs_tol=1e-10)

    def test_cosmic_web_classification(self):
        res = run_node_eval("""
        import { EigenSystem3D } from './src/coordinates/scientific_types.js';
        const eKnot = new EigenSystem3D([3.0, 2.0, 1.0], [[1,0,0],[0,1,0],[0,0,1]]);
        const eFilament = new EigenSystem3D([3.0, 2.0, -1.0], [[1,0,0],[0,1,0],[0,0,1]]);
        const eSheet = new EigenSystem3D([3.0, -1.0, -2.0], [[1,0,0],[0,1,0],[0,0,1]]);
        const eVoid = new EigenSystem3D([-1.0, -2.0, -3.0], [[1,0,0],[0,1,0],[0,0,1]]);
        console.log(JSON.stringify({
          knot: eKnot.classifyWebStructure(0.0),
          filament: eFilament.classifyWebStructure(0.0),
          sheet: eSheet.classifyWebStructure(0.0),
          void: eVoid.classifyWebStructure(0.0)
        }));
        """)
        assert res["knot"]["category"] == "knot"
        assert res["filament"]["category"] == "filament"
        assert res["sheet"]["category"] == "sheet"
        assert res["void"]["category"] == "void"


class TestStreamlineAndManifest:
    """7. StreamlinePath trajectories and RunManifest cryptographic provenance."""

    def test_streamline_arclength_and_downsampling(self):
        res = run_node_eval("""
        import { StreamlinePath, SupergalacticPosition, VelocityVector, BasinID } from './src/coordinates/scientific_types.js';
        const positions = [];
        const velocities = [];
        // Helix trajectory: x = t, y = sin(t), z = cos(t)
        for (let t = 0; t <= 10; t += 0.5) {
          positions.push(new SupergalacticPosition(t, Math.sin(t), Math.cos(t)));
          velocities.push(new VelocityVector(100, 50, -50));
        }
        const stream = new StreamlinePath(positions, velocities, [], { basinAssignment: new BasinID(1) });
        const len = stream.totalArclength();
        const bbox = stream.boundingBox();
        const downsampled = stream.downsample(5);
        console.log(JSON.stringify({
          pointCount: stream.length,
          arclength: len.value,
          bbox,
          downsampledCount: downsampled.length,
          basin: stream.basinAssignment.name
        }));
        """)
        assert res["pointCount"] == 21
        assert res["arclength"] > 10.0
        assert res["downsampledCount"] == 5
        assert res["basin"] == "Laniakea Supercluster"
        assert res["bbox"]["min"][0] == 0.0
        assert res["bbox"]["max"][0] == 10.0

    def test_run_manifest_prov_and_sha256(self):
        res = run_node_eval("""
        import { RunManifest } from './src/coordinates/scientific_types.js';
        const manifest = new RunManifest({
          runId: 'test-run-123',
          algorithmName: 'CosmicFlows.RK45Integrator',
          parameters: { stepSize: 0.01, tolerance: 1e-6 },
          inputDatasetHash: 'input-hash-abc',
          outputFieldHash: 'output-hash-xyz'
        });
        const isValid = manifest.verifyIntegrity();
        const prov = manifest.toW3CProvJSONLD();
        console.log(JSON.stringify({
          isValid,
          runId: manifest.runId,
          manifestChecksum: manifest.manifestChecksum,
          provType: prov['@type'],
          agent: prov['prov:wasAssociatedWith']['@id']
        }));
        """)
        assert res["isValid"] is True
        assert res["runId"] == "test-run-123"
        assert len(res["manifestChecksum"]) == 64
        assert res["provType"] == "prov:Activity"
        assert "CosmicFlows.RK45Integrator" in res["agent"]


# ============================================================================
# TIER 2: 3D TRICUBIC POLYNOMIAL FIELD INTERPOLATOR
# ============================================================================

class TestTricubicInterpolator:
    """8. 3D 64-point local Hermite & Catmull-Rom spline interpolation, gradients, Hessians, and boundaries."""

    def test_tricubic_exact_node_interpolation(self):
        """Rule: At exact mesh nodes, interpolated value must match grid data to machine precision."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const data = new Float64Array(grid.totalCells);
        // Populate with arbitrary random-like values
        for (let i = 0; i < grid.totalCells; i++) {
          data[i] = Math.sin(i * 0.3) * 100.0;
        }
        const interp = new TricubicInterpolator(grid);
        let maxErr = 0.0;
        for (let iz = 1; iz < 7; iz++) {
          for (let iy = 1; iy < 7; iy++) {
            for (let ix = 1; ix < 7; ix++) {
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              const val = interp.interpolate(data, x, y, z);
              const exact = data[grid.getLinearIndex(ix, iy, iz)];
              maxErr = Math.max(maxErr, Math.abs(val - exact));
            }
          }
        }
        console.log(JSON.stringify({ maxErr }));
        """)
        assert res["maxErr"] < 1e-12

    def test_tricubic_exact_linear_field_reproduction(self):
        """Rule: Tricubic interpolator reproduces linear 3D field f(x,y,z) = a*x + b*y + c*z + d exactly everywhere."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const data = new Float64Array(grid.totalCells);
        const a = 2.5, b = -3.0, c = 4.2, d = 15.0;
        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              data[grid.getLinearIndex(ix, iy, iz)] = a * x + b * y + c * z + d;
            }
          }
        }
        const interp = new TricubicInterpolator(grid);
        let maxValErr = 0.0;
        let maxGradErr = 0.0;
        // Test interior continuous positions
        for (let sz = -30; sz <= 30; sz += 11.5) {
          for (let sy = -30; sy <= 30; sy += 12.3) {
            for (let sx = -30; sx <= 30; sx += 9.7) {
              const val = interp.interpolate(data, sx, sy, sz);
              const exactVal = a * sx + b * sy + c * sz + d;
              maxValErr = Math.max(maxValErr, Math.abs(val - exactVal));

              const grad = interp.interpolateGradient(data, sx, sy, sz);
              maxGradErr = Math.max(maxGradErr, Math.abs(grad[0] - a));
              maxGradErr = Math.max(maxGradErr, Math.abs(grad[1] - b));
              maxGradErr = Math.max(maxGradErr, Math.abs(grad[2] - c));
            }
          }
        }
        console.log(JSON.stringify({ maxValErr, maxGradErr }));
        """)
        assert res["maxValErr"] < 1e-10
        assert res["maxGradErr"] < 1e-10

    def test_tricubic_gradient_and_hessian_accuracy(self):
        """Quadratic field f(x,y,z) = x^2 + 2*y^2 + 3*z^2 + x*y -> grad = [2x+y, 4y+x, 6z], Hess = [[2,1,0],[1,4,0],[0,0,6]]."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 32, ny: 32, nz: 32, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const data = new Float64Array(grid.totalCells);
        for (let iz = 0; iz < 32; iz++) {
          for (let iy = 0; iy < 32; iy++) {
            for (let ix = 0; ix < 32; ix++) {
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              data[grid.getLinearIndex(ix, iy, iz)] = x * x + 2.0 * y * y + 3.0 * z * z + x * y;
            }
          }
        }
        const interp = new TricubicInterpolator(grid);
        const qx = 12.34, qy = -8.76, qz = 5.43;
        const grad = interp.interpolateGradient(data, qx, qy, qz);
        const { hessian, laplacian } = interp.interpolateHessian(data, qx, qy, qz);

        const exactGrad = [2.0 * qx + qy, 4.0 * qy + qx, 6.0 * qz];
        const exactLaplacian = 2.0 + 4.0 + 6.0; // 12.0

        console.log(JSON.stringify({
          gradErrX: Math.abs(grad[0] - exactGrad[0]),
          gradErrY: Math.abs(grad[1] - exactGrad[1]),
          gradErrZ: Math.abs(grad[2] - exactGrad[2]),
          lapErr: Math.abs(laplacian - exactLaplacian),
          h00: hessian[0],
          h11: hessian[4],
          h22: hessian[8],
          h01: hessian[1]
        }));
        """)
        assert res["gradErrX"] < 1e-4
        assert res["gradErrY"] < 1e-4
        assert res["gradErrZ"] < 1e-4
        assert res["lapErr"] < 1e-2
        assert math.isclose(res["h00"], 2.0, abs_tol=1e-2)
        assert math.isclose(res["h11"], 4.0, abs_tol=1e-2)
        assert math.isclose(res["h22"], 6.0, abs_tol=1e-2)
        assert math.isclose(res["h01"], 1.0, abs_tol=1e-2)

    def test_tricubic_vector_and_jacobian_interpolation(self):
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const vx = new Float64Array(grid.totalCells);
        const vy = new Float64Array(grid.totalCells);
        const vz = new Float64Array(grid.totalCells);

        // Hubble flow v = 10 * [x, y, z] -> Jacobian = diag(10, 10, 10), div = 30, curl = 0
        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vx[idx] = 10.0 * x;
              vy[idx] = 10.0 * y;
              vz[idx] = 10.0 * z;
            }
          }
        }

        const interp = new TricubicInterpolator(grid);
        const J = interp.interpolateJacobian(vx, vy, vz, 15.0, -25.0, 30.0);
        const vort = J.vorticityVector();
        console.log(JSON.stringify({
          div: J.divergence(),
          det: J.determinant(),
          vortMag: Math.sqrt(vort[0]*vort[0] + vort[1]*vort[1] + vort[2]*vort[2]),
          j00: J.get(0, 0),
          j11: J.get(1, 1),
          j22: J.get(2, 2)
        }));
        """)
        assert math.isclose(res["div"], 30.0, rel_tol=1e-6)
        assert math.isclose(res["det"], 1000.0, rel_tol=1e-6)
        assert res["vortMag"] < 1e-10
        assert math.isclose(res["j00"], 10.0, rel_tol=1e-6)
        assert math.isclose(res["j11"], 10.0, rel_tol=1e-6)
        assert math.isclose(res["j22"], 10.0, rel_tol=1e-6)

    def test_tricubic_batch_evaluation(self):
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const buf = new Float64Array(grid.totalCells);
        for (let i = 0; i < grid.totalCells; i++) buf[i] = i * 2.0;

        const coords = new Float64Array([
          0, 0, 0,
          10, 20, 30,
          -25, 15, -5
        ]);

        const interp = new TricubicInterpolator(grid);
        const batchRes = interp.interpolateScalarBatch(buf, coords);
        const single0 = interp.interpolate(buf, 0, 0, 0);
        const single1 = interp.interpolate(buf, 10, 20, 30);
        const single2 = interp.interpolate(buf, -25, 15, -5);

        console.log(JSON.stringify({
          batch: Array.from(batchRes),
          singles: [single0, single1, single2]
        }));
        """)
        assert math.isclose(res["batch"][0], res["singles"][0], abs_tol=1e-12)
        assert math.isclose(res["batch"][1], res["singles"][1], abs_tol=1e-12)
        assert math.isclose(res["batch"][2], res["singles"][2], abs_tol=1e-12)

    def test_tricubic_polynomial_64_coefficients(self):
        """Test explicit 64-coefficient solver against continuous interpolation."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [0, 0, 0], boxSize: [7, 7, 7] });
        const data = new Float64Array(grid.totalCells);
        for (let i = 0; i < grid.totalCells; i++) data[i] = Math.cos(i);

        const interp = new TricubicInterpolator(grid);
        const coeffs = interp.getPolynomialCoefficients(data, 2, 2, 2);
        const u = 0.4, v = 0.6, w = 0.2;
        const valFromCoeffs = TricubicInterpolator.evaluateFromCoefficients(coeffs, u, v, w);
        const valFromInterp = interp.interpolate(data, 2 + u, 2 + v, 2 + w);

        console.log(JSON.stringify({
          valFromCoeffs,
          valFromInterp,
          diff: Math.abs(valFromCoeffs - valFromInterp)
        }));
        """)
        assert res["diff"] < 1e-12


# ============================================================================
# TIER 3: FIELD DIFFERENTIAL OPERATORS & KINETIC SANITY MAPS
# ============================================================================

class TestFieldDifferentialOperators:
    """9. Divergence, curl, strain, shear, Okubo-Weiss, enstrophy, and helicity operators."""

    def test_divergence_of_expansion_flow(self):
        """Pure isotropic Hubble expansion: v(x) = H * x -> div(v) = 3*H."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const H = 20.0;
        const vx = new Float64Array(grid.totalCells);
        const vy = new Float64Array(grid.totalCells);
        const vz = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vx[idx] = H * x;
              vy[idx] = H * y;
              vz[idx] = H * z;
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const divGrid = ops.computeDivergenceGrid(vx, vy, vz);
        let maxErr = 0.0;
        // Check interior cells
        for (let iz = 2; iz < 14; iz++) {
          for (let iy = 2; iy < 14; iy++) {
            for (let ix = 2; ix < 14; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              maxErr = Math.max(maxErr, Math.abs(divGrid[idx] - 3.0 * H));
            }
          }
        }
        console.log(JSON.stringify({ maxErr, expectedDiv: 3.0 * H }));
        """)
        assert res["maxErr"] < 1e-10

    def test_curl_of_rigid_rotation_vortex(self):
        """Rigid rotation around Z axis: v = Omega x r = [-Omega*y, Omega*x, 0] -> curl(v) = [0, 0, 2*Omega]."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const Omega = 15.0;
        const vx = new Float64Array(grid.totalCells);
        const vy = new Float64Array(grid.totalCells);
        const vz = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vx[idx] = -Omega * y;
              vy[idx] = Omega * x;
              vz[idx] = 0.0;
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const { wx, wy, wz, magnitude } = ops.computeVorticityGrid(vx, vy, vz);
        let maxErr = 0.0;
        for (let iz = 2; iz < 14; iz++) {
          for (let iy = 2; iy < 14; iy++) {
            for (let ix = 2; ix < 14; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              maxErr = Math.max(maxErr, Math.abs(wx[idx] - 0.0));
              maxErr = Math.max(maxErr, Math.abs(wy[idx] - 0.0));
              maxErr = Math.max(maxErr, Math.abs(wz[idx] - 2.0 * Omega));
            }
          }
        }
        console.log(JSON.stringify({ maxErr, expectedWz: 2.0 * Omega }));
        """)
        assert res["maxErr"] < 1e-10

    def test_irrotational_potential_flow_zero_curl(self):
        """Rule: For potential flow v = -grad(Phi), curl(v) must be zero everywhere."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        // Phi = x^2 + 2*y^2 + 3*z^2 -> v = -[2x, 4y, 6z]
        const vx = new Float64Array(grid.totalCells);
        const vy = new Float64Array(grid.totalCells);
        const vz = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vx[idx] = -2.0 * x;
              vy[idx] = -4.0 * y;
              vz[idx] = -6.0 * z;
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const { magnitude } = ops.computeVorticityGrid(vx, vy, vz);
        let maxVort = 0.0;
        for (let iz = 2; iz < 14; iz++) {
          for (let iy = 2; iy < 14; iy++) {
            for (let ix = 2; ix < 14; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              maxVort = Math.max(maxVort, magnitude[idx]);
            }
          }
        }
        console.log(JSON.stringify({ maxVort }));
        """)
        assert res["maxVort"] < 1e-10

    def test_okubo_weiss_vortex_vs_shear_dominance(self):
        """Okubo-Weiss parameter Q: positive in vortex cores, negative in deformation/strain regions."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        
        // 1. Vortex flow: v = [-Omega*y, Omega*x, 0] -> ||omega|| = 2*Omega, S = 0 -> Q = +Omega^2 > 0
        const vxVort = new Float64Array(grid.totalCells);
        const vyVort = new Float64Array(grid.totalCells);
        const vzVort = new Float64Array(grid.totalCells);
        const Omega = 10.0;
        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vxVort[idx] = -Omega * y;
              vyVort[idx] = Omega * x;
              vzVort[idx] = 0.0;
            }
          }
        }

        // 2. Pure hyperbolic strain flow: v = [S*x, -S*y, 0] -> omega = 0, ||S|| > 0 -> Q < 0
        const vxStrain = new Float64Array(grid.totalCells);
        const vyStrain = new Float64Array(grid.totalCells);
        const vzStrain = new Float64Array(grid.totalCells);
        const S0 = 10.0;
        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              vxStrain[idx] = S0 * x;
              vyStrain[idx] = -S0 * y;
              vzStrain[idx] = 0.0;
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const qVort = ops.computeOkuboWeissGrid(vxVort, vyVort, vzVort);
        const qStrain = ops.computeOkuboWeissGrid(vxStrain, vyStrain, vzStrain);

        const midIdx = grid.getLinearIndex(8, 8, 8);
        console.log(JSON.stringify({
          qVortexMid: qVort[midIdx],
          qStrainMid: qStrain[midIdx]
        }));
        """)
        assert res["qVortexMid"] > 0.0
        assert res["qStrainMid"] < 0.0
        assert math.isclose(res["qVortexMid"], 100.0, rel_tol=1e-6)

    def test_kinetic_energy_and_enstrophy_maps(self):
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-10, -10, -10], boxSize: [20, 20, 20] });
        const total = grid.totalCells;
        const vx = new Float64Array(total).fill(100.0);
        const vy = new Float64Array(total).fill(200.0);
        const vz = new Float64Array(total).fill(200.0); // ||v||^2 = 10000 + 40000 + 40000 = 90000 (km/s)^2
        const density = new Float64Array(total).fill(0.5); // 1 + delta = 1.5

        const ops = new FieldDifferentialOperators(grid);
        const ke = ops.computeKineticEnergyGrid(vx, vy, vz, density, 0.31);
        const enstrophy = ops.computeEnstrophyGrid(vx, vy, vz);

        console.log(JSON.stringify({
          keMid: ke[0],
          enstrophyMid: enstrophy[0]
        }));
        """)
        # rho_c0 = 2.77536627e11, rho_bar = 0.31 * 2.77536627e11 = 8.603635437e10
        # rho = 1.5 * rho_bar = 1.29054531555e11
        # E_k = 0.5 * rho * 90000 = 5.80745392e15
        assert res["keMid"] > 5e15
        assert res["enstrophyMid"] == 0.0  # Constant velocity field has zero enstrophy

    def test_helicity_on_beltrami_vortex(self):
        """Beltrami flow has collinear velocity and vorticity -> h = v . omega != 0."""
        res = run_node_eval("""
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { FieldDifferentialOperators } from './src/fields/field_differential_operators.js';
        const grid = new GridIndexer({ nx: 16, ny: 16, nz: 16, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
        const k = 0.1;
        const vx = new Float64Array(grid.totalCells);
        const vy = new Float64Array(grid.totalCells);
        const vz = new Float64Array(grid.totalCells);

        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
              // Beltrami Arnold-Beltrami-Childress (ABC) field
              vx[idx] = Math.sin(k * z) + Math.cos(k * y);
              vy[idx] = Math.sin(k * x) + Math.cos(k * z);
              vz[idx] = Math.sin(k * y) + Math.cos(k * x);
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const helicity = ops.computeHelicityGrid(vx, vy, vz);
        let nonZeroCount = 0;
        for (let i = 0; i < grid.totalCells; i++) {
          if (Math.abs(helicity[i]) > 1e-4) nonZeroCount++;
        }
        console.log(JSON.stringify({ nonZeroCount, total: grid.totalCells }));
        """)
        assert res["nonZeroCount"] > 0.8 * res["total"]


class TestKineticSanitySuiteAndLinearContinuity:
    """10. Full diagnostic kinetic sanity maps dossier, linear continuity verification, and provenance."""

    def test_generate_kinetic_sanity_maps_synthetic_cosmological_field(self):
        """Verifies linear continuity relation div(v) ~= -H0*f*delta on synthetic cosmological field."""
        res = run_node_eval("""
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
        const continuityCoeff = cosmo.continuityCoefficient; // a * H0 * f ~ 39.09 km/s/(Mpc/h)

        // Generate synthetic gravitational collapse field:
        // Phi = -A * cos(kx*x) * cos(ky*y) * cos(kz*z)
        // delta = lap(Phi) / (H0*f)^2 ... or directly:
        // v = -grad(Phi) / (H0*f)
        const k = 2.0 * Math.PI / 100.0;
        const A = 500.0;

        const vx = new Float64Array(total);
        const vy = new Float64Array(total);
        const vz = new Float64Array(total);
        const density = new Float64Array(total);

        for (let iz = 0; iz < 16; iz++) {
          for (let iy = 0; iy < 16; iy++) {
            for (let ix = 0; ix < 16; ix++) {
              const idx = grid.getLinearIndex(ix, iy, iz);
              const [x, y, z] = grid.getNodeCoord(ix, iy, iz);

              // v = -grad(Phi)
              // vx = - A * k * sin(kx) * cos(ky) * cos(kz)
              // vy = - A * k * cos(kx) * sin(ky) * cos(kz)
              // vz = - A * k * cos(kx) * cos(ky) * sin(kz)
              vx[idx] = -A * k * Math.sin(k * x) * Math.cos(k * y) * Math.cos(k * z);
              vy[idx] = -A * k * Math.cos(k * x) * Math.sin(k * y) * Math.cos(k * z);
              vz[idx] = -A * k * Math.cos(k * x) * Math.cos(k * y) * Math.sin(k * z);

              // div(v) = -3 * A * k^2 * cos(kx) * cos(ky) * cos(kz)
              // In linear theory, delta = -div(v) / h0f = (3 * A * k^2 / h0f) * cos(kx) * cos(ky) * cos(kz)
              density[idx] = (3.0 * A * k * k / continuityCoeff) * Math.cos(k * x) * Math.cos(k * y) * Math.cos(k * z);
            }
          }
        }

        const ops = new FieldDifferentialOperators(grid);
        const dossier = ops.generateKineticSanityMaps({ vx, vy, vz, density }, { cosmology: cosmo });

        console.log(JSON.stringify({
          isLinearConsistent: dossier.continuityAudit.isLinearConsistent,
          relResidual: dossier.continuityAudit.relativeResidual,
          vortToDivRatio: dossier.flowIrrotationality.vorticityToDivergenceRatio,
          isPotentialDominant: dossier.flowIrrotationality.isPotentialFlowDominant,
          fractions: dossier.cosmicWeb.volumeFractions,
          manifestType: dossier.manifest.type,
          checksumLen: dossier.manifest.manifestChecksum.length
        }));
        """)
        assert res["isLinearConsistent"] is True
        assert res["relResidual"] < 0.05  # < 5% residual
        assert res["vortToDivRatio"] < 1e-4  # Purely irrotational potential flow
        assert res["isPotentialDominant"] is True
        assert res["manifestType"] == "RunManifest"
        assert res["checksumLen"] == 64
        # Check that web fractions sum to ~100%
        fractions = res["fractions"]
        total_pct = fractions["voidPct"] + fractions["sheetPct"] + fractions["filamentPct"] + fractions["knotPct"]
        assert math.isclose(total_pct, 100.0, rel_tol=1e-5)

class TestTricubicC1FaceCrossingAndConvergence:
    """Tests rigorous C^1 face-crossing continuity and analytical derivative agreement for Tricubic Catmull-Rom."""

    def test_c1_face_crossing_continuity(self):
        """Verify f(x^-) ~= f(x^+) and grad(f)(x^-) ~= grad(f)(x^+) across voxel boundaries."""
        res = run_node_eval("""
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { TricubicInterpolator } from './src/interpolation/tricubic_interpolator.js';

        const grid = new GridIndexer({
          nx: 8, ny: 8, nz: 8,
          origin: [-20, -20, -20],
          boxSize: [40, 40, 40],
          boundaryMode: BoundaryMode.CLAMP
        });
        const total = grid.totalCells;
        const data = new Float64Array(total);

        for (let i = 0; i < total; i++) {
          const [ix, iy, iz] = grid.get3DIndices(i);
          const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
          data[i] = Math.sin(0.15 * x) * Math.cos(0.15 * y) + 0.1 * z;
        }

        const tri = new TricubicInterpolator(grid, data);

        // Test across internal cell boundary at x = 0.0 (between cell ix=3 and ix=4)
        const eps = 1e-5;
        const y = 2.5, z = -1.5;

        const valMinus = tri.interpolate(data, 0.0 - eps, y, z);
        const valPlus = tri.interpolate(data, 0.0 + eps, y, z);

        const gradMinus = tri.interpolateGradient(data, 0.0 - eps, y, z);
        const gradPlus = tri.interpolateGradient(data, 0.0 + eps, y, z);

        console.log(JSON.stringify({
          valJump: Math.abs(valPlus - valMinus),
          gradXJump: Math.abs(gradPlus[0] - gradMinus[0]),
          gradYJump: Math.abs(gradPlus[1] - gradMinus[1]),
          gradZJump: Math.abs(gradPlus[2] - gradMinus[2])
        }));
        """)
        assert res["valJump"] < 1e-4, f"C0 jump across face: {res['valJump']}"
        assert res["gradXJump"] < 1e-3, f"C1 normal derivative jump across face: {res['gradXJump']}"
        assert res["gradYJump"] < 1e-3, f"C1 tangential derivative jump across face: {res['gradYJump']}"
        assert res["gradZJump"] < 1e-3, f"C1 tangential derivative jump across face: {res['gradZJump']}"
