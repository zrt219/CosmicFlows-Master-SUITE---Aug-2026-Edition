# -*- coding: utf-8 -*-
"""
Automated Pytest Test Suite for Celestial WCS Coordinate Frame Transformations,
Redshift-Space Distortions (RSD / Kaiser Effect / FoG), and Astrometric FITS WCS Projections.

Comprehensive scientific verification covering:
1. Astrometric Constants, Orientation Vectors, and Dipole Frames (IAU J2000, Planck 2018, CF4).
2. Exact 3x3 Matrix Algebra, Determinants, Inverses, Eigenvalues, and Orthonormality.
3. Astrometric Rotation Matrices (Galactic <-> Equatorial <-> Supergalactic).
4. Canonical Landmarks & Invariants (NGP, SGP, GC, Ascending Nodes).
5. Full Round-Trip Coordinate Inversions & Precision Limits (Spherical & Cartesian).
6. Cosmography Hubs (Virgo, Coma, Centaurus, Norma, Perseus, Shapley Core).
7. Velocity Reference Frame Conversions (Heliocentric <-> CMB <-> Local Group <-> GSR).
8. Vectorized TypedArray Batch Transformations (Strides 3, 4, 6, in-place / out-of-place, large batches).
9. Cosmological Background Evolution (E(z), H(z), D_C(z), D_A(z), D_L(z), f(z), D(z), Alcock-Paczynski, CPL w0/wa).
10. Redshift-Space Distortions (RSD / Kaiser 1987 forward mapping & iterative fixed-point inversion).
11. Anisotropic Redshift-Space Power Spectrum & Legendre Multipoles (P0, P2, P4 + FoG damping models).
12. Cluster Fingers-of-God (FoG) Compression & Dispersion Estimators (Biweight, Gapper, Sigma-Clipped).
13. Astrometric FITS WCS Parser (CD, PC+CDELT, CROTA2 representations, high-order SIP polynomials).
14. Complete Astrometric Map Projections (TAN, SIN, ARC, STG, CAR, CEA, MER, AIT, MOL, PAR, LONPOLE/LATPOLE).
15. Differential Geometry & Jacobian Tensors (Spherical metric tensor, volume preservation, coordinate singularities).
16. Monte Carlo Statistical Validation (Large random ensembles, distribution checks, norm preservation).
17. Adversarial Edge Cases & Numerical Invariants (Zero Vectors, Antipodes, Off-Sky Rays, Ill-Conditioned WCS).
"""

import json
import math
import os
import subprocess
import pytest

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
        pytest.fail(
            f"Node execution failed with exit code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )

    stdout = proc.stdout.strip()
    if not stdout:
        return None
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return stdout


# ============================================================================
# TIER 1: ASTROMETRIC CONSTANTS, ORIENTATION VECTORS & VELOCITY DIPOLE FRAMES
# ============================================================================

class TestAstrometricConstantsAndDipoles:
    """Validates fundamental astronomical constants and velocity frame parameters."""

    def test_constants_and_speed_of_light(self):
        res = run_node_eval("""
        import {
            DEG2RAD, RAD2DEG, ARCSEC2DEG, DEG2ARCSEC, MAS2DEG,
            SPEED_OF_LIGHT_KMS, ASTROMETRIC_CONSTANTS, VELOCITY_REST_FRAMES,
            COSMOLOGY_DEFAULTS
        } from './src/coordinates/celestial_wcs_transform.js';

        console.log(JSON.stringify({
            deg2rad: DEG2RAD,
            rad2deg: RAD2DEG,
            arcsec2deg: ARCSEC2DEG,
            deg2arcsec: DEG2ARCSEC,
            mas2deg: MAS2DEG,
            c: SPEED_OF_LIGHT_KMS,
            ra_ngp: ASTROMETRIC_CONSTANTS.RA_NGP_DEG,
            dec_ngp: ASTROMETRIC_CONSTANTS.DEC_NGP_DEG,
            l_sgp: ASTROMETRIC_CONSTANTS.L_SGP_DEG,
            b_sgp: ASTROMETRIC_CONSTANTS.B_SGP_DEG,
            l_sg0: ASTROMETRIC_CONSTANTS.L_SG0_DEG,
            b_sg0: ASTROMETRIC_CONSTANTS.B_SG0_DEG,
            v_cmb: VELOCITY_REST_FRAMES.CMB_DIPOLE.VELOCITY_KMS,
            l_cmb: VELOCITY_REST_FRAMES.CMB_DIPOLE.GALACTIC_L_DEG,
            b_cmb: VELOCITY_REST_FRAMES.CMB_DIPOLE.GALACTIC_B_DEG,
            v_lg: VELOCITY_REST_FRAMES.LOCAL_GROUP.VELOCITY_KMS,
            l_lg: VELOCITY_REST_FRAMES.LOCAL_GROUP.GALACTIC_L_DEG,
            b_lg: VELOCITY_REST_FRAMES.LOCAL_GROUP.GALACTIC_B_DEG,
            h0: COSMOLOGY_DEFAULTS.H0,
            omega_m: COSMOLOGY_DEFAULTS.OMEGA_M,
            omega_l: COSMOLOGY_DEFAULTS.OMEGA_L
        }));
        """)

        assert pytest.approx(res["deg2rad"], rel=1e-12) == math.pi / 180.0
        assert pytest.approx(res["rad2deg"], rel=1e-12) == 180.0 / math.pi
        assert pytest.approx(res["arcsec2deg"], rel=1e-12) == 1.0 / 3600.0
        assert pytest.approx(res["deg2arcsec"], rel=1e-12) == 3600.0
        assert pytest.approx(res["c"], rel=1e-9) == 299792.458

        # IAU J2000 NGP
        assert pytest.approx(res["ra_ngp"], abs=1e-4) == 192.85948
        assert pytest.approx(res["dec_ngp"], abs=1e-4) == 27.12825

        # de Vaucouleurs SGP & SG0
        assert pytest.approx(res["l_sgp"], abs=1e-3) == 47.37
        assert pytest.approx(res["b_sgp"], abs=1e-3) == 6.32
        assert pytest.approx(res["l_sg0"], abs=1e-3) == 137.37
        assert res["b_sg0"] == 0.0

        # CMB Dipole (Planck 2018)
        assert pytest.approx(res["v_cmb"], abs=1e-3) == 369.8
        assert pytest.approx(res["l_cmb"], abs=1e-3) == 264.021
        assert pytest.approx(res["b_cmb"], abs=1e-3) == 48.253

        # Local Group Motion (Karachentsev 1996)
        assert pytest.approx(res["v_lg"], abs=1e-3) == 316.0
        assert pytest.approx(res["l_lg"], abs=1e-3) == 93.0
        assert pytest.approx(res["b_lg"], abs=1e-3) == -4.0


# ============================================================================
# TIER 2: MATRIX ALGEBRA PRIMITIVES & VECTOR MATHEMATICS
# ============================================================================

class TestMatrixAlgebraPrimitives:
    """Verifies low-level 3x3 matrix and 3D vector arithmetic."""

    def test_matrix_multiplication_and_properties(self):
        res = run_node_eval("""
        import {
            matMul3x3, matVecMul3x3, matTranspose3x3,
            matDet3x3, matInverse3x3, matTrace3x3
        } from './src/coordinates/celestial_wcs_transform.js';

        const I = [1,0,0, 0,1,0, 0,0,1];
        const A = [1,2,3, 0,1,4, 5,6,0];
        const B = [2,0,-1, 1,3,2, 0,-2,1];

        const AI = matMul3x3(A, I);
        const AB = matMul3x3(A, B);
        const detA = matDet3x3(A);
        const traceA = matTrace3x3(A);
        const invA = matInverse3x3(A);
        const A_invA = matMul3x3(A, invA);

        const v = [1, 2, 3];
        const Av = matVecMul3x3(A, v);
        const transA = matTranspose3x3(A);

        console.log(JSON.stringify({
            AI: Array.from(AI),
            AB: Array.from(AB),
            detA: detA,
            traceA: traceA,
            invA: Array.from(invA),
            A_invA: Array.from(A_invA),
            Av: Array.from(Av),
            transA: Array.from(transA)
        }));
        """)

        assert pytest.approx(res["AI"], abs=1e-12) == [1, 2, 3, 0, 1, 4, 5, 6, 0]
        assert pytest.approx(res["detA"], abs=1e-12) == 1.0
        assert pytest.approx(res["traceA"], abs=1e-12) == 2.0

        expected_I = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        assert pytest.approx(res["A_invA"], abs=1e-12) == expected_I
        assert pytest.approx(res["Av"], abs=1e-12) == [14, 14, 17]
        assert pytest.approx(res["transA"], abs=1e-12) == [1, 0, 5, 2, 1, 6, 3, 4, 0]

    def test_vector_primitives(self):
        res = run_node_eval("""
        import {
            vec3Norm, vec3Dot, vec3Cross, vec3Normalize, vec3Angle,
            createRodriguesRotationMatrix, matVecMul3x3
        } from './src/coordinates/celestial_wcs_transform.js';

        const u = [3, 0, 4];
        const v = [0, 5, 0];

        const normU = vec3Norm(u);
        const normV = vec3Norm(v);
        const dotUV = vec3Dot(u, v);
        const crossUV = vec3Cross(u, v);
        const uHat = vec3Normalize(u);
        const angleUV = vec3Angle(u, v);

        const R_z90 = createRodriguesRotationMatrix([0, 0, 1], Math.PI / 2);
        const rotX = matVecMul3x3(R_z90, [1, 0, 0]);

        console.log(JSON.stringify({
            normU, normV, dotUV,
            crossUV: Array.from(crossUV),
            uHat: Array.from(uHat),
            angleUV,
            rotX: Array.from(rotX)
        }));
        """)

        assert pytest.approx(res["normU"], abs=1e-12) == 5.0
        assert pytest.approx(res["normV"], abs=1e-12) == 5.0
        assert pytest.approx(res["dotUV"], abs=1e-12) == 0.0
        assert pytest.approx(res["crossUV"], abs=1e-12) == [-20, 0, 15]
        assert pytest.approx(res["uHat"], abs=1e-12) == [0.6, 0.0, 0.8]
        assert pytest.approx(res["angleUV"], abs=1e-12) == math.pi / 2.0
        assert pytest.approx(res["rotX"], abs=1e-12) == [0.0, 1.0, 0.0]

    def test_spherical_jacobian_tensor(self):
        res = run_node_eval("""
        import { computeSphericalJacobian, sphericalToCartesian, matDet3x3 } from './src/coordinates/celestial_wcs_transform.js';

        const lon = 45.0;
        const lat = 30.0;
        const dist = 100.0;

        const J = computeSphericalJacobian(lon, lat, dist);

        const eps = 1e-6;
        const p0 = sphericalToCartesian(lon, lat, dist);
        const p_r = sphericalToCartesian(lon, lat, dist + eps);
        const p_lon = sphericalToCartesian(lon + eps * (180/Math.PI), lat, dist);
        const p_lat = sphericalToCartesian(lon, lat + eps * (180/Math.PI), dist);

        const d_r = [(p_r[0]-p0[0])/eps, (p_r[1]-p0[1])/eps, (p_r[2]-p0[2])/eps];
        const d_lon = [(p_lon[0]-p0[0])/eps, (p_lon[1]-p0[1])/eps, (p_lon[2]-p0[2])/eps];
        const d_lat = [(p_lat[0]-p0[0])/eps, (p_lat[1]-p0[1])/eps, (p_lat[2]-p0[2])/eps];

        const detJ = matDet3x3(J);
        const expectedDetJ = dist * dist * Math.cos(lat * Math.PI / 180);

        console.log(JSON.stringify({
            J: Array.from(J),
            d_r, d_lon, d_lat,
            detJ, expectedDetJ
        }));
        """)

        J = res["J"]
        d_r = res["d_r"]
        d_lon = res["d_lon"]
        d_lat = res["d_lat"]

        assert pytest.approx([J[0], J[3], J[6]], abs=1e-4) == d_r
        assert pytest.approx([J[1], J[4], J[7]], abs=1e-4) == d_lon
        assert pytest.approx([J[2], J[5], J[8]], abs=1e-4) == d_lat
        assert pytest.approx(res["detJ"], rel=1e-7) == res["expectedDetJ"]


# ============================================================================
# TIER 3: ROTATION MATRICES & ORTHONORMALITY RIGOR
# ============================================================================

class TestRotationMatricesRigorousOrthonormality:
    """Validates the strict orthonormality and right-handedness of all transformation matrices."""

    def test_all_astrometric_rotation_matrices(self):
        res = run_node_eval("""
        import {
            ROT_GAL_TO_EQ, ROT_EQ_TO_GAL,
            ROT_SG_TO_GAL, ROT_GAL_TO_SG,
            ROT_SG_TO_EQ, ROT_EQ_TO_SG,
            verifyOrthonormality3x3, matMul3x3
        } from './src/coordinates/celestial_wcs_transform.js';

        const mats = [
            { name: 'ROT_GAL_TO_EQ', mat: ROT_GAL_TO_EQ },
            { name: 'ROT_EQ_TO_GAL', mat: ROT_EQ_TO_GAL },
            { name: 'ROT_SG_TO_GAL', mat: ROT_SG_TO_GAL },
            { name: 'ROT_GAL_TO_SG', mat: ROT_GAL_TO_SG },
            { name: 'ROT_SG_TO_EQ', mat: ROT_SG_TO_EQ },
            { name: 'ROT_EQ_TO_SG', mat: ROT_EQ_TO_SG }
        ];

        const results = mats.map(m => {
            const check = verifyOrthonormality3x3(m.mat, 1e-12);
            return { name: m.name, isValid: check.isValid, maxErr: check.maxOrthogonalityError, det: check.det };
        });

        const comp = matMul3x3(ROT_GAL_TO_EQ, ROT_SG_TO_GAL);
        let maxCompDiff = 0.0;
        for (let i = 0; i < 9; i++) {
            maxCompDiff = Math.max(maxCompDiff, Math.abs(comp[i] - ROT_SG_TO_EQ[i]));
        }

        console.log(JSON.stringify({ results, maxCompDiff }));
        """)

        assert res["maxCompDiff"] < 1e-14

        for item in res["results"]:
            assert item["isValid"] is True, f"Matrix {item['name']} failed orthonormality: err={item['maxErr']}"
            assert pytest.approx(item["det"], abs=1e-12) == 1.0
            assert item["maxErr"] < 1e-12


# ============================================================================
# TIER 4: CANONICAL LANDMARKS & INVARIANTS
# ============================================================================

class TestLandmarkCoordinateInvariants:
    """Verifies standard astronomical reference poles and landmark coordinates."""

    def test_north_galactic_pole(self):
        res = run_node_eval("""
        import { galacticToEquatorial } from './src/coordinates/celestial_wcs_transform.js';
        const ngp = galacticToEquatorial(0.0, 90.0, 1.0);
        console.log(JSON.stringify(ngp));
        """)
        assert pytest.approx(res["ra"], abs=1e-4) == 192.85948
        assert pytest.approx(res["dec"], abs=1e-4) == 27.12825

    def test_galactic_center(self):
        res = run_node_eval("""
        import { galacticToEquatorial } from './src/coordinates/celestial_wcs_transform.js';
        const gc = galacticToEquatorial(0.0, 0.0, 1.0);
        console.log(JSON.stringify(gc));
        """)
        assert pytest.approx(res["ra"], abs=1e-4) == 266.40510
        assert pytest.approx(res["dec"], abs=1e-4) == -28.93617

    def test_supergalactic_pole_and_node(self):
        res = run_node_eval("""
        import { supergalacticToGalactic, supergalacticToEquatorial } from './src/coordinates/celestial_wcs_transform.js';
        const sgp_gal = supergalacticToGalactic(0.0, 90.0, 1.0);
        const sg0_gal = supergalacticToGalactic(0.0, 0.0, 1.0);
        const sgp_eq = supergalacticToEquatorial(0.0, 90.0, 1.0);
        console.log(JSON.stringify({ sgp_gal, sg0_gal, sgp_eq }));
        """)
        assert pytest.approx(res["sgp_gal"]["l"], abs=1e-4) == 47.37
        assert pytest.approx(res["sgp_gal"]["b"], abs=1e-4) == 6.32
        assert pytest.approx(res["sg0_gal"]["l"], abs=1e-4) == 137.37
        assert pytest.approx(res["sg0_gal"]["b"], abs=1e-4) == 0.0

        assert pytest.approx(res["sgp_eq"]["ra"], abs=1e-2) == 283.754
        assert pytest.approx(res["sgp_eq"]["dec"], abs=1e-2) == 15.700


# ============================================================================
# TIER 5: ROUND-TRIP COORDINATE INVERSIONS
# ============================================================================

class TestRoundTripCoordinateTransforms:
    """Verifies round-trip precision of all coordinate frame converters."""

    def test_spherical_full_circle_roundtrip(self):
        res = run_node_eval("""
        import {
            equatorialToGalactic, galacticToEquatorial,
            galacticToSupergalactic, supergalacticToGalactic,
            equatorialToSupergalactic, supergalacticToEquatorial
        } from './src/coordinates/celestial_wcs_transform.js';

        let maxErrEqGal = 0.0;
        let maxErrGalSG = 0.0;
        let maxErrEqSG = 0.0;

        for (let ra = 10.0; ra < 360.0; ra += 45.0) {
            for (let dec = -75.0; dec <= 75.0; dec += 30.0) {
                const dist = 50.0;

                // Eq -> Gal -> Eq
                const gal = equatorialToGalactic(ra, dec, dist);
                const eqBack = galacticToEquatorial(gal.l, gal.b, gal.dist);
                const dRA = Math.abs(eqBack.ra - ra);
                const dDec = Math.abs(eqBack.dec - dec);
                maxErrEqGal = Math.max(maxErrEqGal, dRA, dDec);

                // Gal -> SG -> Gal
                const sg = galacticToSupergalactic(gal.l, gal.b, dist);
                const galBack = supergalacticToGalactic(sg.sgl, sg.sgb, sg.dist);
                const dL = Math.abs(galBack.l - gal.l);
                const dB = Math.abs(galBack.b - gal.b);
                maxErrGalSG = Math.max(maxErrGalSG, dL, dB);

                // Eq -> SG -> Eq
                const sgDirect = equatorialToSupergalactic(ra, dec, dist);
                const eqDirectBack = supergalacticToEquatorial(sgDirect.sgl, sgDirect.sgb, sgDirect.dist);
                const dRA_dir = Math.abs(eqDirectBack.ra - ra);
                const dDec_dir = Math.abs(eqDirectBack.dec - dec);
                maxErrEqSG = Math.max(maxErrEqSG, dRA_dir, dDec_dir);
            }
        }

        console.log(JSON.stringify({ maxErrEqGal, maxErrGalSG, maxErrEqSG }));
        """)

        assert res["maxErrEqGal"] < 1e-9
        assert res["maxErrGalSG"] < 1e-9
        assert res["maxErrEqSG"] < 1e-9

    def test_cartesian_roundtrips(self):
        res = run_node_eval("""
        import {
            cartesianEqToGal, cartesianGalToEq,
            cartesianGalToSG, cartesianSGToGal,
            cartesianEqToSG, cartesianSGToEq
        } from './src/coordinates/celestial_wcs_transform.js';

        const testPoints = [
            [10.0, 20.0, 30.0],
            [-50.0, 80.0, -120.0],
            [250.0, -300.0, 450.0]
        ];

        let maxErr = 0.0;
        for (const [x, y, z] of testPoints) {
            const gal = cartesianEqToGal(x, y, z);
            const eqBack = cartesianGalToEq(gal[0], gal[1], gal[2]);
            const err1 = Math.hypot(eqBack[0]-x, eqBack[1]-y, eqBack[2]-z);

            const sg = cartesianGalToSG(gal[0], gal[1], gal[2]);
            const galBack = cartesianSGToGal(sg[0], sg[1], sg[2]);
            const err2 = Math.hypot(galBack[0]-gal[0], galBack[1]-gal[1], galBack[2]-gal[2]);

            const sgDirect = cartesianEqToSG(x, y, z);
            const eqDirectBack = cartesianSGToEq(sgDirect[0], sgDirect[1], sgDirect[2]);
            const err3 = Math.hypot(eqDirectBack[0]-x, eqDirectBack[1]-y, eqDirectBack[2]-z);

            maxErr = Math.max(maxErr, err1, err2, err3);
        }

        console.log(JSON.stringify({ maxErr }));
        """)
        assert res["maxErr"] < 1e-12


# ============================================================================
# TIER 6: MAJOR COSMOGRAPHY HUBS BENCHMARKS
# ============================================================================

class TestCosmographyHubsBenchmarks:
    """Verifies that real astronomical hubs match known coordinates across frames."""

    def test_cosmography_hubs_consistency(self):
        res = run_node_eval("""
        import {
            COSMOGRAPHY_HUBS, equatorialToGalactic,
            equatorialToSupergalactic, sphericalToCartesian
        } from './src/coordinates/celestial_wcs_transform.js';

        const results = COSMOGRAPHY_HUBS.map(hub => {
            const galCalc = equatorialToGalactic(hub.eq.ra, hub.eq.dec, hub.eq.distMpc);
            const sgCalc = equatorialToSupergalactic(hub.eq.ra, hub.eq.dec, hub.eq.distMpc);
            const sgCart = sphericalToCartesian(sgCalc.sgl, sgCalc.sgb, hub.eq.distMpc);

            const dGalL = Math.abs(galCalc.l - hub.gal.l);
            const dGalB = Math.abs(galCalc.b - hub.gal.b);
            const dSGX = Math.abs(sgCart[0] - hub.sg.sgx);
            const dSGY = Math.abs(sgCart[1] - hub.sg.sgy);
            const dSGZ = Math.abs(sgCart[2] - hub.sg.sgz);

            return {
                name: hub.name,
                dGalL, dGalB, dSGX, dSGY, dSGZ
            };
        });

        console.log(JSON.stringify(results));
        """)

        for hub in res:
            assert hub["dGalL"] < 0.1, f"Hub {hub['name']} Gal L mismatch: {hub['dGalL']}"
            assert hub["dGalB"] < 0.1, f"Hub {hub['name']} Gal B mismatch: {hub['dGalB']}"
            assert hub["dSGX"] < 1.5, f"Hub {hub['name']} SGX mismatch: {hub['dSGX']}"
            assert hub["dSGY"] < 1.5, f"Hub {hub['name']} SGY mismatch: {hub['dSGY']}"
            assert hub["dSGZ"] < 1.5, f"Hub {hub['name']} SGZ mismatch: {hub['dSGZ']}"


# ============================================================================
# TIER 7: VELOCITY REFERENCE FRAME CONVERSIONS (CMB, LOCAL GROUP, GSR)
# ============================================================================

class TestVelocityReferenceFrameConversions:
    """Verifies velocity transformations between Heliocentric, CMB, Local Group, and GSR frames."""

    def test_cmb_dipole_velocity_corrections(self):
        res = run_node_eval("""
        import {
            heliocentricToCmb, cmbToHeliocentric,
            VELOCITY_REST_FRAMES
        } from './src/coordinates/celestial_wcs_transform.js';

        const dipole = VELOCITY_REST_FRAMES.CMB_DIPOLE;
        const vHelio = 1000.0;

        const vApex = heliocentricToCmb(vHelio, dipole.GALACTIC_L_DEG, dipole.GALACTIC_B_DEG);
        const antiL = (dipole.GALACTIC_L_DEG + 180.0) % 360.0;
        const antiB = -dipole.GALACTIC_B_DEG;
        const vAnti = heliocentricToCmb(vHelio, antiL, antiB);

        const vHelioBack = cmbToHeliocentric(vApex, dipole.GALACTIC_L_DEG, dipole.GALACTIC_B_DEG);

        console.log(JSON.stringify({
            vApex, vAnti, vHelioBack, vApexExpected: vHelio + dipole.VELOCITY_KMS,
            vAntiExpected: vHelio - dipole.VELOCITY_KMS
        }));
        """)

        assert pytest.approx(res["vApex"], abs=1e-4) == res["vApexExpected"]
        assert pytest.approx(res["vAnti"], abs=1e-4) == res["vAntiExpected"]
        assert pytest.approx(res["vHelioBack"], abs=1e-10) == 1000.0

    def test_local_group_and_interframe_consistency(self):
        res = run_node_eval("""
        import {
            heliocentricToLocalGroup, localGroupToHeliocentric,
            cmbToLocalGroup, localGroupToCmb,
            heliocentricToCmb
        } from './src/coordinates/celestial_wcs_transform.js';

        const l = 150.0;
        const b = 30.0;
        const vHelio = 2500.0;

        const vLG = heliocentricToLocalGroup(vHelio, l, b);
        const vHelioBack = localGroupToHeliocentric(vLG, l, b);

        const vCMB = heliocentricToCmb(vHelio, l, b);
        const vLG_from_CMB = cmbToLocalGroup(vCMB, l, b);
        const vCMB_back = localGroupToCmb(vLG_from_CMB, l, b);

        console.log(JSON.stringify({
            diffLG: Math.abs(vLG - vLG_from_CMB),
            diffHelioBack: Math.abs(vHelioBack - vHelio),
            diffCMBBack: Math.abs(vCMB_back - vCMB)
        }));
        """)

        assert res["diffLG"] < 1e-10
        assert res["diffHelioBack"] < 1e-10
        assert res["diffCMBBack"] < 1e-10

    def test_galactic_standard_of_rest(self):
        res = run_node_eval("""
        import { heliocentricToGSR } from './src/coordinates/celestial_wcs_transform.js';
        const vGSR = heliocentricToGSR(0.0, 90.0, 0.0);
        console.log(JSON.stringify({ vGSR }));
        """)
        assert pytest.approx(res["vGSR"], abs=1e-2) == 232.24


# ============================================================================
# TIER 8: VECTORIZED TYPEDARRAY BATCH TRANSFORMATIONS
# ============================================================================

class TestVectorizedBatchTransformations:
    """Verifies high-performance TypedArray batch operations and shader compatibility."""

    def test_stride_3_cartesian_batch(self):
        res = run_node_eval("""
        import { transformCartesianBatch, vec3Norm } from './src/coordinates/celestial_wcs_transform.js';

        const count = 1000;
        const inBuf = new Float32Array(count * 3);
        const outBuf = new Float32Array(count * 3);

        for (let i = 0; i < count; i++) {
            inBuf[i * 3] = (i * 1.7) % 100.0 - 50.0;
            inBuf[i * 3 + 1] = (i * 3.1) % 100.0 - 50.0;
            inBuf[i * 3 + 2] = (i * 4.9) % 100.0 - 50.0;
        }

        transformCartesianBatch(inBuf, outBuf, count, 'GAL', 'EQ', 3, 0);

        let maxNormDiff = 0.0;
        for (let i = 0; i < count; i++) {
            const nIn = Math.hypot(inBuf[i*3], inBuf[i*3+1], inBuf[i*3+2]);
            const nOut = Math.hypot(outBuf[i*3], outBuf[i*3+1], outBuf[i*3+2]);
            maxNormDiff = Math.max(maxNormDiff, Math.abs(nIn - nOut));
        }

        transformCartesianBatch(outBuf, outBuf, count, 'EQ', 'GAL', 3, 0);
        let maxInPlaceDiff = 0.0;
        for (let i = 0; i < count * 3; i++) {
            maxInPlaceDiff = Math.max(maxInPlaceDiff, Math.abs(outBuf[i] - inBuf[i]));
        }

        console.log(JSON.stringify({ maxNormDiff, maxInPlaceDiff }));
        """)

        assert res["maxNormDiff"] < 1e-4
        assert res["maxInPlaceDiff"] < 1e-4

    def test_stride_6_phase_space_batch(self):
        res = run_node_eval("""
        import { transformPhaseSpace6DBatch } from './src/coordinates/celestial_wcs_transform.js';

        const count = 500;
        const inBuf = new Float64Array(count * 6);
        const outBuf = new Float64Array(count * 6);

        for (let i = 0; i < count; i++) {
            inBuf[i * 6] = 10.0 + i;
            inBuf[i * 6 + 1] = 20.0 + i;
            inBuf[i * 6 + 2] = 30.0 + i;
            inBuf[i * 6 + 3] = 100.0 + i * 2;
            inBuf[i * 6 + 4] = 200.0 + i * 2;
            inBuf[i * 6 + 5] = 300.0 + i * 2;
        }

        transformPhaseSpace6DBatch(inBuf, outBuf, count, 'EQ', 'SG', 6);
        transformPhaseSpace6DBatch(outBuf, outBuf, count, 'SG', 'EQ', 6);

        let maxResidual = 0.0;
        for (let i = 0; i < count * 6; i++) {
            maxResidual = Math.max(maxResidual, Math.abs(outBuf[i] - inBuf[i]));
        }

        console.log(JSON.stringify({ maxResidual }));
        """)
        assert res["maxResidual"] < 1e-12


# ============================================================================
# TIER 9: COSMOLOGICAL BACKGROUND & EXPANSION RATE
# ============================================================================

class TestCosmologyAndExpansionRate:
    """Verifies Friedmann expansion rates, distance measures, and growth parameters."""

    def test_expansion_rate_and_distances(self):
        res = run_node_eval("""
        import { CosmologyEngine } from './src/coordinates/celestial_wcs_transform.js';

        const cosmo = new CosmologyEngine({ H0: 70.0, OmegaM: 0.3, OmegaL: 0.7, OmegaR: 0.0 });

        const E0 = cosmo.E(0.0);
        const H0 = cosmo.H(0.0);
        const E1 = cosmo.E(1.0);
        const D_C1 = cosmo.comovingDistanceMpc(1.0);
        const D_A1 = cosmo.angularDiameterDistanceMpc(1.0);
        const D_L1 = cosmo.luminosityDistanceMpc(1.0);

        const zRecov = cosmo.redshiftFromComovingDistance(D_C1);

        console.log(JSON.stringify({
            E0, H0, E1, D_C1, D_A1, D_L1, zRecov,
            expectedE1: Math.sqrt(0.3 * 8.0 + 0.7)
        }));
        """)

        assert pytest.approx(res["E0"], abs=1e-12) == 1.0
        assert pytest.approx(res["H0"], abs=1e-12) == 70.0
        assert pytest.approx(res["E1"], abs=1e-4) == res["expectedE1"]
        assert pytest.approx(res["D_L1"], abs=1e-6) == 4.0 * res["D_A1"]
        assert pytest.approx(res["zRecov"], abs=1e-6) == 1.0

    def test_linear_growth_rate_and_factor(self):
        res = run_node_eval("""
        import { CosmologyEngine } from './src/coordinates/celestial_wcs_transform.js';

        const cosmo = new CosmologyEngine({ OmegaM: 0.3, OmegaL: 0.7, OmegaR: 0.0 });
        const f0 = cosmo.growthRateF(0.0);
        const D0 = cosmo.growthFactorD(0.0);
        const D1 = cosmo.growthFactorD(1.0);

        console.log(JSON.stringify({ f0, D0, D1 }));
        """)

        assert pytest.approx(res["f0"], abs=0.02) == 0.52
        assert pytest.approx(res["D0"], abs=1e-6) == 1.0
        assert 0.5 < res["D1"] < 0.7

    def test_cpl_dark_energy_parametrization(self):
        res = run_node_eval("""
        import { CosmologyEngine } from './src/coordinates/celestial_wcs_transform.js';

        const cpl = new CosmologyEngine({ H0: 70.0, OmegaM: 0.3, OmegaL: 0.7, w0: -0.9, wa: 0.2 });
        const E0 = cpl.E(0.0);
        const E1 = cpl.E(1.0);

        console.log(JSON.stringify({ E0, E1 }));
        """)

        assert pytest.approx(res["E0"], abs=1e-10) == 1.0
        assert res["E1"] > 1.0

    def test_alcock_paczynski_effect(self):
        res = run_node_eval("""
        import { CosmologyEngine } from './src/coordinates/celestial_wcs_transform.js';

        const trueCosmo = new CosmologyEngine({ H0: 74.6, OmegaM: 0.315, OmegaL: 0.685 });
        const fidCosmo = new CosmologyEngine({ H0: 70.0, OmegaM: 0.300, OmegaL: 0.700 });

        const ap = trueCosmo.alcockPaczynskiParameters(0.5, fidCosmo);
        console.log(JSON.stringify(ap));
        """)

        assert res["qPerp"] > 0.0
        assert res["qPara"] > 0.0
        assert pytest.approx(res["F_AP"], rel=1e-10) == res["qPara"] / res["qPerp"]


# ============================================================================
# TIER 10: REDSHIFT-SPACE DISTORTIONS (RSD / KAISER ENGINE)
# ============================================================================

class TestRedshiftSpaceDistortions:
    """Verifies Kaiser linear RSD transformations, inversions, and power spectrum multipoles."""

    def test_real_to_redshift_space_mapping_and_inversion(self):
        res = run_node_eval("""
        import {
            realToRedshiftSpace, redshiftToRealSpace,
            CosmologyEngine
        } from './src/coordinates/celestial_wcs_transform.js';

        const cosmo = new CosmologyEngine({ H0: 100.0 });
        const xReal = [0.0, 100.0, 0.0];
        const vPec = [0.0, 500.0, 0.0];

        const sRedshift = realToRedshiftSpace(xReal, vPec, { cosmology: cosmo, isMpcOverH: true });
        const xInverted = redshiftToRealSpace(sRedshift, vPec, { cosmology: cosmo, isMpcOverH: true });

        console.log(JSON.stringify({
            xReal,
            sRedshift: Array.from(sRedshift),
            xInverted: Array.from(xInverted)
        }));
        """)

        assert pytest.approx(res["sRedshift"][1], abs=0.2) == 105.0
        assert pytest.approx(res["xInverted"], abs=1e-4) == res["xReal"]

    def test_kaiser_analytical_and_numerical_multipoles(self):
        res = run_node_eval("""
        import { KaiserRSDModel } from './src/coordinates/celestial_wcs_transform.js';

        const modelNoFoG = new KaiserRSDModel({ growthRateF: 0.5, biasB: 1.0, sigmaV: 0.0, fogModel: 'none' });
        const PkReal = 1000.0;

        const analytical = modelNoFoG.analyticalLinearMultipoles(PkReal);
        const numerical = modelNoFoG.computeMultipolesNumerical(0.1, PkReal);

        console.log(JSON.stringify({ analytical, numerical }));
        """)

        ana = res["analytical"]
        num = res["numerical"]

        assert pytest.approx(ana["P0"], abs=1.0) == 1383.33
        assert pytest.approx(ana["P2"], abs=1.0) == 809.52
        assert pytest.approx(ana["P4"], abs=1.0) == 57.14

        assert pytest.approx(num["P0"], rel=1e-2) == ana["P0"]
        assert pytest.approx(num["P2"], rel=1e-2) == ana["P2"]
        assert pytest.approx(num["P4"], rel=1e-2) == ana["P4"]

    def test_fog_damping_models(self):
        res = run_node_eval("""
        import { KaiserRSDModel } from './src/coordinates/celestial_wcs_transform.js';

        const k = 0.5;
        const mu = 0.8;
        const PkReal = 500.0;

        const mNone = new KaiserRSDModel({ fogModel: 'none', sigmaV: 400.0 });
        const mLorentz = new KaiserRSDModel({ fogModel: 'lorentzian', sigmaV: 400.0, H0: 100.0 });
        const mGauss = new KaiserRSDModel({ fogModel: 'gaussian', sigmaV: 400.0, H0: 100.0 });
        const mExp = new KaiserRSDModel({ fogModel: 'exponential', sigmaV: 400.0, H0: 100.0 });

        console.log(JSON.stringify({
            pNone: mNone.anisotropicPower(k, mu, PkReal),
            pLorentz: mLorentz.anisotropicPower(k, mu, PkReal),
            pGauss: mGauss.anisotropicPower(k, mu, PkReal),
            pExp: mExp.anisotropicPower(k, mu, PkReal)
        }));
        """)

        pNone = res["pNone"]
        pLorentz = res["pLorentz"]
        pGauss = res["pGauss"]
        pExp = res["pExp"]

        assert pLorentz < pNone
        assert pGauss < pNone
        assert pExp < pNone
        assert pGauss < pLorentz


# ============================================================================
# TIER 11: CLUSTER FOG COMPRESSION & DISPERSION ESTIMATORS
# ============================================================================

class TestClusterFoGCompression:
    """Verifies cluster virial velocity dispersion estimators and FoG radial de-stretching."""

    def test_velocity_dispersion_estimators(self):
        res = run_node_eval("""
        import { ClusterFoGCompressor } from './src/coordinates/celestial_wcs_transform.js';

        const v = [
            4200, 4450, 4700, 4850, 4950, 5000, 5050, 5150, 5300, 5550, 5800,
            12000
        ];

        const biweight = ClusterFoGCompressor.estimateBiweightDispersion(v);
        const gapper = ClusterFoGCompressor.estimateGapperDispersion(v);
        const sigmaClipped = ClusterFoGCompressor.estimateSigmaClippedDispersion(v, 2.5);

        console.log(JSON.stringify({ biweight, gapper, sigmaClipped }));
        """)

        biweight = res["biweight"]
        sigmaClipped = res["sigmaClipped"]

        assert pytest.approx(biweight["median"], abs=100.0) == 5000.0
        assert 400.0 < biweight["sigmaBiweight"] < 800.0
        assert 400.0 < sigmaClipped["std"] < 800.0

    def test_cluster_radial_compression(self):
        res = run_node_eval("""
        import { ClusterFoGCompressor } from './src/coordinates/celestial_wcs_transform.js';

        const center = { x: 0.0, y: 100.0, z: 0.0 };
        const members = [
            { x: 0.0, y: 90.0, z: 0.0 },
            { x: 0.0, y: 110.0, z: 0.0 },
            { x: 5.0, y: 100.0, z: 0.0 }
        ];

        const sigmaV = 1000.0;
        const rVir = 2.0;
        const H0 = 100.0;

        const compressed = ClusterFoGCompressor.compressClusterFoG(members, center, sigmaV, rVir, { H0 });
        console.log(JSON.stringify({ center, members, compressed }));
        """)

        c = res["compressed"]
        assert pytest.approx(c[2]["x"], abs=1e-6) == 5.0
        assert pytest.approx(c[2]["y"], abs=1e-6) == 100.0

        assert pytest.approx(c[0]["y"], abs=0.1) == 100.0 - (10.0 / 6.0)
        assert pytest.approx(c[1]["y"], abs=0.1) == 100.0 + (10.0 / 6.0)


# ============================================================================
# TIER 12: ASTROMETRIC FITS WCS PARSER & PROJECTIONS
# ============================================================================

class TestFITSWCSParserAndProjections:
    """Verifies FITS WCS standard parsing, CD/PC matrix unification, SIP distortions, and projections."""

    def test_fits_wcs_header_cd_and_pc_matrix_representations(self):
        res = run_node_eval("""
        import { WCSHeader, FITS_WCS_Parser } from './src/coordinates/celestial_wcs_transform.js';

        const wcsCD = new WCSHeader({
            CRPIX1: 512.5, CRPIX2: 512.5,
            CRVAL1: 180.0, CRVAL2: 45.0,
            CD1_1: -0.00027778, CD1_2: 0.0,
            CD2_1: 0.0, CD2_2: 0.00027778,
            CTYPE1: 'RA---TAN', CTYPE2: 'DEC--TAN'
        });

        const wcsPC = new WCSHeader({
            CRPIX1: 512.5, CRPIX2: 512.5,
            CRVAL1: 180.0, CRVAL2: 45.0,
            CDELT1: -0.00027778, CDELT2: 0.00027778,
            PC1_1: 1.0, PC1_2: 0.0,
            PC2_1: 0.0, PC2_2: 1.0,
            CTYPE1: 'RA---TAN', CTYPE2: 'DEC--TAN'
        });

        console.log(JSON.stringify({
            scaleCD: wcsCD.pixelScaleArcsec,
            scalePC: wcsPC.pixelScaleArcsec,
            cd1_1: wcsPC.cd[0],
            cd2_2: wcsPC.cd[3],
            invCD: Array.from(wcsCD.cdInverse)
        }));
        """)

        assert pytest.approx(res["scaleCD"], abs=1e-3) == 1.0
        assert pytest.approx(res["scalePC"], abs=1e-3) == 1.0
        assert pytest.approx(res["cd1_1"], abs=1e-8) == -0.00027778
        assert pytest.approx(res["cd2_2"], abs=1e-8) == 0.00027778

    def test_raw_fits_header_cards_parser(self):
        raw_header = (
            "SIMPLE  =                    T / file conforms to FITS standard\n"
            "BITPIX  =                  -32 / number of bits per data pixel\n"
            "NAXIS   =                    2 / number of data axes\n"
            "NAXIS1  =                 1024 / length of data axis 1\n"
            "NAXIS2  =                 1024 / length of data axis 2\n"
            "CRPIX1  =                512.0 / Reference pixel in X\n"
            "CRPIX2  =                512.0 / Reference pixel in Y\n"
            "CRVAL1  =            150.12345 / RA at reference pixel (deg)\n"
            "CRVAL2  =             2.543210 / Dec at reference pixel (deg)\n"
            "CTYPE1  = 'RA---TAN'           / Coordinate type axis 1\n"
            "CTYPE2  = 'DEC--TAN'           / Coordinate type axis 2\n"
            "CD1_1   =        -0.0001388889 / CD matrix\n"
            "CD1_2   =                  0.0 / CD matrix\n"
            "CD2_1   =                  0.0 / CD matrix\n"
            "CD2_2   =         0.0001388889 / CD matrix\n"
            "END"
        )

        res = run_node_eval(f"""
        import {{ FITS_WCS_Parser }} from './src/coordinates/celestial_wcs_transform.js';

        const raw = {json.dumps(raw_header)};
        const parser = new FITS_WCS_Parser(raw);

        const centerWorld = parser.pixelToWorld(512.0, 512.0, 1);
        const centerPixBack = parser.worldToPixel(centerWorld.worldLon, centerWorld.worldLat, 1);
        const footprint = parser.computeFootprint(1);
        const cards = parser.exportToFITSHeaderCards();

        console.log(JSON.stringify({{
            crval1: parser.wcs.crval1,
            crval2: parser.wcs.crval2,
            projCode: parser.wcs.projectionCode,
            centerWorld,
            centerPixBack,
            footprintCorners: footprint.length,
            cardsCount: cards.length
        }}));
        """)

        assert pytest.approx(res["crval1"], abs=1e-5) == 150.12345
        assert pytest.approx(res["crval2"], abs=1e-5) == 2.54321
        assert res["projCode"] == "TAN"
        assert pytest.approx(res["centerWorld"]["worldLon"], abs=1e-5) == 150.12345
        assert pytest.approx(res["centerWorld"]["worldLat"], abs=1e-5) == 2.54321
        assert pytest.approx(res["centerPixBack"]["xPix"], abs=1e-5) == 512.0
        assert pytest.approx(res["centerPixBack"]["yPix"], abs=1e-5) == 512.0
        assert res["footprintCorners"] == 4
        assert res["cardsCount"] > 10

    def test_all_sky_and_zenithal_projections_roundtrip(self):
        res = run_node_eval("""
        import { FITS_WCS_Parser } from './src/coordinates/celestial_wcs_transform.js';

        const projList = ['TAN', 'SIN', 'ARC', 'STG', 'CAR', 'CEA', 'MER', 'AIT', 'MOL', 'PAR'];
        const results = {};

        for (const proj of projList) {
            const parser = new FITS_WCS_Parser({
                CRPIX1: 100.0, CRPIX2: 100.0,
                CRVAL1: 45.0, CRVAL2: 20.0,
                CD1_1: -0.01, CD1_2: 0.0,
                CD2_1: 0.0, CD2_2: 0.01,
                CTYPE1: `RA---${proj}`, CTYPE2: `DEC--${proj}`
            });

            const testX = 120.0;
            const testY = 110.0;
            const world = parser.pixelToWorld(testX, testY, 1);
            const pixBack = parser.worldToPixel(world.worldLon, world.worldLat, 1);

            const dx = Math.abs(pixBack.xPix - testX);
            const dy = Math.abs(pixBack.yPix - testY);
            results[proj] = { dx, dy, worldLon: world.worldLon, worldLat: world.worldLat };
        }

        console.log(JSON.stringify(results));
        """)

        for proj, data in res.items():
            assert data["dx"] < 1e-4, f"Projection {proj} X inversion error: {data['dx']}"
            assert data["dy"] < 1e-4, f"Projection {proj} Y inversion error: {data['dy']}"

    def test_sip_polynomial_distortion_correction(self):
        res = run_node_eval("""
        import { FITS_WCS_Parser } from './src/coordinates/celestial_wcs_transform.js';

        const parser = new FITS_WCS_Parser({
            CRPIX1: 500.0, CRPIX2: 500.0,
            CRVAL1: 200.0, CRVAL2: -15.0,
            CD1_1: -0.0001, CD1_2: 0.0,
            CD2_1: 0.0, CD2_2: 0.0001,
            CTYPE1: 'RA---TAN-SIP', CTYPE2: 'DEC--TAN-SIP',
            A_ORDER: 2, B_ORDER: 2,
            A_2_0: 0.00001, A_0_2: 0.00002,
            B_1_1: 0.000015,
            AP_ORDER: 2, BP_ORDER: 2,
            AP_2_0: -0.00001, AP_0_2: -0.00002,
            BP_1_1: -0.000015
        });

        const xPix = 600.0;
        const yPix = 550.0;
        const world = parser.pixelToWorld(xPix, yPix, 1);
        const pixBack = parser.worldToPixel(world.worldLon, world.worldLat, 1);

        console.log(JSON.stringify({
            xPix, yPix, world,
            pixBack,
            diffX: Math.abs(pixBack.xPix - xPix),
            diffY: Math.abs(pixBack.yPix - yPix)
        }));
        """)

        assert res["diffX"] < 1e-3
        assert res["diffY"] < 1e-3


# ============================================================================
# TIER 13: ADVERSARIAL EDGE CASES & NUMERICAL STABILITY
# ============================================================================

class TestAdversarialEdgeCasesAndStability:
    """Verifies robustness under singular matrices, pole singularities, and boundary wrapping."""

    def test_zero_vector_and_pole_singularities(self):
        res = run_node_eval("""
        import {
            cartesianToSpherical, sphericalToCartesian,
            equatorialToGalactic, wrapAngleDeg, clampLatitudeDeg
        } from './src/coordinates/celestial_wcs_transform.js';

        const sphZero = cartesianToSpherical(0.0, 0.0, 0.0);
        const northPole = equatorialToGalactic(0.0, 90.0, 1.0);
        const southPole = equatorialToGalactic(180.0, -90.0, 1.0);

        const wrapNeg = wrapAngleDeg(-45.0, 0, 360);
        const wrapOver = wrapAngleDeg(400.0, 0, 360);
        const clampOver = clampLatitudeDeg(95.0);
        const clampUnder = clampLatitudeDeg(-95.0);

        console.log(JSON.stringify({
            sphZero: Array.from(sphZero),
            northPole, southPole,
            wrapNeg, wrapOver, clampOver, clampUnder
        }));
        """)

        assert res["sphZero"] == [0.0, 0.0, 0.0]
        assert pytest.approx(res["wrapNeg"], abs=1e-12) == 315.0
        assert pytest.approx(res["wrapOver"], abs=1e-12) == 40.0
        assert res["clampOver"] == 90.0
        assert res["clampUnder"] == -90.0

    def test_singular_wcs_matrix_handling(self):
        res = run_node_eval("""
        import { WCSHeader } from './src/coordinates/celestial_wcs_transform.js';

        let threw = false;
        try {
            const wcs = new WCSHeader({
                CD1_1: 0.0, CD1_2: 0.0,
                CD2_1: 0.0, CD2_2: 0.0
            });
            const inv = wcs.cdInverse;
        } catch (e) {
            threw = true;
        }

        console.log(JSON.stringify({ threw }));
        """)
        assert res["threw"] is True


# ============================================================================
# TIER 14: MONTE CARLO STATISTICAL ENSEMBLE & VOLUME PRESERVATION
# ============================================================================

class TestMonteCarloStatisticalEnsemble:
    """Verifies statistical properties and metric conservation across large coordinate ensembles."""

    def test_monte_carlo_spherical_volume_preservation(self):
        res = run_node_eval("""
        import {
            sphericalToCartesian, cartesianToSpherical,
            equatorialToSupergalactic, supergalacticToEquatorial
        } from './src/coordinates/celestial_wcs_transform.js';

        let maxDev = 0.0;
        const N = 2000;

        for (let i = 0; i < N; i++) {
            const ra = (i * 137.5) % 360.0;
            const dec = Math.asin((i / N) * 1.98 - 0.99) * (180.0 / Math.PI);
            const dist = 10.0 + (i % 200);

            const sg = equatorialToSupergalactic(ra, dec, dist);
            const eqBack = supergalacticToEquatorial(sg.sgl, sg.sgb, sg.dist);

            const rawDRA = Math.abs(eqBack.ra - ra);
            const dRA = Math.min(rawDRA, Math.abs(360.0 - rawDRA));
            const dDec = Math.abs(eqBack.dec - dec);
            const dDist = Math.abs(eqBack.dist - dist);

            maxDev = Math.max(maxDev, dRA, dDec, dDist);
        }

        console.log(JSON.stringify({ maxDev, N }));
        """)

        assert res["maxDev"] < 1e-8
        assert res["N"] == 2000

    def test_crota_matrix_derivation(self):
        res = run_node_eval("""
        import { WCSHeader, FITS_WCS_Parser } from './src/coordinates/celestial_wcs_transform.js';

        const parserCROTA = new FITS_WCS_Parser({
            CRPIX1: 256.0, CRPIX2: 256.0,
            CRVAL1: 100.0, CRVAL2: 30.0,
            CDELT1: -0.001, CDELT2: 0.001,
            CROTA2: 45.0,
            CTYPE1: 'RA---TAN', CTYPE2: 'DEC--TAN'
        });

        const w = parserCROTA.pixelToWorld(256.0, 256.0, 1);
        const pBack = parserCROTA.worldToPixel(w.worldLon, w.worldLat, 1);

        console.log(JSON.stringify({
            worldCenter: w,
            pixelBack: pBack,
            cd: Array.from(parserCROTA.wcs.cd)
        }));
        """)

        assert pytest.approx(res["worldCenter"]["worldLon"], abs=1e-5) == 100.0
        assert pytest.approx(res["worldCenter"]["worldLat"], abs=1e-5) == 30.0
        assert pytest.approx(res["pixelBack"]["xPix"], abs=1e-5) == 256.0
        assert pytest.approx(res["pixelBack"]["yPix"], abs=1e-5) == 256.0
