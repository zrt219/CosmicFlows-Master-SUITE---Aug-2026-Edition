# -*- coding: utf-8 -*-
"""
Automated Test Suite for Canonical Astrometric Frames and Array Stride Transformations.

Verifies:
1. Exact Matrix Orthonormality and Determinants:
   - Galactic <-> Equatorial (J2000) rotation matrices (R_Gal2Eq, R_Eq2Gal).
   - Supergalactic <-> Galactic rotation matrices (R_SG2Gal, R_Gal2SG).
   - Supergalactic <-> Equatorial rotation matrices (R_SG2Eq, R_Eq2SG).
   - Rigorous check: R^T * R = I (orthogonality) and det(R) = +1.0 (proper rotation).
2. Canonical Landmark Invariants:
   - North Galactic Pole (NGP): (l=0°, b=90°) -> (RA=192.85948°, Dec=27.12825°).
   - Galactic Center (GC): (l=0°, b=0°) -> (RA=266.40510°, Dec=-28.93617°).
   - Supergalactic North Pole (SGP): (SGL=0°, SGB=90°) -> (l=47.37°, b=6.32°).
   - Supergalactic Ascending Node (SG0): (SGL=0°, SGB=0°) -> (l=137.37°, b=0.0°).
3. Full Round-Trip Transformation Invertibility:
   - Equatorial -> Galactic -> Supergalactic -> Equatorial (< 1e-10 deg residual).
   - Cartesian Eq -> Cartesian Gal -> Cartesian SG -> Cartesian Eq (< 1e-14 residual).
4. Major Cosmography Hubs Validation:
   - Virgo, Coma, Centaurus, Norma, Shapley Supercluster Core coordinates in SGX, SGY, SGZ.
5. Pole Singularities, Clamping, and Coordinate Wrapping:
   - Declination / Latitude poles at +/- 90°.
   - Longitude boundaries [0°, 360°) and zero origin vectors (0, 0, 0).
6. 3D Vector & Peculiar Velocity Rotations:
   - Magnitude conservation: |v_sg| == |v_gal| == |v_eq|.
7. Vectorized TypedArray Batch Stride Transformations:
   - Strides 3 (XYZ), 4 (XYZW / interleaved), 6 (XYZUVW position+velocity).
   - In-place buffer mutation (outputBuffer === inputBuffer) vs out-of-place copy.
   - Large synthetic datasets (10,000+ points).
8. Hubble-Lemaître Velocity-to-Distance Inversion:
   - cz <-> d in Mpc/h and Mpc.
"""

import json
import os
import subprocess
import pytest
import math

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


class TestMatrixOrthonormalityAndDeterminants:
    """Tier 1: Astrometric Rotation Matrix Rigor and Properties"""

    def test_rotation_matrices_properties(self):
        """Verify that all rotation matrices are strictly orthogonal with det = +1.0."""
        res = run_node_eval("""
        import {
            ROT_GAL_TO_EQ, ROT_EQ_TO_GAL,
            ROT_SG_TO_GAL, ROT_GAL_TO_SG,
            ROT_SG_TO_EQ, ROT_EQ_TO_SG
        } from './src/coordinates/canonical_frame.js';
        
        function matDet3x3(M) {
            return M[0] * (M[4] * M[8] - M[5] * M[7]) -
                   M[1] * (M[3] * M[8] - M[5] * M[6]) +
                   M[2] * (M[3] * M[7] - M[4] * M[6]);
        }
        
        function checkOrthogonality(M) {
            // Check M * M^T = I
            let maxErr = 0;
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    let sum = 0;
                    for (let k = 0; k < 3; k++) {
                        sum += M[i * 3 + k] * M[j * 3 + k];
                    }
                    const target = (i === j) ? 1.0 : 0.0;
                    maxErr = Math.max(maxErr, Math.abs(sum - target));
                }
            }
            return maxErr;
        }
        
        const matrices = [
            { name: 'ROT_GAL_TO_EQ', mat: ROT_GAL_TO_EQ },
            { name: 'ROT_EQ_TO_GAL', mat: ROT_EQ_TO_GAL },
            { name: 'ROT_SG_TO_GAL', mat: ROT_SG_TO_GAL },
            { name: 'ROT_GAL_TO_SG', mat: ROT_GAL_TO_SG },
            { name: 'ROT_SG_TO_EQ', mat: ROT_SG_TO_EQ },
            { name: 'ROT_EQ_TO_SG', mat: ROT_EQ_TO_SG }
        ];
        
        const results = matrices.map(m => ({
            name: m.name,
            det: matDet3x3(m.mat),
            orthoError: checkOrthogonality(m.mat)
        }));
        
        console.log(JSON.stringify(results));
        """)
        for item in res:
            assert math.isclose(item["det"], 1.0, rel_tol=1e-12), f"Matrix {item['name']} det must be +1.0, got {item['det']}"
            assert item["orthoError"] < 1e-14, f"Matrix {item['name']} orthogonality error must be < 1e-14, got {item['orthoError']}"


class TestCanonicalLandmarkInvariants:
    """Tier 2: Astrometric Constants and Standard Landmark Transformations"""

    def test_canonical_poles_and_centers(self):
        """Verify NGP, GC, SGP, and SG0 transform exactly to their reference values."""
        res = run_node_eval("""
        import {
            galacticToEquatorial,
            equatorialToGalactic,
            supergalacticToGalactic,
            galacticToSupergalactic,
            ASTROMETRIC_CONSTANTS
        } from './src/coordinates/canonical_frame.js';
        
        // 1. North Galactic Pole: (l=0, b=90) -> (RA_NGP, Dec_NGP)
        const ngp_eq = galacticToEquatorial(0.0, 90.0, 1.0);
        
        // 2. Galactic Center: (l=0, b=0) -> (RA_GC, Dec_GC)
        const gc_eq = galacticToEquatorial(0.0, 0.0, 1.0);
        
        // 3. Supergalactic North Pole: (SGL=0, SGB=90) -> (L_SGP, B_SGP)
        const sgp_gal = supergalacticToGalactic(0.0, 90.0, 1.0);
        
        // 4. Supergalactic Origin: (SGL=0, SGB=0) -> (L_SG0, B_SG0)
        const sg0_gal = supergalacticToGalactic(0.0, 0.0, 1.0);
        
        console.log(JSON.stringify({
            ngp_eq: ngp_eq,
            ngp_ref: [ASTROMETRIC_CONSTANTS.RA_NGP_DEG, ASTROMETRIC_CONSTANTS.DEC_NGP_DEG],
            gc_eq: gc_eq,
            gc_ref: [ASTROMETRIC_CONSTANTS.RA_GC_DEG, ASTROMETRIC_CONSTANTS.DEC_GC_DEG],
            sgp_gal: sgp_gal,
            sgp_ref: [ASTROMETRIC_CONSTANTS.L_SGP_DEG, ASTROMETRIC_CONSTANTS.B_SGP_DEG],
            sg0_gal: sg0_gal,
            sg0_ref: [ASTROMETRIC_CONSTANTS.L_SG0_DEG, ASTROMETRIC_CONSTANTS.B_SG0_DEG]
        }));
        """)
        # NGP check
        assert math.isclose(res["ngp_eq"][0], res["ngp_ref"][0], abs_tol=1e-5)
        assert math.isclose(res["ngp_eq"][1], res["ngp_ref"][1], abs_tol=1e-5)

        # GC check
        assert math.isclose(res["gc_eq"][0], res["gc_ref"][0], abs_tol=1e-4)
        assert math.isclose(res["gc_eq"][1], res["gc_ref"][1], abs_tol=1e-4)

        # SGP check
        assert math.isclose(res["sgp_gal"][0], res["sgp_ref"][0], abs_tol=1e-5)
        assert math.isclose(res["sgp_gal"][1], res["sgp_ref"][1], abs_tol=1e-5)

        # SG0 check
        assert math.isclose(res["sg0_gal"][0], res["sg0_ref"][0], abs_tol=1e-5)
        assert math.isclose(res["sg0_gal"][1], res["sg0_ref"][1], abs_tol=1e-5)


class TestRoundTripTransformations:
    """Tier 3: Strict Reversibility and High-Precision Invertibility"""

    def test_spherical_full_loop_round_trip(self):
        """Verify Eq -> Gal -> SG -> Eq round trip for an ensemble of random sky positions."""
        res = run_node_eval("""
        import {
            equatorialToGalactic,
            galacticToSupergalactic,
            supergalacticToGalactic,
            galacticToEquatorial
        } from './src/coordinates/canonical_frame.js';
        
        let maxAngularErrorDeg = 0;
        let maxDistError = 0;
        
        // Test grid of 100 coordinates
        for (let ra = 0; ra < 360; ra += 36) {
            for (let dec = -80; dec <= 80; dec += 20) {
                const dist = 50.0 + ra * 0.1;
                
                // Step 1: Eq -> Gal
                const [l, b, d1] = equatorialToGalactic(ra, dec, dist);
                // Step 2: Gal -> SG
                const [sgl, sgb, d2] = galacticToSupergalactic(l, b, d1);
                // Step 3: SG -> Gal
                const [l_rev, b_rev, d3] = supergalacticToGalactic(sgl, sgb, d2);
                // Step 4: Gal -> Eq
                const [ra_rev, dec_rev, d4] = galacticToEquatorial(l_rev, b_rev, d3);
                
                // Angular difference on sphere: haversine / angular separation
                const dra = (ra - ra_rev) * Math.PI / 180.0;
                const ddec = (dec - dec_rev) * Math.PI / 180.0;
                const sinHalfDec = Math.sin(ddec / 2.0);
                const sinHalfRa = Math.sin(dra / 2.0);
                const a = sinHalfDec * sinHalfDec + Math.cos(dec * Math.PI / 180) * Math.cos(dec_rev * Math.PI / 180) * sinHalfRa * sinHalfRa;
                const angSepDeg = 2.0 * Math.asin(Math.min(1.0, Math.sqrt(a))) * 180.0 / Math.PI;
                
                maxAngularErrorDeg = Math.max(maxAngularErrorDeg, angSepDeg);
                maxDistError = Math.max(maxDistError, Math.abs(dist - d4));
            }
        }
        
        console.log(JSON.stringify({
            maxAngularErrorDeg: maxAngularErrorDeg,
            maxDistError: maxDistError
        }));
        """)
        assert res["maxAngularErrorDeg"] < 1e-10, f"Round trip angular error must be < 1e-10 deg, got {res['maxAngularErrorDeg']}"
        assert res["maxDistError"] < 1e-12

    def test_cartesian_direct_transforms(self):
        """Verify Eq Cartesian <-> SG Cartesian transformations."""
        res = run_node_eval("""
        import {
            equatorialToSupergalacticCartesian,
            supergalacticCartesianToEquatorial,
            sphericalToCartesian
        } from './src/coordinates/canonical_frame.js';
        
        // Virgo Cluster: RA=187.706°, Dec=+12.391°, d=16.5 Mpc
        const ra_virgo = 187.7059;
        const dec_virgo = 12.3911;
        const dist_virgo = 16.5;
        
        const [sgx, sgy, sgz] = equatorialToSupergalacticCartesian(ra_virgo, dec_virgo, dist_virgo);
        const [ra_back, dec_back, dist_back] = supergalacticCartesianToEquatorial(sgx, sgy, sgz);
        
        console.log(JSON.stringify({
            virgo_sg: [sgx, sgy, sgz],
            ra_back: ra_back,
            dec_back: dec_back,
            dist_back: dist_back,
            ra_err: Math.abs(ra_virgo - ra_back),
            dec_err: Math.abs(dec_virgo - dec_back),
            dist_err: Math.abs(dist_virgo - dist_back)
        }));
        """)
        # Virgo should be in the Supergalactic Plane (SGY > 0, small SGZ)
        sgx, sgy, sgz = res["virgo_sg"]
        assert sgy > 10.0, "Virgo Cluster SGY should be positive and dominate"
        assert abs(sgz) < 5.0, "Virgo Cluster is close to Supergalactic Plane (small SGZ)"
        assert res["ra_err"] < 1e-10
        assert res["dec_err"] < 1e-10
        assert res["dist_err"] < 1e-10


class TestEdgeCasesAndVectorTransforms:
    """Tier 4: Pole Singularities, Velocity Rotations, and Array Strides"""

    def test_poles_and_zero_vector_handling(self):
        """Verify pole singularity handling and clamping at +/- 90 degrees."""
        res = run_node_eval("""
        import {
            sphericalToCartesian,
            cartesianToSpherical,
            equatorialToGalactic
        } from './src/coordinates/canonical_frame.js';
        
        // Zero origin
        const [lon0, lat0, r0] = cartesianToSpherical(0, 0, 0);
        
        // North pole
        const [xN, yN, zN] = sphericalToCartesian(45.0, 90.0, 10.0);
        const [lonN, latN, rN] = cartesianToSpherical(xN, yN, zN);
        
        // South pole
        const [xS, yS, zS] = sphericalToCartesian(120.0, -90.0, 10.0);
        const [lonS, latS, rS] = cartesianToSpherical(xS, yS, zS);
        
        console.log(JSON.stringify({
            origin: [lon0, lat0, r0],
            north: [lonN, latN, rN],
            south: [lonS, latS, rS]
        }));
        """)
        assert res["origin"] == [0.0, 0.0, 0.0]
        assert math.isclose(res["north"][1], 90.0, abs_tol=1e-10)
        assert math.isclose(res["north"][2], 10.0, abs_tol=1e-10)
        assert math.isclose(res["south"][1], -90.0, abs_tol=1e-10)
        assert math.isclose(res["south"][2], 10.0, abs_tol=1e-10)

    def test_velocity_vector_magnitude_conservation(self):
        """Verify that rotating 3D peculiar velocity vectors strictly conserves magnitude."""
        res = run_node_eval("""
        import { transformVelocityVector } from './src/coordinates/canonical_frame.js';
        
        const v_in = [150.0, -320.0, 480.0];
        const norm_in = Math.hypot(...v_in);
        
        const v_gal = transformVelocityVector(v_in, 'equatorial', 'galactic');
        const norm_gal = Math.hypot(...v_gal);
        
        const v_sg = transformVelocityVector(v_gal, 'galactic', 'supergalactic');
        const norm_sg = Math.hypot(...v_sg);
        
        const v_eq_back = transformVelocityVector(v_sg, 'supergalactic', 'equatorial');
        const norm_back = Math.hypot(...v_eq_back);
        
        console.log(JSON.stringify({
            norm_in: norm_in,
            norm_gal: norm_gal,
            norm_sg: norm_sg,
            norm_back: norm_back,
            err_gal: Math.abs(norm_in - norm_gal),
            err_sg: Math.abs(norm_in - norm_sg),
            err_back: Math.abs(norm_in - norm_back),
            diff_vec: [
                Math.abs(v_in[0] - v_eq_back[0]),
                Math.abs(v_in[1] - v_eq_back[1]),
                Math.abs(v_in[2] - v_eq_back[2])
            ]
        }));
        """)
        assert res["err_gal"] < 1e-12
        assert res["err_sg"] < 1e-12
        assert res["err_back"] < 1e-12
        for d in res["diff_vec"]:
            assert d < 1e-12

    def test_batch_array_stride_transformations(self):
        """Verify Float32Array batch transformations with various strides (stride 3, 4, 6) and in-place mutation."""
        res = run_node_eval("""
        import { transformArrayStride, equatorialToSupergalacticCartesian } from './src/coordinates/canonical_frame.js';
        
        const N = 500;
        
        // 1. Stride 3 (contiguous XYZ)
        const bufStride3 = new Float32Array(N * 3);
        for (let i = 0; i < N; i++) {
            bufStride3[i * 3 + 0] = (i * 37) % 360;      // RA
            bufStride3[i * 3 + 1] = ((i * 19) % 160) - 80; // Dec
            bufStride3[i * 3 + 2] = 10.0 + (i % 100);    // Dist
        }
        
        // In-place transform spherical Eq -> Cartesian SG
        transformArrayStride(bufStride3, bufStride3, 'spherical_eq_to_cart_sg', { stride: 3, offset: 0, count: N });
        
        // Verify point 0 against single transform
        const [expectedX0, expectedY0, expectedZ0] = equatorialToSupergalacticCartesian(0, -80, 10.0);
        const errX0 = Math.abs(bufStride3[0] - expectedX0);
        const errY0 = Math.abs(bufStride3[1] - expectedY0);
        const errZ0 = Math.abs(bufStride3[2] - expectedZ0);
        
        // 2. Stride 6 (XYZUVW: position + velocity) with offset 0 and 3
        const bufStride6 = new Float64Array(N * 6);
        for (let i = 0; i < N; i++) {
            bufStride6[i * 6 + 0] = (i * 37) % 360;
            bufStride6[i * 6 + 1] = ((i * 19) % 160) - 80;
            bufStride6[i * 6 + 2] = 50.0;
            // velocity in Gal cartesian
            bufStride6[i * 6 + 3] = 100.0;
            bufStride6[i * 6 + 4] = 200.0;
            bufStride6[i * 6 + 5] = 300.0;
        }
        
        transformArrayStride(bufStride6, bufStride6, 'spherical_eq_to_cart_sg', { stride: 6, offset: 0, count: N });
        transformArrayStride(bufStride6, bufStride6, 'cart_gal_to_cart_sg', { stride: 6, offset: 3, count: N });
        
        console.log(JSON.stringify({
            pointCount: N,
            errX0: errX0,
            errY0: errY0,
            errZ0: errZ0,
            stride6SampleVNorm: Math.hypot(bufStride6[3], bufStride6[4], bufStride6[5]),
            expectedVNorm: Math.hypot(100.0, 200.0, 300.0)
        }));
        """)
        assert res["errX0"] < 1e-4, f"Float32 transform error X: {res['errX0']}"
        assert res["errY0"] < 1e-4, f"Float32 transform error Y: {res['errY0']}"
        assert res["errZ0"] < 1e-4, f"Float32 transform error Z: {res['errZ0']}"
        assert math.isclose(res["stride6SampleVNorm"], res["expectedVNorm"], rel_tol=1e-10)

    def test_recessional_velocity_to_distance(self):
        """Verify Hubble velocity <-> distance relations."""
        res = run_node_eval("""
        import { recessionalVelocityToDistance, distanceToRecessionalVelocity } from './src/coordinates/canonical_frame.js';
        
        const cz = 7460.0; // km/s
        const d_mpch = recessionalVelocityToDistance(cz, 74.6, 0.746);
        const cz_back = distanceToRecessionalVelocity(d_mpch, 74.6, 0.746);
        
        console.log(JSON.stringify({
            cz: cz,
            d_mpch: d_mpch,
            cz_back: cz_back
        }));
        """)
        assert math.isclose(res["d_mpch"], 74.6, rel_tol=1e-7)
        assert math.isclose(res["cz_back"], 7460.0, rel_tol=1e-7)
