# -*- coding: utf-8 -*-
"""
Automated Comprehensive Test Suite for Multi-Band Tully-Fisher Relation (TFR)
Extragalactic Distance Calibrator, HI Linewidth De-projection Engine, Dust Extinctions,
Metric Distance & Mpc/h Derivations, Radial Peculiar Velocities, and Malmquist Bias Compensators.

Complies with Courtois et al. 2023 (ApJ) / Cosmicflows-4 (CF4) Standards.

Verifies:
1. Cosmological & Astrometric Constants, Dipole Vectors & Units (SPEED_OF_LIGHT, H0, h, CMB_DIPOLE, LG_APEX).
2. Multi-Band Photometric Calibrations (WISE W1, W2, Spitzer [3.6], 2MASS Ks, SDSS i, r, Johnson B, Cousins I).
3. Morphological Flattening & Axial Ratio Inversion (q_0(T), cos^2 i deprojection, edge-on/face-on boundaries).
4. 21cm HI Profile Width De-projection & Kinematics (W50, W20, Tully-Fouque 1985, (1+z)^-1 redshift dilation, inclination clamping i >= 45°).
5. Galactic Foreground (SFD98 / Schlafly & Finkbeiner 2011) & Internal Dust Extinction Corrections (A_gal, A_int, K-corrections).
6. Astrometric Reference Frame Velocity Transformations (Heliocentric to CMB and Local Group).
7. Complete Tully-Fisher Distance Solver (mu_0, d in Mpc, d_h in Mpc/h, v_Hubble, v_pec, v_pec_rel).
8. Analytical Gaussian Error Propagation & Covariances (sigma_mu, sigma_d, sigma_d_h, sigma_vpec).
9. Homogeneous & Inhomogeneous Malmquist Bias (IMB) Compensations (Delta mu_hom, Delta d_IMB).
10. Baryonic Tully-Fisher Relation (BTFR) & Mass Synthesis (M_*, M_HI, M_gas, M_bary, V_rot^4 scaling).
11. Group & Cluster Multi-Galaxy Hierarchical Averaging Engine (Inverse-variance <mu>, Biweight C_BI/S_BI, Chauvenet rejection).
12. Cosmicflows-4 Benchmark Spirals (NGC 4501, NGC 1365, NGC 7331, NGC 2841, M31, M33, Circinus, Centaurus A).
13. High-Throughput Columnar Store (TypedArray Float32Array buffers, dynamic resizing, bulk statistics, CSV export).
14. IVOA / W3C PROV-O JSON-LD Provenance & Reproducibility Packager.
15. Adversarial Edge Cases, Numerical Invariants & Boundary Guards.
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
        pytest.fail(f"Node execution failed with exit code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")

    stdout = proc.stdout.strip()
    if not stdout:
        return None
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return stdout


class TestCosmologicalConstantsAndUnits:
    """Tier 1: Cosmological Constants, Dipole Motion Vectors, and Solar Band Magnitudes."""

    def test_speed_of_light_and_hubble_constants(self):
        """Verify speed of light, CF4 default H0 (74.6), and reduced Hubble parameter h (0.746)."""
        res = run_node_eval("""
        import { SPEED_OF_LIGHT_KM_S, DEFAULT_H0, DEFAULT_SIGMA_H0, DEFAULT_LITTLE_H, LN10_DIV_5, MALMQUIST_HOMOGENEOUS_COEFF } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            c: SPEED_OF_LIGHT_KM_S,
            H0: DEFAULT_H0,
            sigmaH0: DEFAULT_SIGMA_H0,
            h: DEFAULT_LITTLE_H,
            ln10Div5: LN10_DIV_5,
            malmquistHom: MALMQUIST_HOMOGENEOUS_COEFF
        }));
        """)
        assert abs(res["c"] - 299792.458) < 1e-3
        assert abs(res["H0"] - 74.6) < 1e-5
        assert abs(res["sigmaH0"] - 0.8) < 1e-5
        assert abs(res["h"] - 0.746) < 1e-5
        assert abs(res["ln10Div5"] - (math.log(10.0) / 5.0)) < 1e-12
        assert abs(res["malmquistHom"] - (-3.0 * math.log(10.0) / 5.0)) < 1e-12

    def test_cmb_and_local_group_dipole_constants(self):
        """Verify Planck 2018 / Fixsen CMB dipole apex and Local Group barycentric motion vectors."""
        res = run_node_eval("""
        import { CMB_DIPOLE, LOCAL_GROUP_APEX } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            cmb: CMB_DIPOLE,
            lg: LOCAL_GROUP_APEX
        }));
        """)
        # CMB Dipole: v = 371 km/s, l = 264.14°, b = 48.26°
        assert abs(res["cmb"]["vApex"] - 371.0) < 1e-3
        assert abs(res["cmb"]["lDeg"] - 264.14) < 1e-2
        assert abs(res["cmb"]["bDeg"] - 48.26) < 1e-2
        assert abs(res["cmb"]["raDeg"] - 168.01) < 1e-2
        assert abs(res["cmb"]["decDeg"] - (-6.98)) < 1e-2
        # Local Group: v = 316 km/s, l = 93.0°, b = -4.0°
        assert abs(res["lg"]["vApex"] - 316.0) < 1e-3
        assert abs(res["lg"]["lDeg"] - 93.0) < 1e-2
        assert abs(res["lg"]["bDeg"] - (-4.0)) < 1e-2

    def test_solar_absolute_magnitudes_spectrum(self):
        """Verify absolute solar magnitudes M_sun across NIR, mid-IR, and optical bandpasses."""
        res = run_node_eval("""
        import { SOLAR_ABSOLUTE_MAGNITUDES } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify(SOLAR_ABSOLUTE_MAGNITUDES));
        """)
        assert abs(res["W1"] - 3.24) < 0.01
        assert abs(res["W2"] - 3.27) < 0.01
        assert abs(res["SPITZER_36"] - 3.24) < 0.01
        assert abs(res["SPITZER_45"] - 3.27) < 0.01
        assert abs(res["Ks"] - 3.28) < 0.01
        assert abs(res["H"] - 3.32) < 0.01
        assert abs(res["J"] - 3.64) < 0.01
        assert abs(res["i"] - 4.53) < 0.01
        assert abs(res["r"] - 4.65) < 0.01
        assert abs(res["g"] - 5.12) < 0.01
        assert abs(res["z"] - 4.50) < 0.01
        assert abs(res["B"] - 5.48) < 0.01
        assert abs(res["V"] - 4.83) < 0.01
        assert abs(res["R"] - 4.42) < 0.01
        assert abs(res["I"] - 4.08) < 0.01

    def test_galactic_extinction_coefficients(self):
        """Verify Schlafly & Finkbeiner (2011) / SFD98 extinction reddening ratios R_lambda."""
        res = run_node_eval("""
        import { GALACTIC_EXTINCTION_COEFFICIENTS } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify(GALACTIC_EXTINCTION_COEFFICIENTS));
        """)
        assert res["W1"] < res["Ks"] < res["i"] < res["r"] < res["B"]
        assert abs(res["W1"] - 0.180) < 1e-3
        assert abs(res["W2"] - 0.160) < 1e-3
        assert abs(res["SPITZER_36"] - 0.170) < 1e-3
        assert abs(res["SPITZER_45"] - 0.150) < 1e-3
        assert abs(res["Ks"] - 0.367) < 1e-3
        assert abs(res["H"] - 0.574) < 1e-3
        assert abs(res["J"] - 0.902) < 1e-3
        assert abs(res["i"] - 1.684) < 1e-3
        assert abs(res["r"] - 2.285) < 1e-3
        assert abs(res["g"] - 3.303) < 1e-3
        assert abs(res["z"] - 1.263) < 1e-3
        assert abs(res["B"] - 4.315) < 1e-3
        assert abs(res["V"] - 3.100) < 1e-3
        assert abs(res["R"] - 2.673) < 1e-3
        assert abs(res["I"] - 1.940) < 1e-3


class TestMultiBandPhotometricCalibrations:
    """Tier 2: Multi-Band Photometric Tully-Fisher Calibrations (Direct, Inverse, and Covariances)."""

    def test_calibrations_dictionary_structure(self):
        """Check all 8 standard photometric calibrations (W1, W2, Spitzer 3.6, Ks, i, r, B, I_cousins)."""
        res = run_node_eval("""
        import { TULLY_FISHER_CALIBRATIONS } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            keys: Object.keys(TULLY_FISHER_CALIBRATIONS),
            w1: TULLY_FISHER_CALIBRATIONS.W1,
            w2: TULLY_FISHER_CALIBRATIONS.W2,
            spitzer: TULLY_FISHER_CALIBRATIONS.SPITZER_36,
            ks: TULLY_FISHER_CALIBRATIONS.Ks,
            iBand: TULLY_FISHER_CALIBRATIONS.i
        }));
        """)
        assert "W1" in res["keys"]
        assert "W2" in res["keys"]
        assert "SPITZER_36" in res["keys"]
        assert "Ks" in res["keys"]
        assert "i" in res["keys"]
        assert "r" in res["keys"]
        assert "B" in res["keys"]
        assert "I_cousins" in res["keys"]

        # Check WISE W1: slope = 9.75, zeroPoint = -21.84, pivot = 2.50
        assert abs(res["w1"]["slope"] - 9.75) < 1e-3
        assert abs(res["w1"]["zeroPoint"] - (-21.84)) < 1e-3
        assert abs(res["w1"]["pivotLogW"] - 2.50) < 1e-3
        assert abs(res["w1"]["intrinsicScatter"] - 0.35) < 1e-3

        # Check inverse slope consistency
        assert abs(res["w1"]["inverseSlope"] - (1.0 / 9.75)) < 1e-4
        assert abs(res["ks"]["inverseSlope"] - (1.0 / res["ks"]["slope"])) < 1e-4
        assert abs(res["iBand"]["inverseSlope"] - (1.0 / res["iBand"]["slope"])) < 1e-4

    def test_infrared_vs_optical_slope_steepness(self):
        """Verify that NIR / mid-IR bands (W1, Ks) have steeper slopes than optical bands (B, i)."""
        res = run_node_eval("""
        import { TULLY_FISHER_CALIBRATIONS } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            slopeW1: TULLY_FISHER_CALIBRATIONS.W1.slope,
            slopeKs: TULLY_FISHER_CALIBRATIONS.Ks.slope,
            slopeI: TULLY_FISHER_CALIBRATIONS.i.slope,
            slopeB: TULLY_FISHER_CALIBRATIONS.B.slope
        }));
        """)
        assert res["slopeW1"] > res["slopeKs"] > res["slopeI"] > res["slopeB"]

    def test_spitzer_and_wise_agreement(self):
        """Verify Spitzer 3.6um and WISE W1 (3.4um) calibrations exhibit tight physical concordance."""
        res = run_node_eval("""
        import { TULLY_FISHER_CALIBRATIONS } from './src/data/tfr_multiband_calibrator.js';
        const w1 = TULLY_FISHER_CALIBRATIONS.W1;
        const sp36 = TULLY_FISHER_CALIBRATIONS.SPITZER_36;
        console.log(JSON.stringify({
            deltaSlope: Math.abs(w1.slope - sp36.slope),
            deltaZeroPoint: Math.abs(w1.zeroPoint - sp36.zeroPoint)
        }));
        """)
        assert res["deltaSlope"] < 0.10
        assert res["deltaZeroPoint"] < 0.10


class TestMorphologicalFlatteningAndInclination:
    """Tier 3: Intrinsic Disc Flattening q_0(T) and Axial Ratio (b/a) Inversion."""

    def test_intrinsic_flattening_across_t_types(self):
        """Verify q_0 varies smoothly from 0.20 (Sa, S0) down to 0.13 (Sd, Irr)."""
        res = run_node_eval("""
        import { getIntrinsicFlattening } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            t_neg2: getIntrinsicFlattening(-2),
            t_0: getIntrinsicFlattening(0),
            t_1: getIntrinsicFlattening(1),
            t_3: getIntrinsicFlattening(3),
            t_5: getIntrinsicFlattening(5),
            t_7: getIntrinsicFlattening(7),
            t_10: getIntrinsicFlattening(10)
        }));
        """)
        assert abs(res["t_neg2"] - 0.20) < 1e-3
        assert abs(res["t_0"] - 0.20) < 1e-3
        assert res["t_1"] < 0.20
        assert res["t_3"] < res["t_1"]
        assert abs(res["t_7"] - 0.13) < 1e-3
        assert abs(res["t_10"] - 0.13) < 1e-3

    def test_axial_ratio_to_inclination_closed_form(self):
        """Verify cos^2(i) = ((b/a)^2 - q_0^2) / (1 - q_0^2) inversion."""
        res = run_node_eval("""
        import { computeInclinationFromAxisRatio } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            faceOn: computeInclinationFromAxisRatio(1.0, 4),
            edgeOn: computeInclinationFromAxisRatio(0.10, 4),
            intermediate: computeInclinationFromAxisRatio(0.50, 4)
        }));
        """)
        # Face on: b/a = 1.0 => i = 0°
        assert abs(res["faceOn"]["inclinationDeg"] - 0.0) < 1e-3
        assert abs(res["faceOn"]["cosI"] - 1.0) < 1e-3
        assert abs(res["faceOn"]["sinI"] - 0.0) < 1e-3
        assert res["faceOn"]["isFaceOn"] is True

        # Edge on: b/a <= q0 => i = 90°
        assert abs(res["edgeOn"]["inclinationDeg"] - 90.0) < 1e-3
        assert abs(res["edgeOn"]["cosI"] - 0.0) < 1e-3
        assert abs(res["edgeOn"]["sinI"] - 1.0) < 1e-3
        assert res["edgeOn"]["isEdgeOn"] is True

        # Intermediate: b/a = 0.50, q0 = 0.20 - 0.011667*4 = 0.15333
        assert 59.0 < res["intermediate"]["inclinationDeg"] < 63.0

    def test_axial_ratio_out_of_bounds_handling(self):
        """Verify handling of invalid axis ratios (e.g. b/a > 1.0 or b/a < 0.0)."""
        res = run_node_eval("""
        import { computeInclinationFromAxisRatio } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            tooBig: computeInclinationFromAxisRatio(1.5, 4),
            tooSmall: computeInclinationFromAxisRatio(-0.5, 4)
        }));
        """)
        assert res["tooBig"]["inclinationDeg"] == 0.0
        assert res["tooSmall"]["inclinationDeg"] == 90.0


class TestHIProfileWidthDeprojection:
    """Tier 4: 21cm HI Linewidth De-projection, Turbulent Motion & Inclination Clamping."""

    def test_standard_w50_deprojection_edge_on(self):
        """Verify W50 de-projection at i = 90° with standard turbulent offset (6.5 km/s)."""
        res = run_node_eval("""
        import { correctHIProfileWidth } from './src/data/tfr_multiband_calibrator.js';
        const sol = correctHIProfileWidth(450.0, {
            inclinationDeg: 90.0,
            redshiftZ: 0.0,
            profileType: 'W50'
        });
        console.log(JSON.stringify(sol));
        """)
        # wRest = 450.0, wDeTurb = 450.0 - 6.5 = 443.5 km/s
        # sin(90) = 1.0 => wMx = 443.5 km/s
        # vRot = 443.5 / 2 = 221.75 km/s
        assert abs(res["wRest"] - 450.0) < 1e-3
        assert abs(res["wDeTurb"] - 443.5) < 1e-3
        assert abs(res["wMx"] - 443.5) < 1e-3
        assert abs(res["vRot"] - 221.75) < 1e-3
        assert abs(res["logWMx"] - math.log10(443.5)) < 1e-5
        assert res["isInclinationClamped"] is False
        assert res["isValidCalibrationQuality"] is True
        assert res["qualityFlag"] == "EXCELLENT"

    def test_cosmological_redshift_dilation_stretch(self):
        """Verify (1 + z)^-1 cosmological time dilation correction on spectral linewidth."""
        res = run_node_eval("""
        import { correctHIProfileWidth } from './src/data/tfr_multiband_calibrator.js';
        const sol = correctHIProfileWidth(440.0, {
            inclinationDeg: 90.0,
            redshiftZ: 0.10,
            profileType: 'W50',
            turbulentVelocityDispersion: 0.0
        });
        console.log(JSON.stringify(sol));
        """)
        # wRest = 440.0 / 1.10 = 400.0 km/s
        assert abs(res["wRest"] - 400.0) < 1e-3
        assert abs(res["wMx"] - 400.0) < 1e-3

    def test_instrumental_channel_broadening_deconvolution(self):
        """Verify instrumental channel width deconvolution Delta v_inst = 2 * Delta v_chan * eta."""
        res = run_node_eval("""
        import { correctHIProfileWidth } from './src/data/tfr_multiband_calibrator.js';
        const sol = correctHIProfileWidth(400.0, {
            inclinationDeg: 90.0,
            channelWidthKmS: 10.0,
            instrumentalEta: 0.5,
            turbulentVelocityDispersion: 0.0
        });
        console.log(JSON.stringify(sol));
        """)
        # Delta v_inst = 2 * 10 * 0.5 = 10 km/s => wInst = 390 km/s
        assert abs(res["wInst"] - 390.0) < 1e-3
        assert abs(res["wMx"] - 390.0) < 1e-3

    def test_tully_fouque_1985_non_linear_w20_deprojection(self):
        """Verify Tully-Fouque (1985) non-linear quadratic formulation for W20 profile."""
        res = run_node_eval("""
        import { correctHIProfileWidth } from './src/data/tfr_multiband_calibrator.js';
        const sol = correctHIProfileWidth(480.0, {
            inclinationDeg: 90.0,
            profileType: 'W20'
        });
        console.log(JSON.stringify(sol));
        """)
        # Should execute quadratic formula with W_t = 38.0 km/s, W_c = 120.0 km/s
        assert res["wMx"] > 0
        assert res["wMx"] < 480.0
        assert math.isfinite(res["logWMx"])

    def test_inclination_clamping_at_45_degrees(self):
        """Verify strict Courtois et al. 2023 inclination clamping for i < 45°."""
        res = run_node_eval("""
        import { correctHIProfileWidth } from './src/data/tfr_multiband_calibrator.js';
        const solClamped = correctHIProfileWidth(300.0, {
            inclinationDeg: 35.0,
            minInclinationDeg: 45.0,
            allowExtrapolatedInclination: false
        });
        const solExtrap = correctHIProfileWidth(300.0, {
            inclinationDeg: 35.0,
            minInclinationDeg: 45.0,
            allowExtrapolatedInclination: true
        });
        const solUnreliable = correctHIProfileWidth(300.0, {
            inclinationDeg: 20.0,
            minInclinationDeg: 45.0,
            allowExtrapolatedInclination: false
        });
        console.log(JSON.stringify({ clamped: solClamped, extrap: solExtrap, unreliable: solUnreliable }));
        """)
        # Clamped: appliedInclination = 45.0°, sin(45°) = 0.7071
        assert abs(res["clamped"]["appliedInclinationDeg"] - 45.0) < 1e-3
        assert res["clamped"]["isInclinationClamped"] is True
        assert res["clamped"]["isValidCalibrationQuality"] is False
        assert res["clamped"]["qualityFlag"] == "CLAMPED_BELOW_MIN"
        assert abs(res["clamped"]["sinI"] - math.sin(math.radians(45.0))) < 1e-4

        # Extrapolated: appliedInclination = 35.0°
        assert abs(res["extrap"]["appliedInclinationDeg"] - 35.0) < 1e-3
        assert abs(res["extrap"]["sinI"] - math.sin(math.radians(35.0))) < 1e-4

        # Unreliable: raw inclination < 30°
        assert res["unreliable"]["qualityFlag"] == "UNRELIABLE"


class TestDustExtinctionAndCorrections:
    """Tier 5: Galactic Foreground, Internal Dust Extinction & Cosmological K-corrections."""

    def test_galactic_extinction_scaling(self):
        """Verify A_gal = R_band * E(B-V) across bands."""
        res = run_node_eval("""
        import { computeGalacticExtinction } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            aGal_W1: computeGalacticExtinction(0.10, 'W1'),
            aGal_Ks: computeGalacticExtinction(0.10, 'Ks'),
            aGal_i: computeGalacticExtinction(0.10, 'i'),
            aGal_B: computeGalacticExtinction(0.10, 'B')
        }));
        """)
        assert abs(res["aGal_W1"] - 0.018) < 1e-4
        assert abs(res["aGal_Ks"] - 0.0367) < 1e-4
        assert abs(res["aGal_i"] - 0.1684) < 1e-4
        assert abs(res["aGal_B"] - 0.4315) < 1e-4

    def test_internal_dust_extinction_inclination_dependence(self):
        """Verify internal dust extinction increases with inclination (log10(1/cos i))."""
        res = run_node_eval("""
        import { computeInternalDustExtinction } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            aInt_W1_i0: computeInternalDustExtinction(0.0, 2.50, 'W1'),
            aInt_W1_i60: computeInternalDustExtinction(60.0, 2.50, 'W1'),
            aInt_W1_i80: computeInternalDustExtinction(80.0, 2.50, 'W1'),
            aInt_i_i60: computeInternalDustExtinction(60.0, 2.50, 'i'),
            aInt_i_i80: computeInternalDustExtinction(80.0, 2.50, 'i')
        }));
        """)
        # Face on: log10(1 / cos(0)) = log10(1) = 0.0 => A_int = 0.0
        assert abs(res["aInt_W1_i0"] - 0.0) < 1e-6
        # At i=60°: cos(60) = 0.5 => log10(2) = 0.30103. W1 gamma = 0.08 => A_int = 0.08 * 0.30103 = 0.02408 mag
        assert abs(res["aInt_W1_i60"] - 0.02408) < 1e-3
        # Optical i-band has much higher internal extinction than W1
        assert res["aInt_i_i60"] > 5.0 * res["aInt_W1_i60"]
        assert res["aInt_i_i80"] > res["aInt_i_i60"]

    def test_k_correction_linear_scaling(self):
        """Verify K-correction for low redshifts."""
        res = run_node_eval("""
        import { computeKCorrection } from './src/data/tfr_multiband_calibrator.js';
        console.log(JSON.stringify({
            k_w1_z01: computeKCorrection(0.01, 'W1'),
            k_w1_z05: computeKCorrection(0.05, 'W1'),
            k_i_z01: computeKCorrection(0.01, 'i')
        }));
        """)
        # W1: k1 = 0.20 => K(0.01) = 0.002, K(0.05) = 0.010
        assert abs(res["k_w1_z01"] - 0.002) < 1e-5
        assert abs(res["k_w1_z05"] - 0.010) < 1e-5
        assert res["k_i_z01"] > res["k_w1_z01"]


class TestAstrometricVelocityFrames:
    """Tier 6: Heliocentric to CMB and Local Group Velocity Reference Frames."""

    def test_cmb_velocity_apex_boost_and_anti_apex(self):
        """Verify velocity projection along CMB dipole apex (l=264.14°, b=48.26°)."""
        res = run_node_eval("""
        import { heliocentricToCMBVelocity, CMB_DIPOLE } from './src/data/tfr_multiband_calibrator.js';
        const vApex = heliocentricToCMBVelocity(1000.0, CMB_DIPOLE.lDeg, CMB_DIPOLE.bDeg);
        const vAntiApex = heliocentricToCMBVelocity(1000.0, CMB_DIPOLE.lDeg + 180.0, -CMB_DIPOLE.bDeg);
        console.log(JSON.stringify({ vApex, vAntiApex }));
        """)
        # Apex boost: 1000 + 371 = 1371 km/s
        assert abs(res["vApex"] - 1371.0) < 1e-2
        # Anti-apex reduction: 1000 - 371 = 629 km/s
        assert abs(res["vAntiApex"] - 629.0) < 1e-2

    def test_local_group_velocity_transformation(self):
        """Verify Local Group barycentre apex frame transformation."""
        res = run_node_eval("""
        import { heliocentricToLocalGroupVelocity, LOCAL_GROUP_APEX } from './src/data/tfr_multiband_calibrator.js';
        const vLGApex = heliocentricToLocalGroupVelocity(500.0, LOCAL_GROUP_APEX.lDeg, LOCAL_GROUP_APEX.bDeg);
        console.log(JSON.stringify({ vLGApex }));
        """)
        # Apex boost: 500 + 316 = 816 km/s
        assert abs(res["vLGApex"] - 816.0) < 1e-2


class TestTullyFisherSolverAndMetricDistances:
    """Tier 7: Full Multi-Band TFR Solver, Distance Modulus, Mpc, and Peculiar Velocities."""

    def test_canonical_wise_w1_calibrator_solution(self):
        """Verify end-to-end analytical math for standard benchmark input in W1."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({
            wObs: 450.0,
            inclinationDeg: 90.0,
            mApp: 8.12,
            band: 'W1',
            aExt: 0.15,
            czCMB: 528.0,
            H0: 74.6
        });
        console.log(JSON.stringify(sol));
        """)
        # 1. z = 528 / 299792.458 = 0.0017612
        # 2. wRest = 450 / 1.0017612 = 449.2088 km/s
        # 3. wMx = 449.2088 - 6.5 = 442.7088 km/s
        # 4. log10(442.7088) = 2.646123
        # 5. mAbs = -9.75 * (2.646123 - 2.5) - 21.84 = -23.2647 mag
        # 6. mu = 7.9696 - (-23.2647) = 31.234 mag
        # 7. distMpc = 10^((31.234 - 25)/5) = 17.65 Mpc
        # 8. vPec = 528 - 17.65 * 74.6 = -789.2 km/s (Infall)
        assert abs(res["wCorr"] - 442.71) < 0.2
        assert abs(res["vRot"] - 221.35) < 0.2
        assert abs(res["mu"] - 31.10) < 0.05
        assert abs(res["distMpc"] - 16.57) < 0.2
        assert abs(res["distMpcPerH"] - 12.36) < 0.2
        assert res["vPec"] < 0  # Infall confirmation

    def test_multi_band_solution_consistency(self):
        """Verify solver across W1, W2, Spitzer 3.6, Ks, i, and B bands for the same galaxy."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const bands = ['W1', 'W2', 'SPITZER_36', 'Ks', 'i', 'B'];
        const results = {};
        for (const b of bands) {
            results[b] = solveTullyFisher({
                wObs: 400.0,
                inclinationDeg: 80.0,
                mApp: 10.0,
                band: b,
                czCMB: 1500.0,
                H0: 75.0
            });
        }
        console.log(JSON.stringify(results));
        """)
        for b in ['W1', 'W2', 'SPITZER_36', 'Ks', 'i', 'B']:
            sol = res[b]
            assert sol["distMpc"] > 0
            assert sol["distMpcPerH"] > 0
            assert math.isfinite(sol["vPec"])
            assert sol["sigmaMu"] > 0
            assert sol["sigmaDistMpc"] > 0
            assert sol["sigmaVPec"] > 0

    def test_custom_calibration_override(self):
        """Verify custom calibration object overrides default tables."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const customCal = {
            name: 'Custom Synthetic Calibrator',
            slope: 10.0,
            zeroPoint: -22.0,
            pivotLogW: 2.50,
            intrinsicScatter: 0.20,
            sigmaZeroPoint: 0.02,
            internalExtinctionGamma0: 0.0,
            internalExtinctionGamma1: 0.0,
            kCorrCoeff1: 0.0
        };
        const sol = solveTullyFisher({
            wObs: 316.227766, // log10(W) = 2.50 exactly
            inclinationDeg: 90.0,
            mApp: 10.0,
            aExt: 0.0,
            customCalibration: customCal,
            turbulentVelocityDispersion: 0.0,
            redshiftZ: 0.0
        });
        console.log(JSON.stringify(sol));
        """)
        # At log10(W) = 2.50, mAbs = -22.0 exactly
        # mu = 10.0 - (-22.0) = 32.0
        # distMpc = 10^((32 - 25)/5) = 10^1.4 = 25.11886 Mpc
        assert abs(res["mAbs"] - (-22.0)) < 1e-3
        assert abs(res["mu"] - 32.0) < 1e-3
        assert abs(res["distMpc"] - 25.11886) < 1e-3


class TestGaussianErrorPropagation:
    """Tier 8: Analytical Gaussian Error Propagation & Covariance Derivations."""

    def test_distance_modulus_and_metric_distance_uncertainty(self):
        """Verify sigma_d / d = (ln(10) / 5) * sigma_mu relationship."""
        res = run_node_eval("""
        import { solveTullyFisher, LN10_DIV_5 } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({
            wObs: 400.0,
            wObsErr: 15.0,
            mApp: 9.0,
            mAppErr: 0.08,
            band: 'W1',
            H0: 74.6
        });
        console.log(JSON.stringify({
            sigmaMu: sol.sigmaMu,
            distMpc: sol.distMpc,
            sigmaDistMpc: sol.sigmaDistMpc,
            expectedSigmaDist: sol.distMpc * LN10_DIV_5 * sol.sigmaMu
        }));
        """)
        assert abs(res["sigmaDistMpc"] - res["expectedSigmaDist"]) < 1e-6
        assert res["sigmaMu"] > 0.35  # Must exceed intrinsic scatter floor

    def test_peculiar_velocity_uncertainty_propagation(self):
        """Verify sigma_vpec^2 = sigma_cz^2 + (H0 * sigma_d)^2 + (d * sigma_H0)^2."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({
            wObs: 450.0,
            mApp: 8.5,
            czCMB: 2000.0,
            czErr: 20.0,
            H0: 74.6,
            sigmaH0: 1.0
        });
        const term1 = 20.0 * 20.0;
        const term2 = Math.pow(74.6 * sol.sigmaDistMpc, 2);
        const term3 = Math.pow(sol.distMpc * 1.0, 2);
        const expectedSigmaVPec = Math.sqrt(term1 + term2 + term3);
        console.log(JSON.stringify({
            sigmaVPec: sol.sigmaVPec,
            expectedSigmaVPec
        }));
        """)
        assert abs(res["sigmaVPec"] - res["expectedSigmaVPec"]) < 1e-4


class TestMalmquistBiasCompensators:
    """Tier 9: Homogeneous and Inhomogeneous Malmquist Bias (IMB) Compensations."""

    def test_homogeneous_malmquist_bias_formula(self):
        """Verify Delta mu_hom = - 1.381551 * sigma_mu^2 analytic formula."""
        res = run_node_eval("""
        import { solveTullyFisher, MALMQUIST_HOMOGENEOUS_COEFF } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({ wObs: 400.0, mApp: 9.0 });
        const expectedDeltaMu = MALMQUIST_HOMOGENEOUS_COEFF * (sol.sigmaMu * sol.sigmaMu);
        console.log(JSON.stringify({
            actualDeltaMu: sol.malmquistHomogeneousDeltaMu,
            expectedDeltaMu,
            distRaw: sol.distMpc,
            distHom: sol.distMpcMalmquistHomogeneous
        }));
        """)
        assert abs(res["actualDeltaMu"] - res["expectedDeltaMu"]) < 1e-6
        # Distance after homogeneous Malmquist correction must be slightly smaller than raw distance
        assert res["distHom"] < res["distRaw"]

    def test_inhomogeneous_malmquist_bias_positive_and_negative_density_gradients(self):
        """Verify Delta d_IMB = - sigma_d^2 * d(ln n)/dr pulls distance towards high density."""
        res = run_node_eval("""
        import { computeInhomogeneousMalmquistBias } from './src/data/tfr_multiband_calibrator.js';
        // Approaching cluster: d(ln n)/dr > 0 => correction should decrease distance
        const resInfall = computeInhomogeneousMalmquistBias(20.0, 2.0, 1500.0, 74.6, {
            densityGradientFunc: (r) => 0.10 // +0.10 Mpc^-1
        });
        // Past cluster peak: d(ln n)/dr < 0 => correction should increase distance
        const resOutfall = computeInhomogeneousMalmquistBias(20.0, 2.0, 1500.0, 74.6, {
            densityGradientFunc: (r) => -0.10 // -0.10 Mpc^-1
        });
        console.log(JSON.stringify({ infall: resInfall, outfall: resOutfall }));
        """)
        # Delta d = - (2.0)^2 * 0.10 = -0.40 Mpc
        assert abs(res["infall"]["deltaDistMpc"] - (-0.40)) < 1e-5
        assert abs(res["infall"]["correctedDistMpc"] - 19.60) < 1e-5
        # Delta d = - (2.0)^2 * (-0.10) = +0.40 Mpc
        assert abs(res["outfall"]["deltaDistMpc"] - 0.40) < 1e-5
        assert abs(res["outfall"]["correctedDistMpc"] - 20.40) < 1e-5

    def test_inhomogeneous_malmquist_bias_power_law_model(self):
        """Verify power-law background model n(r) ~ r^-gamma."""
        res = run_node_eval("""
        import { computeInhomogeneousMalmquistBias } from './src/data/tfr_multiband_calibrator.js';
        const resPower = computeInhomogeneousMalmquistBias(50.0, 5.0, 3700.0, 74.6, {
            densityPowerLawGamma: 1.5 // d(ln n)/dr = -1.5 / 50 = -0.03 Mpc^-1
        });
        console.log(JSON.stringify(resPower));
        """)
        # Delta d = - (5.0)^2 * (-0.03) = +0.75 Mpc
        assert abs(res["deltaDistMpc"] - 0.75) < 1e-3
        assert abs(res["correctedDistMpc"] - 50.75) < 1e-3


class TestBaryonicTullyFisherRelation:
    """Tier 10: Baryonic Tully-Fisher Relation (BTFR) & Component Mass Synthesis."""

    def test_btfr_power_law_scaling_exponent_4(self):
        """Verify M_bary = 47.0 * V_rot^4.0 power-law scaling."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol100 = solveTullyFisher({ wObs: 206.5, inclinationDeg: 90.0, redshiftZ: 0.0 }); // vRot = 100 km/s
        const sol200 = solveTullyFisher({ wObs: 406.5, inclinationDeg: 90.0, redshiftZ: 0.0 }); // vRot = 200 km/s
        console.log(JSON.stringify({
            mb100: sol100.mBaryonic,
            mb200: sol200.mBaryonic,
            ratio: sol200.mBaryonic / sol100.mBaryonic
        }));
        """)
        # Doubling V_rot should increase M_bary by exactly 2^4 = 16
        assert abs(res["mb100"] - 4.7e9) < 1e7
        assert abs(res["mb200"] - 7.52e10) < 1e8
        assert abs(res["ratio"] - 16.0) < 1e-3

    def test_hi_gas_mass_synthesis_from_flux(self):
        """Verify M_HI = 2.356e5 * d_Mpc^2 * S_HI and M_gas = 1.33 * M_HI."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({
            wObs: 450.0,
            mApp: 8.0,
            sHI: 100.0 // Jy * km/s
        });
        const expectedMHI = 2.356e5 * (sol.distMpc * sol.distMpc) * 100.0;
        const expectedMGas = 1.33 * expectedMHI;
        console.log(JSON.stringify({
            mGas: sol.mGas,
            expectedMGas,
            mStar: sol.mStar,
            mBaryonic: sol.mBaryonic
        }));
        """)
        assert abs(res["mGas"] - res["expectedMGas"]) < 1e-2
        assert res["mStar"] > 0
        assert res["mBaryonic"] > 0


class TestGroupClusterHierarchicalAverager:
    """Tier 11: Multi-Galaxy Cluster/Group Distance Aggregation & Outlier Rejection."""

    def test_synthetic_virgo_cluster_group_aggregation(self):
        """Verify weighted mean modulus and biweight location across 8 member galaxies."""
        res = run_node_eval("""
        import { aggregateGroupTFR } from './src/data/tfr_multiband_calibrator.js';
        // 7 true Virgo members (mu ~ 31.0) + 1 extreme foreground outlier (mu ~ 27.0)
        const members = [
            { wObs: 450, mApp: 7.73, band: 'W1', czCMB: 1200 },
            { wObs: 420, mApp: 7.95, band: 'W1', czCMB: 1150 },
            { wObs: 480, mApp: 7.50, band: 'W1', czCMB: 1300 },
            { wObs: 400, mApp: 8.10, band: 'W1', czCMB: 1100 },
            { wObs: 460, mApp: 7.65, band: 'W1', czCMB: 1250 },
            { wObs: 380, mApp: 8.30, band: 'W1', czCMB: 1180 },
            { wObs: 440, mApp: 7.80, band: 'W1', czCMB: 1220 },
            { wObs: 450, mApp: 3.73, band: 'W1', czCMB: 200 }   // Outlier: mu = 27.0
        ];
        const groupSummary = aggregateGroupTFR('Virgo Cluster Core Test', members, { H0: 74.6 });
        console.log(JSON.stringify(groupSummary));
        """)
        # Outlier should be rejected (accepted: 7, rejected: 1)
        assert res["memberCount"] == 8
        assert res["acceptedCount"] == 7
        assert res["rejectedCount"] == 1
        # Mean modulus should be around 30.8 - 31.0 mag
        assert abs(res["weightedMeanMu"] - 30.80) < 0.2
        assert 13.0 < res["meanDistMpc"] < 17.0
        assert res["reducedChiSquare"] < 2.0
        assert res["sigmaWeightedMeanMu"] < 0.20


class TestBenchmarkGalaxyPresets:
    """Tier 12: Cosmicflows-4 Primary Calibration Galaxy Presets."""

    def test_benchmark_presets_availability_and_properties(self):
        """Verify all 8 CF4 benchmark calibration spirals exist and solve accurately."""
        res = run_node_eval("""
        import { TULLY_FISHER_BENCHMARK_PRESETS, solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const solved = {};
        for (const [key, preset] of Object.entries(TULLY_FISHER_BENCHMARK_PRESETS)) {
            solved[key] = solveTullyFisher(preset);
        }
        console.log(JSON.stringify({
            keys: Object.keys(TULLY_FISHER_BENCHMARK_PRESETS),
            solved
        }));
        """)
        required_keys = ["M31", "M33", "NGC4501", "NGC1365", "NGC7331", "NGC2841", "CIRCINUS", "NGC5128"]
        for k in required_keys:
            assert k in res["keys"]
            sol = res["solved"][k]
            assert sol["distMpc"] > 0
            assert sol["wCorr"] > 0
            assert math.isfinite(sol["vPec"])

        # Specific physical checks:
        # M31 / Andromeda: nearby Local Group, distance ~ 0.7 - 0.9 Mpc
        assert 0.5 < res["solved"]["M31"]["distMpc"] < 1.2
        # NGC 4501 in Virgo: distance ~ 14 - 19 Mpc
        assert 13.0 < res["solved"]["NGC4501"]["distMpc"] < 20.0
        # NGC 1365 in Fornax: distance ~ 16 - 22 Mpc
        assert 15.0 < res["solved"]["NGC1365"]["distMpc"] < 23.0


class TestColumnarStoreAndBatchProcessor:
    """Tier 13: High-Density TypedArray Columnar Store, Bulk Statistics & CSV Export."""

    def test_columnar_store_push_growth_and_statistics(self):
        """Verify TFRColumnarStore capacity growth and statistical aggregations."""
        res = run_node_eval("""
        import { TFRColumnarStore, solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const store = new TFRColumnarStore(4); // Start small to trigger _grow
        for (let i = 0; i < 20; i++) {
            const sol = solveTullyFisher({
                wObs: 400 + i * 5,
                mApp: 8.0 + i * 0.1,
                czCMB: 1000 + i * 50
            });
            store.push(sol, i * 1.5, i * 2.0, i * 0.5);
        }
        const stats = store.getStatistics();
        console.log(JSON.stringify({
            count: store.count,
            capacity: store.capacity,
            stats
        }));
        """)
        assert res["count"] == 20
        assert res["capacity"] >= 20
        assert res["stats"]["sampleSize"] == 20
        assert res["stats"]["meanDistMpc"] > 0
        assert math.isfinite(res["stats"]["meanVPec"])

    def test_batch_processor_and_csv_generation(self):
        """Verify batchSolveTullyFisher and CSV serialization."""
        res = run_node_eval("""
        import { batchSolveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const catalog = [
            { wObs: 450, mApp: 8.12, raDeg: 180.0, decDeg: 10.0, czCMB: 1200 },
            { wObs: 380, mApp: 9.00, raDeg: 45.0, decDeg: -30.0, czCMB: 1600 }
        ];
        const store = batchSolveTullyFisher(catalog);
        const csv = store.toCSV();
        console.log(JSON.stringify({
            count: store.count,
            csvLines: csv.split('\\n')
        }));
        """)
        assert res["count"] == 2
        assert len(res["csvLines"]) == 3  # Header + 2 rows
        assert "dist_mpc" in res["csvLines"][0]
        assert "v_pec_kms" in res["csvLines"][0]


class TestProvenanceAndReproducibility:
    """Tier 14: IVOA / W3C PROV-O JSON-LD Metadata Generator."""

    def test_jsonld_provenance_schema(self):
        """Verify structured JSON-LD entity generation."""
        res = run_node_eval("""
        import { solveTullyFisher, generateTFRProvenanceJSONLD } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({ wObs: 450, mApp: 8.12 });
        const jsonld = generateTFRProvenanceJSONLD(sol, { catalogVersion: 'CF4-v1.0' });
        console.log(JSON.stringify(jsonld));
        """)
        assert res["@type"] == "Entity"
        assert "prov:wasGeneratedBy" in res
        assert "cf4:distanceMpc" in res
        assert "cf4:distanceModulus" in res
        assert "cf4:peculiarVelocityKmS" in res
        assert res["cf4:metadata"]["catalogVersion"] == "CF4-v1.0"


class TestAdversarialAndBoundaryConditions:
    """Tier 15: Extreme Values, Zero Division Protection, and Cosmological Invariants."""

    def test_extreme_low_and_high_linewidths(self):
        """Verify numerical stability for dwarf galaxies (W=30 km/s) and giant spirals (W=1000 km/s)."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const solDwarf = solveTullyFisher({ wObs: 30.0, mApp: 15.0 });
        const solGiant = solveTullyFisher({ wObs: 1000.0, mApp: 5.0 });
        console.log(JSON.stringify({ dwarf: solDwarf, giant: solGiant }));
        """)
        assert res["dwarf"]["distMpc"] > 0
        assert math.isfinite(res["dwarf"]["mu"])
        assert res["giant"]["distMpc"] > 0
        assert math.isfinite(res["giant"]["mu"])
        # Giant spiral must be significantly brighter in absolute magnitude than dwarf
        assert res["giant"]["mAbs"] < res["dwarf"]["mAbs"]

    def test_face_on_inclination_clamping_safety(self):
        """Verify face-on galaxy with i = 5° does not cause 1/sin(0) infinity or NaN."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({ wObs: 300.0, inclinationDeg: 5.0 });
        console.log(JSON.stringify(sol));
        """)
        assert res["isInclinationClamped"] is True
        assert res["appliedInclinationDeg"] == 45.0
        assert math.isfinite(res["wCorr"])
        assert math.isfinite(res["distMpc"])
        assert res["distMpc"] > 0

    def test_negative_recession_velocity_blueshift(self):
        """Verify blueshifted approaching galaxies (e.g. M31 cz = -300 km/s) yield negative peculiar velocity."""
        res = run_node_eval("""
        import { solveTullyFisher } from './src/data/tfr_multiband_calibrator.js';
        const sol = solveTullyFisher({
            wObs: 510.0,
            mApp: 0.98,
            czCMB: -150.0,
            H0: 74.6
        });
        console.log(JSON.stringify(sol));
        """)
        assert res["vPec"] < 0
        assert math.isfinite(res["vPecRelativistic"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
