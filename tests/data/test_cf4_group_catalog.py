# -*- coding: utf-8 -*-
"""
Automated Test Suite for Cosmicflows-4 Grouped Galaxy Catalog (CF4gp),
Hierarchical Friends-of-Friends (FoF), Projected Phase-Space Kinematic Engine,
and Multi-Estimator Virial Mass Suite.

Verifies:
1. Cosmological Constants, Units & Coordinate Conversions:
   - Universal Gravitational Constant G_cosmo in Mpc*(km/s)^2 / M_sun.
   - Critical density rho_crit(z) at z=0 and higher redshifts.
   - Hubble expansion parameter H(z) for flat LambdaCDM.
   - Great-circle Vincenty angular separation and projected physical separations.
2. Master Cosmological Structures Reference Catalog:
   - Verification of 11 named structures (Laniakea Core, Coma, Virgo, Centaurus, Norma, Shapley Core, etc.).
   - Cross-matching algorithm in 3D Supergalactic Cartesian and velocity space.
3. Robust Statistical Estimators (Beers, Flynn & Gebhardt 1990):
   - Median, MAD, sample standard deviation.
   - Biweight Location (C_BI) vs sample mean with extreme outlier rejection.
   - Biweight Scale (S_BI) vs sample standard deviation.
   - Beers Gapper Scale Estimator for small N (N <= 15).
   - Measurement error deconvolution and cosmological (1+z)^-1 redshift correction.
4. Projected Harmonic Radius & Radii Metrics (Limber & Mathews 1960):
   - Pairwise projected separation matrix.
   - Harmonic radius R_H on regular polygons against closed analytic formulas.
   - Projected moments (R_mean, R_rms, R_med, R_20, R_80, concentration index).
   - Virial radius R_200 and R_500 against critical overdensity criteria.
5. Multi-Estimator Virial Mass Suite (Heisler et al. 1985, Bahcall & Tremaine 1981):
   - Standard Virial Mass: M_vir = (3 * pi / 2) * (R_H * sigma_v^2 / G).
   - Projected Mass Estimator M_proj (isotropic, radial, circular orbits).
   - Median Mass Estimator M_med and Average Mass Estimator M_avg.
   - Bahcall & Tremaine central dominant point-mass estimator.
   - Surface pressure correction factor (The & White 1986).
   - Non-parametric Bootstrap Resampling (B >= 100) for asymmetric confidence intervals.
6. Group Luminosity, Schechter Completeness, and Morphology:
   - Lanczos log-gamma and standard gamma functions.
   - Upper incomplete gamma Gamma(a, x) via series and continued fractions.
   - Schechter (1976) faint-end completeness integration C(L_lim).
   - Total group luminosity in K-band and B-band.
   - de Vaucouleurs T-type morphological fractions (f_early, f_late, f_irr).
   - Mass-to-light ratios (M_vir / L_K, M_vir / L_B).
   - Dressler (1980) local projected surface density Sigma_10.
7. Hierarchical Friends-of-Friends (FoF) & Membership Probabilities:
   - DisjointSet (Union-Find) with path compression and rank optimization.
   - Adaptive linking lengths D_link(cz) and V_link(cz) (Tully 2015, Kourkchi & Tully 2017).
   - Projected phase-space membership probabilities P_mem.
   - Multi-group percolation and cluster recovery on synthetic galaxy fields.
8. CF4 Grouped Catalog Columnar Store & Ingestion:
   - High-density TypedArray columnar store and dynamic resizing.
   - Group aggregation from member records.
   - 3D Uniform Spatial Hash Grid and cone queries.
   - Parsers for CSV, TSV, IP2I ASCII, and VizieR tabular data.
   - GeoJSON and CSV export validation.
9. Scientific Invariants, Boundary Conditions, and Adversarial Edge Cases:
   - Single galaxy groups (N = 1).
   - Binary pairs (N = 2) with exact pairwise velocity dispersion.
   - Rich clusters (N = 500+).
   - Zero velocity dispersion and identical coordinates.
   - Outlier resistance and numerical stability.
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


class TestCosmologicalConstantsAndCoordinates:
    """Tier 1: Cosmological Constants, Units, and Spherical / Cartesian Astrometry"""

    def test_gravitational_constant_and_critical_density(self):
        """Verify cosmological gravitational constant G and critical density rho_crit(z=0)."""
        res = run_node_eval("""
        import { G_COSMO, RHO_CRIT_0, DEFAULT_H0, criticalDensityAtZ } from './src/data/cf4_group_catalog.js';

        // SI calculation: G = 6.67430e-11 m^3 / (kg * s^2)
        // 1 Mpc = 3.085677581e22 m, 1 M_sun = 1.98847e30 kg
        // G_cosmo in Mpc * (km/s)^2 / M_sun
        const G_SI = 6.67430e-11;
        const MPC_M = 3.085677581e22;
        const MSUN_KG = 1.98847e30;
        const KM_M = 1.0e3;

        const G_expected = (G_SI * MSUN_KG) / (MPC_M * KM_M * KM_M);
        const rho_crit_expected = (3.0 * DEFAULT_H0 * DEFAULT_H0) / (8.0 * Math.PI * G_COSMO);
        const rho_crit_z0 = criticalDensityAtZ(0.0, DEFAULT_H0);

        console.log(JSON.stringify({
            G_COSMO,
            G_expected,
            RHO_CRIT_0,
            rho_crit_expected,
            rho_crit_z0
        }));
        """)
        assert math.isclose(res["G_COSMO"], res["G_expected"], rel_tol=1e-4)
        assert math.isclose(res["RHO_CRIT_0"], res["rho_crit_expected"], rel_tol=1e-7)
        assert math.isclose(res["RHO_CRIT_0"], res["rho_crit_z0"], rel_tol=1e-7)
        # Check standard cosmological range: ~1.54e11 M_sun/Mpc^3 for H0=74.6
        assert 1.4e11 < res["RHO_CRIT_0"] < 1.7e11

    def test_hubble_parameter_and_redshift_evolution(self):
        """Verify H(z) evolution in flat LambdaCDM: H(z) = H0 * sqrt(Omega_m * (1+z)^3 + Omega_L)."""
        res = run_node_eval("""
        import { hubbleParameter, criticalDensityAtZ, DEFAULT_H0, DEFAULT_OMEGA_M, DEFAULT_OMEGA_LAMBDA } from './src/data/cf4_group_catalog.js';

        const redshifts = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0];
        const results = redshifts.map(z => {
            const hz = hubbleParameter(z, DEFAULT_H0, DEFAULT_OMEGA_M, DEFAULT_OMEGA_LAMBDA);
            const expectedHz = DEFAULT_H0 * Math.sqrt(DEFAULT_OMEGA_M * Math.pow(1.0 + z, 3) + DEFAULT_OMEGA_LAMBDA);
            const rhoCritZ = criticalDensityAtZ(z, DEFAULT_H0);
            return { z, hz, expectedHz, rhoCritZ };
        });

        console.log(JSON.stringify(results));
        """)
        for item in res:
            assert math.isclose(item["hz"], item["expectedHz"], rel_tol=1e-7)
            assert item["rhoCritZ"] > 0
            if item["z"] > 0:
                assert item["hz"] > res[0]["hz"]
                assert item["rhoCritZ"] > res[0]["rhoCritZ"]

    def test_great_circle_vincenty_angular_separation(self):
        """Verify Great-Circle angular separation formula across poles, equator, and known pairs."""
        res = run_node_eval("""
        import { angularSeparationDeg, angularSeparationRad, DEG_TO_RAD } from './src/data/cf4_group_catalog.js';

        // 1. Same point: sep = 0
        const dZero = angularSeparationDeg(120.0, 30.0, 120.0, 30.0);

        // 2. North pole to South pole: sep = 180 deg
        const dPoles = angularSeparationDeg(0.0, 90.0, 180.0, -90.0);

        // 3. Equator quarter-circle: (0, 0) to (90, 0) => sep = 90 deg
        const dEquatorQuarter = angularSeparationDeg(0.0, 0.0, 90.0, 0.0);

        // 4. Triangle on sphere: (0, 0) to (0, 60) => sep = 60 deg
        const dMeridian = angularSeparationDeg(0.0, 0.0, 0.0, 60.0);

        // 5. Small separation limit
        const dSmall = angularSeparationDeg(187.7, 12.3, 187.7 + 0.01, 12.3 + 0.01);

        console.log(JSON.stringify({
            dZero,
            dPoles,
            dEquatorQuarter,
            dMeridian,
            dSmall
        }));
        """)
        assert math.isclose(res["dZero"], 0.0, abs_tol=1e-7)
        assert math.isclose(res["dPoles"], 180.0, rel_tol=1e-6)
        assert math.isclose(res["dEquatorQuarter"], 90.0, rel_tol=1e-6)
        assert math.isclose(res["dMeridian"], 60.0, rel_tol=1e-6)
        assert 0.01 < res["dSmall"] < 0.02

    def test_projected_physical_separation(self):
        """Verify projected transverse physical separation R_proj in Mpc at metric distance d."""
        res = run_node_eval("""
        import { projectedPhysicalSeparationMpc, angularSeparationRad } from './src/data/cf4_group_catalog.js';

        // Virgo cluster scale: distance ~ 16.5 Mpc, angular separation 1 degree
        const dGroup = 16.5; // Mpc
        const sep1Deg = projectedPhysicalSeparationMpc(187.0, 12.0, 188.0, 12.0, dGroup);
        const thetaRad = angularSeparationRad(187.0, 12.0, 188.0, 12.0);
        const expectedMpc = 2.0 * dGroup * Math.sin(thetaRad / 2.0);

        // Coma cluster scale: distance ~ 100 Mpc, angular separation 2 degrees
        const sepComa = projectedPhysicalSeparationMpc(194.0, 27.0, 196.0, 27.0, 100.0);

        console.log(JSON.stringify({
            sep1Deg,
            expectedMpc,
            sepComa
        }));
        """)
        assert math.isclose(res["sep1Deg"], res["expectedMpc"], rel_tol=1e-7)
        assert 0.2 < res["sep1Deg"] < 0.4  # ~0.28 Mpc at 16.5 Mpc for ~1 deg
        assert 2.0 < res["sepComa"] < 4.0


class TestMasterCosmologicalStructures:
    """Tier 2: Verification and Cross-Matching of Named Master Cosmological Structures"""

    def test_master_structures_registry(self):
        """Verify all 11 canonical master structures are present with valid coordinates and properties."""
        res = run_node_eval("""
        import { MASTER_STRUCTURES } from './src/data/cf4_group_catalog.js';

        const structures = MASTER_STRUCTURES.map(s => ({
            id: s.id,
            name: s.name,
            sgx: s.sgx,
            sgy: s.sgy,
            sgz: s.sgz,
            cz: s.cz,
            rMatch: s.rMatch,
            vMatch: s.vMatch,
            mVirNominal: s.mVirNominal
        }));

        console.log(JSON.stringify(structures));
        """)
        assert len(res) == 11
        ids = [s["id"] for s in res]
        assert "LANIAKEA_CORE" in ids
        assert "COMA_CLUSTER" in ids
        assert "VIRGO_CLUSTER" in ids
        assert "CENTAURUS_CLUSTER" in ids
        assert "NORMA_CLUSTER" in ids
        assert "SHAPLEY_CORE" in ids
        assert "PERSEUS_PISCES" in ids
        assert "FORNAX_CLUSTER" in ids
        assert "LOCAL_GROUP" in ids

        # Verify Coma properties
        coma = next(s for s in res if s["id"] == "COMA_CLUSTER")
        assert math.isclose(coma["cz"], 6925.0, abs_tol=100.0)
        assert coma["mVirNominal"] >= 1.0e15

        # Verify Virgo properties
        virgo = next(s for s in res if s["id"] == "VIRGO_CLUSTER")
        assert math.isclose(virgo["cz"], 1035.0, abs_tol=100.0)
        assert coma["sgx"] < 0

    def test_cross_match_master_structures_exact_and_perturbed(self):
        """Verify cross-matching algorithm correctly identifies structures when given exact and perturbed coordinates."""
        res = run_node_eval("""
        import { crossMatchMasterStructure, MASTER_STRUCTURES } from './src/data/cf4_group_catalog.js';

        // 1. Exact match with Coma Cluster
        const comaExact = crossMatchMasterStructure(-7.1, 95.8, 6.8, 6925.0, 194.953, 27.981);

        // 2. Perturbed match with Virgo Cluster (offset by 2 Mpc/h and 150 km/s)
        const virgoPerturbed = crossMatchMasterStructure(-3.2 + 2.0, 16.2 - 1.5, -0.7 + 0.5, 1035.0 + 150.0);

        // 3. Exact match with Shapley Supercluster Core
        const shapleyExact = crossMatchMasterStructure(-135.0, 72.0, -38.0, 14500.0);

        // 4. Far outlier in the Great Cosmic Void: no match
        const voidPoint = crossMatchMasterStructure(150.0, 150.0, 150.0, 25000.0);

        console.log(JSON.stringify({
            comaExact,
            virgoPerturbed,
            shapleyExact,
            voidPoint
        }));
        """)
        assert res["comaExact"]["matched"] is True
        assert res["comaExact"]["structure"]["id"] == "COMA_CLUSTER"
        assert res["comaExact"]["matchScore"] > 0.95

        assert res["virgoPerturbed"]["matched"] is True
        assert res["virgoPerturbed"]["structure"]["id"] == "VIRGO_CLUSTER"
        assert res["virgoPerturbed"]["matchScore"] > 0.50

        assert res["shapleyExact"]["matched"] is True
        assert res["shapleyExact"]["structure"]["id"] == "SHAPLEY_CORE"

        assert res["voidPoint"]["matched"] is False
        assert res["voidPoint"]["structure"] is None


class TestRobustStatisticalEstimators:
    """Tier 3: Beers et al. (1990) Biweight Location & Scale, Gapper Scale, and Dispersion Correction"""

    def test_median_and_mad_computation(self):
        """Verify sample median and MAD on known discrete arrays and symmetric distributions."""
        res = run_node_eval("""
        import { calculateMedian, calculateMAD } from './src/data/cf4_group_catalog.js';

        const oddArr = [10.0, 20.0, 30.0, 40.0, 50.0];
        const evenArr = [10.0, 20.0, 30.0, 40.0];
        const symmetricArr = [1000, 1100, 1200, 1300, 1400];

        const medOdd = calculateMedian(oddArr);
        const medEven = calculateMedian(evenArr);
        const madOdd = calculateMAD(oddArr, medOdd);
        const madSymmetric = calculateMAD(symmetricArr);

        console.log(JSON.stringify({
            medOdd,
            medEven,
            madOdd,
            madSymmetric
        }));
        """)
        assert math.isclose(res["medOdd"], 30.0, abs_tol=1e-7)
        assert math.isclose(res["medEven"], 25.0, abs_tol=1e-7)
        assert math.isclose(res["madOdd"], 10.0, abs_tol=1e-7)
        assert math.isclose(res["madSymmetric"], 100.0, abs_tol=1e-7)

    def test_biweight_location_outlier_resistance(self):
        """Verify Beers Biweight Location (C_BI) rejects extreme velocity interlopers compared to standard mean."""
        res = run_node_eval("""
        import { biweightLocation } from './src/data/cf4_group_catalog.js';

        // Group of 10 galaxies centered at 5000 km/s with sigma ~ 300 km/s
        const cleanVelocities = [4600, 4750, 4900, 4950, 5000, 5050, 5100, 5200, 5350, 5100];
        // Add 1 extreme background interloper at 25,000 km/s
        const contaminatedVelocities = [...cleanVelocities, 25000];

        const meanClean = cleanVelocities.reduce((a, b) => a + b, 0) / cleanVelocities.length;
        const meanContaminated = contaminatedVelocities.reduce((a, b) => a + b, 0) / contaminatedVelocities.length;

        const biweightClean = biweightLocation(cleanVelocities);
        const biweightContaminated = biweightLocation(contaminatedVelocities);

        console.log(JSON.stringify({
            meanClean,
            meanContaminated,
            biweightClean,
            biweightContaminated
        }));
        """)
        assert math.isclose(res["meanClean"], 5000.0, abs_tol=50.0)
        assert res["meanContaminated"] > 6800.0  # Sample mean pulled up significantly by interloper
        # Biweight location should reject the 25,000 km/s interloper and stay near 5000 km/s
        assert math.isclose(res["biweightContaminated"], 5000.0, abs_tol=100.0)

    def test_biweight_and_gapper_scale_estimators(self):
        """Verify Beers Biweight Scale (S_BI) and Gapper Estimator on Gaussian samples and small N."""
        res = run_node_eval("""
        import { biweightScale, gapperScale, sampleStandardDeviation } from './src/data/cf4_group_catalog.js';

        // Sample of 12 velocities drawn from Normal(3000, 500)
        const velocities = [
            2350, 2540, 2710, 2850, 2980, 3020,
            3080, 3190, 3320, 3480, 3650, 3830
        ];

        const stdDev = sampleStandardDeviation(velocities);
        const biweight = biweightScale(velocities);
        const gapper = gapperScale(velocities);

        // Binary pair: N=2
        const pair = [3000, 3400];
        const pairBiweight = biweightScale(pair);
        const pairGapper = gapperScale(pair);
        const pairExpected = Math.abs(3400 - 3000) / Math.SQRT2;

        console.log(JSON.stringify({
            stdDev,
            biweight,
            gapper,
            pairBiweight,
            pairGapper,
            pairExpected
        }));
        """)
        assert 400.0 < res["stdDev"] < 600.0
        assert 400.0 < res["biweight"] < 600.0
        assert 400.0 < res["gapper"] < 600.0
        assert math.isclose(res["pairBiweight"], res["pairExpected"], rel_tol=1e-5)
        assert math.isclose(res["pairGapper"], res["pairExpected"], rel_tol=1e-5)

    def test_deconvolve_velocity_dispersion(self):
        """Verify measurement error deconvolution: sigma_{v,corr}^2 = max(0, sigma_v^2 - <sigma_err^2>) / (1+z)^2."""
        res = run_node_eval("""
        import { deconvolveVelocityDispersion } from './src/data/cf4_group_catalog.js';

        const sigmaV = 500.0; // km/s
        const errors = [50.0, 60.0, 40.0, 70.0, 50.0]; // mean error ~ 54 km/s
        const z = 0.05; // redshift

        const sigmaCorr = deconvolveVelocityDispersion(sigmaV, errors, z);

        const meanErrSq = (50*50 + 60*60 + 40*40 + 70*70 + 50*50) / 5;
        const expected = Math.sqrt(Math.max(0, sigmaV * sigmaV - meanErrSq)) / (1.0 + z);

        // Case where error exceeds dispersion
        const smallSigma = 30.0;
        const largeErrors = [50.0, 60.0, 70.0];
        const sigmaClamped = deconvolveVelocityDispersion(smallSigma, largeErrors, z);

        console.log(JSON.stringify({
            sigmaCorr,
            expected,
            sigmaClamped
        }));
        """)
        assert math.isclose(res["sigmaCorr"], res["expected"], rel_tol=1e-6)
        assert res["sigmaClamped"] == 0.0


class TestProjectedHarmonicRadiusAndRadiiMetrics:
    """Tier 4: Pairwise Projected Harmonic Radius R_H and Projected Radial Moments"""

    def test_harmonic_radius_equilateral_triangle(self):
        """Verify projected harmonic radius on equilateral triangle with exact analytic geometry."""
        res = run_node_eval("""
        import { calculateProjectedHarmonicRadius } from './src/data/cf4_group_catalog.js';

        // 3 galaxies forming equilateral triangle with side L = 1.0 Mpc at distance d = 50 Mpc
        // L = 2 * d * sin(theta / 2) => theta = 2 * asin(L / (2*d))
        const dGroup = 50.0;
        const L = 1.0; // Mpc
        const thetaRad = 2.0 * Math.asin(L / (2.0 * dGroup));
        const thetaDeg = thetaRad * (180.0 / Math.PI);

        const members = [
            { ra: 0.0, dec: thetaDeg },
            { ra: -thetaDeg * Math.cos(Math.PI / 6), dec: -thetaDeg * Math.sin(Math.PI / 6) },
            { ra: thetaDeg * Math.cos(Math.PI / 6), dec: -thetaDeg * Math.sin(Math.PI / 6) }
        ];

        const rH = calculateProjectedHarmonicRadius(members, dGroup, 0.0);
        // Exact R_H for 3 points with circumradius R_c = L has side s = sqrt(3)*L: R_H = (pi / 2) * sqrt(3) * L
        const expectedRH = (Math.PI / 2.0) * Math.sqrt(3.0) * L;

        console.log(JSON.stringify({
            rH,
            expectedRH
        }));
        """)
        # Due to slight spherical distortion on small angles, should match within 0.1%
        assert math.isclose(res["rH"], res["expectedRH"], rel_tol=1e-3)

    def test_projected_radii_moments_and_concentration(self):
        """Verify R_mean, R_rms, R_med, and concentration index C_proj = R_80 / R_20."""
        res = run_node_eval("""
        import { calculateProjectedRadiiMoments } from './src/data/cf4_group_catalog.js';

        const dGroup = 20.0; // Mpc
        const members = [
            { ra: 10.0, dec: 0.0 },
            { ra: 10.1, dec: 0.0 },
            { ra: 10.2, dec: 0.0 },
            { ra: 10.5, dec: 0.0 },
            { ra: 11.0, dec: 0.0 }
        ];

        const moments = calculateProjectedRadiiMoments(members, 10.0, 0.0, dGroup);

        console.log(JSON.stringify({
            rMean: moments.rMean,
            rRms: moments.rRms,
            rMed: moments.rMed,
            rMax: moments.rMax,
            concentration: moments.concentration,
            nRadii: moments.projectedRadii.length
        }));
        """)
        assert res["nRadii"] == 5
        assert res["rMean"] > 0
        assert res["rRms"] >= res["rMean"]
        assert res["rMax"] >= res["rMed"]
        assert res["concentration"] > 1.0

    def test_virial_radius_r200_and_r500(self):
        """Verify R_200 = sqrt(3) * sigma_v / (10 * H(z)) and R_500 ~ 0.66 * R_200."""
        res = run_node_eval("""
        import { calculateR200Mpc, calculateR500Mpc, DEFAULT_H0 } from './src/data/cf4_group_catalog.js';

        const sigmaV = 1000.0; // km/s (rich cluster like Coma)
        const z = 0.023; // Coma redshift
        const r200 = calculateR200Mpc(sigmaV, z, DEFAULT_H0);
        const r500 = calculateR500Mpc(r200);

        console.log(JSON.stringify({
            r200,
            r500
        }));
        """)
        # For sigma_v = 1000 km/s, H0 = 74.6: R_200 ~ 2.3 Mpc
        assert 2.0 < res["r200"] < 2.6
        assert math.isclose(res["r500"], res["r200"] * 0.66, rel_tol=1e-5)


class TestVirialMassAndMultiEstimators:
    """Tier 5: Virial Mass, Heisler et al. Estimators, and Bootstrap Uncertainty"""

    def test_standard_virial_mass_synthetic_cluster(self):
        """Verify standard virial mass M_vir = (3 * pi / 2) * (R_H * sigma_v^2 / G)."""
        res = run_node_eval("""
        import { calculateStandardVirialMass, G_COSMO } from './src/data/cf4_group_catalog.js';

        const rH = 1.5; // Mpc
        const sigmaV = 800.0; // km/s
        const mVir = calculateStandardVirialMass(rH, sigmaV);

        const expectedMVir = (3.0 * Math.PI / 2.0) * (rH * sigmaV * sigmaV) / G_COSMO;

        console.log(JSON.stringify({
            mVir,
            expectedMVir,
            mVirSolarOrder: Math.log10(mVir)
        }));
        """)
        assert math.isclose(res["mVir"], res["expectedMVir"], rel_tol=1e-7)
        # Coma-like rich cluster mass should be ~10^15 M_sun
        assert 14.8 < res["mVirSolarOrder"] < 15.3

    def test_heisler_tremaine_bahcall_multi_estimators(self):
        """Verify Heisler et al. (1985) Projected, Median, and Average mass estimators consistency."""
        res = run_node_eval("""
        import {
            calculateProjectedMassEstimator,
            calculateMedianMassEstimator,
            calculateAverageMassEstimator,
            calculateStandardVirialMass,
            calculateProjectedHarmonicRadius,
            biweightScale
        } from './src/data/cf4_group_catalog.js';

        // Synthetic 6-galaxy virialized group
        const dGroup = 30.0; // Mpc
        const members = [
            { ra: 150.0, dec: 20.0, cz: 3800.0, vRec: 3800.0 },
            { ra: 150.2, dec: 20.1, cz: 4200.0, vRec: 4200.0 },
            { ra: 149.9, dec: 20.3, cz: 3950.0, vRec: 3950.0 },
            { ra: 150.3, dec: 19.8, cz: 4100.0, vRec: 4100.0 },
            { ra: 149.8, dec: 19.9, cz: 3750.0, vRec: 3750.0 },
            { ra: 150.1, dec: 20.2, cz: 4300.0, vRec: 4300.0 }
        ];

        const vMean = 4016.67;
        const raCentroid = 150.05;
        const decCentroid = 20.05;

        const mProj = calculateProjectedMassEstimator(members, raCentroid, decCentroid, vMean, dGroup, 'isotropic');
        const mMed = calculateMedianMassEstimator(members, dGroup);
        const mAvg = calculateAverageMassEstimator(members, dGroup);

        const rH = calculateProjectedHarmonicRadius(members, dGroup);
        const sigV = biweightScale(members.map(m => m.vRec));
        const mVir = calculateStandardVirialMass(rH, sigV);

        console.log(JSON.stringify({
            mVir,
            mProj,
            mMed,
            mAvg
        }));
        """)
        m_vir = res["mVir"]
        m_proj = res["mProj"]
        m_med = res["mMed"]
        m_avg = res["mAvg"]

        assert m_vir > 0
        assert m_proj > 0
        assert m_med > 0
        assert m_avg > 0
        # All 4 estimators on a virialized system should agree within a factor of ~2-3
        assert 0.3 < (m_proj / m_vir) < 3.0
        assert 0.3 < (m_med / m_vir) < 3.0
        assert 0.3 < (m_avg / m_vir) < 3.0

    def test_bahcall_tremaine_central_dominant_estimator(self):
        """Verify Bahcall & Tremaine (1981) central-mass estimator for satellite systems."""
        res = run_node_eval("""
        import { calculateBahcallTremaineMassEstimator } from './src/data/cf4_group_catalog.js';

        const dGroup = 16.5; // Mpc
        const satellites = [
            { ra: 187.5, dec: 12.5, cz: 1150.0, vRec: 1150.0 },
            { ra: 187.9, dec: 12.2, cz: 950.0, vRec: 950.0 },
            { ra: 188.0, dec: 12.6, cz: 1250.0, vRec: 1250.0 }
        ];

        const mCentral = calculateBahcallTremaineMassEstimator(satellites, 187.7, 12.4, 1050.0, dGroup);

        console.log(JSON.stringify({ mCentral }));
        """)
        assert res["mCentral"] > 0
        assert 1.0e12 < res["mCentral"] < 1.0e15

    def test_surface_pressure_correction(self):
        """Verify The & White (1986) surface correction factor reduces mass by ~18%."""
        res = run_node_eval("""
        import { applySurfacePressureCorrection } from './src/data/cf4_group_catalog.js';

        const mVirRaw = 1.0e15;
        const mVirCorr = applySurfacePressureCorrection(mVirRaw, 0.82);

        console.log(JSON.stringify({ mVirRaw, mVirCorr }));
        """)
        assert math.isclose(res["mVirCorr"], 0.82e15, rel_tol=1e-7)

    def test_bootstrap_virial_mass_errors(self):
        """Verify non-parametric bootstrap resampling generates valid confidence intervals."""
        res = run_node_eval("""
        import { bootstrapVirialMassErrors } from './src/data/cf4_group_catalog.js';

        const dGroup = 25.0; // Mpc
        const members = [
            { ra: 120.0, dec: 30.0, cz: 2500.0, vRec: 2500.0 },
            { ra: 120.2, dec: 30.1, cz: 2700.0, vRec: 2700.0 },
            { ra: 119.8, dec: 29.9, cz: 2400.0, vRec: 2400.0 },
            { ra: 120.1, dec: 30.3, cz: 2650.0, vRec: 2650.0 },
            { ra: 119.9, dec: 30.2, cz: 2550.0, vRec: 2550.0 }
        ];

        const boot = bootstrapVirialMassErrors(members, dGroup, 150, 12345);

        console.log(JSON.stringify(boot));
        """)
        assert res["mVirMedian"] > 0
        assert res["mVirLow"] <= res["mVirMedian"] <= res["mVirHigh"]
        assert res["sigmaVLow"] <= res["sigmaVMedian"] <= res["sigmaVHigh"]


class TestLuminosityCompletenessAndMorphology:
    """Tier 6: Incomplete Gamma, Schechter Completeness, and Morphological Segregation"""

    def test_lanczos_log_gamma_and_standard_gamma(self):
        """Verify Gamma(x) for integer factorials: Gamma(n) = (n-1)!."""
        res = run_node_eval("""
        import { standardGamma, logGamma } from './src/data/cf4_group_catalog.js';

        const g1 = standardGamma(1.0); // 0! = 1
        const g2 = standardGamma(2.0); // 1! = 1
        const g3 = standardGamma(3.0); // 2! = 2
        const g4 = standardGamma(4.0); // 3! = 6
        const g5 = standardGamma(5.0); // 4! = 24
        const gHalf = standardGamma(0.5); // sqrt(pi) ~ 1.77245385

        console.log(JSON.stringify({
            g1, g2, g3, g4, g5, gHalf,
            expectedGHalf: Math.sqrt(Math.PI)
        }));
        """)
        assert math.isclose(res["g1"], 1.0, rel_tol=1e-6)
        assert math.isclose(res["g2"], 1.0, rel_tol=1e-6)
        assert math.isclose(res["g3"], 2.0, rel_tol=1e-6)
        assert math.isclose(res["g4"], 6.0, rel_tol=1e-6)
        assert math.isclose(res["g5"], 24.0, rel_tol=1e-6)
        assert math.isclose(res["gHalf"], res["expectedGHalf"], rel_tol=1e-6)

    def test_upper_incomplete_gamma(self):
        """Verify upper incomplete gamma function Gamma(a, x) against asymptotic and special cases."""
        res = run_node_eval("""
        import { upperIncompleteGamma, standardGamma } from './src/data/cf4_group_catalog.js';

        // 1. Gamma(a, 0) = Gamma(a)
        const g0 = upperIncompleteGamma(2.5, 0.0);
        const gExpected = standardGamma(2.5);

        // 2. Gamma(1, x) = exp(-x)
        const xVal = 2.0;
        const g1x = upperIncompleteGamma(1.0, xVal);
        const g1Expected = Math.exp(-xVal);

        console.log(JSON.stringify({
            g0,
            gExpected,
            g1x,
            g1Expected
        }));
        """)
        assert math.isclose(res["g0"], res["gExpected"], rel_tol=1e-6)
        assert math.isclose(res["g1x"], res["g1Expected"], rel_tol=1e-6)

    def test_schechter_completeness_correction(self):
        """Verify Schechter completeness fraction C(L_lim) decreases with distance."""
        res = run_node_eval("""
        import { schechterCompletenessCorrection } from './src/data/cf4_group_catalog.js';

        const appMagLim = 11.75; // 2MASS K-band limit
        const c10Mpc = schechterCompletenessCorrection(appMagLim, 10.0);
        const c50Mpc = schechterCompletenessCorrection(appMagLim, 50.0);
        const c200Mpc = schechterCompletenessCorrection(appMagLim, 200.0);

        console.log(JSON.stringify({
            c10Mpc,
            c50Mpc,
            c200Mpc
        }));
        """)
        assert 0.0 < res["c200Mpc"] <= res["c50Mpc"] <= res["c10Mpc"] <= 1.0

    def test_group_photometry_and_morphology_analysis(self):
        """Verify total luminosity integration, morphological fractions, and M/L ratios."""
        res = run_node_eval("""
        import { analyzeGroupLuminosityAndMorphology } from './src/data/cf4_group_catalog.js';

        const dGroup = 20.0; // Mpc
        const mVir = 2.0e14; // M_sun
        const members = [
            { kMag: 8.5, bMag: 11.0, tType: -2.0 }, // Early-type (E)
            { kMag: 9.0, bMag: 11.5, tType: -1.0 }, // Early-type (S0)
            { kMag: 9.5, bMag: 12.0, tType: 3.0 },  // Late-type (Sb)
            { kMag: 10.0, bMag: 12.5, tType: 5.0 }, // Late-type (Sc)
            { kMag: 11.0, bMag: 13.5, tType: 10.0 } // Irregular (Irr)
        ];

        const analysis = analyzeGroupLuminosityAndMorphology(members, dGroup, mVir);

        console.log(JSON.stringify(analysis));
        """)
        assert res["lObsK"] > 0
        assert res["lTotK"] >= res["lObsK"]
        assert math.isclose(res["fEarly"], 0.4, abs_tol=1e-5)  # 2 out of 5
        assert math.isclose(res["fLate"], 0.4, abs_tol=1e-5)   # 2 out of 5
        assert math.isclose(res["fIrr"], 0.2, abs_tol=1e-5)    # 1 out of 5
        assert res["mToLK"] > 0
        assert res["mToLB"] > 0

    def test_dressler_surface_density(self):
        """Verify Dressler (1980) local projected surface density Sigma_10 computation."""
        res = run_node_eval("""
        import { calculateLocalSurfaceDensity } from './src/data/cf4_group_catalog.js';

        const dGroup = 15.0; // Mpc
        const galaxies = [];
        // Create 15 galaxies around (180, 10)
        for (let i = 0; i < 15; i++) {
            galaxies.push({
                ra: 180.0 + (i % 4) * 0.05,
                dec: 10.0 + Math.floor(i / 4) * 0.05
            });
        }

        const sigma10 = calculateLocalSurfaceDensity(180.0, 10.0, galaxies, dGroup, 5);

        console.log(JSON.stringify({ sigma10 }));
        """)
        assert res["sigma10"] > 0


class TestHierarchicalFriendsOfFriendsAndMembership:
    """Tier 7: Hierarchical FoF Group Finding and Projected Phase-Space Membership"""

    def test_disjoint_set_union_find(self):
        """Verify DisjointSet path compression and rank union."""
        res = run_node_eval("""
        import { DisjointSet } from './src/data/cf4_group_catalog.js';

        const ds = new DisjointSet(10);
        ds.union(0, 1);
        ds.union(1, 2);
        ds.union(3, 4);

        const root0 = ds.find(0);
        const root2 = ds.find(2);
        const root3 = ds.find(3);
        const root4 = ds.find(4);
        const root8 = ds.find(8);

        console.log(JSON.stringify({
            same0_2: root0 === root2,
            same3_4: root3 === root4,
            diff0_3: root0 !== root3,
            root8Self: root8 === 8
        }));
        """)
        assert res["same0_2"] is True
        assert res["same3_4"] is True
        assert res["diff0_3"] is True
        assert res["root8Self"] is True

    def test_adaptive_fof_linking_lengths(self):
        """Verify Tully (2015) adaptive linking lengths increase with recession velocity."""
        res = run_node_eval("""
        import { getAdaptiveFoFLinkingLengths } from './src/data/cf4_group_catalog.js';

        const linkZ0 = getAdaptiveFoFLinkingLengths(0.0);
        const linkZ3000 = getAdaptiveFoFLinkingLengths(3000.0);
        const linkZ10000 = getAdaptiveFoFLinkingLengths(10000.0);

        console.log(JSON.stringify({
            linkZ0,
            linkZ3000,
            linkZ10000
        }));
        """)
        assert res["linkZ0"]["dLinkMpc"] < res["linkZ3000"]["dLinkMpc"] < res["linkZ10000"]["dLinkMpc"]
        assert res["linkZ0"]["vLinkKmS"] < res["linkZ3000"]["vLinkKmS"] < res["linkZ10000"]["vLinkKmS"]

    def test_membership_probability_radial_and_velocity_falloff(self):
        """Verify membership probability P_mem decreases with projected radius and velocity offset."""
        res = run_node_eval("""
        import { calculateMembershipProbability } from './src/data/cf4_group_catalog.js';

        const r200 = 1.5; // Mpc
        const sigmaV = 500.0; // km/s

        // Center galaxy
        const pCenter = calculateMembershipProbability(0.0, 0.0, r200, sigmaV);

        // Core member at 0.5 * R_200 and 0.5 * sigma_v
        const pCore = calculateMembershipProbability(0.75, 250.0, r200, sigmaV);

        // Virial boundary member at R_200 and sigma_v
        const pVirial = calculateMembershipProbability(1.5, 500.0, r200, sigmaV);

        // Outlier at 3 * R_200 and 3 * sigma_v
        const pOutlier = calculateMembershipProbability(4.5, 1500.0, r200, sigmaV);

        console.log(JSON.stringify({
            pCenter,
            pCore,
            pVirial,
            pOutlier
        }));
        """)
        assert math.isclose(res["pCenter"], 1.0, rel_tol=1e-5)
        assert res["pCore"] > res["pVirial"] > res["pOutlier"]
        assert res["pOutlier"] < 0.01

    def test_run_hierarchical_fof_synthetic_cluster_field(self):
        """Verify FoF algorithm separates two distinct clusters in projected phase space."""
        res = run_node_eval("""
        import { runHierarchicalFoF } from './src/data/cf4_group_catalog.js';

        // Cluster 1 at (10.0, 20.0) with cz ~ 3000 km/s
        const cluster1 = [
            { pgc: 1, ra: 10.0, dec: 20.0, cz: 3000.0 },
            { pgc: 2, ra: 10.1, dec: 20.05, cz: 3050.0 },
            { pgc: 3, ra: 9.95, dec: 19.95, cz: 2950.0 }
        ];

        // Cluster 2 at (50.0, -30.0) with cz ~ 8000 km/s
        const cluster2 = [
            { pgc: 4, ra: 50.0, dec: -30.0, cz: 8000.0 },
            { pgc: 5, ra: 50.05, dec: -29.95, cz: 8050.0 },
            { pgc: 6, ra: 49.95, dec: -30.05, cz: 7950.0 }
        ];

        const allGalaxies = [...cluster1, ...cluster2];
        const groupsMap = runHierarchicalFoF(allGalaxies, { d0: 0.5, v0: 400.0, minMembers: 2 });

        const groupsArray = [];
        for (const [root, members] of groupsMap.entries()) {
            groupsArray.push(members.map(m => m.pgc));
        }

        console.log(JSON.stringify(groupsArray));
        """)
        assert len(res) == 2
        flat_pgcs = [sorted(g) for g in res]
        assert [1, 2, 3] in flat_pgcs
        assert [4, 5, 6] in flat_pgcs


class TestCF4GroupedCatalogColumnarStoreAndParsers:
    """Tier 8: Columnar Store, Aggregation Engine, Spatial Index, and File Ingestion"""

    def test_columnar_store_ingestion_and_spatial_query(self):
        """Verify TypedArray columnar store, spatial hash grid indexing, and 3D cone queries."""
        res = run_node_eval("""
        import { CF4GroupedCatalog } from './src/data/cf4_group_catalog.js';

        const catalog = new CF4GroupedCatalog(100);

        // Add 10 galaxies
        for (let i = 0; i < 10; i++) {
            catalog.addGalaxy({
                pgc: 1000 + i,
                groupId: i < 5 ? 101 : 102,
                ra: 180.0 + i * 0.1,
                dec: 10.0 + i * 0.1,
                cz: 3000.0 + i * 20.0,
                kMag: 10.0 + i * 0.2
            });
        }

        catalog.aggregateGroups();

        const g1005 = catalog.getGalaxyByPGC(1005);
        const group101 = catalog.groups.get(101);
        const group102 = catalog.groups.get(102);

        // 3D query near first group
        const coneHits = catalog.queryCone3D(catalog.sgx[0], catalog.sgy[0], catalog.sgz[0], 5.0);

        console.log(JSON.stringify({
            catalogLength: catalog.length,
            g1005_pgc: g1005 ? g1005.pgc : null,
            group101_nMembers: group101 ? group101.nMembers : null,
            group102_nMembers: group102 ? group102.nMembers : null,
            coneHitsCount: coneHits.length
        }));
        """)
        assert res["catalogLength"] == 10
        assert res["g1005_pgc"] == 1005
        assert res["group101_nMembers"] == 5
        assert res["group102_nMembers"] == 5
        assert res["coneHitsCount"] >= 1

    def test_parse_cf4_group_catalog_csv(self):
        """Verify CSV text parsing and automatic group aggregation."""
        res = run_node_eval("""
        import { parseCF4GroupCatalogText } from './src/data/cf4_group_catalog.js';

        const csvText = [
            'PGC,GroupID,RA,Dec,cz,vErr,Dist,BMag,KMag,TType',
            '101,501,187.70,12.39,1035.0,25.0,16.5,10.2,7.5,-2.0',
            '102,501,187.75,12.45,1120.0,30.0,16.8,11.5,8.8,3.0',
            '103,501,187.65,12.30,980.0,20.0,16.2,12.0,9.2,5.0',
            '201,602,194.95,27.98,6925.0,40.0,95.0,11.0,8.0,-1.0',
            '202,602,195.05,28.05,7150.0,35.0,98.0,12.5,9.5,-2.0'
        ].join('\\n');

        const catalog = parseCF4GroupCatalogText(csvText);
        const g501 = catalog.groups.get(501);
        const g602 = catalog.groups.get(602);

        console.log(JSON.stringify({
            totalGalaxies: catalog.length,
            nGroups: catalog.groups.size,
            g501_n: g501 ? g501.nMembers : 0,
            g501_mVir: g501 ? g501.mVir : 0,
            g501_matchedStructure: g501 && g501.matchedStructure ? g501.matchedStructure.name : null,
            g602_n: g602 ? g602.nMembers : 0,
            g602_matchedStructure: g602 && g602.matchedStructure ? g602.matchedStructure.name : null
        }));
        """)
        assert res["totalGalaxies"] == 5
        assert res["nGroups"] == 2
        assert res["g501_n"] == 3
        assert res["g501_mVir"] > 0
        # Group 501 is at Virgo coordinates
        assert "Virgo" in str(res["g501_matchedStructure"])
        # Group 602 is at Coma coordinates
        assert "Coma" in str(res["g602_matchedStructure"])

    def test_export_to_geojson_and_csv(self):
        """Verify catalog exports valid GeoJSON FeatureCollection and CSV formatting."""
        res = run_node_eval("""
        import { CF4GroupedCatalog } from './src/data/cf4_group_catalog.js';

        const catalog = new CF4GroupedCatalog(10);
        catalog.addGalaxy({ pgc: 1, groupId: 10, ra: 180.0, dec: 10.0, cz: 2000.0 });
        catalog.addGalaxy({ pgc: 2, groupId: 10, ra: 180.1, dec: 10.1, cz: 2100.0 });
        catalog.aggregateGroups();

        const geojson = catalog.exportGroupsToGeoJSON();
        const csv = catalog.exportToCSV();

        console.log(JSON.stringify({
            geoJsonType: geojson.type,
            nFeatures: geojson.features.length,
            firstFeatureProps: geojson.features[0].properties.groupId,
            csvHeaderPresent: csv.startsWith('PGC,GroupID,NEST,RA,Dec')
        }));
        """)
        assert res["geoJsonType"] == "FeatureCollection"
        assert res["nFeatures"] == 1
        assert res["firstFeatureProps"] == 10
        assert res["csvHeaderPresent"] is True


class TestScientificInvariantsAndEdgeCases:
    """Tier 9: Scientific Invariants, Singletons, Binary Pairs, Large Clusters, and Error Handling"""

    def test_single_galaxy_group_n1(self):
        """Verify N=1 single galaxy group handles metrics safely without crashing or NaN."""
        res = run_node_eval("""
        import { CF4Group } from './src/data/cf4_group_catalog.js';

        const singleMember = [{ pgc: 777, ra: 120.0, dec: 45.0, cz: 3500.0, dist: 46.9 }];
        const group = new CF4Group(777, singleMember);

        console.log(JSON.stringify({
            nMembers: group.nMembers,
            meanVRec: group.meanVRec,
            sigmaVCorr: group.sigmaVCorr,
            mVir: group.mVir,
            rHarmonicMpc: group.rHarmonicMpc,
            isNanMeanV: isNaN(group.meanVRec),
            isNanMVir: isNaN(group.mVir)
        }));
        """)
        assert res["nMembers"] == 1
        assert res["meanVRec"] == 3500.0
        assert res["sigmaVCorr"] == 0.0
        assert res["mVir"] == 0.0
        assert res["isNanMeanV"] is False
        assert res["isNanMVir"] is False

    def test_binary_galaxy_pair_n2(self):
        """Verify N=2 galaxy pair calculates exact pairwise velocity dispersion and harmonic radius."""
        res = run_node_eval("""
        import { CF4Group } from './src/data/cf4_group_catalog.js';

        const pairMembers = [
            { pgc: 1, ra: 150.0, dec: 20.0, cz: 4000.0, dist: 50.0 },
            { pgc: 2, ra: 150.1, dec: 20.0, cz: 4300.0, dist: 50.0 }
        ];

        const group = new CF4Group(101, pairMembers, { computeBootstrap: false });

        // Expected pairwise velocity dispersion: |4300 - 4000| / sqrt(2) ~ 212.13 km/s
        const expectedSigma = 300.0 / Math.SQRT2;

        console.log(JSON.stringify({
            nMembers: group.nMembers,
            meanVRec: group.meanVRec,
            sigmaVRaw: group.sigmaVRaw,
            sigmaVGapper: group.sigmaVGapper,
            expectedSigma,
            rHarmonicMpc: group.rHarmonicMpc,
            mVir: group.mVir
        }));
        """)
        assert res["nMembers"] == 2
        assert res["meanVRec"] == 4150.0
        assert math.isclose(res["sigmaVGapper"], res["expectedSigma"], rel_tol=1e-5)
        assert res["rHarmonicMpc"] > 0
        assert res["mVir"] > 0

    def test_large_galaxy_cluster_n500(self):
        """Verify numerical stability and performance on large rich cluster (N = 500 members)."""
        res = run_node_eval("""
        import { CF4Group } from './src/data/cf4_group_catalog.js';

        const members = [];
        const n = 500;
        const dGroup = 100.0; // Mpc

        // Deterministic pseudo-random distribution
        for (let i = 0; i < n; i++) {
            const angle = (i * 137.5) * (Math.PI / 180.0);
            const radiusDeg = Math.sqrt(i / n) * 1.5; // up to 1.5 deg radius
            const ra = 195.0 + radiusDeg * Math.cos(angle);
            const dec = 28.0 + radiusDeg * Math.sin(angle);
            // Gaussian-like velocity distribution via Box-Muller approx
            const u1 = ((i * 17 + 5) % 1000) / 1000.0 + 1.0e-5;
            const u2 = ((i * 31 + 13) % 1000) / 1000.0;
            const gVel = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
            const cz = 7000.0 + gVel * 850.0;

            members.push({ pgc: i + 1, ra, dec, cz, dist: dGroup });
        }

        const group = new CF4Group(999, members, { computeBootstrap: false });

        console.log(JSON.stringify({
            nMembers: group.nMembers,
            meanVRec: group.meanVRec,
            biweightVRec: group.biweightVRec,
            sigmaVBiweight: group.sigmaVBiweight,
            rHarmonicMpc: group.rHarmonicMpc,
            r200Mpc: group.r200Mpc,
            mVir: group.mVir,
            mVirLog10: Math.log10(group.mVir)
        }));
        """)
        assert res["nMembers"] == 500
        assert math.isclose(res["biweightVRec"], 7000.0, abs_tol=100.0)
        assert 700.0 < res["sigmaVBiweight"] < 1000.0
        assert 1.0 < res["rHarmonicMpc"] < 3.0
        assert 1.5 < res["r200Mpc"] < 3.0
        # Cluster mass ~ 10^15 M_sun
        assert 14.8 < res["mVirLog10"] < 15.5

    def test_zero_velocity_dispersion_and_softening_behavior(self):
        """Verify identical velocities (sigma = 0) and identical coordinates handle softening cleanly."""
        res = run_node_eval("""
        import { CF4Group } from './src/data/cf4_group_catalog.js';

        // 3 galaxies at identical coordinates and identical velocity
        const identicalMembers = [
            { pgc: 1, ra: 100.0, dec: 20.0, cz: 5000.0, dist: 50.0 },
            { pgc: 2, ra: 100.0, dec: 20.0, cz: 5000.0, dist: 50.0 },
            { pgc: 3, ra: 100.0, dec: 20.0, cz: 5000.0, dist: 50.0 }
        ];

        const group = new CF4Group(888, identicalMembers, { computeBootstrap: false });

        console.log(JSON.stringify({
            sigmaVRaw: group.sigmaVRaw,
            sigmaVCorr: group.sigmaVCorr,
            rHarmonicMpc: group.rHarmonicMpc,
            mVir: group.mVir,
            isFiniteMVir: Number.isFinite(group.mVir)
        }));
        """)
        assert res["sigmaVRaw"] == 0.0
        assert res["sigmaVCorr"] == 0.0
        assert res["mVir"] == 0.0
        assert res["isFiniteMVir"] is True
        assert res["rHarmonicMpc"] > 0  # Softened separation


class TestExtendedAdversarialAndExhaustiveStress:
    """Tier 10: Extended Multi-Band Completeness, Orbit Variations, All 11 Master Structures, and Fuzzing"""

    def test_all_eleven_master_structures_individual_matching(self):
        """Verify each of the 11 named master structures is uniquely and correctly matched."""
        res = run_node_eval("""
        import { crossMatchMasterStructure, MASTER_STRUCTURES } from './src/data/cf4_group_catalog.js';

        const matchResults = MASTER_STRUCTURES.map(s => {
            const match = crossMatchMasterStructure(s.sgx, s.sgy, s.sgz, s.cz, s.ra, s.dec);
            return {
                id: s.id,
                matched: match.matched,
                matchedId: match.structure ? match.structure.id : null,
                score: match.matchScore
            };
        });

        console.log(JSON.stringify(matchResults));
        """)
        assert len(res) == 11
        for item in res:
            assert item["matched"] is True, f"Failed to match structure {item['id']}"
            assert item["matchedId"] == item["id"], f"Structure mismatch: expected {item['id']} got {item['matchedId']}"
            assert item["score"] > 0.90

    def test_orbit_anisotropy_variations(self):
        """Verify isotropic vs radial vs circular orbit mass estimator scalings: M_proj(radial) = 2 * M_proj(isotropic) = 4 * M_proj(circular)."""
        res = run_node_eval("""
        import { calculateProjectedMassEstimator } from './src/data/cf4_group_catalog.js';

        const dGroup = 30.0;
        const members = [
            { ra: 100.0, dec: 10.0, cz: 2800.0 },
            { ra: 100.2, dec: 10.1, cz: 3200.0 },
            { ra: 99.8, dec: 9.9, cz: 3000.0 }
        ];

        const mIso = calculateProjectedMassEstimator(members, 100.0, 10.0, 3000.0, dGroup, 'isotropic');
        const mRad = calculateProjectedMassEstimator(members, 100.0, 10.0, 3000.0, dGroup, 'radial');
        const mCirc = calculateProjectedMassEstimator(members, 100.0, 10.0, 3000.0, dGroup, 'circular');

        console.log(JSON.stringify({
            mIso,
            mRad,
            mCirc,
            ratioRadIso: mRad / mIso,
            ratioCircIso: mCirc / mIso
        }));
        """)
        assert math.isclose(res["ratioRadIso"], 2.0, rel_tol=1e-5)
        assert math.isclose(res["ratioCircIso"], 0.5, rel_tol=1e-5)

    def test_schechter_integration_parameter_sweep(self):
        """Verify Schechter completeness integral monotonicity across a grid of alpha and M* parameters."""
        res = run_node_eval("""
        import { schechterCompletenessCorrection } from './src/data/cf4_group_catalog.js';

        const alphas = [-1.5, -1.25, -1.05, -0.8];
        const mStars = [-24.0, -23.5, -22.5, -21.0];
        const distances = [10.0, 50.0, 100.0, 200.0];

        const results = [];
        for (const alpha of alphas) {
            for (const mStar of mStars) {
                for (const d of distances) {
                    const c = schechterCompletenessCorrection(11.75, d, mStar, alpha);
                    results.push({ alpha, mStar, d, c, isValid: c > 0 && c <= 1.0 });
                }
            }
        }

        console.log(JSON.stringify(results));
        """)
        assert len(res) == 64
        for item in res:
            assert item["isValid"] is True
            assert 0.0 < item["c"] <= 1.0

    def test_tsv_and_ip2i_ascii_table_parsing(self):
        """Verify tab-delimited and IP2I ASCII tables with header comments are cleanly parsed."""
        res = run_node_eval("""
        import { parseCF4GroupCatalogText } from './src/data/cf4_group_catalog.js';

        const tsvText = [
            '# Cosmicflows-4 Grouped Catalog IP2I release',
            '# Column definitions: PGC, GroupID, RA, Dec, cz, Dist',
            'PGC\\tGroupID\\tRA\\tDec\\tcz\\tDist',
            '1001\\t9901\\t150.25\\t35.10\\t4500.0\\t60.0',
            '1002\\t9901\\t150.30\\t35.15\\t4620.0\\t61.2',
            '1003\\t9901\\t150.20\\t35.05\\t4480.0\\t59.5'
        ].join('\\n');

        const cat = parseCF4GroupCatalogText(tsvText);
        const group = cat.groups.get(9901);

        console.log(JSON.stringify({
            nGalaxies: cat.length,
            nGroups: cat.groups.size,
            groupMembers: group ? group.nMembers : 0,
            groupMeanCz: group ? group.meanVRec : 0,
            groupMVir: group ? group.mVir : 0
        }));
        """)
        assert res["nGalaxies"] == 3
        assert res["nGroups"] == 1
        assert res["groupMembers"] == 3
        assert math.isclose(res["groupMeanCz"], 4533.33, abs_tol=1.0)
        assert res["groupMVir"] > 0

    def test_high_density_catalog_mass_aggregation_benchmark(self):
        """Verify performance and memory integrity on 1,000 groups (3,000 galaxies)."""
        res = run_node_eval("""
        import { CF4GroupedCatalog } from './src/data/cf4_group_catalog.js';

        const cat = new CF4GroupedCatalog(4000);
        const nGroups = 1000;
        const membersPerGroup = 3;

        for (let g = 1; g <= nGroups; g++) {
            const baseRa = (g * 0.35) % 360;
            const baseDec = ((g * 0.17) % 160) - 80;
            const baseCz = 1000.0 + (g * 15.0);

            for (let m = 0; m < membersPerGroup; m++) {
                cat.addGalaxy({
                    pgc: g * 10 + m,
                    groupId: g,
                    ra: baseRa + m * 0.02,
                    dec: baseDec + m * 0.02,
                    cz: baseCz + m * 50.0,
                    kMag: 10.0 + m * 0.5
                });
            }
        }

        const startTime = Date.now();
        cat.aggregateGroups({ computeBootstrap: false });
        const elapsedMs = Date.now() - startTime;

        let totalVirialMass = 0.0;
        for (const grp of cat.groups.values()) {
            totalVirialMass += grp.mVir;
        }

        console.log(JSON.stringify({
            totalGalaxies: cat.length,
            totalGroups: cat.groups.size,
            totalVirialMass,
            elapsedMs,
            isFast: elapsedMs < 2000
        }));
        """)
        assert res["totalGalaxies"] == 3000
        assert res["totalGroups"] == 1000
        assert res["totalVirialMass"] > 0
        assert res["isFast"] is True

    def test_extreme_coordinate_projections_poles_and_antimeridian(self):
        """Verify coordinate transforms and angular separations at North Pole, South Pole, and RA=0/360 boundary."""
        res = run_node_eval("""
        import { angularSeparationDeg, projectedPhysicalSeparationMpc, CF4Group } from './src/data/cf4_group_catalog.js';

        // Pair spanning RA=0/360 boundary
        const sepAntimeridian = angularSeparationDeg(359.9, 10.0, 0.1, 10.0);

        // Pair near North Celestial Pole
        const sepPole = angularSeparationDeg(10.0, 89.9, 190.0, 89.9);

        // Group across RA=0/360
        const polarMembers = [
            { pgc: 1, ra: 359.9, dec: 0.0, cz: 2000.0, dist: 25.0 },
            { pgc: 2, ra: 0.1, dec: 0.0, cz: 2100.0, dist: 25.0 }
        ];
        const group = new CF4Group(1, polarMembers, { computeBootstrap: false });

        console.log(JSON.stringify({
            sepAntimeridian,
            sepPole,
            groupRH: group.rHarmonicMpc,
            groupMVir: group.mVir
        }));
        """)
        expectedSep = 0.2 * math.cos(math.radians(10.0))
        assert math.isclose(res["sepAntimeridian"], expectedSep, rel_tol=1e-3)
        assert res["sepPole"] < 0.5
        assert res["groupRH"] > 0
        assert res["groupMVir"] > 0

    def test_morphology_classification_edge_cases(self):
        """Verify de Vaucouleurs T-type classification across extreme T values (-5, 0, 1, 9, 10, 99)."""
        res = run_node_eval("""
        import { analyzeGroupLuminosityAndMorphology } from './src/data/cf4_group_catalog.js';

        const members = [
            { tType: -5.0, kMag: 9.0 }, // Pure Elliptical (E0-E6)
            { tType: 0.0, kMag: 9.5 },  // Lenticular / S0
            { tType: 1.0, kMag: 10.0 }, // Sa
            { tType: 9.0, kMag: 10.5 }, // Sm
            { tType: 10.0, kMag: 11.0 },// Irr
            { tType: 11.0, kMag: 11.5 },// Dwarf Irr
            { tType: -99.0, kMag: 12.0 } // Unclassified
        ];

        const analysis = analyzeGroupLuminosityAndMorphology(members, 20.0, 1.0e14);

        console.log(JSON.stringify(analysis));
        """)
        # 2 early (-5, 0), 2 late (1, 9), 2 irr (10, 11), 1 unclassified
        assert math.isclose(analysis_early := res["fEarly"], 2.0 / 6.0, rel_tol=1e-5)
        assert math.isclose(analysis_late := res["fLate"], 2.0 / 6.0, rel_tol=1e-5)
        assert math.isclose(analysis_irr := res["fIrr"], 2.0 / 6.0, rel_tol=1e-5)


class TestHierarchicalGroupPercolationAndDendrogram:
    """Tier 11: Hierarchical Percolation, Dendrogram Branching, and Multi-Scale FoF"""

    def test_multi_scale_fof_cluster_substructure(self):
        """Verify FoF with tight vs loose linking lengths resolves cluster cores and extended halos."""
        res = run_node_eval("""
        import { runHierarchicalFoF } from './src/data/cf4_group_catalog.js';

        // Coma-like cluster with tight core (10 galaxies) and loose outer envelope (5 galaxies)
        const galaxies = [];
        // Core: radius < 0.2 Mpc, cz = 6900 +- 150 km/s
        for (let i = 0; i < 10; i++) {
            galaxies.push({
                pgc: i + 1,
                ra: 194.95 + (i % 3 - 1) * 0.05,
                dec: 27.98 + (Math.floor(i / 3) - 1) * 0.05,
                cz: 6900.0 + (i * 25.0 - 100.0),
                dist: 95.0
            });
        }
        // Outer envelope: radius ~ 1.0 Mpc, cz = 7200 +- 400 km/s
        for (let i = 10; i < 15; i++) {
            galaxies.push({
                pgc: i + 1,
                ra: 194.95 + (i - 12) * 0.4,
                dec: 27.98 + (i - 12) * 0.3,
                cz: 7200.0 + (i * 60.0 - 150.0),
                dist: 95.0
            });
        }

        // Tight linking: should find core only
        const tightGroups = runHierarchicalFoF(galaxies, { d0: 0.15, v0: 200.0, minMembers: 3 });
        // Loose linking: should merge envelope and core
        const looseGroups = runHierarchicalFoF(galaxies, { d0: 0.80, v0: 800.0, minMembers: 3 });

        console.log(JSON.stringify({
            nTightGroups: tightGroups.size,
            nLooseGroups: looseGroups.size,
            tightCoreCount: tightGroups.values().next().value ? tightGroups.values().next().value.length : 0,
            looseClusterCount: looseGroups.values().next().value ? looseGroups.values().next().value.length : 0
        }));
        """)
        assert res["nTightGroups"] >= 1
        assert res["nLooseGroups"] == 1
        assert res["tightCoreCount"] >= 8
        assert res["looseClusterCount"] == 15


class TestMultiBandLuminosityFunctionAndFaintEndLimits:
    """Tier 12: Multi-band Completeness, Surface Brightness, and Morphological Segregation Invariants"""

    def test_multi_band_luminosity_ratio_invariants(self):
        """Verify K-band and B-band luminosities obey empirical color (B - K) cosmological bounds."""
        res = run_node_eval("""
        import { analyzeGroupLuminosityAndMorphology, M_SUN_K, M_SUN_B } from './src/data/cf4_group_catalog.js';

        const dGroup = 25.0; // Mpc
        // Group of elliptical galaxies: typical (B - K) ~ 4.0 - 4.5 mag
        const earlyMembers = [
            { kMag: 8.0, bMag: 12.2, tType: -3.0 },
            { kMag: 8.5, bMag: 12.7, tType: -2.0 },
            { kMag: 9.0, bMag: 13.2, tType: -1.0 }
        ];

        const analysis = analyzeGroupLuminosityAndMorphology(earlyMembers, dGroup, 5.0e14);

        // L_K / L_B ratio in solar units
        const lRatio = analysis.lTotK / analysis.lTotB;

        console.log(JSON.stringify({
            lTotK: analysis.lTotK,
            lTotB: analysis.lTotB,
            lRatio,
            mToLK: analysis.mToLK,
            mToLB: analysis.mToLB
        }));
        """)
        assert res["lTotK"] > res["lTotB"]
        # Typical early-type galaxy L_K / L_B is ~ 3 to 10 in solar units
        assert 2.0 < res["lRatio"] < 15.0
        assert res["mToLK"] > 0
        assert res["mToLB"] > res["mToLK"]


class TestStatisticalUncertaintiesAndBootstrapConvergence:
    """Tier 13: Bootstrap Resampling Convergence and Asymmetric Confidence Intervals"""

    def test_bootstrap_error_sample_size_scaling(self):
        """Verify bootstrap confidence interval width contracts with larger sample size N."""
        res = run_node_eval("""
        import { bootstrapVirialMassErrors } from './src/data/cf4_group_catalog.js';

        const dGroup = 50.0;
        // Small group: N = 4
        const smallMembers = [];
        for (let i = 0; i < 4; i++) {
            smallMembers.push({ ra: 100.0 + i * 0.1, dec: 10.0 + i * 0.1, cz: 3000.0 + (i * 100 - 150) });
        }

        // Large group: N = 25
        const largeMembers = [];
        for (let i = 0; i < 25; i++) {
            largeMembers.push({ ra: 100.0 + (i % 5) * 0.05, dec: 10.0 + Math.floor(i / 5) * 0.05, cz: 3000.0 + ((i * 37) % 300 - 150) });
        }

        const bootSmall = bootstrapVirialMassErrors(smallMembers, dGroup, 100, 42);
        const bootLarge = bootstrapVirialMassErrors(largeMembers, dGroup, 100, 42);

        const relWidthSmall = (bootSmall.mVirHigh - bootSmall.mVirLow) / bootSmall.mVirMedian;
        const relWidthLarge = (bootLarge.mVirHigh - bootLarge.mVirLow) / bootLarge.mVirMedian;

        console.log(JSON.stringify({
            relWidthSmall,
            relWidthLarge,
            largeIsTighter: relWidthLarge < relWidthSmall
        }));
        """)
        assert res["relWidthSmall"] > 0
        assert res["relWidthLarge"] > 0
        assert res["largeIsTighter"] is True


class TestAdversarialSentinelsAndCorruptInputRecovery:
    """Tier 14: Robust Parsing with Extreme Sentinels, Nulls, and Corrupt Lines"""

    def test_sentinel_replacement_and_missing_field_defaults(self):
        """Verify parser ignores corrupt lines and sentinels (-99.99, 99.99, NaN, null)."""
        res = run_node_eval("""
        import { parseCF4GroupCatalogText } from './src/data/cf4_group_catalog.js';

        const dirtyCsv = [
            '# Corrupted Header and metadata',
            'PGC,GroupID,RA,Dec,cz,vErr,Dist,BMag,KMag,TType',
            '101,901,180.0,10.0,2000.0,-99.99,null,99.99,NaN,-99.0',
            '102,901,180.1,10.1,2100.0,25.0,27.5,12.0,9.0,3.0',
            'CORRUPT LINE WITH TEXT ONLY NO NUMBERS',
            '103,901,180.2,10.2,2050.0,20.0,28.0,13.0,10.0,5.0',
            ',,,,',
            '104,901,NaN,NaN,2050.0,20.0,28.0,13.0,10.0,5.0'
        ].join('\\n');

        const cat = parseCF4GroupCatalogText(dirtyCsv);
        const group = cat.groups.get(901);

        console.log(JSON.stringify({
            totalGalaxies: cat.length,
            groupNMembers: group ? group.nMembers : 0,
            groupMeanCz: group ? group.meanVRec : 0,
            groupMVir: group ? group.mVir : 0
        }));
        """)
        # Only rows 101, 102, 103 are valid (row 104 has NaN RA/Dec, corrupt line is skipped, empty line skipped)
        assert res["totalGalaxies"] == 3
        assert res["groupNMembers"] == 3
        assert math.isclose(res["groupMeanCz"], 2050.0, abs_tol=1.0)
        assert res["groupMVir"] > 0

    def test_schechter_faint_end_slope_divergence_guard(self):
        """Verify numerical stability of upper incomplete gamma and Schechter integral near critical bounds."""
        res = run_node_eval("""
        import { upperIncompleteGamma, schechterCompletenessCorrection } from './src/data/cf4_group_catalog.js';

        // Extreme faint limit (very distant or very shallow survey)
        const cDistant = schechterCompletenessCorrection(8.0, 500.0, -23.5, -1.05);
        // Very deep survey (completeness -> 1.0)
        const cDeep = schechterCompletenessCorrection(25.0, 5.0, -23.5, -1.05);

        // Near-zero x evaluation
        const gNearZero = upperIncompleteGamma(0.95, 1.0e-8);

        console.log(JSON.stringify({
            cDistant,
            cDeep,
            gNearZero,
            isFiniteDistant: Number.isFinite(cDistant),
            isFiniteDeep: Number.isFinite(cDeep),
            isFiniteGamma: Number.isFinite(gNearZero)
        }));
        """)
        assert res["isFiniteDistant"] is True
        assert res["isFiniteDeep"] is True
        assert res["isFiniteGamma"] is True
        assert 0.0 < res["cDistant"] <= 1.0
        assert math.isclose(res["cDeep"], 1.0, rel_tol=1e-3)

    def test_spatial_hash_nearest_neighbor_ordering(self):
        """Verify 3D Spatial Hash Grid queries return results strictly sorted in ascending Euclidean distance."""
        res = run_node_eval("""
        import { SpatialHashGrid3D, euclideanDistance3D } from './src/data/cf4_group_catalog.js';

        const grid = new SpatialHashGrid3D(10.0);
        // Insert points at concentric distances: 2, 5, 8, 12, 18 Mpc/h from (10, 20, 30)
        const radii = [18.0, 5.0, 2.0, 12.0, 8.0];
        radii.forEach((r, idx) => {
            grid.insert({ id: idx, targetR: r }, 10.0 + r, 20.0, 30.0);
        });

        const hits = grid.queryRadius(10.0, 20.0, 30.0, 20.0);
        const distances = hits.map(h => h.distance);

        let isAscending = true;
        for (let i = 1; i < distances.length; i++) {
            if (distances[i] < distances[i - 1]) {
                isAscending = false;
                break;
            }
        }

        console.log(JSON.stringify({
            nHits: hits.length,
            distances,
            isAscending
        }));
        """)
        assert res["nHits"] == 5
        assert res["isAscending"] is True
        assert math.isclose(res["distances"][0], 2.0, abs_tol=1e-5)
        assert math.isclose(res["distances"][-1], 18.0, abs_tol=1e-5)



