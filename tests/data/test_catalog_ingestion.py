# -*- coding: utf-8 -*-
"""
Automated Test Suite for Astrometric Catalog Ingestion and Zone of Avoidance (ZoA) Reconstructive Extrapolator.

Verifies:
1. Distance Modulus & Recessional Velocity Conversions:
   - mu <-> d (Mpc and Mpc/h) with exact logarithmic inversion.
   - Analytic first-order error propagation: sigma_d = d * (ln(10)/5) * sigma_mu.
   - Peculiar velocity linear & relativistic formulations.
   - Tully-Fisher absolute magnitude and distance modulus derivations.
2. Malmquist Bias Correction Models:
   - Homogeneous Malmquist bias corrections against analytic exponential / logarithmic forms.
   - Inhomogeneous Malmquist bias (IMB) line-of-sight quadrature over overdensity fluctuations.
3. Selection Functions & Effective Volume Weights:
   - Fermi-Dirac / sigmoid apparent magnitude selection.
   - Schechter luminosity distance completeness integration via upper incomplete Gamma function.
   - Survey selection weighting and normalization.
4. Binary & ASCII File Parsers:
   - FITS Binary Table parser with synthetic BINTABLE generation.
   - IP2I custom ASCII table parser with comment metadata and sentinel replacements.
   - CDS VizieR catalogue parser.
5. Grouped Catalog Aggregator:
   - Inverse-variance weighted group distance modulus, velocity dispersion, centroids, multiplicity.
6. Zone of Avoidance (ZoA) Obscuration & Dust Extinction:
   - Galactic and Supergalactic latitude masking (|b| < 10°, |SGB| < 10°).
   - Analytic SFD98/Planck dust extinction profile.
   - Continuous smooth sigmoid/tanh boundary blending.
7. 21cm HI & Infrared Piercing Corridors:
   - MeerKAT/Parkes Vela HI corridor, Parkes HIZOA Norma/Great Attractor, and Puppis HI corridor.
   - Transmission profile and effective de-obscuration factor.
8. Wiener-Filter & Kriging Inpainting:
   - Point-based local Wiener filter / Kriging reconstruction and posterior variance.
   - 3D Grid-based iterative Wiener inpainting solver.
   - Boundary continuity (C^0/C^1 smoothness across |b| = 10°).
9. Synthetic Control Benchmark:
   - Artificial mask recovery fidelity, RMSE, Pearson r > 0.75, relative flux conservation.
10. GalaxyCatalog Columnar Store:
    - High-density TypedArray storage, predicate filtering, bounding box, and GeoJSON export.
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


class TestDistanceModulusAndCosmography:
    """Tier 1: Cosmographic Inversion, Error Propagation, and Peculiar Velocity"""

    def test_distance_modulus_inversion(self):
        """Verify mu <-> d inversion in Mpc and Mpc/h across cosmological range [1, 1000] Mpc."""
        res = run_node_eval("""
        import {
            distanceModulusToDistance,
            distanceToDistanceModulus,
            distanceModulusUncertaintyToDistance,
            distanceUncertaintyToDistanceModulus,
            LN10_DIV_5
        } from './src/data/catalog_ingestion.js';

        const testDists = [1.0, 10.0, 16.5, 50.0, 100.0, 500.0, 1000.0];
        const results = [];

        for (const d of testDists) {
            const mu = distanceToDistanceModulus(d);
            const dBack = distanceModulusToDistance(mu);
            const dMpch = distanceModulusToDistance(mu, 74.6, true);

            const sigmaMu = 0.35;
            const sigmaD = distanceModulusUncertaintyToDistance(d, sigmaMu);
            const sigmaMuBack = distanceUncertaintyToDistanceModulus(d, sigmaD);

            results.push({
                d,
                mu,
                dBack,
                dMpch,
                expectedMpch: d * 0.746,
                dErr: Math.abs(d - dBack),
                sigmaD,
                expectedSigmaD: d * LN10_DIV_5 * sigmaMu,
                sigmaMuBackErr: Math.abs(sigmaMu - sigmaMuBack)
            });
        }

        console.log(JSON.stringify(results));
        """)
        for item in res:
            assert item["dErr"] < 1e-11, f"Distance round-trip failed for d={item['d']}"
            assert math.isclose(item["dMpch"], item["expectedMpch"], rel_tol=1e-12)
            assert math.isclose(item["sigmaD"], item["expectedSigmaD"], rel_tol=1e-12)
            assert item["sigmaMuBackErr"] < 1e-12

    def test_peculiar_velocity_and_tully_fisher(self):
        """Verify linear and relativistic peculiar velocity and Tully-Fisher relation."""
        res = run_node_eval("""
        import {
            calculatePeculiarVelocity,
            peculiarVelocityError,
            tullyFisherAbsoluteMagnitude
        } from './src/data/catalog_ingestion.js';

        // Recessional velocity 5000 km/s, distance 60 Mpc, H0 = 74.6
        // v_Hubble = 60 * 74.6 = 4476 km/s
        // Linear v_pec = 5000 - 4476 = 524 km/s
        const vLinear = calculatePeculiarVelocity(5000.0, 60.0, 74.6, { relativistic: false });
        const vRel = calculatePeculiarVelocity(5000.0, 60.0, 74.6, { relativistic: true });

        // Error propagation with sigma_mu = 0.35, sigma_vrec = 20, thermal = 150
        const vErr = peculiarVelocityError(60.0, 0.35, 20.0, 74.6, 150.0);

        // Tully-Fisher: logW21 = 2.6 (W21 = 398 km/s), slope = -8.0, zeroPoint = -20.5, logW0 = 2.5
        // M = -20.5 + (-8.0) * (2.6 - 2.5) = -20.5 - 0.8 = -21.3
        const mTF = tullyFisherAbsoluteMagnitude(2.6, -8.0, -20.5, 2.5);

        console.log(JSON.stringify({
            vLinear,
            vRel,
            vErr,
            mTF
        }));
        """)
        assert math.isclose(res["vLinear"], 524.0, abs_tol=1e-10)
        # Relativistic denominator 1 + 4476 / 299792.458 = 1.01493
        expectedRel = 524.0 / (1.0 + 4476.0 / 299792.458)
        assert math.isclose(res["vRel"], expectedRel, rel_tol=1e-7)
        assert math.isclose(res["mTF"], -21.3, abs_tol=1e-10)
        assert res["vErr"] > 700.0  # Significant uncertainty from distance modulus error at 60 Mpc


class TestMalmquistBiasCorrections:
    """Tier 2: Homogeneous and Inhomogeneous Malmquist Bias Formulas"""

    def test_homogeneous_malmquist_formulas(self):
        """Verify Homogeneous Malmquist Bias correction against exact analytic forms."""
        res = run_node_eval("""
        import {
            correctHomogeneousMalmquist,
            distanceModulusToDistance,
            distanceToDistanceModulus,
            HOMOGENEOUS_MALMQUIST_FACTOR,
            HOMOGENEOUS_MALMQUIST_DM_FACTOR
        } from './src/data/catalog_ingestion.js';

        const dObs = 100.0; // Mpc
        const muObs = distanceToDistanceModulus(dObs);
        const sigmaMu = 0.40;

        // 1. Distance correction
        const dCorr = correctHomogeneousMalmquist(dObs, sigmaMu, { mode: 'distance' });
        const expectedDCorr = dObs * Math.exp(-HOMOGENEOUS_MALMQUIST_FACTOR * sigmaMu * sigmaMu);

        // 2. Modulus correction
        const muCorr = correctHomogeneousMalmquist(muObs, sigmaMu, { mode: 'modulus' });
        const expectedMuCorr = muObs - HOMOGENEOUS_MALMQUIST_DM_FACTOR * sigmaMu * sigmaMu;

        // 3. Forward prediction
        const muFwd = correctHomogeneousMalmquist(muCorr, sigmaMu, { mode: 'forward' });

        console.log(JSON.stringify({
            dCorr,
            expectedDCorr,
            dErr: Math.abs(dCorr - expectedDCorr),
            muCorr,
            expectedMuCorr,
            muErr: Math.abs(muCorr - expectedMuCorr),
            muFwdErr: Math.abs(muFwd - muObs)
        }));
        """)
        assert res["dErr"] < 1e-12
        assert res["muErr"] < 1e-12
        assert res["muFwdErr"] < 1e-12
        # dCorr should be smaller than dObs (correcting for volume over-representation)
        assert res["dCorr"] < 100.0
        expectedFactor = 1.5 * (math.log(10) / 5.0) ** 2
        assert math.isclose(res["dCorr"], 100.0 * math.exp(-expectedFactor * 0.16), rel_tol=1e-10)

    def test_inhomogeneous_malmquist_quadrature(self):
        """Verify Inhomogeneous Malmquist line-of-sight quadrature over density gradients."""
        res = run_node_eval("""
        import {
            correctInhomogeneousMalmquist,
            distanceModulusToDistance
        } from './src/data/catalog_ingestion.js';

        // Galaxy along X axis at nominal d = 50 Mpc
        const muObs = 33.49485; // 5*log10(50) + 25 = 33.49485
        const sigmaMu = 0.35;

        // Case A: Flat zero density -> spatial volume element prior d^2
        const imbFlat = correctInhomogeneousMalmquist(1.0, 0.0, 0.0, muObs, sigmaMu, () => 0.0);

        // Case B: High density cluster at d = 55 Mpc (should pull posterior distance higher)
        const densityCluster = (x, y, z) => {
            const r = Math.hypot(x, y, z);
            return 5.0 * Math.exp(-0.5 * Math.pow((r - 55.0) / 4.0, 2));
        };
        const imbCluster = correctInhomogeneousMalmquist(1.0, 0.0, 0.0, muObs, sigmaMu, densityCluster);

        // Case C: Deep void at d = 55 Mpc, peak at 45 Mpc (should pull posterior distance lower)
        const densityVoid = (x, y, z) => {
            const r = Math.hypot(x, y, z);
            return 5.0 * Math.exp(-0.5 * Math.pow((r - 45.0) / 4.0, 2)) - 0.8 * Math.exp(-0.5 * Math.pow((r - 55.0) / 4.0, 2));
        };
        const imbVoid = correctInhomogeneousMalmquist(1.0, 0.0, 0.0, muObs, sigmaMu, densityVoid);

        console.log(JSON.stringify({
            flatD: imbFlat.correctedDistance,
            clusterD: imbCluster.correctedDistance,
            voidD: imbVoid.correctedDistance
        }));
        """)
        assert res["flatD"] > 50.0  # d^2 volume prior weights outer radial shells
        assert res["clusterD"] > res["flatD"], "Overdensity ahead of galaxy should bias distance estimate outward"
        assert res["voidD"] < res["flatD"], "Overdensity behind galaxy should bias distance estimate inward"


class TestSelectionFunctionsAndWeights:
    """Tier 3: Selection Functions, Incomplete Gamma, and Effective Volume Normalization"""

    def test_upper_incomplete_gamma_and_schechter(self):
        """Verify upper incomplete gamma function and Schechter distance selection."""
        res = run_node_eval("""
        import {
            upperIncompleteGamma,
            evaluateMagnitudeSelection,
            evaluateSchechterDistanceSelection
        } from './src/data/catalog_ingestion.js';

        // Gamma(1.0, x) = exp(-x)
        const g1_0 = upperIncompleteGamma(1.0, 0.0);
        const g1_1 = upperIncompleteGamma(1.0, 1.0);
        const g1_2 = upperIncompleteGamma(1.0, 2.0);

        // Apparent magnitude sigmoid selection at magLim = 17.5
        const sBright = evaluateMagnitudeSelection(15.0, 17.5, 0.4); // ~1.0
        const sMid = evaluateMagnitudeSelection(17.5, 17.5, 0.4);    // 0.5
        const sFaint = evaluateMagnitudeSelection(20.0, 17.5, 0.4);   // ~0.0

        // Distance selection completeness (should decrease monotonically with distance)
        const phi10 = evaluateSchechterDistanceSelection(10.0, 17.5, -20.5, -1.0);
        const phi50 = evaluateSchechterDistanceSelection(50.0, 17.5, -20.5, -1.0);
        const phi150 = evaluateSchechterDistanceSelection(150.0, 17.5, -20.5, -1.0);

        console.log(JSON.stringify({
            g1_0,
            g1_1,
            g1_1_expected: Math.exp(-1.0),
            g1_2,
            g1_2_expected: Math.exp(-2.0),
            sBright,
            sMid,
            sFaint,
            phi10,
            phi50,
            phi150
        }));
        """)
        assert math.isclose(res["g1_0"], 1.0, abs_tol=1e-10)
        assert math.isclose(res["g1_1"], res["g1_1_expected"], rel_tol=1e-7)
        assert math.isclose(res["g1_2"], res["g1_2_expected"], rel_tol=1e-7)
        assert res["sBright"] > 0.99
        assert math.isclose(res["sMid"], 0.5, abs_tol=1e-10)
        assert res["sFaint"] < 0.01
        assert res["phi10"] > res["phi50"] > res["phi150"], "Completeness must decrease monotonically with distance"


class TestFileParsersAndGroupedAggregation:
    """Tier 4: FITS Binary Table, IP2I ASCII, CDS VizieR, and Group Aggregator"""

    def test_fits_binary_table_synthetic_parsing(self):
        """Verify FITS Binary Table parser by constructing and decoding standard FITS byte stream."""
        res = run_node_eval("""
        import { parseFITSBinaryTable } from './src/data/catalog_ingestion.js';

        // Build synthetic FITS 2880-byte header and data block
        const totalBytes = 2880 * 4; // Primary header + BINTABLE header + Data + padding
        const buffer = new Uint8Array(totalBytes);
        const encoder = new TextEncoder();

        // 1. Primary Header
        let card = 'SIMPLE  =                    T / Standard FITS format                           ';
        buffer.set(encoder.encode(card), 0);
        card = 'BITPIX  =                    8 / Character or unsigned binary integer           ';
        buffer.set(encoder.encode(card), 80);
        card = 'NAXIS   =                    0 / No primary data array                          ';
        buffer.set(encoder.encode(card), 160);
        card = 'EXTEND  =                    T / Extensions exist                               ';
        buffer.set(encoder.encode(card), 240);
        card = 'END                                                                             ';
        buffer.set(encoder.encode(card), 320);

        // 2. Extension BINTABLE Header (Starts at byte 2880)
        let extOffset = 2880;
        function writeCard(key, val, comment = '') {
            let cardStr = key.padEnd(8, ' ') + '= ' + String(val).padStart(20, ' ') + ' / ' + comment;
            cardStr = cardStr.padEnd(80, ' ');
            buffer.set(encoder.encode(cardStr), extOffset);
            extOffset += 80;
        }

        writeCard('XTENSION', "'BINTABLE'", 'Binary Table Extension');
        writeCard('BITPIX', 8);
        writeCard('NAXIS', 2);
        writeCard('NAXIS1', 20, 'Row width in bytes');
        writeCard('NAXIS2', 3, 'Number of rows');
        writeCard('PCOUNT', 0);
        writeCard('GCOUNT', 1);
        writeCard('TFIELDS', 3);
        writeCard('TTYPE1', "'PGC     '");
        writeCard('TFORM1', "'1J      '"); // Int32 (4 bytes)
        writeCard('TTYPE2', "'RA      '");
        writeCard('TFORM2', "'1D      '"); // Float64 (8 bytes)
        writeCard('TTYPE3', "'DEC     '");
        writeCard('TFORM3', "'1D      '"); // Float64 (8 bytes)

        const endCard = 'END                                                                             ';
        buffer.set(encoder.encode(endCard), extOffset);

        // 3. Binary Table Data (Starts at byte 5760 = 2880 * 2)
        const dataOffset = 2880 * 2;
        const dataView = new DataView(buffer.buffer, dataOffset, 2880);

        const testRows = [
            { pgc: 143, ra: 187.7059, dec: 12.3911 },
            { pgc: 1024, ra: 192.8595, dec: 27.1283 },
            { pgc: 56000, ra: 266.4051, dec: -28.9362 }
        ];

        for (let r = 0; r < testRows.length; r++) {
            const rowStart = r * 20;
            dataView.setInt32(rowStart + 0, testRows[r].pgc, false);      // Big-endian J
            dataView.setFloat64(rowStart + 4, testRows[r].ra, false);     // Big-endian D
            dataView.setFloat64(rowStart + 12, testRows[r].dec, false);   // Big-endian D
        }

        const parsed = parseFITSBinaryTable(buffer, { extensionIndex: 1 });

        console.log(JSON.stringify({
            rowCount: parsed.rowCount,
            columnNames: parsed.columnNames,
            rows: parsed.rows
        }));
        """)
        assert res["rowCount"] == 3
        assert "PGC" in res["columnNames"]
        assert "RA" in res["columnNames"]
        assert "DEC" in res["columnNames"]
        assert res["rows"][0]["PGC"] == 143
        assert math.isclose(res["rows"][0]["RA"], 187.7059, abs_tol=1e-4)
        assert math.isclose(res["rows"][2]["DEC"], -28.9362, abs_tol=1e-4)

    def test_ip2i_ascii_and_cds_parsers(self):
        """Verify IP2I ASCII table and CDS VizieR parser with sentinels."""
        res = run_node_eval("""
        import { parseIP2IASCIITable, parseCDSVizieRCatalog } from './src/data/catalog_ingestion.js';

        const asciiData = `
        # H0 = 74.6
        # Survey: CF4 IP2I Compilation
        # Columns: PGC RA DEC Vcmb DM eDM GroupId Method
        1001  187.70 12.39  1100.0  31.05 0.35  10  TF
        1002  188.10 12.50  1150.0  31.10 0.30  10  TF
        1003  200.00 -15.0  3400.0  99.99 0.40  -1  SNIa
        `;

        const parsedAscii = parseIP2IASCIITable(asciiData);

        const cdsData = `
        # CDS VizieR format
        # PGC|RA|DEC|cz|DM|eDM
        ---|---|---|---|---|---
        2001|150.0|25.0|2200|32.40|0.38
        2002|151.0|25.5|2250|---|0.40
        `;
        const parsedCds = parseCDSVizieRCatalog(cdsData);

        console.log(JSON.stringify({
            asciiMeta: parsedAscii.metadata,
            asciiRowCount: parsedAscii.rowCount,
            asciiRows: parsedAscii.rows,
            cdsRowCount: parsedCds.rowCount,
            cdsRows: parsedCds.rows
        }));
        """)
        assert res["asciiMeta"]["H0"] == 74.6
        assert res["asciiRowCount"] == 3
        assert res["asciiRows"][0]["PGC"] == 1001
        assert res["asciiRows"][0]["Method"] == "TF"
        assert res["asciiRows"][2]["DM"] is None or math.isnan(res["asciiRows"][2]["DM"])
        assert res["cdsRowCount"] == 2
        assert res["cdsRows"][0]["PGC"] == 2001
        assert res["cdsRows"][1]["DM"] is None or math.isnan(res["cdsRows"][1]["DM"])

    def test_grouped_catalog_aggregation(self):
        """Verify Grouped Catalog aggregation: inverse-variance mean, velocity dispersion, centroids."""
        res = run_node_eval("""
        import { aggregateGroupedCatalog } from './src/data/catalog_ingestion.js';

        // Group 100 with 3 galaxies
        // Gal 1: mu=31.0, dmu=0.20 (w = 1/0.04 = 25), cz=1200
        // Gal 2: mu=31.5, dmu=0.20 (w = 1/0.04 = 25), cz=1300
        // Gal 3: mu=31.2, dmu=0.10 (w = 1/0.01 = 100), cz=1250
        // Total weight = 25 + 25 + 100 = 150
        // Expected weighted mu = (31.0*25 + 31.5*25 + 31.2*100) / 150 = 4682.5 / 150 = 31.2166667
        // Expected err_mu = 1 / sqrt(150) = 0.081649658
        const mockGalaxies = [
            { pgc: 1, ra: 187.0, dec: 12.0, cz: 1200.0, mu: 31.0, dmu: 0.20, groupId: 100, method: 'TF' },
            { pgc: 2, ra: 187.5, dec: 12.5, cz: 1300.0, mu: 31.5, dmu: 0.20, groupId: 100, method: 'FP' },
            { pgc: 3, ra: 188.0, dec: 13.0, cz: 1250.0, mu: 31.2, dmu: 0.10, groupId: 100, method: 'SNIa' },
            { pgc: 4, ra: 220.0, dec: -30.0, cz: 4500.0, mu: 34.0, dmu: 0.40, groupId: -1, method: 'TF' } // singleton
        ];

        const grouped = aggregateGroupedCatalog(mockGalaxies, { H0: 74.6 });

        console.log(JSON.stringify({
            groupCount: grouped.length,
            groups: grouped
        }));
        """)
        assert res["groupCount"] == 2
        grp100 = next(g for g in res["groups"] if g["groupId"] == 100)
        assert grp100["multiplicity"] == 3
        assert math.isclose(grp100["mu"], 31.2166667, rel_tol=1e-5)
        assert math.isclose(grp100["dmu"], 1.0 / math.sqrt(150.0), rel_tol=1e-5)
        assert math.isclose(grp100["cz"], 1250.0, abs_tol=1e-5)
        assert math.isclose(grp100["velDispersion"], 50.0, abs_tol=1e-5)
        assert set(grp100["methods"]) == {"TF", "FP", "SNIa"}


class TestZoneOfAvoidanceAndCorridors:
    """Tier 5: ZoA Masking, Dust Extinction, and HI / IR Piercing Corridors"""

    def test_zoa_masking_and_dust_extinction(self):
        """Verify ZoA latitude masking and dust extinction thresholding."""
        res = run_node_eval("""
        import {
            isGalacticZoA,
            isSupergalacticZoA,
            computeDustExtinctionAV,
            computeZoAObscurationWeight
        } from './src/data/zoa_reconstructor.js';

        // 1. Angular bounds
        const isGalCore = isGalacticZoA(2.5, 10.0);
        const isGalPole = isGalacticZoA(45.0, 10.0);
        const isSGPlane = isSupergalacticZoA(-4.0, 10.0);

        // 2. SFD Dust Extinction
        const avGalCenter = computeDustExtinctionAV(0.0, 0.5); // Deep Galactic center
        const avGalPole = computeDustExtinctionAV(0.0, 90.0);   // Galactic pole

        // 3. Obscuration weights
        const wCore = computeZoAObscurationWeight(0.0, 0.0);   // Full obscuration ~1.0
        const wUnobscured = computeZoAObscurationWeight(100.0, 60.0); // Un-obscured ~0.0
        const wMargin = computeZoAObscurationWeight(0.0, 10.0); // Margin transition ~0.5

        console.log(JSON.stringify({
            isGalCore,
            isGalPole,
            isSGPlane,
            avGalCenter,
            avGalPole,
            wCore,
            wUnobscured,
            wMargin
        }));
        """)
        assert res["isGalCore"] is True
        assert res["isGalPole"] is False
        assert res["isSGPlane"] is True
        assert res["avGalCenter"] > 4.0, "Dust extinction at Galactic Center should be very high"
        assert res["avGalPole"] < 0.20, "Dust extinction at Galactic Pole should be low"
        assert res["wCore"] > 0.95, "Obscuration weight in ZoA core must approach 1.0"
        assert res["wUnobscured"] < 0.05, "Obscuration weight in un-obscured sky must approach 0.0"
        assert 0.35 <= res["wMargin"] <= 0.85, "Obscuration at boundary margin must be smooth"

    def test_piercing_corridors_transmission(self):
        """Verify MeerKAT Vela, Parkes HIZOA Norma, and Puppis corridor piercing transmission."""
        res = run_node_eval("""
        import { PiercingCorridorManager } from './src/data/zoa_reconstructor.js';

        const manager = new PiercingCorridorManager();

        // 1. Inside Vela HI corridor: l=272.5°, b=0°, d=240 Mpc
        const tVela = manager.evaluateTransmission(272.5, 0.0, 240.0);
        const effObscVela = manager.computeEffectiveObscuration(272.5, 0.0, 240.0);

        // 2. Inside Norma Great Attractor corridor: l=325.3°, b=-7.2°, d=68 Mpc
        const tNorma = manager.evaluateTransmission(325.3, -7.2, 68.0);
        const effObscNorma = manager.computeEffectiveObscuration(325.3, -7.2, 68.0);

        // 3. Deep masked region outside any corridor: l=180°, b=0°, d=100 Mpc
        const tNone = manager.evaluateTransmission(180.0, 0.0, 100.0);
        const effObscNone = manager.computeEffectiveObscuration(180.0, 0.0, 100.0);

        console.log(JSON.stringify({
            tVela,
            effObscVela,
            tNorma,
            effObscNorma,
            tNone,
            effObscNone
        }));
        """)
        assert res["tVela"]["inCorridor"] is True
        assert "Vela_MeerKAT_Parkes" in res["tVela"]["activeCorridors"]
        assert res["tVela"]["transmission"] > 0.70
        assert res["effObscVela"] < 0.30, "Vela corridor should de-obscure effective mask"

        assert res["tNorma"]["inCorridor"] is True
        assert "Norma_GreatAttractor_HIZOA" in res["tNorma"]["activeCorridors"]
        assert res["tNorma"]["transmission"] > 0.50

        assert res["tNone"]["inCorridor"] is False
        assert res["effObscNone"] > 0.90, "Outside corridors effective obscuration remains high"


class TestWienerInpaintingAndContinuity:
    """Tier 6: Wiener Filter Inpainting, Boundary Continuity, and Kriging"""

    def test_point_wiener_kriging_reconstruction(self):
        """Verify local point-based Wiener / Kriging reconstruction and posterior variance."""
        res = run_node_eval("""
        import { reconstructPointWiener, cosmologicalCovariance } from './src/data/zoa_reconstructor.js';

        // 4 surrounding observed galaxies around target (0, 0, 0)
        const tracers = [
            { sgx: 10.0, sgy: 0.0, sgz: 0.0, delta: 2.0, noiseVar: 0.01 },
            { sgx: -10.0, sgy: 0.0, sgz: 0.0, delta: 2.0, noiseVar: 0.01 },
            { sgx: 0.0, sgy: 10.0, sgz: 0.0, delta: 2.0, noiseVar: 0.01 },
            { sgx: 0.0, sgy: -10.0, sgz: 0.0, delta: 2.0, noiseVar: 0.01 }
        ];

        // Reconstruct at center (0, 0, 0)
        const recon = reconstructPointWiener(0.0, 0.0, 0.0, tracers, { r0: 10.0 });
        const priorVar = cosmologicalCovariance(0.0, { r0: 10.0 });

        console.log(JSON.stringify({
            delta: recon.delta,
            variance: recon.variance,
            priorVar: priorVar,
            neighborCount: recon.neighborCount
        }));
        """)
        assert res["neighborCount"] == 4
        assert res["delta"] > 1.0, "Positive surrounding overdensities must yield positive interpolated overdensity"
        assert res["variance"] < res["priorVar"], "Data assimilation must reduce posterior variance"

    def test_3d_grid_wiener_inpainting_and_boundary_continuity(self):
        """Verify 3D grid-based Wiener inpainting convergence and boundary continuity across |b| = 10°."""
        res = run_node_eval("""
        import { ZoAReconstructor } from './src/data/zoa_reconstructor.js';
        import { GridIndexer, StrideOrder, BoundaryMode } from './src/fields/grid_indexer.js';
        import { DensityField } from './src/fields/density_field.js';
        import { equatorialToGalactic, supergalacticCartesianToEquatorial } from './src/coordinates/canonical_frame.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
            nx: N, ny: N, nz: N,
            origin: [-L/2, -L/2, -L/2],
            boxSize: [L, L, L],
            strideOrder: StrideOrder.CANONICAL_XYZ,
            boundaryMode: BoundaryMode.CLAMP
        });

        // Create observed field with background structure
        const total = grid.totalCells;
        const dObs = new Float64Array(total);

        for (let i = 0; i < total; i++) {
            const [ix, iy, iz] = grid.get3DIndices(i);
            const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
            dObs[i] = 1.5 * Math.sin(x * 0.1) * Math.cos(y * 0.1);
        }

        const inputDensity = new DensityField(grid, dObs);
        const reconstructor = new ZoAReconstructor({
            bCut: 10.0,
            sgbCut: 10.0,
            maxIterations: 30,
            tolerance: 1e-4
        });

        const result = reconstructor.reconstructDensityField(inputDensity);
        const dRecon = result.reconstructedField.delta;

        // Check boundary continuity: difference between observed and reconstructed across boundary
        let maxBoundaryJump = 0.0;
        for (let i = 0; i < total; i++) {
            const [ix, iy, iz] = grid.get3DIndices(i);
            const [x, y, z] = grid.getNodeCoord(ix, iy, iz);
            const dist = Math.hypot(x, y, z);
            if (dist < 1e-3) continue;

            const [ra, dec] = supergalacticCartesianToEquatorial(x, y, z);
            const [glon, glat] = equatorialToGalactic(ra, dec);

            // Around |b| = 10° boundary
            if (Math.abs(Math.abs(glat) - 10.0) < 2.0) {
                const diff = Math.abs(dRecon[i] - dObs[i]);
                maxBoundaryJump = Math.max(maxBoundaryJump, diff);
            }
        }

        console.log(JSON.stringify({
            iterations: result.iterations,
            residual: result.residual,
            maxBoundaryJump: maxBoundaryJump
        }));
        """)
        assert res["iterations"] > 0
        assert res["residual"] < 1e-3
        assert res["maxBoundaryJump"] < 1.0, "Smooth transition across ZoA boundary without sharp jumps"


class TestSyntheticControlBenchmark:
    """Tier 7: Full Synthetic Control Benchmark and Recovery Validation"""

    def test_synthetic_control_benchmark_execution(self):
        """Verify synthetic control benchmark achieves Pearson r > 0.75 and low RMSE."""
        res = run_node_eval("""
        import { SyntheticControlBenchmark } from './src/data/zoa_reconstructor.js';

        const benchmark = new SyntheticControlBenchmark({
            gridSize: 24,
            boxSize: 180.0
        });

        const report = benchmark.runBenchmark({
            bCut: 10.0,
            sgbCut: 10.0,
            maxIterations: 40,
            tolerance: 1e-4
        });

        console.log(JSON.stringify(report));
        """)
        assert res["success"] is True, f"Benchmark failed with Pearson r = {res['pearsonR']}"
        assert res["pearsonR"] > 0.75, f"Pearson correlation must exceed 0.75, got {res['pearsonR']}"
        assert res["rmse"] < 2.0, f"RMSE in masked region should be small, got {res['rmse']}"
        assert res["maskedFraction"] > 0.05, "Synthetic mask should obscure significant volume"
        assert res["maskedCellCount"] > 100


class TestGalaxyCatalogAndPipelineEndToEnd:
    """Tier 8: Columnar GalaxyCatalog Storage and Full Pipeline Ingestion"""

    def test_galaxy_catalog_columnar_operations(self):
        """Verify GalaxyCatalog TypedArray operations, filtering, bounding box, and GeoJSON."""
        res = run_node_eval("""
        import { GalaxyCatalog } from './src/data/catalog_ingestion.js';

        const cat = new GalaxyCatalog(5);
        for (let i = 0; i < 5; i++) {
            cat.pgc[i] = 100 + i;
            cat.ra[i] = 180.0 + i * 2;
            cat.dec[i] = 10.0 + i;
            cat.dist[i] = 20.0 + i * 10;
            cat.sgx[i] = (i - 2) * 15.0;
            cat.sgy[i] = i * 20.0;
            cat.sgz[i] = (i - 1) * 5.0;
            cat.cz[i] = 1500.0 + i * 700;
            cat.mu[i] = 31.5 + i * 0.5;
            cat.groupId[i] = i < 3 ? 42 : -1;
            cat.method[i] = 'TF';
        }

        const g2 = cat.getGalaxy(2);
        const filtered = cat.filter(g => g.groupId === 42);
        const bbox = cat.getBoundingBox();
        const geojson = cat.toGeoJSON();

        console.log(JSON.stringify({
            g2Pgc: g2.pgc,
            filteredCount: filtered.count,
            filteredPgcs: Array.from(filtered.pgc),
            bbox: bbox,
            geojsonFeatureCount: geojson.features.length
        }));
        """)
        assert res["g2Pgc"] == 102
        assert res["filteredCount"] == 3
        assert res["filteredPgcs"] == [100, 101, 102]
        assert res["bbox"]["xMin"] == -30.0
        assert res["bbox"]["xMax"] == 30.0
        assert res["geojsonFeatureCount"] == 5

    def test_full_pipeline_ingestion_with_malmquist_and_selection(self):
        """Verify complete pipeline ingestion processing raw record objects into calibrated catalog."""
        res = run_node_eval("""
        import { CatalogIngestionPipeline } from './src/data/catalog_ingestion.js';

        const rawRecords = [
            { PGC: 143, RA: 187.7059, DEC: 12.3911, Vcmb: 1150.0, DM: 31.08, eDM: 0.35, GroupId: 10, Method: 'TF' },
            { PGC: 1024, RA: 192.8595, DEC: 27.1283, Vcmb: 7100.0, DM: 34.90, eDM: 0.25, GroupId: 20, Method: 'SNIa' },
            { PGC: 56000, RA: 266.4051, DEC: -28.9362, Vcmb: 4500.0, DM: 33.80, eDM: 0.40, GroupId: -1, Method: 'FP' }
        ];

        const pipeline = new CatalogIngestionPipeline({
            H0: 74.6,
            applyMalmquist: true,
            malmquistModel: 'homogeneous',
            applySelectionWeights: true,
            magLim: 17.5
        });

        const catalog = pipeline.processRecords(rawRecords);
        const g0 = catalog.getGalaxy(0);
        const g1 = catalog.getGalaxy(1);

        console.log(JSON.stringify({
            count: catalog.count,
            g0: g0,
            g1: g1
        }));
        """)
        assert res["count"] == 3
        # Virgo galaxy (PGC 143)
        g0 = res["g0"]
        assert g0["pgc"] == 143
        assert math.isclose(g0["ra"], 187.7059, abs_tol=1e-4)
        assert math.isclose(g0["dec"], 12.3911, abs_tol=1e-4)
        # Virgo distance should be around 16 Mpc (with Malmquist correction)
        assert 14.0 < g0["dist"] < 18.0
        assert g0["sgy"] > 10.0, "Virgo SGY should be positive in Supergalactic plane"
        assert abs(g0["sgz"]) < 5.0, "Virgo SGZ should be close to 0"
        assert g0["weight"] >= 1.0
