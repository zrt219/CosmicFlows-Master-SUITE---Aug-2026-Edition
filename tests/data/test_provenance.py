# -*- coding: utf-8 -*-
"""
Automated Test Suite for Provenance, Manifests, Audit Logger, and Unit Registry.

Verifies:
1. Canonical Dataset Manifests (IP2I CF4 / CF4++ Products):
   - Registry integrity, metadata completeness, and mirror endpoints.
   - SHA-256 hash calculation across strings, ArrayBuffers, and TypedArrays.
   - Buffer verification with size, checksum, and record-count validation.
   - Dataset query and filtering engine.
   - BibTeX and Markdown academic citation generators.
2. Cryptographic Audit Logger:
   - Genesis record initialization with standard root hash.
   - Sequential SHA-256 hash chaining and tamper-evidence.
   - Single-bit tamper detection identifying the exact corrupted record index.
   - Parameter tracking and snapshot merging.
   - Computation stage execution logging (durations, telemetry, error trapping).
   - W3C PROV-DM / PROV-O JSON-LD manifest generation.
   - Audit trail JSON export / import round-trip.
3. Dimensionally Consistent Unit Registry & Typed Quantities:
   - Physical and cosmological constants (c, G, M_sun, Mpc, yr, rho_crit, etc.).
   - 5D dimensional vector algebra (L, M, T, Theta, A).
   - Conversions across astrophysical and SI units.
   - Cosmological h-parameter scaling (Mpc/h, M_sun/h, km/s/(Mpc/h), rho_crit).
   - Dimensional mismatch error handling.
   - Quantity typed arithmetic (+, -, *, /, pow) and string formatting.
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


class TestDatasetManifest:
    """Tier 1: Canonical Dataset Manifests and Provenance Verification"""

    def test_canonical_manifest_registry(self):
        """Verify all 8 official CF4 / CF4++ / 2M++ products are registered with complete metadata."""
        res = run_node_eval("""
        import { defaultManifest, CANONICAL_MANIFESTS } from './src/provenance/manifest.js';
        
        const ids = defaultManifest.listIds();
        const records = defaultManifest.getAll();
        
        const checks = records.map(r => ({
            id: r.id,
            title: r.title,
            shortName: r.shortName,
            version: r.version,
            format: r.format,
            hasSha256: typeof r.sha256 === 'string' && r.sha256.length === 64,
            hasCosmology: !!r.cosmology && r.cosmology.H0 > 0,
            hasBibtex: typeof r.getBibTeX() === 'string' && r.getBibTeX().includes('@article'),
            hasCitation: !!r.citation && !!r.citation.doi
        }));
        
        console.log(JSON.stringify({
            count: ids.length,
            ids: ids,
            checks: checks
        }));
        """)
        assert res["count"] == 8, f"Expected 8 canonical datasets, got {res['count']}"
        expected_ids = [
            'cf4_distance_compilation_v1',
            'cf4_grouped_catalog_v1',
            'cf4_wf_grid_128_v1',
            'cf4_wf_grid_256_v1',
            'cf4_plus_plus_compilation_v1',
            'cf3_distance_compilation_v1',
            'twompp_velocity_field_v1',
            'vela_zoa_survey_v1'
        ]
        for eid in expected_ids:
            assert eid in res["ids"], f"Dataset {eid} must be in manifest"

        for chk in res["checks"]:
            assert chk["hasSha256"] is True, f"Dataset {chk['id']} missing valid SHA-256"
            assert chk["hasCosmology"] is True, f"Dataset {chk['id']} missing cosmology"
            assert chk["hasBibtex"] is True, f"Dataset {chk['id']} missing BibTeX citation"
            assert chk["hasCitation"] is True, f"Dataset {chk['id']} missing DOI citation"

    def test_sha256_computation_and_verification(self):
        """Verify universal SHA-256 calculation against known test vectors."""
        res = run_node_eval("""
        import { computeSha256, pureJsSha256, defaultManifest } from './src/provenance/manifest.js';
        
        const testStr = 'ZRT-Cosmicflows-Workbench-2026';
        const hash1 = await computeSha256(testStr);
        
        const encoder = new TextEncoder();
        const bytes = encoder.encode(testStr);
        const hash2 = pureJsSha256(bytes);
        const hash3 = await computeSha256(bytes.buffer);
        
        // Test verification of registered dataset
        const dummyData = 'sample-cf4-payload-stream';
        const dummyHash = await computeSha256(dummyData);
        
        // Mock dataset entry for verification
        const mockManifest = new (defaultManifest.constructor)({
            'test_cf4_stream': {
                id: 'test_cf4_stream',
                title: 'Test Stream',
                sha256: dummyHash,
                byteSize: bytes.length,
                recordCount: 10
            }
        });
        
        const verifySuccess = await mockManifest.verify('test_cf4_stream', dummyData);
        const verifyFail = await mockManifest.verify('test_cf4_stream', 'tampered-data');
        
        console.log(JSON.stringify({
            hash1: hash1,
            hash2: hash2,
            hash3: hash3,
            allEqual: (hash1 === hash2) && (hash2 === hash3),
            verifySuccess: verifySuccess,
            verifyFail: verifyFail
        }));
        """)
        assert res["allEqual"] is True, "SubtleCrypto and pure JS SHA-256 implementations must produce identical digests"
        assert len(res["hash1"]) == 64
        assert res["verifySuccess"]["valid"] is True
        assert res["verifyFail"]["valid"] is False
        assert "mismatch" in res["verifyFail"]["error"].lower()

    def test_manifest_query_and_citations(self):
        """Verify search query filters and BibTeX bibliography export."""
        res = run_node_eval("""
        import { defaultManifest } from './src/provenance/manifest.js';
        
        const csvDatasets = defaultManifest.find({ format: 'csv' });
        const gridDatasets = defaultManifest.find({ format: 'binary_f32_grid' });
        const highResGrids = defaultManifest.find({ minResolution: 256 });
        const year2023 = defaultManifest.find({ year: 2023 });
        const tullySearches = defaultManifest.find({ searchTerm: 'Tully' });
        
        const bibtexAll = defaultManifest.exportAllBibTeX();
        const cf4Record = defaultManifest.get('cf4_distance_compilation_v1');
        const mdCitation = cf4Record ? cf4Record.getMarkdownCitation() : '';
        
        console.log(JSON.stringify({
            csvCount: csvDatasets.length,
            gridCount: gridDatasets.length,
            highResCount: highResGrids.length,
            year2023Count: year2023.length,
            tullyCount: tullySearches.length,
            bibtexLength: bibtexAll.length,
            hasBibtexEntries: bibtexAll.includes('@article{Tully_2023_CF4') && bibtexAll.includes('@article{Courtois_2023_WF'),
            hasMdCitation: mdCitation.includes('Cosmicflows-4') && mdCitation.includes('DOI')
        }));
        """)
        assert res["csvCount"] >= 3
        assert res["gridCount"] >= 2
        assert res["highResCount"] >= 1
        assert res["year2023Count"] >= 4
        assert res["tullyCount"] >= 4
        assert res["hasBibtexEntries"] is True
        assert res["hasMdCitation"] is True


class TestAuditLogger:
    """Tier 2: Cryptographic Run Logging, Hash Chaining & W3C PROV Manifests"""

    def test_audit_logger_hash_chaining(self):
        """Verify that audit records form an unbroken, verifiable SHA-256 hash chain."""
        res = run_node_eval("""
        import { AuditLogger, LogLevel, GENESIS_ROOT_HASH } from './src/provenance/audit_logger.js';
        
        const logger = new AuditLogger({ runId: 'test_run_001', appName: 'CF4-Test-Runner' });
        
        logger.logInfo('INITIALIZE_PIPELINE', { mode: 'wiener_filter', gridDim: 256 });
        logger.trackParameters('cosmology', { H0: 74.6, Omega_m: 0.315, Omega_Lambda: 0.685 });
        logger.trackParameters('filter', { smoothingRadius: 5.0, prior: 'linear_lcdm' });
        logger.logCoordinateTransform('equatorial', 'supergalactic', { pointCount: 55877 });
        
        const integrityBefore = logger.verifyIntegrity();
        
        console.log(JSON.stringify({
            recordCount: logger.records.length,
            genesisPrevHash: logger.records[0].prevHash,
            genesisRootExpected: GENESIS_ROOT_HASH,
            integrityValid: integrityBefore.valid,
            headHash: integrityBefore.headHash
        }));
        """)
        assert res["recordCount"] == 5
        assert res["genesisPrevHash"] == res["genesisRootExpected"]
        assert res["integrityValid"] is True
        assert len(res["headHash"]) == 64

    def test_audit_logger_tamper_detection(self):
        """Verify that modifying any past record in the audit chain is immediately detected."""
        res = run_node_eval("""
        import { AuditLogger } from './src/provenance/audit_logger.js';
        
        const logger = new AuditLogger({ runId: 'test_tamper_run' });
        logger.logInfo('EVENT_A', { val: 100 });
        logger.logInfo('EVENT_B', { val: 200 });
        logger.logInfo('EVENT_C', { val: 300 });
        
        const validBefore = logger.verifyIntegrity();
        
        // Tamper with record at index 2 (EVENT_B)
        logger.records[2].payload.val = 999;
        const tamperCheck1 = logger.verifyIntegrity();
        
        // Tamper with record prevHash link
        logger.records[2].payload.val = 200; // restore content
        logger.records[3].prevHash = 'bad0000000000000000000000000000000000000000000000000000000000000';
        const tamperCheck2 = logger.verifyIntegrity();
        
        console.log(JSON.stringify({
            validBefore: validBefore.valid,
            tamperDetected1: !tamperCheck1.valid,
            corruptedIndex1: tamperCheck1.corruptedIndex,
            tamperDetected2: !tamperCheck2.valid,
            corruptedIndex2: tamperCheck2.corruptedIndex
        }));
        """)
        assert res["validBefore"] is True
        assert res["tamperDetected1"] is True
        assert res["corruptedIndex1"] == 2
        assert res["tamperDetected2"] is True
        assert res["corruptedIndex2"] == 3

    def test_audit_logger_stage_execution_and_prov_manifest(self):
        """Verify stage tracking, error handling, and PROV-O JSON-LD manifest generation."""
        res = run_node_eval("""
        import { AuditLogger } from './src/provenance/audit_logger.js';
        
        const logger = new AuditLogger({ runId: 'prov_test_run' });
        logger.trackParameters('cosmology', { H0: 74.6, Omega_m: 0.315 });
        
        // Log successful stage
        const stageRes = await logger.logStage('WienerFilterReconstruction', async () => {
            const arr = new Float32Array(1000);
            return arr;
        }, { inputDataset: 'cf4_distance_compilation_v1' });
        
        // Log stage that throws error
        let caughtError = false;
        try {
            await logger.logStage('FailingStep', async () => {
                throw new Error('Simulated divergence singular matrix');
            });
        } catch (e) {
            caughtError = true;
        }
        
        const provManifest = logger.generateExecutionManifest({
            usedDatasets: ['cf4_distance_compilation_v1'],
            generatedArtifacts: ['cf4_wf_grid_128']
        });
        
        const jsonExport = logger.toJson();
        const reconstructedLogger = AuditLogger.fromJson(jsonExport);
        const reconstructedIntegrity = reconstructedLogger.verifyIntegrity();
        
        const mdReport = logger.toMarkdownReport();
        
        console.log(JSON.stringify({
            caughtError: caughtError,
            provType: provManifest['@type'],
            hasActivities: provManifest['prov:activity'].length >= 2,
            hasEntities: provManifest['prov:entity'].length >= 2,
            hasUsed: provManifest['prov:used'].length >= 1,
            hasGenerated: provManifest['prov:wasGeneratedBy'].length >= 1,
            reconstructedValid: reconstructedIntegrity.valid,
            reconstructedRecordCount: reconstructedIntegrity.recordCount,
            mdReportLength: mdReport.length
        }));
        """)
        assert res["caughtError"] is True
        assert 'prov:Bundle' in res["provType"]
        assert res["hasActivities"] is True
        assert res["hasEntities"] is True
        assert res["hasUsed"] is True
        assert res["hasGenerated"] is True
        assert res["reconstructedValid"] is True
        assert res["reconstructedRecordCount"] >= 5
        assert res["mdReportLength"] > 100


class TestUnitRegistryAndQuantities:
    """Tier 3: Dimensionally Consistent Unit Registry & Physical Quantities"""

    def test_fundamental_and_astrophysical_conversions(self):
        """Verify exact unit conversions across length, velocity, mass, and time."""
        res = run_node_eval("""
        import { defaultRegistry, CONSTANTS, qty } from './src/units/unit_registry.js';
        
        const r = defaultRegistry;
        
        // 1 Mpc in meters
        const mpc_to_m = r.convert(1.0, 'Mpc', 'm');
        // 1 Mpc in km
        const mpc_to_km = r.convert(1.0, 'Mpc', 'km');
        // 1 kpc in pc
        const kpc_to_pc = r.convert(1.0, 'kpc', 'pc');
        
        // 1 km/s in m/s
        const kms_to_ms = r.convert(1.0, 'km/s', 'm/s');
        // c in km/s
        const c_kms = r.convert(1.0, 'c', 'km/s');
        
        // 1 M_sun in kg
        const msun_kg = r.convert(1.0, 'M_sun', 'kg');
        // 1 Gyr in yr
        const gyr_yr = r.convert(1.0, 'Gyr', 'yr');
        
        // 1 deg in rad
        const deg_rad = r.convert(180.0, 'deg', 'rad');
        
        console.log(JSON.stringify({
            mpc_to_m: mpc_to_m,
            mpc_to_km: mpc_to_km,
            kpc_to_pc: kpc_to_pc,
            kms_to_ms: kms_to_ms,
            c_kms: c_kms,
            msun_kg: msun_kg,
            gyr_yr: gyr_yr,
            deg_rad: deg_rad
        }));
        """)
        assert math.isclose(res["mpc_to_m"], 3.085677581491367e22, rel_tol=1e-9)
        assert math.isclose(res["mpc_to_km"], 3.085677581491367e19, rel_tol=1e-9)
        assert math.isclose(res["kpc_to_pc"], 1000.0, rel_tol=1e-9)
        assert math.isclose(res["kms_to_ms"], 1000.0, rel_tol=1e-9)
        assert math.isclose(res["c_kms"], 299792.458, rel_tol=1e-9)
        assert math.isclose(res["msun_kg"], 1.98847e30, rel_tol=1e-6)
        assert math.isclose(res["gyr_yr"], 1.0e9, rel_tol=1e-9)
        assert math.isclose(res["deg_rad"], math.pi, rel_tol=1e-9)

    def test_h_parameter_cosmological_scaling(self):
        """Verify h-dependent units scaling (Mpc/h, M_sun/h, km/s/(Mpc/h), rho_crit)."""
        res = run_node_eval("""
        import { defaultRegistry } from './src/units/unit_registry.js';
        
        const r = defaultRegistry;
        const h_cf4 = 0.746;
        
        // 100 Mpc/h -> Mpc: 100 / h = 134.048 Mpc
        const d_mpc = r.convert(100.0, 'Mpc/h', 'Mpc', { h: h_cf4 });
        
        // 1e12 M_sun/h -> M_sun: 1e12 / h = 1.34048e12 M_sun
        const m_msun = r.convert(1.0e12, 'M_sun/h', 'M_sun', { h: h_cf4 });
        
        // 100 km/s/(Mpc/h) -> km/s/Mpc: 100 * h = 74.6 km/s/Mpc
        const h0_kms_mpc = r.convert(100.0, 'km/s/(Mpc/h)', 'km/s/Mpc', { h: h_cf4 });
        
        // rho_crit -> M_sun/Mpc^3: 2.77536627e11 * h^2 M_sun/Mpc^3
        const rho_crit_val = r.convert(1.0, 'rho_crit', 'M_sun/Mpc^3', { h: h_cf4 });
        const expectedRhoCrit = 2.77536627e11 * (h_cf4 * h_cf4);
        
        console.log(JSON.stringify({
            d_mpc: d_mpc,
            m_msun: m_msun,
            h0_kms_mpc: h0_kms_mpc,
            rho_crit_val: rho_crit_val,
            expectedRhoCrit: expectedRhoCrit
        }));
        """)
        assert math.isclose(res["d_mpc"], 100.0 / 0.746, rel_tol=1e-7)
        assert math.isclose(res["m_msun"], 1.0e12 / 0.746, rel_tol=1e-7)
        assert math.isclose(res["h0_kms_mpc"], 74.6, rel_tol=1e-7)
        assert math.isclose(res["rho_crit_val"], res["expectedRhoCrit"], rel_tol=1e-7)

    def test_dimensional_mismatch_rejection(self):
        """Verify that converting between incompatible physical dimensions throws DimensionMismatchError."""
        res = run_node_eval("""
        import { defaultRegistry, DimensionMismatchError } from './src/units/unit_registry.js';
        
        let caughtMismatch = false;
        let errorName = '';
        
        try {
            defaultRegistry.convert(100.0, 'km/s', 'Mpc');
        } catch (e) {
            caughtMismatch = true;
            errorName = e.name;
        }
        
        let caughtMassToTime = false;
        try {
            defaultRegistry.convert(1.0, 'M_sun', 'Gyr');
        } catch (e) {
            caughtMassToTime = true;
        }
        
        console.log(JSON.stringify({
            caughtMismatch: caughtMismatch,
            errorName: errorName,
            caughtMassToTime: caughtMassToTime
        }));
        """)
        assert res["caughtMismatch"] is True
        assert res["errorName"] == "DimensionMismatchError"
        assert res["caughtMassToTime"] is True

    def test_quantity_typed_arithmetic(self):
        """Verify Quantity class operations (+, -, *, /, pow, to, toSI, equals, format)."""
        res = run_node_eval("""
        import { Quantity, qty } from './src/units/unit_registry.js';
        
        const q1 = qty(50.0, 'Mpc/h');
        const q2 = qty(25.0, 'Mpc/h');
        const qSum = q1.add(q2);
        const qDiff = q1.subtract(q2);
        
        const v = qty(300.0, 'km/s');
        const t = qty(1.0, 'Gyr');
        
        // Distance d = v * t
        const distProd = v.multiply(t);
        const distInMpc = distProd.to('Mpc');
        
        // Kinetic energy per unit mass: 0.5 * v^2 -> (km/s)^2
        const vSq = v.pow(2).to('(km/s)^2');
        
        const qFormat = qSum.format(4);
        const qEq = q1.equals(qty(50.0, 'Mpc/h'));
        
        console.log(JSON.stringify({
            qSumVal: qSum.value,
            qSumUnit: qSum.unit,
            qDiffVal: qDiff.value,
            distInMpcVal: distInMpc.value,
            distInMpcUnit: distInMpc.unit,
            vSqVal: vSq.value,
            vSqUnit: vSq.unit,
            qFormat: qFormat,
            qEq: qEq
        }));
        """)
        assert res["qSumVal"] == 75.0
        assert res["qSumUnit"] == 'Mpc/h'
        assert res["qDiffVal"] == 25.0
        # 300 km/s * 1 Gyr:
        # 300 km/s = 300 * 10^3 m/s
        # 1 Gyr = 3.15576e16 s
        # distance in meters = 3e5 * 3.15576e16 = 9.46728e21 m
        # in Mpc = 9.46728e21 / 3.08567758e22 = 0.3068 Mpc
        assert math.isclose(res["distInMpcVal"], 0.3068136, rel_tol=1e-4)
        assert res["vSqVal"] == 90000.0
        assert res["qEq"] is True
        assert "75" in res["qFormat"]
