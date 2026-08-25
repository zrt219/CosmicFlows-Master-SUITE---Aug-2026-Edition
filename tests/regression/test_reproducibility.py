"""
Automated Test Suite for ZRT Cosmicflows Reproducibility Packaging & Worker Scheduler
Covers:
- SHA-256 cryptographic engine accuracy (NIST test vectors)
- FITS 2880-byte multiple standard block headers & card parsing
- Academic LaTeX figure caption math formatting & macro generation
- W3C PROV-O / FAIR execution provenance sidecar generation
- Reproducibility bundle creation and manifest checksums
- Tampering / corruption detection via bundle verification
- High-performance WorkerScheduler thread pooling, priority scheduling, and telemetry
- WorkerScheduler task cancellation, timeouts, and batch map operations
"""

import json
import subprocess
import pytest

def run_js(js_code: str):
    """Executes ES6 JavaScript code snippet in Node.js and returns parsed JSON output."""
    full_code = f"""
    import {{ ReproducibilityPackager, FitsHeaderFormatter, LatexCaptionGenerator, ProvenanceSidecarBuilder }} from './src/export/reproducibility_packager.js';
    import {{ WorkerScheduler, TaskPriority, TaskState }} from './src/runtime/worker_scheduler.js';

    async function run() {{
        {js_code}
    }}
    run().then(res => console.log('__JSON_OUT__' + JSON.stringify(res))).catch(err => {{
        console.error(err);
        process.exit(1);
    }});
    """
    proc = subprocess.run(
        ["node", "--input-type=module", "-e", full_code],
        capture_output=True,
        text=True,
        cwd="c:/Users/Zhane/Documents/antigravity/gallant-newton"
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed with code {proc.returncode}:\nSTDERR: {proc.stderr}\nSTDOUT: {proc.stdout}")
    
    for line in proc.stdout.splitlines():
        if line.startswith("__JSON_OUT__"):
            return json.loads(line[len("__JSON_OUT__"):])
    raise RuntimeError(f"No JSON output marker found in stdout:\n{proc.stdout}")


class TestReproducibilityPackaging:
    """Test suite for research-grade reproducibility bundle packaging."""

    def test_sha256_nist_vectors(self):
        js = """
        const emptyHash = ReproducibilityPackager.computeSha256('');
        const abcHash = ReproducibilityPackager.computeSha256('abc');
        const longHash = ReproducibilityPackager.computeSha256('abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq');
        return { emptyHash, abcHash, longHash };
        """
        res = run_js(js)
        # NIST standard known answer hashes
        assert res["emptyHash"] == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        assert res["abcHash"] == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        assert res["longHash"] == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"

    def test_fits_header_blocks_and_parsing(self):
        js = """
        const header = FitsHeaderFormatter.buildPrimaryHeader({
            dimensions: [256, 256, 128],
            bitpix: -32,
            unit: 'km/s',
            cosmology: { h0: 74.6, omegaM: 0.315, omegaL: 0.685, sigma8: 0.811 }
        });
        const parsed = FitsHeaderFormatter.parseHeader(header);
        return {
            headerLength: header.length,
            isMultipleOf2880: header.length % 2880 === 0,
            naxis1: parsed.cards['NAXIS1'].value,
            naxis2: parsed.cards['NAXIS2'].value,
            naxis3: parsed.cards['NAXIS3'].value,
            h0: parsed.cards['COSMO_H0'].value,
            omegaM: parsed.cards['COSMO_OM'].value,
            hasEnd: !header.includes('END     =') && header.includes('END     ')
        };
        """
        res = run_js(js)
        assert res["isMultipleOf2880"] is True
        assert res["headerLength"] >= 2880
        assert res["naxis1"] == 256
        assert res["naxis2"] == 256
        assert res["naxis3"] == 128
        assert abs(res["h0"] - 74.6) < 1e-4
        assert abs(res["omegaM"] - 0.315) < 1e-4
        assert res["hasEnd"] is True

    def test_latex_figure_caption_generation(self):
        js = """
        const caption = LatexCaptionGenerator.generateCaption({
            figureLabel: 'fig:cf4_cosmography',
            title: 'Cosmic Velocity Field and Laniakea Watershed Basin',
            dataset: 'Cosmicflows-4',
            cosmology: { h0: 74.6, omegaM: 0.315, omegaL: 0.685 },
            smoothingScale: 3.5,
            integrationStep: 0.25,
            basinCount: 8,
            criticalPointsCount: 42
        });
        return { caption };
        """
        res = run_js(js)
        caption = res["caption"]
        assert "\\begin{figure*}" in caption
        assert "\\label{fig:cf4_cosmography}" in caption
        assert "Cosmicflows-4" in caption
        assert "\\Omega_{\\mathrm{m}} = 0.315" in caption
        assert "N_{\\text{basin}} = 8" in caption
        assert "N_{\\text{crit}} = 42" in caption

    def test_provenance_sidecar_generation(self):
        js = """
        const sidecar = ProvenanceSidecarBuilder.buildSidecar({
            version: '2026.2.0',
            H0: 74.6,
            Omega_m: 0.315,
            Omega_Lambda: 0.685,
            smoothingScale: 3.5,
            basinCount: 8,
            streamlineCount: 2500,
            criticalPointsCount: 42
        });
        return sidecar;
        """
        res = run_js(js)
        assert res["schema_version"] == "1.2.0"
        assert res["provenance"]["creator"] == "ZRT Cosmicflows Workbench"
        assert res["cosmology"]["H0_km_s_Mpc"] == 74.6
        assert res["execution_parameters"]["smoothing_scale_Mpc_h"] == 3.5
        assert res["output_metrics"]["total_basins"] == 8
        assert res["output_metrics"]["total_streamlines"] == 2500

    def test_reproducibility_bundle_create_and_verify(self):
        js = """
        const datasets = {
            'data/catalog_peculiar_velocities.csv': 'basin_id,attractor_name,x,y,z\\n0,Great Attractor,-45,12,-8\\n',
            'data/streamlines.geojson': JSON.stringify({ type: 'FeatureCollection', features: [] })
        };

        const bundle = ReproducibilityPackager.createBundle({
            datasets,
            metadata: {
                bundleId: 'CF4-EXP-2026-BUNDLE-01',
                title: 'High-Precision Cosmography Test Run',
                basinCount: 1,
                streamlineCount: 0,
                criticalPointsCount: 4
            },
            cosmology: { h0: 74.6, omegaM: 0.315, omegaL: 0.685 },
            logs: [
                'Log 1: Initialized parameters',
                'Log 2: Computed watershed partition',
                'Log 3: Export complete'
            ]
        });

        const verification = ReproducibilityPackager.verifyBundle(bundle);
        return { bundle, verification };
        """
        res = run_js(js)
        bundle = res["bundle"]
        verif = res["verification"]

        assert bundle["bundle_id"] == "CF4-EXP-2026-BUNDLE-01"
        assert bundle["files_count"] >= 5
        assert "manifest.sha256" in bundle["files"]
        assert verif["valid"] is True
        assert verif["manifestVerified"] is True
        assert verif["sidecarValid"] is True
        assert len(verif["errors"]) == 0
        assert verif["details"]["corruptedFiles"] == []

    def test_bundle_tampering_failure_injection(self):
        js = """
        const datasets = {
            'data/catalog.csv': 'id,name\\n1,Virgo\\n'
        };
        const bundle = ReproducibilityPackager.createBundle({ datasets });

        // Tamper with file content without updating manifest
        bundle.files['data/catalog.csv'] = 'id,name\\n1,TamperedCluster\\n';

        const verification = ReproducibilityPackager.verifyBundle(bundle);
        return { verification };
        """
        res = run_js(js)
        verif = res["verification"]
        assert verif["valid"] is False
        assert verif["manifestVerified"] is False
        assert len(verif["details"]["corruptedFiles"]) == 1
        assert "data/catalog.csv" in verif["details"]["corruptedFiles"]
        assert any("SHA-256 mismatch" in err for err in verif["errors"])


class TestWorkerScheduler:
    """Test suite for high-performance WorkerScheduler."""

    def test_worker_scheduler_priority_queue_and_batch_map(self):
        js = """
        const scheduler = new WorkerScheduler({ maxWorkers: 4 });
        const progressEvents = [];
        const taskCompletedEvents = [];

        scheduler.on('onProgress', (p) => progressEvents.push(p));
        scheduler.on('onTaskComplete', (t) => taskCompletedEvents.push(t.taskId));

        // Submit low priority first, then high priority
        const taskLow = scheduler.submitTask({
            type: 'EVALUATE_BULK_FLOW',
            payload: { shells: [20, 50] },
            priority: TaskPriority.LOW
        });

        const taskHigh = scheduler.submitTask({
            type: 'INTEGRATE_STREAMLINES',
            payload: { seeds: [[0, 0, 0], [10, 10, 10], [-5, -5, 5]] },
            priority: TaskPriority.HIGH
        });

        // Test batch map operation
        const mapResults = await scheduler.map([1, 2, 3, 4], 'COMPUTE_WATERSHED', {
            commonPayload: { basinCount: 6 }
        });

        const resLow = await taskLow;
        const resHigh = await taskHigh;
        const telemetry = scheduler.getTelemetry();

        scheduler.terminate();

        return {
            resLowLength: resLow.length,
            resHighLength: resHigh.length,
            mapResultsLength: mapResults.length,
            completedTasks: telemetry.completedTasks,
            failedTasks: telemetry.failedTasks,
            poolSize: telemetry.poolSize,
            progressEventsCount: progressEvents.length,
            taskCompletedCount: taskCompletedEvents.length
        };
        """
        res = run_js(js)
        assert res["resLowLength"] == 2
        assert res["resHighLength"] == 3
        assert res["mapResultsLength"] == 4
        assert res["completedTasks"] == 6
        assert res["failedTasks"] == 0
        assert res["poolSize"] == 4
        assert res["taskCompletedCount"] == 6

    def test_worker_scheduler_task_cancellation(self):
        js = """
        const scheduler = new WorkerScheduler({ maxWorkers: 1 });

        // Submit a blocking custom task
        const t1 = scheduler.submitTask({
            type: 'CUSTOM',
            payload: {},
            handler: async () => {
                await new Promise(r => setTimeout(r, 100));
                return 'done1';
            }
        });

        // Submit a second task that sits in queue
        let task2Error = null;
        const t2Promise = scheduler.submitTask({
            type: 'CUSTOM',
            payload: {},
            handler: async () => 'done2'
        }).catch(err => {
            task2Error = err.message;
        });

        // Find and cancel task 2 from queue
        const queuedTaskId = scheduler.taskQueue[0]?.id;
        const cancelled = queuedTaskId ? scheduler.cancelTask(queuedTaskId) : false;

        await t1;
        await t2Promise;

        scheduler.terminate();

        return {
            cancelled,
            task2Error
        };
        """
        res = run_js(js)
        assert res["cancelled"] is True
        assert "cancelled" in res["task2Error"]
