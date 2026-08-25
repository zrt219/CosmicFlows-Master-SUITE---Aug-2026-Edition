"""
Automated Test Suite for ZRT Cosmicflows Data Exporters and Acceptance Gates
Covers:
- Gate A (Data Integrity Gate: coordinates, speed of light, NaN/Inf, grid buffers)
- Gate B (Units & Dimensional Consistency Gate: H0, Omega_m, scale factor, redshift)
- Gate C (Field Divergence & Curl Gate: potential flow, vorticity fraction)
- Gate D (Integration Accuracy Gate: RK45 truncation error, reversibility)
- Gate E (Watershed Segmentation Consistency Gate: basin volumes, conservation)
- Gate F (Critical Point Roots & Morse Index Gate: residuals, Hessian nonsingularity, Euler sum)
- Gate G (Eigenvalue & Tensor Orthonormality Gate: Jacobi diagonalization, orthonormality)
- Gate H (Uncertainty & Covariance Gate: positive semi-definiteness, PSD, correlation bounds)
- Gate I (Export & Serialization Gate: FITS blocks, manifest hashes, GeoJSON)
- AcceptanceGatekeeper (full sequential pipeline validation and failure rejection)
- Basin Catalog Exporter (CSV and Astropy ECSV format, roundtrip parsing)
- Streamline Exporter (GeoJSON FeatureCollection, ZRT-BIN compact binary serialization/deserialization)
- Topology Graph Exporter (GraphML XML specification, node-link JSON)
- Bulk Flow Exporter (Multi-shell CSV, IVOA VOTable 1.4 XML)
"""

import json
import subprocess
import pytest
import xml.etree.ElementTree as ET
import csv
import io

def run_js(js_code: str):
    """Executes ES6 JavaScript code snippet in Node.js and returns parsed JSON output."""
    full_code = f"""
    import {{ BasinCatalogExporter, StreamlineExporter, TopologyGraphExporter, BulkFlowExporter }} from './src/export/data_exporters.js';
    import {{ AcceptanceGatekeeper, GateStatus, GateA_DataIntegrity, GateB_UnitsConsistency, GateC_FieldDivergenceAndCurl, GateD_IntegrationAccuracy, GateE_WatershedSegmentation, GateF_CriticalPointRoots, GateG_EigenvalueOrthonormality, GateH_UncertaintyCovariance, GateI_ExportSerialization }} from './src/validation/acceptance_gates.js';

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


class TestValidationGates:
    """Test suite for Validation Gates A through I and Unified Gatekeeper."""

    def test_gate_a_valid_data(self):
        js = """
        const payload = {
            records: [
                { id: 1, position: [10.5, -20.2, 30.1], vMag: 450, distance: 35.0 },
                { id: 2, position: [0.0, 15.0, -10.0], vMag: 1200, distance: 18.0 }
            ],
            grid: {
                dimensions: [2, 2, 2],
                data: new Float32Array([1, 2, 3, 4, 5, 6, 7, 8])
            }
        };
        return GateA_DataIntegrity.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["status"] == "PASSED"
        assert len(res["errors"]) == 0
        assert res["metrics"]["recordCount"] == 2
        assert res["metrics"]["speedOfLightViolations"] == 0

    def test_gate_a_failure_injection_nan_and_c_violation(self):
        js = """
        const payload = {
            records: [
                { id: 1, position: [NaN, 0, 0], vMag: 400000, distance: -5.0 }
            ]
        };
        return GateA_DataIntegrity.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is False
        assert res["status"] == "FAILED"
        assert len(res["errors"]) >= 2
        assert any("NaN" in e for e in res["errors"])
        assert any("speed of light" in e for e in res["errors"])

    def test_gate_b_units_and_cosmology(self):
        js = """
        const payload = {
            cosmology: { h0: 74.6, omegaM: 0.315, omegaL: 0.685 },
            scaleFactor: 1.0,
            redshift: 0.0
        };
        return GateB_UnitsConsistency.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["metrics"]["h0"] == 74.6

    def test_gate_b_failure_unphysical_params(self):
        js = """
        const payload = {
            cosmology: { h0: 15.0, omegaM: -0.5, omegaL: 5.0 },
            scaleFactor: -0.2,
            redshift: -2.0
        };
        return GateB_UnitsConsistency.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is False
        assert len(res["errors"]) >= 4

    def test_gate_c_field_curl_and_divergence(self):
        js = """
        const validPayload = { divNorm: 10.0, curlNorm: 0.5, maxVorticityRatio: 0.2 };
        const invalidPayload = { divNorm: 2.0, curlNorm: 8.0, maxVorticityRatio: 0.2 };
        return {
            valid: GateC_FieldDivergenceAndCurl.evaluate(validPayload),
            invalid: GateC_FieldDivergenceAndCurl.evaluate(invalidPayload)
        };
        """
        res = run_js(js)
        assert res["valid"]["passed"] is True
        assert res["invalid"]["passed"] is False
        assert "vorticity fraction" in res["invalid"]["errors"][0]

    def test_gate_d_integration_accuracy(self):
        js = """
        const payload = {
            truncationError: 1e-6,
            tolerance: 1e-4,
            reversibilityResidual: 1e-5,
            reversibilityTol: 1e-3,
            stepSizes: [0.25, 0.25, 0.20, 0.15]
        };
        return GateD_IntegrationAccuracy.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is True

    def test_gate_e_watershed_volume_conservation(self):
        js = """
        const payload = {
            basins: [
                { id: 0, volume: 400.0 },
                { id: 1, volume: 600.0 }
            ],
            totalVolume: 1000.0,
            strictVolumeConservation: true
        };
        return GateE_WatershedSegmentation.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["metrics"]["basinCount"] == 2
        assert res["metrics"]["volumeConservationRelDiff"] == 0.0

    def test_gate_f_critical_point_roots_and_morse_census(self):
        js = """
        const payload = {
            criticalPoints: [
                { id: 0, name: 'Sink', morseIndex: 0, residual: 1e-5, hessian: [[2,0,0],[0,3,0],[0,0,1]] },
                { id: 1, name: 'Saddle1', morseIndex: 1, residual: 2e-5, hessian: [[-1,0,0],[0,2,0],[0,0,2]] }
            ],
            rootTolerance: 1e-3
        };
        return GateF_CriticalPointRoots.evaluate(payload);
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["metrics"]["criticalPointsCount"] == 2
        assert res["metrics"]["morseCensus"]["0"] == 1
        assert res["metrics"]["morseCensus"]["1"] == 1

    def test_gate_g_eigenvalue_and_tensor_orthonormality(self):
        js = """
        const symmetricTensor = [
            [4.0, 1.0, -2.0],
            [1.0, 5.0, 0.0],
            [-2.0, 0.0, 3.0]
        ];
        return GateG_EigenvalueOrthonormality.evaluate({ tensor: symmetricTensor });
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["metrics"]["isSymmetric"] is True
        assert res["metrics"]["diagonalizationConverged"] is True
        assert len(res["metrics"]["eigenvalues"]) == 3
        # Check ascending order
        eigs = res["metrics"]["eigenvalues"]
        assert eigs[0] <= eigs[1] <= eigs[2]

    def test_gate_h_covariance_matrix_psd(self):
        js = """
        const validCov = [
            [25.0, 5.0, 0.0],
            [5.0, 16.0, 2.0],
            [0.0, 2.0, 9.0]
        ];
        const invalidCov = [
            [1.0, 10.0, 0.0],
            [10.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ];
        return {
            valid: GateH_UncertaintyCovariance.evaluate({ covarianceMatrix: validCov }),
            invalid: GateH_UncertaintyCovariance.evaluate({ covarianceMatrix: invalidCov })
        };
        """
        res = run_js(js)
        assert res["valid"]["passed"] is True
        assert res["invalid"]["passed"] is False

    def test_gate_i_export_and_fits_blocks(self):
        js = """
        const fitsHeader = ('SIMPLE  =                    T / Conforms                                       ' +
                           'END                                                                             ').padEnd(2880, ' ');
        const manifest = {
            'data/catalog.csv': 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'
        };
        const geoJson = {
            type: 'FeatureCollection',
            features: []
        };
        return GateI_ExportSerialization.evaluate({ fitsHeader, manifest, geoJson });
        """
        res = run_js(js)
        assert res["passed"] is True
        assert res["metrics"]["fitsHeaderBytes"] == 2880

    def test_gatekeeper_full_audit(self):
        js = """
        const context = {
            gateA: { records: [{ id: 1, position: [0,0,0], vMag: 100 }] },
            gateB: { cosmology: { h0: 74.6, omegaM: 0.315, omegaL: 0.685 } },
            gateC: { divNorm: 1.0, curlNorm: 0.01 },
            gateD: { truncationError: 1e-6 },
            gateE: { basins: [{ id: 0, volume: 100 }] },
            gateF: { criticalPoints: [{ id: 0, morseIndex: 0, residual: 1e-6 }] },
            gateG: { tensor: [[2,0,0],[0,2,0],[0,0,2]] },
            gateH: { covarianceMatrix: [[4,0,0],[0,4,0],[0,0,4]] },
            gateI: { manifest: { 'test.csv': 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855' } }
        };
        return AcceptanceGatekeeper.evaluateAllGates(context);
        """
        res = run_js(js)
        assert res["status"] == "VALIDATED"
        assert res["allPassed"] is True
        assert res["passedCount"] == 9
        assert res["failedCount"] == 0


class TestDataExporters:
    """Test suite for research-grade data format converters."""

    def test_basin_catalog_csv_and_ecsv_roundtrip(self):
        js = """
        const sampleBasins = [
            {
                id: 0,
                name: 'Laniakea Great Attractor',
                position: [-45.2, 12.8, -8.4],
                velocity: [-350.0, 180.0, -120.0],
                volume: 125000.5,
                massProxy: 4.8,
                galaxyCount: 8420,
                percolationThreshold: -1.45,
                stabilityIndex: 0.985,
                boundingBox: [-120, 40, -60, 80, -70, 50]
            },
            {
                id: 1,
                name: 'Perseus-Pisces Supercluster',
                position: [55.0, -22.1, -15.3],
                velocity: [280.0, -90.0, 40.0],
                volume: 98000.2,
                massProxy: 3.9,
                galaxyCount: 6150,
                percolationThreshold: -1.20,
                stabilityIndex: 0.940,
                boundingBox: [10, 110, -70, 30, -60, 30]
            }
        ];
        const plainCsv = BasinCatalogExporter.exportToCsv(sampleBasins, { includeEcsvHeader: false });
        const ecsv = BasinCatalogExporter.exportToCsv(sampleBasins, { includeEcsvHeader: true });
        const parsed = BasinCatalogExporter.parseCsv(plainCsv);
        return { plainCsv, ecsv, parsed };
        """
        res = run_js(js)
        assert "Laniakea Great Attractor" in res["plainCsv"]
        assert "%ECSV 1.0" in res["ecsv"]
        assert len(res["parsed"]) == 2
        assert res["parsed"][0]["name"] == "Laniakea Great Attractor"
        assert abs(res["parsed"][0]["position"][0] - (-45.2)) < 1e-4

        # Validate CSV using python csv reader
        reader = csv.reader(io.StringIO(res["plainCsv"]))
        rows = list(reader)
        assert len(rows) == 3  # Header + 2 data rows
        assert rows[0][1] == "attractor_name"

    def test_streamlines_geojson_export(self):
        js = """
        const streamlines = [
            {
                id: 101,
                basinId: 0,
                attractorName: 'Great Attractor',
                points: [
                    [0.0, 0.0, 0.0, 50.0, 10.0, -20.0],
                    [5.0, 2.0, -1.0, 80.0, 15.0, -30.0],
                    [10.0, 4.0, -2.0, 120.0, 20.0, -40.0]
                ]
            }
        ];
        return StreamlineExporter.exportToGeoJson(streamlines);
        """
        res = run_js(js)
        assert res["type"] == "FeatureCollection"
        assert len(res["features"]) == 1
        feature = res["features"][0]
        assert feature["geometry"]["type"] == "LineString"
        assert len(feature["geometry"]["coordinates"]) == 3
        assert feature["properties"]["attractor_name"] == "Great Attractor"
        assert feature["properties"]["step_count"] == 3
        assert feature["properties"]["arc_length_mpc_h"] > 0

    def test_streamlines_compact_binary_roundtrip(self):
        js = """
        const streamlines = [
            {
                basinId: 3,
                arcLength: 42.5,
                points: [
                    [1.0, 2.0, 3.0, 100.0, -50.0, 25.0, 0.0],
                    [1.5, 2.5, 3.5, 110.0, -45.0, 30.0, 0.25],
                    [2.0, 3.0, 4.0, 120.0, -40.0, 35.0, 0.50]
                ]
            },
            {
                basinId: 5,
                arcLength: 18.2,
                points: [
                    [-10.0, -20.0, -30.0, -200.0, 150.0, -80.0, 0.0],
                    [-11.0, -21.0, -31.0, -210.0, 155.0, -85.0, 0.25]
                ]
            }
        ];

        const buffer = StreamlineExporter.exportToCompactBinary(streamlines);
        const decoded = StreamlineExporter.importFromCompactBinary(buffer);
        return {
            byteLength: buffer.byteLength,
            version: decoded.version,
            totalPoints: decoded.totalPoints,
            bounds: decoded.bounds,
            streamlinesCount: decoded.streamlines.length,
            firstLinePointsCount: decoded.streamlines[0].points.length,
            firstLineFirstPoint: decoded.streamlines[0].points[0],
            secondLineBasinId: decoded.streamlines[1].basinId
        };
        """
        res = run_js(js)
        assert res["version"] == 1
        assert res["totalPoints"] == 5
        assert res["streamlinesCount"] == 2
        assert res["firstLinePointsCount"] == 3
        assert res["secondLineBasinId"] == 5
        # Verify first vertex values
        assert abs(res["firstLineFirstPoint"][0] - 1.0) < 1e-4
        assert abs(res["firstLineFirstPoint"][3] - 100.0) < 1e-4

    def test_topology_graphml_xml_validity(self):
        js = """
        const graphData = {
            criticalPoints: [
                { id: 'node_sink_1', name: 'Shapley Core Sink', morseIndex: 0, position: [-130.0, 30.0, -20.0], potential: -54000.0, eigenvalues: [1.2, 2.4, 3.1] },
                { id: 'node_saddle_1', name: 'Hydra-Centaurus 1-Saddle', morseIndex: 1, position: [-50.0, 10.0, -5.0], potential: -32000.0, eigenvalues: [-0.8, 1.5, 2.2] }
            ],
            separatrices: [
                { id: 'edge_1', source: 'node_saddle_1', target: 'node_sink_1', type: 'stable_manifold', weight: 4.5, length: 82.5, flux: 120.0 }
            ]
        };
        return {
            graphml: TopologyGraphExporter.exportToGraphML(graphData),
            nodeLink: TopologyGraphExporter.exportToNodeLink(graphData)
        };
        """
        res = run_js(js)
        graphml_str = res["graphml"]
        assert "<?xml version" in graphml_str
        assert "<graphml" in graphml_str
        assert "Shapley Core Sink" in graphml_str

        # Parse XML with python ElementTree to verify strict XML compliance
        root = ET.fromstring(graphml_str)
        assert root.tag.endswith("graphml")

        # Test Node-Link JSON
        node_link = res["nodeLink"]
        assert len(node_link["nodes"]) == 2
        assert len(node_link["links"]) == 1
        assert node_link["links"][0]["source"] == "node_saddle_1"
        assert node_link["links"][0]["target"] == "node_sink_1"

    def test_bulk_flow_multishell_and_votable(self):
        js = """
        const shells = [
            { id: 0, rMin: 0, rMax: 50, rEff: 25.0, vBulk: [120.0, -180.0, 240.0], galacticL: 290.4, galacticB: 15.2, error: 18.5, sampleSize: 1200 },
            { id: 1, rMin: 50, rMax: 100, rEff: 75.0, vBulk: [140.0, -210.0, 270.0], galacticL: 288.1, galacticB: 12.8, error: 22.0, sampleSize: 3400 }
        ];
        return {
            csv: BulkFlowExporter.exportToCsv(shells),
            votable: BulkFlowExporter.exportToVOTable(shells)
        };
        """
        res = run_js(js)
        assert "shell_id" in res["csv"]
        assert "vx_bulk_kms" in res["csv"]
        assert "290.400" in res["csv"]

        votable_str = res["votable"]
        assert "<VOTABLE" in votable_str
        assert "<TABLE name=\"MultiShellBulkFlows\">" in votable_str
        root = ET.fromstring(votable_str)
        assert root.tag.endswith("VOTABLE")
