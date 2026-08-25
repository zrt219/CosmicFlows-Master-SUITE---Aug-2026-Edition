# -*- coding: utf-8 -*-
"""
Automated Pytest Suite for W3C PROV-JSONLD Exporter and Lineage Manifests.
"""

import pytest
import json
from tests.utils import run_node_snippet

class TestProvenanceJsonLdExporter:
    def test_exporter_instantiation(self):
        res = run_node_snippet("""
        import { ProvenanceJsonLdExporter } from './src/export/provenance_jsonld_exporter.js';
        const exporter = new ProvenanceJsonLdExporter();
        console.log(JSON.stringify({
            hasExporter: !!exporter,
            version: exporter.version,
            hasContext: !!exporter.defaultContext['prov']
        }));
        """)
        assert res["hasExporter"] is True
        assert res["version"] == "1.0.0"
        assert res["hasContext"] is True

    def test_jsonld_bundle_serialization(self):
        res = run_node_snippet("""
        import { LineageGraph } from './src/provenance/lineage_graph.js';
        import { ProvenanceJsonLdExporter } from './src/export/provenance_jsonld_exporter.js';
        
        const g = new LineageGraph();
        g.addEntity('cf4_velocity_grid', {
            sha256: '9f86d081884c7d659a2feaa0c55ad015a3bf4f1b2b0b822cd15d6c15b0f00a08',
            dimensions: '128x128x128',
            units: 'km/s'
        });
        g.recordActivity('watershed_segmentation', {
            algorithm: 'Cash-Karp-RK45-Watershed',
            velocityScale: 52.0,
            smoothingMpc: 5.0
        });

        const exporter = new ProvenanceJsonLdExporter();
        const doc = exporter.exportJsonLd(g, { title: 'CF4 Basin Lineage', author: 'ZRT Tester' });
        console.log(JSON.stringify({
            contextProv: doc['@context']['prov'],
            contextZrt: doc['@context']['zrt'],
            title: doc['dc:title'],
            entityCount: doc.entities.length,
            activityCount: doc.activities.length,
            citationCount: doc.citations.length
        }));
        """)
        assert res["contextProv"] == "http://www.w3.org/ns/prov#"
        assert res["contextZrt"] == "https://zrt.science/ontology/cosmicflows#"
        assert res["entityCount"] == 1
        assert res["activityCount"] == 1
        assert res["citationCount"] == 2

    def test_latex_figure_caption_generation(self):
        res = run_node_snippet("""
        import { ProvenanceJsonLdExporter } from './src/export/provenance_jsonld_exporter.js';
        const exporter = new ProvenanceJsonLdExporter();
        const caption = exporter.generateLatexFigureCaption({
            algorithm: 'Cash-Karp RK45',
            gridResolution: '128^3',
            velocityScale: '52.0'
        });
        console.log(JSON.stringify({
            hasCaption: caption.includes('\\caption{'),
            hasCitation1: caption.includes('Courtois2023_CF4'),
            hasCitation2: caption.includes('Dupuy2023_Watershed')
        }));
        """)
        assert res["hasCaption"] is True
        assert res["hasCitation1"] is True
        assert res["hasCitation2"] is True

    def test_bibtex_generation(self):
        res = run_node_snippet("""
        import { ProvenanceJsonLdExporter } from './src/export/provenance_jsonld_exporter.js';
        const exporter = new ProvenanceJsonLdExporter();
        const bib = exporter.generateBibTeX();
        console.log(JSON.stringify({
            hasCF4: bib.includes('@article{Courtois2023_CF4'),
            hasWatershed: bib.includes('@article{Dupuy2023_Watershed')
        }));
        """)
        assert res["hasCF4"] is True
        assert res["hasWatershed"] is True
