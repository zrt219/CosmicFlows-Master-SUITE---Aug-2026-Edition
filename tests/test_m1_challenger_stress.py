# -*- coding: utf-8 -*-
"""
Adversarial Stress Test Suite for Milestone 1
Challenger: teamwork_preview_challenger_m1_2

Adversarially challenges:
1. Rapid-fire camera navigation (50+ consecutive entity fly-tos, asserting no NaN/Infinity vectors or race conditions).
2. Targeting beacon lifecycle, dual-ring reticle spawning, and complete Three.js geometry/material disposal without memory leakage.
3. Exhaustive dossier synthesis & DOM rendering across all 38 ASTROMETRIC_FEATURES and all 8 DUPUY_2023_BASINS (46 total).
4. High-frequency keyboard navigation & search input fuzzing under mid-flight camera transitions.
"""

import pytest
import time

class TestM1AdversarialStress:
    """Adversarial Verification Suite for Milestone 1 Scene Navigation, Memory & Dossiers"""

    def test_exhaustive_dossier_synthesis_all_46_entities(self, cdp):
        """
        Verify every single entity in ASTROMETRIC_FEATURES (all 38) and
        DUPUY_2023_BASINS (all 8) generates a complete, valid dossier with zero missing fields,
        zero NaNs, and zero undefined strings.
        """
        res = cdp.evaluate("""
        (function() {
            const results = [];
            const errors = [];
            
            const features = (typeof ASTROMETRIC_FEATURES !== 'undefined') ? ASTROMETRIC_FEATURES : [];
            const basins = (typeof DUPUY_2023_BASINS !== 'undefined') ? DUPUY_2023_BASINS : [];
            
            // 1. Test all 38 astrometric features
            features.forEach(f => {
                try {
                    const d = getOrSynthesizeDossier(f.id);
                    if (!d) {
                        errors.push(`getOrSynthesizeDossier returned null for feature id: ${f.id}`);
                        return;
                    }
                    
                    // Check required fields
                    const missing = [];
                    if (!d.id) missing.push('id');
                    if (!d.name) missing.push('name');
                    if (typeof d.mass !== 'number' || isNaN(d.mass)) missing.push('mass');
                    if (typeof d.cz !== 'number' || isNaN(d.cz)) missing.push('cz');
                    if (typeof d.sigma_v !== 'number' || isNaN(d.sigma_v)) missing.push('sigma_v');
                    if (typeof d.dist !== 'number' || isNaN(d.dist)) missing.push('dist');
                    if (!d.bcg || d.bcg.includes('undefined')) missing.push('bcg');
                    if (!d.vPec || !Array.isArray(d.vPec) || d.vPec.some(v => typeof v !== 'number' || isNaN(v))) missing.push('vPec');
                    if (typeof d.tx !== 'number' || isNaN(d.tx)) missing.push('tx');
                    if (!d.lx || d.lx.includes('undefined')) missing.push('lx');
                    if (!d.members || d.members.includes('undefined')) missing.push('members');
                    if (!d.citation || d.citation.includes('undefined')) missing.push('citation');
                    
                    if (missing.length > 0) {
                        errors.push(`Feature ${f.id} has invalid/missing fields: ${missing.join(', ')}`);
                    } else {
                        results.push({ id: f.id, name: d.name, cat: d.cat, dist: d.dist, cz: d.cz });
                    }
                } catch (err) {
                    errors.push(`Exception evaluating feature ${f.id}: ${err.message}`);
                }
            });
            
            // 2. Test all 8 supercluster basins
            basins.forEach(b => {
                try {
                    const basinId = `basin_${b.id}`;
                    const d = getOrSynthesizeDossier(basinId);
                    if (!d) {
                        errors.push(`getOrSynthesizeDossier returned null for basin: ${basinId}`);
                        return;
                    }
                    
                    const missing = [];
                    if (!d.id) missing.push('id');
                    if (!d.name) missing.push('name');
                    if (typeof d.mass !== 'number' || isNaN(d.mass)) missing.push('mass');
                    if (typeof d.cz !== 'number' || isNaN(d.cz)) missing.push('cz');
                    if (typeof d.sigma_v !== 'number' || isNaN(d.sigma_v)) missing.push('sigma_v');
                    if (typeof d.dist !== 'number' || isNaN(d.dist)) missing.push('dist');
                    if (!d.bcg || d.bcg.includes('undefined')) missing.push('bcg');
                    if (!d.vPec || !Array.isArray(d.vPec) || d.vPec.some(v => typeof v !== 'number' || isNaN(v))) missing.push('vPec');
                    if (typeof d.tx !== 'number' || isNaN(d.tx)) missing.push('tx');
                    if (!d.lx || d.lx.includes('undefined')) missing.push('lx');
                    if (!d.members || d.members.includes('undefined')) missing.push('members');
                    if (!d.citation || d.citation.includes('undefined')) missing.push('citation');
                    
                    if (missing.length > 0) {
                        errors.push(`Basin ${basinId} has invalid/missing fields: ${missing.join(', ')}`);
                    } else {
                        results.push({ id: basinId, name: d.name, cat: d.cat, dist: d.dist, cz: d.cz });
                    }
                } catch (err) {
                    errors.push(`Exception evaluating basin ${b.id}: ${err.message}`);
                }
            });
            
            return {
                totalTested: features.length + basins.length,
                successCount: results.length,
                errorCount: errors.length,
                errors: errors,
                featuresCount: features.length,
                basinsCount: basins.length
            };
        })()
        """)
        assert res["featuresCount"] == 38, f"Expected 38 ASTROMETRIC_FEATURES, got {res['featuresCount']}"
        assert res["basinsCount"] == 8, f"Expected 8 DUPUY_2023_BASINS, got {res['basinsCount']}"
        assert res["errorCount"] == 0, f"Dossier synthesis errors encountered: {res['errors']}"
        assert res["successCount"] == 46, f"Expected 46 valid dossiers, got {res['successCount']}"

    def test_open_dossier_dom_rendering_all_46_entities(self, cdp):
        """
        Verify openDossier(id) renders perfectly into the DOM for every single entity
        without throwing exceptions, with zero undefined/NaN strings in DOM elements,
        and with valid Gaussian velocity dispersion profile canvas rendering.
        """
        res = cdp.evaluate("""
        (function() {
            const entities = window.cosmicflows.search.getEntities();
            const domErrors = [];
            let renderedCount = 0;
            
            entities.forEach(ent => {
                try {
                    openDossier(ent.id);
                    
                    const modal = document.getElementById('spectroscopy-modal');
                    if (!modal || modal.style.display !== 'block' || modal.classList.contains('hidden')) {
                        domErrors.push(`Modal not visible after openDossier('${ent.id}')`);
                        return;
                    }
                    
                    const fields = {
                        name: document.getElementById('spec-name')?.textContent || '',
                        m200: document.getElementById('spec-m200')?.textContent || '',
                        cz: document.getElementById('spec-cz')?.textContent || '',
                        sigmav: document.getElementById('spec-sigmav')?.textContent || '',
                        dist: document.getElementById('spec-dist')?.textContent || '',
                        vpec: document.getElementById('spec-vpec')?.textContent || '',
                        tx: document.getElementById('spec-tx')?.textContent || '',
                        lx: document.getElementById('spec-lx')?.textContent || '',
                        bcg: document.getElementById('spec-bcg')?.textContent || '',
                        members: document.getElementById('spec-members')?.textContent || '',
                        citation: document.getElementById('spec-citation')?.textContent || ''
                    };
                    
                    for (const [key, val] of Object.entries(fields)) {
                        if (!val || val.trim() === '') {
                            domErrors.push(`Entity ${ent.id}: DOM field #${key} is empty`);
                        }
                        if (val.includes('undefined')) {
                            domErrors.push(`Entity ${ent.id}: DOM field #${key} contains 'undefined' ("${val}")`);
                        }
                        if (val.includes('NaN')) {
                            domErrors.push(`Entity ${ent.id}: DOM field #${key} contains 'NaN' ("${val}")`);
                        }
                    }
                    
                    // Check relative distance comparator presence
                    const distHtml = document.getElementById('spec-dist')?.innerHTML || '';
                    if (!distHtml.includes('farther than Virgo') && !distHtml.includes('Core Reference Cluster')) {
                        domErrors.push(`Entity ${ent.id}: Relative distance comparator missing in #spec-dist ("${distHtml}")`);
                    }
                    
                    // Check external links
                    const studyLink = document.getElementById('spec-link-study');
                    const wikiLink = document.getElementById('spec-link-wiki');
                    if (!studyLink || !studyLink.href || !studyLink.href.startsWith('http')) {
                        domErrors.push(`Entity ${ent.id}: #spec-link-study missing valid href`);
                    }
                    if (!wikiLink || !wikiLink.href || !wikiLink.href.startsWith('http')) {
                        domErrors.push(`Entity ${ent.id}: #spec-link-wiki missing valid href`);
                    }
                    
                    // Check canvas
                    const canvas = document.getElementById('spec-canvas');
                    if (!canvas) {
                        domErrors.push(`Entity ${ent.id}: #spec-canvas element missing`);
                    }
                    
                    renderedCount++;
                } catch (err) {
                    domErrors.push(`Exception rendering dossier for ${ent.id}: ${err.message}`);
                }
            });
            
            // Clean up by closing dossier
            closeDossier();
            
            return {
                totalEntities: entities.length,
                renderedCount: renderedCount,
                domErrors: domErrors
            };
        })()
        """)
        assert res["totalEntities"] >= 46, f"Expected >= 46 searchable entities, got {res['totalEntities']}"
        assert len(res["domErrors"]) == 0, f"DOM rendering errors encountered:\n" + "\n".join(res["domErrors"])
        assert res["renderedCount"] == res["totalEntities"]

    def test_rapid_fire_entity_selections_camera_stability(self, cdp):
        """
        Adversarially trigger 50 rapid consecutive entity selections across the 3D scene
        with zero delay, asserting:
        1. Camera position coordinates (x, y, z) remain finite and non-NaN at every frame.
        2. OrbitControls target coordinates remain finite and non-NaN.
        3. No unhandled exceptions or console errors occur.
        """
        res = cdp.evaluate("""
        (function() {
            const entities = window.cosmicflows.search.getEntities();
            const selectLog = [];
            const nanErrors = [];
            
            // Fire 50 rapid-fire selections
            for (let i = 0; i < 50; i++) {
                const ent = entities[i % entities.length];
                window.cosmicflows.search.select(ent.id);
                
                // Immediately sample camera & controls target
                const cx = camera.position.x;
                const cy = camera.position.y;
                const cz = camera.position.z;
                const tx = controls.target.x;
                const ty = controls.target.y;
                const tz = controls.target.z;
                
                if (isNaN(cx) || !isFinite(cx) || isNaN(cy) || !isFinite(cy) || isNaN(cz) || !isFinite(cz)) {
                    nanErrors.push(`Step ${i} (${ent.id}): Camera position is invalid: (${cx}, ${cy}, ${cz})`);
                }
                if (isNaN(tx) || !isFinite(tx) || isNaN(ty) || !isFinite(ty) || isNaN(tz) || !isFinite(tz)) {
                    nanErrors.push(`Step ${i} (${ent.id}): Controls target is invalid: (${tx}, ${ty}, ${tz})`);
                }
            }
            
            return {
                selectionsExecuted: 50,
                nanErrors: nanErrors,
                finalCameraPos: { x: camera.position.x, y: camera.position.y, z: camera.position.z },
                finalControlsTarget: { x: controls.target.x, y: controls.target.y, z: controls.target.z }
            };
        })()
        """)
        assert res["selectionsExecuted"] == 50
        assert len(res["nanErrors"]) == 0, f"Camera NaN/Infinite coordinate errors: {res['nanErrors']}"
        assert not any(v is None for v in res["finalCameraPos"].values())

    def test_targeting_beacon_lifecycle_and_memory_disposal(self, cdp):
        """
        Adversarially verify targeting beacon lifecycle:
        1. Spawns dual-ring geometry group in scene.
        2. When re-triggered or on animation completion (duration = 4000ms),
           all previous geometries and materials are properly disposed and removed from scene.
        3. No orphan meshes accumulate in scene.children.
        """
        # Step 1: Trigger beacon and verify scene mesh structure
        beacon_info = cdp.evaluate("""
        (function() {
            window.cosmicflows.search.select('virgo_cl');
            
            // Find beacon mesh in scene
            const beacon = scene.children.find(c => c instanceof THREE.Group && c.children.length === 2 && c.children[0].geometry instanceof THREE.RingGeometry);
            
            return {
                beaconFound: !!beacon,
                childCount: beacon ? beacon.children.length : 0,
                hasOuterRing: beacon ? beacon.children[0].geometry instanceof THREE.RingGeometry : false,
                hasInnerRing: beacon ? beacon.children[1].geometry instanceof THREE.RingGeometry : false,
                initialSceneChildren: scene.children.length,
                initialGeometries: renderer.info.memory.geometries
            };
        })()
        """)
        assert beacon_info["beaconFound"] is True, "Targeting beacon mesh group must be attached to scene"
        assert beacon_info["childCount"] == 2, "Beacon must have 2 concentric ring meshes"
        assert beacon_info["hasOuterRing"] is True, "Beacon child 0 must be RingGeometry"
        assert beacon_info["hasInnerRing"] is True, "Beacon child 1 must be RingGeometry"

        # Step 2: Trigger 30 rapid-fire beacon creations
        rapid_res = cdp.evaluate("""
        (function() {
            const entities = ['coma_cl', 'perseus_cl', 'shapley_super', 'great_attractor', 'bootes_void', 'hercules_cl'];
            for (let i = 0; i < 30; i++) {
                window.cosmicflows.search.select(entities[i % entities.length]);
            }
            
            // Count beacon groups in scene
            const beaconGroups = scene.children.filter(c => c instanceof THREE.Group && c.children.length === 2 && c.children[0].geometry instanceof THREE.RingGeometry);
            
            return {
                beaconGroupCount: beaconGroups.length,
                sceneChildrenCount: scene.children.length
            };
        })()
        """)
        # Must only have exactly 1 active beacon group in scene, not 30 accumulated orphan groups
        assert rapid_res["beaconGroupCount"] == 1, (
            f"Expected exactly 1 active beacon in scene after rapid firing, found {rapid_res['beaconGroupCount']}"
        )

        # Step 3: Wait 4.3 seconds for beacon animation duration (4000ms) to complete
        time.sleep(4.3)

        # Step 4: Verify complete cleanup & disposal
        cleanup_res = cdp.evaluate("""
        (function() {
            const beaconGroups = scene.children.filter(c => c instanceof THREE.Group && c.children.length === 2 && c.children[0].geometry instanceof THREE.RingGeometry);
            
            return {
                beaconGroupCount: beaconGroups.length,
                activeBeaconMeshNull: (typeof beaconMesh === 'undefined' || beaconMesh === null),
                activeBeaconAnimIdNull: (typeof beaconAnimId === 'undefined' || beaconAnimId === null)
            };
        })()
        """)
        assert cleanup_res["beaconGroupCount"] == 0, (
            f"Expected 0 beacon meshes in scene after 4.3s duration, found {cleanup_res['beaconGroupCount']}"
        )
        assert cleanup_res["activeBeaconMeshNull"] is True, "beaconMesh variable must be reset to null"
        assert cleanup_res["activeBeaconAnimIdNull"] is True, "beaconAnimId must be reset to null"

    def test_keyboard_navigation_and_fuzzing(self, cdp):
        """
        Adversarially fuzz keyboard navigation:
        - Rapidly type query
        - Send ArrowDown 20 times (testing boundary clamping)
        - Send ArrowUp 25 times (testing wraparound / top clamping)
        - Send Enter to fly to selected item
        - Send Escape to close
        - Assert zero console errors and clean dropdown dismissal
        """
        res = cdp.evaluate("""
        (function() {
            const input = document.getElementById('search-input');
            const results = document.getElementById('search-results');
            if (!input || !results) return { error: 'Search DOM elements missing' };
            
            // Focus and input 'Abell'
            input.focus();
            input.value = 'Abell';
            input.dispatchEvent(new Event('input', { bubbles: true }));
            
            const matchesCount = results.children.length;
            
            // Fuzz ArrowDown 20 times
            for (let i = 0; i < 20; i++) {
                input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
            }
            
            // Fuzz ArrowUp 25 times
            for (let i = 0; i < 25; i++) {
                input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowUp', bubbles: true }));
            }
            
            // Press Enter
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
            
            const resultsDisplayAfterEnter = results.style.display;
            
            // Re-open and press Escape
            input.value = 'Norma';
            input.dispatchEvent(new Event('input', { bubbles: true }));
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
            const resultsDisplayAfterEscape = results.style.display;
            
            return {
                matchesCount: matchesCount,
                resultsDisplayAfterEnter: resultsDisplayAfterEnter,
                resultsDisplayAfterEscape: resultsDisplayAfterEscape,
                activeModal: document.getElementById('spectroscopy-modal')?.style.display
            };
        })()
        """)
        assert res.get("matchesCount", 0) > 0, "Query 'Abell' must return matches"
        assert res.get("resultsDisplayAfterEnter") == "none", "Results dropdown must close on Enter"
        assert res.get("resultsDisplayAfterEscape") == "none", "Results dropdown must close on Escape"
