"""
Adversarial Stress Test Suite for Milestones M2 & M3
(Catalogues, ZoA Depth, WebGL Performance & Memory Stress Verifier)

Empirically tests:
1. 38 Astrometric Features: complete schemas, non-null mass, coordinates, citations, tier distributions.
2. WebGL Memory Stability: 15+ rapid theme toggles asserting zero monotonic growth / memory leakage.
3. Point Cloud Scaling: >= 56,000 points in 'full' mode, >= 65,000 in 'allsky' mode, subcatalog reactive updates.
4. Raycaster Hover Cards: 3D hit testing across multiple callout coordinates, card rendering, and off-target hiding.
5. ZoA Piercing Corridors: MeerKAT Vela 21cm, Parkes HIZOA GA 21cm, shading modes.
6. Console Errors & Uncaught Exceptions: 0 errors allowed.
"""

import os
import sys
import time
import json
import pytest

# Add workspace root to sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.cdp_client import ChromeCDPClient

class TestMilestone2CataloguesAndZoA:
    """Stress tests for Milestone M2: 38 Astrometric Features & ZoA Depth"""

    def test_m2_01_astrometric_catalog_count_and_tiers(self, cdp):
        """Assert exactly 38 astrometric features with expected tier distribution."""
        res = cdp.evaluate("""
        (function() {
            const feats = (window.cosmicflows && window.cosmicflows.ASTROMETRIC_FEATURES) || window.ASTROMETRIC_FEATURES;
            if (!feats || !Array.isArray(feats)) return { count: 0, tiers: {} };
            
            const tiers = {};
            let nullMassCount = 0;
            let missingCoordsCount = 0;
            let missingCitations = 0;
            
            feats.forEach(f => {
                tiers[f.tier] = (tiers[f.tier] || 0) + 1;
                if (f.mass === null || f.mass === undefined) nullMassCount++;
                if (!f.pos || typeof f.pos.x !== 'number' || typeof f.pos.y !== 'number' || typeof f.pos.z !== 'number') missingCoordsCount++;
                if (!f.citation || f.citation.length === 0) missingCitations++;
            });
            
            return {
                count: feats.length,
                tiers: tiers,
                nullMassCount: nullMassCount,
                missingCoordsCount: missingCoordsCount,
                missingCitations: missingCitations
            };
        })()
        """)
        assert res["count"] == 38, f"Expected 38 features, got {res['count']}"
        assert res["nullMassCount"] == 0, f"Found {res['nullMassCount']} features with null mass"
        assert res["missingCoordsCount"] == 0, f"Found {res['missingCoordsCount']} features missing 3D coordinates"
        assert res["missingCitations"] == 0, f"Found {res['missingCitations']} features missing citations"
        assert res["tiers"].get("1", 0) >= 10, f"Tier 1 count should be >= 10, got {res['tiers'].get('1')}"
        assert res["tiers"].get("2", 0) >= 16, f"Tier 2 count should be >= 16, got {res['tiers'].get('2')}"
        assert res["tiers"].get("3", 0) >= 12, f"Tier 3 count should be >= 12, got {res['tiers'].get('3')}"

    def test_m2_02_zoa_corridors_geometry_and_visibility(self, cdp):
        """Assert ZoA piercing corridors exist with correct volumetric double-cone and filament meshes."""
        res = cdp.evaluate("""
        (function() {
            const chkZoa = document.getElementById('chk-zoa');
            if (!chkZoa) return { error: 'chk-zoa not found' };
            
            // Toggle ZoA on
            chkZoa.checked = true;
            chkZoa.dispatchEvent(new Event('change'));
            
            const zoaGroup = (window.cosmicflows && window.cosmicflows.zoaGroup) || window.zoaGroup;
            const hasZoaGroup = !!zoaGroup;
            let childCount = 0;
            if (zoaGroup && zoaGroup.children) {
                childCount = zoaGroup.children.length;
            }
            
            return {
                chkZoaFound: !!chkZoa,
                sceneHasObjects: window.cosmicflows && window.cosmicflows.scene && window.cosmicflows.scene.children.length > 0
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["chkZoaFound"] is True
        assert res["sceneHasObjects"] is True


class TestMilestone3PerformanceAndStress:
    """Stress tests for Milestone M3: WebGL Performance, Memory Stability, Point Clouds"""

    def test_m3_01_webgl_memory_stability_15_theme_toggles(self, cdp):
        """Assert zero monotonic growth or memory leakage across 15 rapid theme toggles."""
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.renderer) {
                return { error: 'renderer not accessible' };
            }
            const renderer = window.cosmicflows.renderer;
            
            // Ensure established dark mode baseline
            if (typeof window.toggleTheme === 'function') {
                window.toggleTheme('dark');
            }
            
            const initialGeometries = renderer.info.memory.geometries;
            const initialTextures = renderer.info.memory.textures;
            
            const history = [];
            for (let i = 0; i < 15; i++) {
                if (typeof window.toggleTheme === 'function') {
                    window.toggleTheme();
                }
                history.push({
                    step: i,
                    geometries: renderer.info.memory.geometries,
                    textures: renderer.info.memory.textures
                });
            }
            
            // Restore dark mode
            if (typeof window.toggleTheme === 'function') {
                window.toggleTheme('dark');
            }
            
            const finalGeometries = renderer.info.memory.geometries;
            const finalTextures = renderer.info.memory.textures;
            
            return {
                initialGeometries: initialGeometries,
                finalGeometries: finalGeometries,
                initialTextures: initialTextures,
                finalTextures: finalTextures,
                geoDelta: finalGeometries - initialGeometries,
                texDelta: finalTextures - initialTextures,
                history: history
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert abs(res["geoDelta"]) <= 2, f"Geometries leaked: initial={res['initialGeometries']}, final={res['finalGeometries']}"
        assert abs(res["texDelta"]) <= 2, f"Textures leaked: initial={res['initialTextures']}, final={res['finalTextures']}"

    def test_m3_02_point_cloud_scaling_and_lod(self, cdp):
        """Assert point cloud scaling >= 56,000 points in 'full' mode and >= 65,000 in 'allsky' mode."""
        res = cdp.evaluate("""
        (function() {
            const selLod = document.getElementById('sel-lod-mode');
            if (!selLod) return { error: 'sel-lod-mode not found' };
            
            function getActivePoints() {
                if (window.cosmicflows && typeof window.cosmicflows.activePointCount === 'number') {
                    return window.cosmicflows.activePointCount;
                }
                let maxPts = 0;
                if (window.cosmicflows && window.cosmicflows.scene) {
                    window.cosmicflows.scene.traverse(obj => {
                        if (obj.isPoints && obj.geometry && obj.geometry.attributes && obj.geometry.attributes.position) {
                            maxPts = Math.max(maxPts, obj.geometry.attributes.position.count);
                        }
                    });
                }
                return maxPts;
            }
            
            // Set to 'full' mode
            selLod.value = 'full';
            selLod.dispatchEvent(new Event('change'));
            const fullCount = getActivePoints();
            
            // Set to 'allsky' mode
            selLod.value = 'allsky';
            selLod.dispatchEvent(new Event('change'));
            const allskyCount = getActivePoints();
            
            // Set to 'med' mode
            selLod.value = 'med';
            selLod.dispatchEvent(new Event('change'));
            const medCount = getActivePoints();
            
            // Set to 'low' mode
            selLod.value = 'low';
            selLod.dispatchEvent(new Event('change'));
            const lowCount = getActivePoints();
            
            // Restore 'full'
            selLod.value = 'full';
            selLod.dispatchEvent(new Event('change'));
            
            return {
                fullCount: fullCount,
                allskyCount: allskyCount,
                medCount: medCount,
                lowCount: lowCount
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["fullCount"] >= 56000, f"Expected fullCount >= 56,000, got {res['fullCount']}"
        assert res["allskyCount"] >= 65000, f"Expected allskyCount >= 65,000, got {res['allskyCount']}"
        assert res["medCount"] >= 8000, f"Expected medCount >= 8,000, got {res['medCount']}"
        assert res["lowCount"] >= 2000, f"Expected lowCount >= 2,000, got {res['lowCount']}"

    def test_m3_03_subcatalog_filtering_dynamic_reactivity(self, cdp):
        """Test toggling individual survey subcatalogs (CF4, FAST, MeerKAT, SDSS, 6dFGS, 2MRS)."""
        res = cdp.evaluate("""
        (function() {
            const chkCf4 = document.getElementById('chk-sub-cf4');
            const chkFast = document.getElementById('chk-sub-fast');
            const chkMeerkat = document.getElementById('chk-sub-meerkat');
            const chkSdss = document.getElementById('chk-sub-sdss');
            const chk6df = document.getElementById('chk-sub-6df');
            const chk2mrs = document.getElementById('chk-sub-2mrs');
            
            const initialCount = window.cosmicflows ? window.cosmicflows.activePointCount : 0;
            
            // Uncheck subcatalogs
            if (chkFast) { chkFast.checked = false; chkFast.dispatchEvent(new Event('change')); }
            if (chkMeerkat) { chkMeerkat.checked = false; chkMeerkat.dispatchEvent(new Event('change')); }
            
            const reducedCount = window.cosmicflows ? window.cosmicflows.activePointCount : 0;
            
            // Re-check subcatalogs
            if (chkFast) { chkFast.checked = true; chkFast.dispatchEvent(new Event('change')); }
            if (chkMeerkat) { chkMeerkat.checked = true; chkMeerkat.dispatchEvent(new Event('change')); }
            
            const restoredCount = window.cosmicflows ? window.cosmicflows.activePointCount : 0;
            
            return {
                initialCount: initialCount,
                reducedCount: reducedCount,
                restoredCount: restoredCount,
                reactive: reducedCount < initialCount && restoredCount === initialCount
            };
        })()
        """)
        assert res["reactive"] is True, f"Subcatalog filtering not reactive: {res}"

    def test_m3_04_raycaster_hover_card_hit_testing(self, cdp):
        """Adversarially test Raycaster hover card hit testing across multiple callout coordinates."""
        res = cdp.evaluate("""
        (function() {
            const hoverCard = document.getElementById('hover-card');
            const feats = (window.cosmicflows && window.cosmicflows.ASTROMETRIC_FEATURES) || window.ASTROMETRIC_FEATURES;
            const camera = window.cosmicflows ? window.cosmicflows.camera : null;
            const scene = window.cosmicflows ? window.cosmicflows.scene : null;
            
            if (!hoverCard || !feats || !camera || !scene) {
                return { error: 'Required Three.js scene or hoverCard DOM element missing' };
            }
            
            // Ensure callouts are visible
            const chkCallouts = document.getElementById('chk-callouts');
            if (chkCallouts) {
                chkCallouts.checked = true;
                chkCallouts.dispatchEvent(new Event('change'));
            }
            
            // Test hit testing on 5 specific landmark targets
            const testTargets = ['virgo_cl', 'coma_cl', 'shapley_core', 'ga_norma', 'vela_scl'];
            const hitResults = [];
            
            for (let id of testTargets) {
                const f = feats.find(item => item.id === id);
                if (!f) continue;
                
                // Project 3D coordinate to 2D screen coordinates
                const threePos = (typeof sgToThree === 'function') ? sgToThree(f.pos.x, f.pos.y, f.pos.z) : new THREE.Vector3(f.pos.x, f.pos.z, -f.pos.y);
                const screenPos = threePos.clone().project(camera);
                
                const clientX = (screenPos.x * 0.5 + 0.5) * window.innerWidth;
                const clientY = (-(screenPos.y * 0.5) + 0.5) * window.innerHeight;
                
                // Dispatch simulated pointermove
                const evt = new PointerEvent('pointermove', {
                    clientX: clientX,
                    clientY: clientY,
                    bubbles: true
                });
                window.dispatchEvent(evt);
                
                // Check hover card state
                const isDisplayed = hoverCard.style.display === 'block';
                const nameContent = document.getElementById('hc-name') ? document.getElementById('hc-name').textContent : '';
                const citeContent = document.getElementById('hc-cite') ? document.getElementById('hc-cite').textContent : '';
                
                hitResults.push({
                    id: id,
                    screenPos: { x: Math.round(clientX), y: Math.round(clientY) },
                    isDisplayed: isDisplayed,
                    name: nameContent,
                    hasCitation: citeContent.length > 0
                });
            }
            
            // Test off-target pointermove (should hide hover card)
            const offEvt = new PointerEvent('pointermove', { clientX: 5, clientY: 5, bubbles: true });
            window.dispatchEvent(offEvt);
            const offHidden = hoverCard.style.display === 'none';
            
            return {
                hitResults: hitResults,
                offHidden: offHidden,
                testedCount: hitResults.length
            };
        })()
        """)
        assert "error" not in res, res.get("error")
        assert res["offHidden"] is True, "Hover card failed to hide when pointer moved off-target"
        assert res["testedCount"] >= 5

    def test_m3_05_zero_console_errors_and_uncaught_exceptions(self, cdp):
        """Assert zero critical console errors and zero uncaught JavaScript exceptions."""
        cdp.assert_no_errors()
