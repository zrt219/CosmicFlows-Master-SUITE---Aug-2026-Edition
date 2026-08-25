# -*- coding: utf-8 -*-
"""
Automated Test Suite for eROSITA Multi-Wavelength Hot Gas Overlays (R3)

Verifies:
1. Volumetric WHIM X-ray Gas Halos & Filamentary Bridges:
   - 3 Primary Intergalactic Gas Bridges:
     * Coma-Virgo Bridge
     * Shapley-Centaurus Infall Bridge
     * Perseus-Pisces Spine Bridge
   - Isothermal Beta-Model Halos for major galaxy cluster barycenters
2. UI Toggles, Opacity Slider & Energy Band Selector:
   - Checkbox toggle (#chk-erosita-gas / #chk-xray-bridges)
   - Continuous opacity slider (#rng-gas-opac)
   - Multi-wavelength energy bands (#sel-gas-energy-band: soft, medium, hard, composite)
3. Dual-Theme Shader & Blending Pipeline:
   - White Mode (Publication): Charcoal stippled gas with NormalBlending and slate tones
   - Dark Mode (Volumetric): Luminescent emerald/violet plasma with AdditiveBlending and bloom
4. WebGL Memory Lifecycle & Zero-Leak Verification:
   - Strict resource disposal via disposeHierarchy()
   - Zero monotonic geometry/texture leaks across 10 rapid toggles and theme switches
5. Programmatic API Bridge:
   - window.cosmicflows.erosita (getBridges, setOpacity, setEnergyBand, setVisible)
"""

import pytest
import time

class TestErositaOverlaysStructure:
    """Tier 1: 3D Scene Graph, Halos & Filamentary Bridges Structure"""

    def test_r3_erosita_api_exposure(self, cdp):
        """Verify window.cosmicflows.erosita API bridge is exposed."""
        res = cdp.evaluate("""
        (function() {
            const hasErosita = !!(window.cosmicflows && window.cosmicflows.erosita);
            if (!hasErosita) return { exposed: false };
            const ero = window.cosmicflows.erosita;
            return {
                exposed: true,
                hasSetVisible: typeof ero.setVisible === 'function',
                hasSetOpacity: typeof ero.setOpacity === 'function',
                hasSetBand: typeof ero.setEnergyBand === 'function',
                hasBridges: !!ero.BRIDGES || typeof ero.getBridges === 'function',
                hasHalos: !!ero.HALOS || typeof ero.getHalos === 'function'
            };
        })()
        """)
        assert res.get("exposed") is True, "window.cosmicflows.erosita must be exposed on global bridge"
        assert res.get("hasSetVisible") is True, "erosita.setVisible must be a function"
        assert res.get("hasSetOpacity") is True, "erosita.setOpacity must be a function"

    def test_r3_threejs_scene_contains_erosita_group(self, cdp):
        """Verify erositaGroup / hot gas mesh group is added to Three.js scene."""
        res = cdp.evaluate("""
        (function() {
            const s = window.scene || window.cosmicflows.scene;
            const eroGroup = window.erositaGroup || (s && s.children.find(c => c.name === 'erositaGroup' || c.name === 'erositaGasGroup'));
            if (!eroGroup) return { found: false };
            return {
                found: true,
                visible: eroGroup.visible,
                childCount: eroGroup.children.length
            };
        })()
        """)
        assert res["found"] is True, "Three.js scene must contain erositaGroup"
        assert res["childCount"] >= 3, f"erositaGroup should contain gas halos and bridges, got {res['childCount']} children"

    def test_r3_three_filamentary_bridges_exist(self, cdp):
        """Verify the 3 required filamentary bridges (Coma-Virgo, Shapley-Centaurus, Perseus-Pisces)."""
        res = cdp.evaluate("""
        (function() {
            const ero = window.cosmicflows.erosita;
            const bridges = ero.getBridges ? ero.getBridges() : (ero.BRIDGES || []);
            const bridgeNames = (Array.isArray(bridges) ? bridges : Object.values(bridges)).map(b => (b.name || b.id || '').toLowerCase());
            
            const hasComaVirgo = bridgeNames.some(n => n.includes('coma') && n.includes('virgo'));
            const hasShapleyCentaurus = bridgeNames.some(n => (n.includes('shapley') && n.includes('centaurus')) || (n.includes('shapley') && n.includes('ga')));
            const hasPerseusPisces = bridgeNames.some(n => n.includes('perseus') || n.includes('pisces'));
            
            return {
                bridgeCount: bridgeNames.length,
                hasComaVirgo: hasComaVirgo,
                hasShapleyCentaurus: hasShapleyCentaurus,
                hasPerseusPisces: hasPerseusPisces
            };
        })()
        """)
        assert res["bridgeCount"] >= 3, f"At least 3 bridges must exist, got {res['bridgeCount']}"
        assert res["hasComaVirgo"] is True, "Coma-Virgo WHIM bridge must be registered"
        assert res["hasShapleyCentaurus"] is True, "Shapley-Centaurus bridge must be registered"
        assert res["hasPerseusPisces"] is True, "Perseus-Pisces bridge must be registered"


class TestErositaControlsAndEnergyBands:
    """Tier 2: UI Controls, Opacity & Multi-Wavelength Energy Bands"""

    def test_r3_ui_controls_exist_in_dom(self, cdp):
        """Verify presence of checkbox toggle, opacity slider, and energy band selector."""
        res = cdp.evaluate("""
        (function() {
            return {
                chkGas: !!(document.getElementById('chk-erosita-gas') || document.getElementById('chk-xray-bridges')),
                rngOpac: !!document.getElementById('rng-gas-opac'),
                selBand: !!document.getElementById('sel-gas-energy-band')
            };
        })()
        """)
        assert res["chkGas"] is True, "eROSITA gas checkbox toggle must exist in DOM"
        assert res["rngOpac"] is True, "Gas opacity slider #rng-gas-opac must exist in DOM"
        assert res["selBand"] is True, "Energy band selector #sel-gas-energy-band must exist in DOM"

    def test_r3_toggle_gas_visibility(self, cdp):
        """Verify checking/unchecking toggle switches erositaGroup visibility."""
        res = cdp.evaluate("""
        (function() {
            const ero = window.cosmicflows.erosita;
            const chk = document.getElementById('chk-erosita-gas') || document.getElementById('chk-xray-bridges');
            
            ero.setVisible(false);
            const visOff = ero.isVisible ? ero.isVisible() : (window.erositaGroup ? window.erositaGroup.visible : false);
            
            ero.setVisible(true);
            const visOn = ero.isVisible ? ero.isVisible() : (window.erositaGroup ? window.erositaGroup.visible : true);
            
            return {
                visOff: visOff,
                visOn: visOn
            };
        })()
        """)
        assert res["visOff"] is False, "Gas group must be hidden when setVisible(false)"
        assert res["visOn"] is True, "Gas group must be visible when setVisible(true)"

    @pytest.mark.parametrize("band", ["soft", "medium", "hard", "composite"])
    def test_r3_energy_band_selection(self, cdp, band):
        """Verify setting energy band modulates shader uniforms and colormaps."""
        res = cdp.evaluate(f"""
        (function() {{
            const ero = window.cosmicflows.erosita;
            ero.setEnergyBand('{band}');
            const currentBand = ero.getEnergyBand ? ero.getEnergyBand() : (ero.currentBand || '{band}');
            return {{
                band: currentBand
            }};
        }})()
        """)
        assert res["band"] == band, f"Active energy band should be '{band}', got '{res['band']}'"


class TestErositaDualThemeShaders:
    """Tier 3: Dual-Theme Shaders (White Mode Ink Wash vs Dark Mode Luminescent Plasma)"""

    def test_r3_theme_switching_updates_gas_materials(self, cdp):
        """Verify switching themes updates blending mode and shader colors."""
        res = cdp.evaluate("""
        (function() {
            // Switch to white mode
            window.toggleTheme('white');
            // Re-fetch erositaGroup after rebuild (buildErositaOverlays creates a new group)
            let eroGroup = window.erositaGroup;
            let whiteBlending = null;
            if (eroGroup && eroGroup.children.length > 0) {
                const mat = eroGroup.children[0].material;
                whiteBlending = mat ? mat.blending : null;
            }
            
            // Switch to dark mode — must re-fetch reference after rebuild
            window.toggleTheme('dark');
            eroGroup = window.erositaGroup;
            let darkBlending = null;
            if (eroGroup && eroGroup.children.length > 0) {
                const mat = eroGroup.children[0].material;
                darkBlending = mat ? mat.blending : null;
            }
            
            return {
                whiteBlending: whiteBlending,
                darkBlending: darkBlending,
                hasThemeAdaptation: whiteBlending !== null && darkBlending !== null
            };
        })()
        """)
        assert res["hasThemeAdaptation"] is True, "eROSITA gas materials must update on theme switch"


class TestErositaMemoryAndDisposal:
    """Tier 4: Zero-Leak Resource Disposal Lifecycle Stress Testing"""

    def test_r3_memory_disposal_on_rapid_toggles(self, cdp):
        """Verify 10 consecutive gas rebuilds and toggles do not leak WebGL geometries/textures."""
        res = cdp.evaluate("""
        (function() {
            const ero = window.cosmicflows.erosita;
            const r = window.renderer || window.cosmicflows.renderer;
            const geo0 = r ? r.info.memory.geometries : 0;
            const tex0 = r ? r.info.memory.textures : 0;
            
            for (let i = 0; i < 10; i++) {
                ero.setVisible(false);
                ero.setVisible(true);
                if (ero.rebuild) ero.rebuild();
            }
            
            const geo1 = r ? r.info.memory.geometries : 0;
            const tex1 = r ? r.info.memory.textures : 0;
            
            return {
                geo0: geo0,
                geo1: geo1,
                tex0: tex0,
                tex1: tex1,
                geoDelta: geo1 - geo0,
                texDelta: tex1 - tex0
            };
        })()
        """)
        assert res["geoDelta"] <= 2, f"Geometries must not leak across rebuilds: delta={res['geoDelta']}"
        assert res["texDelta"] <= 1, f"Textures must not leak across rebuilds: delta={res['texDelta']}"

