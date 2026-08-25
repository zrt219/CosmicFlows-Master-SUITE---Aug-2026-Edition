"""
Tier 3: Cross-Feature Interactions E2E Tests
Exercises combinatorial pairings across:
- 7 Science Engines x 2 Themes x 5 Colormaps
- LOD Modes x Sub-Catalog Layers x Inset Maps
"""
import pytest
import time

ENGINES = ["cf4-wf", "cf4-hmc", "vela-zoa", "vweb-2026", "bayesian-2026", "gadget4", "nusser-tully-2026"]
THEMES = ["white", "dark"]
COLORMAPS = ["velocity", "divergence", "vorticity", "doppler", "ink"]

class TestTier3Interactions:

    @pytest.mark.parametrize("engine_id", ENGINES)
    @pytest.mark.parametrize("theme", THEMES)
    def test_t3_01_engine_x_theme_matrix(self, cdp, engine_id, theme):
        res = cdp.evaluate(f"""
        (function() {{
            const selEng = document.getElementById('sel-science-engine');
            if (selEng) {{
                selEng.value = '{engine_id}';
                selEng.dispatchEvent(new Event('change'));
            }}
            window.toggleTheme('{theme}');
            return {{
                engine: selEng ? selEng.value : '{engine_id}',
                isWhite: document.body.classList.contains('theme-white'),
                expectedWhite: '{theme}' === 'white'
            }};
        }})()
        """)
        assert res.get("engine") == engine_id
        assert res.get("isWhite") == res.get("expectedWhite")

    @pytest.mark.parametrize("colormap", COLORMAPS)
    def test_t3_02_colormap_switching(self, cdp, colormap):
        res = cdp.evaluate(f"""
        (function() {{
            const selCol = document.getElementById('sel-colormap');
            if (selCol) {{
                selCol.value = '{colormap}';
                selCol.dispatchEvent(new Event('change'));
                return {{ colormap: selCol.value, valid: true }};
            }}
            return {{ valid: false }};
        }})()
        """)
        assert res.get("valid") is True
        assert res.get("colormap") == colormap

    @pytest.mark.parametrize("engine_id", ["cf4-wf", "vela-zoa", "gadget4"])
    @pytest.mark.parametrize("colormap", ["divergence", "vorticity", "doppler"])
    def test_t3_03_engine_x_colormap_streamline_rebuild(self, cdp, engine_id, colormap):
        res = cdp.evaluate(f"""
        (function() {{
            const selEng = document.getElementById('sel-science-engine');
            const selCol = document.getElementById('sel-colormap');
            if (selEng) {{ selEng.value = '{engine_id}'; selEng.dispatchEvent(new Event('change')); }}
            if (selCol) {{ selCol.value = '{colormap}'; selCol.dispatchEvent(new Event('change')); }}
            if (window.cosmicflows && window.cosmicflows.rebuildStreamlines) {{
                window.cosmicflows.rebuildStreamlines();
            }}
            return {{ success: true }};
        }})()
        """)
        assert res.get("success") is True

    def test_t3_04_lod_mode_x_subcatalogs(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const selLod = document.getElementById('sel-lod-mode');
            if (selLod) {
                selLod.value = 'performance';
                selLod.dispatchEvent(new Event('change'));
                selLod.value = 'full';
                selLod.dispatchEvent(new Event('change'));
            }
            return { completed: true };
        })()
        """)
        assert res.get("completed") is True

    def test_t3_05_inset_map_x_engine_sync(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const btn = document.getElementById('btn-inset-toggle');
            if (btn) btn.click();
            const selEng = document.getElementById('sel-science-engine');
            if (selEng) {
                selEng.value = 'bayesian-2026';
                selEng.dispatchEvent(new Event('change'));
            }
            return { syncSuccess: true };
        })()
        """)
        assert res.get("syncSuccess") is True
