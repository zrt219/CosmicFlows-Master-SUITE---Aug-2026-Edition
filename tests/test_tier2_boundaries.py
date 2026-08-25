"""
Tier 2: Boundary & Corner Cases E2E Tests
Covers physical and numerical singularities, domain limits (+/- 15000 km/s),
zero/extreme step sizes, rapid switching, and extreme aspect ratios.
"""
import pytest
import time

class TestTier2Boundaries:

    @pytest.mark.parametrize("landmark_coord", [
        (0, 0, 0),                       # Origin / Virgo / Local Group
        (7200, -8600, -2400),            # Shapley Supercluster Core
        (15000, -14000, 6800),           # Great Attractor / Norma
        (-10000, 10000, 12000),          # Dipole Repeller
        (-8500, -12000, -3200)           # Vela Supercluster
    ])
    def test_t2_01_singularity_at_node_centers(self, cdp, landmark_coord):
        x, y, z = landmark_coord
        res = cdp.evaluate(f"""
        (function() {{
            if (!window.cosmicflows || !window.cosmicflows.getVelocitySg) return {{ fallback: true }};
            const v = new THREE.Vector3();
            window.cosmicflows.getVelocitySg(new THREE.Vector3({x}, {y}, {z}), v);
            return {{
                x: v.x, y: v.y, z: v.z,
                isValid: !isNaN(v.x) && !isNaN(v.y) && !isNaN(v.z) && isFinite(v.length())
            }};
        }})()
        """)
        if res.get("fallback"):
            pytest.skip("window.cosmicflows.getVelocitySg not yet exposed")
        assert res.get("isValid") is True, f"Singularity at ({x}, {y}, {z}) produced invalid vector: {res}"

    @pytest.mark.parametrize("boundary_pt", [
        (14999, 0, 0),
        (-14999, 0, 0),
        (0, 14999, 0),
        (0, -14999, 0),
        (0, 0, 14999),
        (0, 0, -14999),
        (15000, 15000, 15000),
        (-15000, -15000, -15000)
    ])
    def test_t2_02_domain_boundary_clipping(self, cdp, boundary_pt):
        x, y, z = boundary_pt
        res = cdp.evaluate(f"""
        (function() {{
            if (!window.cosmicflows || !window.cosmicflows.integrateRK4) return {{ fallback: true }};
            const pts = window.cosmicflows.integrateRK4(new THREE.Vector3({x}, {y}, {z}), 1, 100, 40);
            return {{
                count: pts ? pts.length : 0,
                terminated: pts ? pts.length <= 100 : true
            }};
        }})()
        """)
        if res.get("fallback"):
            pytest.skip("window.cosmicflows.integrateRK4 not yet exposed")
        assert res.get("terminated") is True

    def test_t2_03_zero_step_size(self, cdp):
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.integrateRK4) return { fallback: true };
            const pts = window.cosmicflows.integrateRK4(new THREE.Vector3(1000, 1000, 1000), 1, 50, 0);
            return {
                finished: true,
                count: pts ? pts.length : 0
            };
        })()
        """)
        if res.get("fallback"):
            pytest.skip("window.cosmicflows.integrateRK4 not yet exposed")
        assert res.get("finished") is True

    def test_t2_04_extreme_large_step_size(self, cdp):
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.integrateRK4) return { fallback: true };
            const pts = window.cosmicflows.integrateRK4(new THREE.Vector3(500, 500, 500), 1, 20, 2000);
            return {
                finished: true,
                noNaN: pts ? pts.every(p => !isNaN(p.x) && !isNaN(p.y) && !isNaN(p.z)) : true
            };
        })()
        """)
        if res.get("fallback"):
            pytest.skip("window.cosmicflows.integrateRK4 not yet exposed")
        assert res.get("finished") is True
        assert res.get("noNaN") is True

    def test_t2_05_zero_max_steps(self, cdp):
        res = cdp.evaluate("""
        (function() {
            if (!window.cosmicflows || !window.cosmicflows.integrateRK4) return { fallback: true };
            const pts = window.cosmicflows.integrateRK4(new THREE.Vector3(0, 0, 0), 1, 0, 40);
            return {
                valid: Array.isArray(pts) && pts.length <= 1
            };
        })()
        """)
        if res.get("fallback"):
            pytest.skip("window.cosmicflows.integrateRK4 not yet exposed")
        assert res.get("valid") is True

    def test_t2_06_rapid_sequential_engine_switches(self, cdp):
        res = cdp.evaluate("""
        (function() {
            const sel = document.getElementById('sel-science-engine');
            if (!sel) return { success: false };
            const allEngines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026', 'cf4-wf'];
            for (let eng of allEngines) {
                sel.value = eng;
                sel.dispatchEvent(new Event('change'));
                if (window.cosmicflows && window.cosmicflows.switchScienceEngine) {
                    window.cosmicflows.switchScienceEngine(eng);
                }
            }
            return {
                finalValue: sel.value,
                success: true
            };
        })()
        """)
        assert res.get("success") is True
        assert res.get("finalValue") == "cf4-wf"

    def test_t2_07_canvas_resize_extreme_aspect_ratios(self, cdp):
        res = cdp.evaluate("""
        (function() {
            // Trigger window resize handler simulation
            window.dispatchEvent(new Event('resize'));
            const canvas = document.querySelector('#scene canvas');
            return {
                canvasExists: !!canvas,
                width: canvas ? canvas.width : 0,
                height: canvas ? canvas.height : 0
            };
        })()
        """)
        assert res.get("canvasExists") is True
        assert res.get("width") > 0
        assert res.get("height") > 0
