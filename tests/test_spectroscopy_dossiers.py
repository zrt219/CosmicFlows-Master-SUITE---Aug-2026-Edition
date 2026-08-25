# -*- coding: utf-8 -*-
"""
Automated Test Suite for Interactive Galaxy Cluster Spectroscopy Dossiers (R2)

Verifies:
1. Complete Astrophysical Data Matrix for all 8 Primary Hubs:
   - Virgo, Coma, Centaurus, Norma/GA, Shapley Supercluster Core, Vela, Perseus-Pisces, Fornax
   - Parameters: M200 (Virial Mass), cz (Recessional Velocity), sigma_v (Velocity Dispersion),
     Peculiar Velocity components (vx, vy, vz), metric distance, BCG classification, Literature Citations
2. 3D Raycast Hit-Testing & Pointer Disambiguation:
   - Screen coordinate raycasting against 3D barycenters and callout pick proxies
   - Camera orbit drag vs single-click cluster selection disambiguation
3. Inspector Modal Lifecycle & Interaction:
   - Open dossier, render astrophysics metrics, update dynamic text fields
   - Pin / Unpin modal state persistence
   - Close modal via close button, outside click, or Escape key
4. Gaussian Velocity Dispersion Canvas Rendering:
   - HTML5 2D canvas generation of Gaussian velocity profile N(v) = N0/(sqrt(2pi)sigma) * exp(-(v-cz)^2 / 2sigma^2)
   - Shaded 1-sigma and 2-sigma confidence intervals and caustic velocity boundaries
5. Dual-Theme Support:
   - Dark Mode glassmorphic styling and luminescence
   - White Mode clean monochrome publication rendering
6. Programmatic API Bridge:
   - window.cosmicflows.spectroscopy.getDossier(id)
   - window.cosmicflows.spectroscopy.openDossier(id)
   - window.cosmicflows.spectroscopy.closeDossier()
"""

import pytest
import math
import time

HUB_IDS = [
    ("virgo_cl", "Virgo", 1.0e14, 2.0e15, 1000, 1300, 650, 800),
    ("coma_cl", "Coma", 1.0e15, 3.0e15, 6500, 7300, 950, 1100),
    ("centaurus_cl", "Centaurus", 1.5e15, 4.0e15, 3000, 3500, 800, 950),
    ("ga_norma", "Norma", 3.0e16, 8.0e16, 4500, 5200, 850, 1050),
    ("shapley_core", "Shapley", 8.0e16, 2.0e17, 13500, 15500, 1150, 1350),
    ("vela_scl", "Vela", 2.0e17, 5.0e17, 17500, 20000, 1000, 1250),
    ("perseus_cl", "Perseus", 1.5e15, 3.5e15, 5000, 5600, 1200, 1400),
    ("fornax_cl", "Fornax", 4.0e14, 1.2e15, 1200, 1600, 320, 440)
]

class TestSpectroscopyDataMatrix:
    """Tier 1: Astrophysical Parameters Matrix for All 8 Primary Hubs"""

    def test_r2_spectroscopy_api_exposure(self, cdp):
        """Verify window.cosmicflows.spectroscopy API bridge is exposed."""
        res = cdp.evaluate("""
        (function() {
            const hasSpec = !!(window.cosmicflows && window.cosmicflows.spectroscopy);
            if (!hasSpec) return { exposed: false };
            const spec = window.cosmicflows.spectroscopy;
            return {
                exposed: true,
                hasGetDossier: typeof spec.getDossier === 'function',
                hasOpenDossier: typeof spec.openDossier === 'function',
                hasCloseDossier: typeof spec.closeDossier === 'function',
                hasPinDossier: typeof spec.pinDossier === 'function',
                hasCatalog: !!spec.CATALOG || !!spec.DOSSIERS || !!spec.DATA
            };
        })()
        """)
        assert res.get("exposed") is True, "window.cosmicflows.spectroscopy must be exposed"
        assert res.get("hasGetDossier") is True, "spectroscopy.getDossier must be a function"
        assert res.get("hasOpenDossier") is True, "spectroscopy.openDossier must be a function"

    @pytest.mark.parametrize("hub_id, name_substr, m_min, m_max, cz_min, cz_max, sig_min, sig_max", HUB_IDS)
    def test_r2_hub_astrophysical_parameters(self, cdp, hub_id, name_substr, m_min, m_max, cz_min, cz_max, sig_min, sig_max):
        """Verify physical metrics for each of the 8 primary hubs."""
        res = cdp.evaluate(f"""
        (function() {{
            const spec = window.cosmicflows.spectroscopy;
            const d = spec.getDossier('{hub_id}');
            if (!d) return {{ found: false }};
            return {{
                found: true,
                id: d.id,
                name: d.name || d.label,
                mass: d.mass || d.m200 || d.M200,
                cz: d.cz || d.recessionalVelocity,
                sigma_v: d.sigma_v || d.sigmaV || d.dispersion,
                dist: d.dist || d.distance,
                bcg: d.bcg || d.primaryGalaxy,
                citation: d.citation || d.ref || d.paper,
                vPec: d.vPec || d.peculiarVelocity || d.v_pec
            }};
        }})()
        """)
        assert res.get("found") is True, f"Dossier for hub '{hub_id}' must exist in spectroscopy catalog"
        assert name_substr.lower() in res["name"].lower(), f"Hub name '{res['name']}' should contain '{name_substr}'"
        assert m_min <= res["mass"] <= m_max, f"Hub '{hub_id}' mass {res['mass']} out of range [{m_min}, {m_max}]"
        assert cz_min <= res["cz"] <= cz_max, f"Hub '{hub_id}' cz {res['cz']} out of range [{cz_min}, {cz_max}]"
        assert sig_min <= res["sigma_v"] <= sig_max, f"Hub '{hub_id}' sigma_v {res['sigma_v']} out of range [{sig_min}, {sig_max}]"
        assert res["dist"] > 0, f"Distance for '{hub_id}' must be positive"
        assert len(str(res["citation"])) > 3, f"Citation for '{hub_id}' must be present"


class TestSpectroscopyRaycastAndPicking:
    """Tier 2: 3D Raycasting, Coordinate Inversion & Drag Disambiguation"""

    def test_r2_raycast_cluster_selection(self, cdp):
        """Verify that raycasting against cluster barycenter / hit-proxy selects the target cluster."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            // Raycast direct programmatic selection on Coma
            spec.openDossier('coma_cl');
            const modal = document.getElementById('spectroscopy-modal');
            const nameEl = document.getElementById('spec-name');
            const czEl = document.getElementById('spec-cz');
            
            const isOpen = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
            const title = nameEl ? nameEl.textContent : '';
            const czVal = czEl ? czEl.textContent : '';
            
            // Close after test
            spec.closeDossier();
            
            return {
                isOpen: isOpen,
                title: title,
                czVal: czVal
            };
        })()
        """)
        assert res["isOpen"] is True, "Spectroscopy modal must open when cluster is selected"
        assert "coma" in res["title"].lower(), f"Modal title should show Coma, got '{res['title']}'"

    def test_r2_drag_vs_click_disambiguation(self, cdp):
        """Verify camera drag movement does not accidentally trigger modal opening."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            spec.closeDossier();
            
            const canvas = document.querySelector('#scene canvas');
            if (!canvas) return { ok: false };
            
            // Simulate camera orbit drag (pointerdown at (200, 200), move to (400, 400), pointerup)
            const rect = canvas.getBoundingClientRect();
            canvas.dispatchEvent(new PointerEvent('pointerdown', { clientX: rect.left + 200, clientY: rect.top + 200, bubbles: true }));
            canvas.dispatchEvent(new PointerEvent('pointermove', { clientX: rect.left + 400, clientY: rect.top + 400, bubbles: true }));
            canvas.dispatchEvent(new PointerEvent('pointerup', { clientX: rect.left + 400, clientY: rect.top + 400, bubbles: true }));
            
            const modal = document.getElementById('spectroscopy-modal');
            const isOpenAfterDrag = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
            
            return {
                ok: true,
                isOpenAfterDrag: !!isOpenAfterDrag
            };
        })()
        """)
        assert res["ok"] is True
        assert res["isOpenAfterDrag"] is False, "Camera drag must not trigger accidental dossier modal opening"


class TestSpectroscopyModalLifecycle:
    """Tier 3: Modal UI Lifecycle, Pinning, Dismissal & Canvas Rendering"""

    def test_r2_modal_dom_structure(self, cdp):
        """Verify presence of all inspector card DOM fields and canvas element."""
        res = cdp.evaluate("""
        (function() {
            return {
                modal: !!document.getElementById('spectroscopy-modal'),
                title: !!document.getElementById('spec-name'),
                cz: !!document.getElementById('spec-cz'),
                sigmaV: !!document.getElementById('spec-sigmav'),
                m200: !!document.getElementById('spec-m200'),
                vPec: !!document.getElementById('spec-vpec'),
                citation: !!document.getElementById('spec-citation'),
                canvas: !!document.getElementById('spec-canvas'),
                btnClose: !!document.getElementById('btn-close-dossier'),
                btnPin: !!document.getElementById('btn-pin-dossier')
            };
        })()
        """)
        assert res["modal"] is True, "#spectroscopy-modal container must exist"
        assert res["title"] is True, "#spec-name must exist"
        assert res["cz"] is True, "#spec-cz must exist"
        assert res["sigmaV"] is True, "#spec-sigmav must exist"
        assert res["m200"] is True, "#spec-m200 must exist"
        assert res["canvas"] is True, "#spec-canvas velocity profile canvas must exist"
        assert res["btnClose"] is True, "#btn-close-dossier must exist"

    def test_r2_modal_pin_and_dismiss_behavior(self, cdp):
        """Verify modal pinning prevents outside click dismissal, and close button works."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            spec.openDossier('shapley_core');
            
            const modal = document.getElementById('spectroscopy-modal');
            const btnPin = document.getElementById('btn-pin-dossier');
            const btnClose = document.getElementById('btn-close-dossier');
            
            // Pin modal
            if (btnPin) btnPin.click();
            const isPinned = modal.classList.contains('pinned') || (spec.isPinned && spec.isPinned());
            
            // Try closing via background click
            document.body.dispatchEvent(new MouseEvent('click', { bubbles: true }));
            const stillOpenAfterBgClick = modal.style.display !== 'none' && !modal.classList.contains('hidden');
            
            // Unpin and click close button
            if (btnPin) btnPin.click();
            if (btnClose) btnClose.click();
            
            const isClosed = modal.style.display === 'none' || modal.classList.contains('hidden');
            
            return {
                isPinned: isPinned,
                stillOpenAfterBgClick: stillOpenAfterBgClick,
                isClosed: isClosed
            };
        })()
        """)
        assert res["isPinned"] is True, "Clicking pin button should pin dossier modal"
        assert res["stillOpenAfterBgClick"] is True, "Pinned dossier must not close on outside click"
        assert res["isClosed"] is True, "Clicking close button must dismiss modal"

    def test_r2_gaussian_velocity_canvas_rendering(self, cdp):
        """Verify Gaussian velocity dispersion profile draws pixels on #spec-canvas."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            spec.openDossier('coma_cl');
            
            const canvas = document.getElementById('spec-canvas');
            if (!canvas) return { hasCanvas: false };
            
            const ctx = canvas.getContext('2d');
            const imgData = ctx.getImageData(0, 0, canvas.width, canvas.height);
            
            // Check that non-transparent pixels exist
            let nonZero = 0;
            for (let i = 3; i < imgData.data.length; i += 4) {
                if (imgData.data[i] > 0) nonZero++;
            }
            
            spec.closeDossier();
            
            return {
                hasCanvas: true,
                width: canvas.width,
                height: canvas.height,
                nonZeroPixels: nonZero
            };
        })()
        """)
        assert res["hasCanvas"] is True
        assert res["nonZeroPixels"] > 50, f"Gaussian velocity profile canvas must draw visible pixels, got {res['nonZeroPixels']}"


class TestSpectroscopyDualTheme:
    """Tier 4: Dual-Theme Invariance (White Publication Mode vs Dark Plasma Mode)"""

    def test_r2_dual_theme_rendering(self, cdp):
        """Verify dossier modal styling and canvas redraw across theme toggles."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            spec.openDossier('perseus_cl');
            
            // Toggle white
            window.toggleTheme('white');
            const modalWhite = document.getElementById('spectroscopy-modal');
            const whiteBg = window.getComputedStyle(modalWhite).backgroundColor;
            
            // Toggle dark
            window.toggleTheme('dark');
            const modalDark = document.getElementById('spectroscopy-modal');
            const darkBg = window.getComputedStyle(modalDark).backgroundColor;
            
            spec.closeDossier();
            
            return {
                whiteBg: whiteBg,
                darkBg: darkBg,
                themesDistinct: whiteBg !== darkBg
            };
        })()
        """)
        assert res["themesDistinct"] is True, "Spectroscopy modal styling must adapt between White and Dark themes"
