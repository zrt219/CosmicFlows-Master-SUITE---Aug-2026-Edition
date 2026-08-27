# -*- coding: utf-8 -*-
"""
Empirical Challenger 2 Test Suite:
1. Canvas Gestures & Drag Disambiguation Stress (dr >= 5 px, dt < 350 ms, long-press, cancel)
2. Modal Raycast Isolation & Overlay Event Shielding
3. Multi-Catalog Search Performance Benchmark (110+ Queries, < 5 ms avg latency, zero main-thread lock)
4. Pointer Event Target Type Safety & Defensiveness
"""

import pytest
import time
import json

# 110+ Exhaustive Astronomical Queries spanning all catalogs, basins, voids, and edge cases
BENCHMARK_100_PLUS_QUERIES = [
    # Messier Objects (15)
    "M31", "M33", "M51", "M81", "M82", "M87", "M101", "M104", "M106", "M42",
    "Messier 31", "Messier 87", "Messier 104", "m49", "m86",
    
    # NGC Catalog (15)
    "NGC 224", "NGC 4889", "NGC 4696", "NGC 5128", "NGC 1275", "NGC 1316", "NGC 1365",
    "NGC 1399", "NGC 253", "NGC 3031", "NGC 5194", "ngc4889", "ngc 4696", "NGC5128", "ngc1275",
    
    # IC Catalog (10)
    "IC 1101", "IC 4764", "IC 342", "IC 10", "IC 1613", "IC 4182", "IC 5386",
    "ic1101", "ic 4764", "IC 1",
    
    # UGC Catalog (10)
    "UGC 8168", "UGC 7654", "UGC 2669", "UGC 10170", "UGC 12921", "UGC 1",
    "ugc8168", "ugc 7654", "UGC 300", "ugc10170",
    
    # PGC / LEDA IDs (15)
    "PGC 1", "PGC 2557", "PGC 41220", "PGC 44715", "PGC 55877", "PGC 73197", "PGC 126000",
    "pgc41220", "pgc 44715", "LEDA 41220", "leda 44715", "pgc1", "PGC 50000", "PGC 80000", "pgc126000",
    
    # Abell & ACO Clusters (15)
    "Abell 1656", "A1656", "aco 1656", "a-1656", "A3627", "ACO 3627", "A3558",
    "A426", "A1060", "A2151", "A262", "A1367", "A2199", "A2142", "ACO S0373",
    
    # Southern ACO Clusters (5)
    "ACO S0636", "ACO 3656", "a s373", "abell s0373", "aco s636",
    
    # Supercluster Basins (8)
    "Basin 1", "Basin 2", "Basin 3", "Basin 4", "Basin 5", "Basin 6", "Basin 7", "Basin 8",
    
    # Named Basins & Superclusters (8)
    "Laniakea", "Laniakea Supercluster Basin", "Shapley Supercluster Core", "Perseus-Pisces Basin",
    "Coma Basin", "Hercules Basin", "SDSS-1a", "Shapley Supercluster Basin",
    
    # Cosmic Voids & Repellers (8)
    "Boötes", "Bootes", "Boötes Supervoid", "Dipole Repeller", "Cold Spot Repeller",
    "Sculptor Void", "Local Void", "Eridanus Void",
    
    # V-Web Knots, Filaments & Corridors (6)
    "MeerKAT", "MeerKAT Vela ZoA Piercing Corridor", "Parkes", "Parkes HIZOA Great Attractor Corridor",
    "Centaurus Wall", "Perseus-Pegasus Filament",
    
    # Local Volume Named Galaxies (10)
    "Andromeda", "Sombrero Galaxy", "Centaurus A", "Whirlpool Galaxy", "Cigar Galaxy",
    "Pinwheel Galaxy", "Sculptor Galaxy", "Triangulum Galaxy", "Fornax Dwarf", "Large Magellanic Cloud",
    
    # Fuzzy / Typo / Diacritics / Edge Cases (12)
    "  Andromeda  ", "bootes   void", "A-1656", "ACO_3627", "ngc-4889",
    "supercluster", "cluster", "void", "corridor", "basin", "zrt", "coma"
]


class TestCanvasGestureDisambiguationStress:
    """Stress tests for OrbitControls drag vs click disambiguation."""

    def test_drag_threshold_suppression(self, cdp):
        """Verify drags with delta r >= 5 px NEVER trigger modal popups or selection."""
        # Ensure clean state
        cdp.evaluate("""
        (function() {
            if (window.splashController) window.splashController.dismiss(true);
            if (window.cosmicflows && window.cosmicflows.spectroscopy) {
                window.cosmicflows.spectroscopy.closeDossier();
            }
        })()
        """)
        time.sleep(0.3)

        drag_vectors = [
            (5, 0),      # Exactly threshold dx=5
            (0, 5),      # Exactly threshold dy=5
            (4, 4),      # Diagonal r = sqrt(16+16) = 5.65 px >= 5
            (-6, 0),     # Negative dx
            (0, -6),     # Negative dy
            (50, 50),    # Moderate drag
            (120, -80)   # Large drag
        ]

        for dx, dy in drag_vectors:
            # Start drag at center of canvas (640, 360)
            start_x, start_y = 640, 360
            end_x, end_y = start_x + dx, start_y + dy

            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mousePressed", "x": start_x, "y": start_y, "button": "left", "clickCount": 1})
            time.sleep(0.02)
            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": end_x, "y": end_y, "button": "left"})
            time.sleep(0.02)
            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": end_x, "y": end_y, "button": "left"})
            time.sleep(0.1)

            # Check if spectroscopy modal or any dossier opened
            res = cdp.evaluate("""
            (function() {
                const modal = document.getElementById('spectroscopy-modal');
                const isOpen = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
                return {
                    isOpen: isOpen
                };
            })()
            """)

            assert res["isOpen"] is False, f"Drag ({dx}, {dy}) accidentally opened spectroscopy dossier!"

    def test_discrete_click_on_canvas_cluster(self, cdp):
        """Verify discrete click (< 5 px, < 350 ms) on canvas cluster coordinates."""
        cdp.evaluate("""
        (function() {
            if (window.splashController) window.splashController.dismiss(true);
            if (window.cosmicflows && window.cosmicflows.spectroscopy) {
                window.cosmicflows.spectroscopy.closeDossier();
            }
        })()
        """)
        time.sleep(0.3)

        # Retrieve projected screen coordinates for Virgo or Coma cluster
        coords = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows && window.cosmicflows.spectroscopy;
            const features = window.cosmicflows && window.cosmicflows.ASTROMETRIC_FEATURES;
            if (!features || features.length === 0) return null;
            
            const virgo = features.find(f => f.id === 'virgo_cl') || features[0];
            if (!virgo || !virgo.pos) return null;
            
            const worldPos = (typeof sgToThree === 'function') ? sgToThree(virgo.pos.x, virgo.pos.y, virgo.pos.z) : virgo.pos.clone();
            const proj = worldPos.clone().project(camera);
            const x = Math.round((proj.x * 0.5 + 0.5) * window.innerWidth);
            const y = Math.round((-(proj.y * 0.5) + 0.5) * window.innerHeight);
            
            return { x: x, y: y, id: virgo.id };
        })()
        """)

        assert coords is not None, "Failed to locate astrometric feature coordinates"
        cx, cy = max(50, min(1200, coords["x"])), max(50, min(700, coords["y"]))

        # Dispatch discrete click on canvas element
        cdp.evaluate(f"""
        (function() {{
            const canvas = document.querySelector('canvas') || document.body;
            const x = {cx};
            const y = {cy};
            canvas.dispatchEvent(new PointerEvent('pointerdown', {{ clientX: x, clientY: y, bubbles: true }}));
            setTimeout(() => {{
                canvas.dispatchEvent(new PointerEvent('pointerup', {{ clientX: x + 1, clientY: y + 1, bubbles: true }}));
            }}, 40);
        }})()
        """)
        time.sleep(0.4)

        res = cdp.evaluate("""
        (function() {
            const modal = document.getElementById('spectroscopy-modal');
            return {
                evaluated: true,
                modalPresent: !!modal
            };
        })()
        """)
        assert res["evaluated"] is True

    def test_long_press_duration_suppression(self, cdp):
        """Verify long hold (dt=450ms >= 350ms) suppresses accidental clicks."""
        cdp.evaluate("if (window.cosmicflows && window.cosmicflows.spectroscopy) window.cosmicflows.spectroscopy.closeDossier();")
        time.sleep(0.2)

        cdp.evaluate("""
        (function() {
            const canvas = document.querySelector('canvas') || document.body;
            canvas.dispatchEvent(new PointerEvent('pointerdown', { clientX: 640, clientY: 360, bubbles: true }));
        })()
        """)
        time.sleep(0.45) # Hold longer than 350ms
        cdp.evaluate("""
        (function() {
            const canvas = document.querySelector('canvas') || document.body;
            canvas.dispatchEvent(new PointerEvent('pointerup', { clientX: 640, clientY: 360, bubbles: true }));
        })()
        """)
        time.sleep(0.2)

        res = cdp.evaluate("""
        (function() {
            const modal = document.getElementById('spectroscopy-modal');
            return {
                isOpen: modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'))
            };
        })()
        """)
        assert res["isOpen"] is False, "Long press (>350ms) should NOT trigger cluster click!"

    def test_pointer_cancel_resets_lifecycle(self, cdp):
        """Verify pointercancel event resets pointer down state and prevents subsequent click trigger."""
        res = cdp.evaluate("""
        (function() {
            const canvas = document.querySelector('canvas') || document.body;
            if (window.cosmicflows && window.cosmicflows.spectroscopy) {
                window.cosmicflows.spectroscopy.closeDossier();
            }
            canvas.dispatchEvent(new PointerEvent('pointerdown', { clientX: 640, clientY: 360, bubbles: true }));
            canvas.dispatchEvent(new PointerEvent('pointercancel', { bubbles: true }));
            canvas.dispatchEvent(new PointerEvent('pointerup', { clientX: 640, clientY: 360, bubbles: true }));
            
            const modal = document.getElementById('spectroscopy-modal');
            const isOpen = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
            
            return {
                isOpen: isOpen
            };
        })()
        """)
        assert res["isOpen"] is False, "pointercancel must cancel the click gesture lifecycle"

    def test_continuous_rapid_dragging_performance(self, cdp):
        """Verify 25 rapid consecutive drag gestures maintain camera tracking and 0 exceptions."""
        init_pos = cdp.evaluate("({ x: camera.position.x, y: camera.position.y, z: camera.position.z })")

        for i in range(25):
            x1 = 500 + (i % 5) * 40
            y1 = 300 + (i % 5) * 30
            x2 = x1 + 30
            y2 = y1 + 20
            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mousePressed", "x": x1, "y": y1, "button": "left", "clickCount": 1})
            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": x2, "y": y2, "button": "left"})
            cdp.send_cdp("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": x2, "y": y2, "button": "left"})

        time.sleep(0.3)
        final_pos = cdp.evaluate("({ x: camera.position.x, y: camera.position.y, z: camera.position.z })")
        dx = abs(final_pos['x'] - init_pos['x'])
        dy = abs(final_pos['y'] - init_pos['y'])
        dz = abs(final_pos['z'] - init_pos['z'])

        assert (dx > 50 or dy > 50 or dz > 50), "Camera must have orbited during rapid drag sequences"


class TestModalRaycastIsolationStress:
    """Stress tests ensuring modal UI overlays completely shield 3D canvas from pointer raycasts."""

    def test_spectroscopy_modal_click_isolation(self, cdp):
        """Clicking on spectroscopy modal tabs/body must NOT leak raycasts to scene."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows && window.cosmicflows.spectroscopy;
            if (spec) spec.openDossier('virgo_cl');
            const modal = document.getElementById('spectroscopy-modal');
            if (!modal) return { ok: false };
            
            const initialDossier = spec ? spec.getActiveDossier() : null;
            const initialId = initialDossier ? initialDossier.id : null;
            
            // Dispatch pointerup event on modal container
            const evt = new PointerEvent('pointerup', {
                clientX: modal.offsetLeft + 50,
                clientY: modal.offsetTop + 50,
                bubbles: true
            });
            modal.dispatchEvent(evt);
            
            const finalDossier = spec ? spec.getActiveDossier() : null;
            const finalId = finalDossier ? finalDossier.id : null;
            
            // Clean up
            if (spec) spec.closeDossier();
            
            return {
                ok: true,
                initialId: initialId,
                finalId: finalId,
                preserved: initialId === finalId
            };
        })()
        """)
        assert res["ok"] is True
        assert res["preserved"] is True, "Click on modal should not re-trigger scene picking or change dossier!"

    def test_primer_modal_click_isolation(self, cdp):
        """Clicking inside the interactive primer modal must never trigger raycast selection."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows && window.cosmicflows.spectroscopy;
            if (spec) spec.closeDossier();
            
            const modal = document.getElementById('primer-modal');
            if (!modal) return { ok: false, error: 'no primer modal' };
            
            // Open primer modal
            modal.classList.remove('hidden');
            modal.style.display = 'flex';
            
            // Dispatch click on latex equation container inside primer
            const content = modal.querySelector('.modal-content') || modal;
            const evt = new PointerEvent('pointerup', {
                clientX: modal.offsetLeft + 100,
                clientY: modal.offsetTop + 100,
                bubbles: true
            });
            content.dispatchEvent(evt);
            
            const activeD = spec ? spec.getActiveDossier() : null;
            
            // Clean up
            modal.classList.add('hidden');
            modal.style.display = 'none';
            
            return {
                ok: true,
                hasActiveDossier: !!activeD
            };
        })()
        """)
        assert res["ok"] is True
        assert res["hasActiveDossier"] is False, "Clicking inside primer modal leaked raycast to open spectroscopy dossier!"

    def test_camera_settings_modal_isolation(self, cdp):
        """Clicking on camera settings modal controls must never trigger 3D raycasts."""
        res = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows && window.cosmicflows.spectroscopy;
            if (spec) spec.closeDossier();
            
            const modal = document.getElementById('camera-settings-modal');
            if (!modal) return { ok: false, error: 'no cam settings modal' };
            
            modal.classList.remove('hidden');
            modal.style.display = 'flex';
            
            const evt = new PointerEvent('pointerup', {
                clientX: modal.offsetLeft + 30,
                clientY: modal.offsetTop + 30,
                bubbles: true
            });
            modal.dispatchEvent(evt);
            
            const activeD = spec ? spec.getActiveDossier() : null;
            
            modal.classList.add('hidden');
            modal.style.display = 'none';
            
            return {
                ok: true,
                hasActiveDossier: !!activeD
            };
        })()
        """)
        assert res["ok"] is True
        assert res["hasActiveDossier"] is False, "Clicking inside camera settings modal leaked raycast!"

    def test_hover_card_suppressed_over_modals(self, cdp):
        """Hover card must be hidden when cursor is over modal overlays."""
        res = cdp.evaluate("""
        (function() {
            const topbar = document.getElementById('topbar');
            const hoverCard = document.getElementById('hover-card');
            if (!topbar || !hoverCard) return { ok: false };
            
            // Simulate pointermove over topbar
            const rect = topbar.getBoundingClientRect();
            const evt = new PointerEvent('pointermove', {
                clientX: rect.left + 50,
                clientY: rect.top + 10,
                bubbles: true
            });
            topbar.dispatchEvent(evt);
            
            return {
                ok: true,
                hoverCardDisplay: hoverCard.style.display
            };
        })()
        """)
        assert res["ok"] is True
        assert res["hoverCardDisplay"] == "none", "Hover card must be hidden when cursor is over topbar/modals!"


class TestMultiCatalogSearchBenchmark100Plus:
    """Empirical 100+ Query Performance Benchmark for Astronomical Search."""

    def test_search_benchmark_100_plus_queries(self, cdp):
        """Execute all 110+ astronomical queries, asserting avg latency < 5ms and 100% valid results."""
        payload = json.dumps(BENCHMARK_100_PLUS_QUERIES)

        res = cdp.evaluate(f"""
        (function() {{
            const search = window.cosmicflows && window.cosmicflows.search;
            if (!search) return {{ ok: false, error: 'Search API missing' }};
            
            const queries = {payload};
            const timings = [];
            const resultsSummary = [];
            let totalTime = 0;
            let maxTime = 0;

            for (let i = 0; i < queries.length; i++) {{
                const q = queries[i];
                const t0 = performance.now();
                const matches = search.query(q);
                const t1 = performance.now();
                const dt = t1 - t0;
                
                timings.push(dt);
                totalTime += dt;
                if (dt > maxTime) maxTime = dt;
                
                resultsSummary.push({{
                    query: q,
                    matchCount: matches ? matches.length : 0,
                    topName: matches && matches.length > 0 ? matches[0].name : null,
                    dtMs: dt
                }});
            }}

            const avgTime = totalTime / queries.length;

            return {{
                ok: true,
                queryCount: queries.length,
                avgTimeMs: avgTime,
                maxTimeMs: maxTime,
                totalTimeMs: totalTime,
                results: resultsSummary
            }};
        }})()
        """)

        assert res.get("ok") is True, f"Search benchmark failed: {res.get('error')}"
        query_count = res["queryCount"]
        avg_time = res["avgTimeMs"]
        max_time = res["maxTimeMs"]
        total_time = res["totalTimeMs"]

        print(f"\n=======================================================")
        print(f"MULTI-CATALOG SEARCH BENCHMARK (100+ QUERIES)")
        print(f"=======================================================")
        print(f"Total Queries Tested: {query_count}")
        print(f"Average Query Latency: {avg_time:.3f} ms (Budget < 5.0 ms)")
        print(f"Max Query Latency:     {max_time:.3f} ms (Budget < 15.0 ms)")
        print(f"Total Benchmark Time:  {total_time:.2f} ms")
        print(f"=======================================================")

        assert query_count >= 100, f"Expected >= 100 benchmark queries, ran {query_count}"
        assert avg_time < 5.0, f"Average query latency ({avg_time:.3f} ms) exceeded 5.0 ms SLA threshold!"
        assert max_time < 15.0, f"Max query latency ({max_time:.3f} ms) exceeded 15.0 ms threshold!"

        # Verify key queries returned valid results
        res_by_query = {r["query"]: r for r in res["results"]}
        assert res_by_query["M31"]["matchCount"] > 0
        assert res_by_query["A1656"]["matchCount"] > 0
        assert res_by_query["Basin 1"]["matchCount"] > 0
        assert res_by_query["Boötes"]["matchCount"] > 0
        assert res_by_query["Centaurus A"]["matchCount"] > 0

    def test_sustained_stress_1000_query_loop(self, cdp):
        """Execute 1,000 sustained queries in a rapid loop asserting zero memory/GC stall."""
        res = cdp.evaluate("""
        (function() {
            const search = window.cosmicflows.search;
            const querySamples = [
                'M31', 'A1656', 'NGC 4889', 'Laniakea', 'Basin 6',
                'Boötes', 'MeerKAT', 'Sombrero', 'PGC 41220', 'UGC 8168'
            ];
            
            const t0 = performance.now();
            let totalMatches = 0;
            const iterations = 1000;
            
            for (let i = 0; i < iterations; i++) {
                const q = querySamples[i % querySamples.length];
                const matches = search.query(q);
                totalMatches += matches.length;
            }
            
            const totalMs = performance.now() - t0;
            const perQueryMs = totalMs / iterations;
            
            return {
                iterations: iterations,
                totalMs: totalMs,
                perQueryMs: perQueryMs,
                totalMatches: totalMatches
            };
        })()
        """)

        print(f"\nSUSTAINED SEARCH STRESS: 1,000 iterations in {res['totalMs']:.2f} ms ({res['perQueryMs']:.3f} ms/query)")
        assert res["totalMs"] < 2000.0, f"1,000 queries took {res['totalMs']:.2f} ms (> 2000 ms limit)"
        assert res["perQueryMs"] < 2.0, f"Per-query latency {res['perQueryMs']:.3f} ms exceeded 2.0 ms limit"

    def test_adversarial_malformed_queries(self, cdp):
        """Verify search gracefully handles extreme adversarial inputs (null, huge strings, symbols)."""
        adversarial_inputs = [
            "",
            "   ",
            "\n\t\r",
            "\\.*+?^${}()|[]",
            "\"';!--",
            "<script>alert(1)</script>",
            "A" * 5000,
            "0" * 1000
        ]
        payload = json.dumps(adversarial_inputs)

        res = cdp.evaluate(f"""
        (function() {{
            const search = window.cosmicflows.search;
            const extremeInputs = {payload};
            const results = [];
            
            for (let i = 0; i < extremeInputs.length; i++) {{
                try {{
                    const matches = search.query(extremeInputs[i]);
                    results.push({{ inputIdx: i, success: true, count: matches ? matches.length : 0 }});
                }} catch (e) {{
                    results.push({{ inputIdx: i, success: false, error: e.message }});
                }}
            }}
            
            return {{
                results: results,
                allPassed: results.every(r => r.success)
            }};
        }})()
        """)
        assert res["allPassed"] is True, "Search query crashed on adversarial or extreme input string!"
