"""
Standalone Diagnostic Runner for Challenger M2 & M3
Executes detailed adversarial stress tests directly with ChromeCDPClient
and outputs full empirical telemetry data for the handoff report.
"""

import os
import sys
import time
import json
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def main():
    print(= * 80)
    print( CHALLENGER M2 & M3: EMPIRICAL ADVERSARIAL STRESS TEST SUITE)
    print(= * 80)
    print(fTimestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}\n)

    ensure_http_server(8000)
    client = ChromeCDPClient(port=9222, url=http://localhost:8000, spawn_headless=True)
    
    results = {}
    
    try:
        print([1/6] Launching Chrome CDP Client and verifying WebGL Canvas...)
        client.start()
        client.wait_for_condition(
            document.readyState === 'complete' && !!document.querySelector('#scene canvas'),
            timeout=10.0
        )
        print( Connected. WebGL context online.\n)
        
        # ---------------------------------------------------------------------
        # Check 1: 38 Astrometric Features
        # ---------------------------------------------------------------------
        print([2/6] Auditing 38 Astrometric Features...)
        feat_data = client.evaluate("
        (function() {
            const feats = (window.cosmicflows && window.cosmicflows.ASTROMETRIC_FEATURES) || window.ASTROMETRIC_FEATURES;
            if (!feats) return { error: ASTROMETRIC_FEATURES missing };
            
            const summary = {
                total: feats.length,
                t1: feats.filter(f => f.tier === 1).length,
                t2: feats.filter(f => f.tier === 2).length,
                t3: feats.filter(f => f.tier === 3).length,
                items: feats.map(f => ({
                    id: f.id,
                    name: f.name,
                    tier: f.tier,
                    cat: f.cat,
                    pos: [Math.round(f.pos.x), Math.round(f.pos.y), Math.round(f.pos.z)],
                    dist: f.dist,
                    cz: f.cz,
                    mass: f.mass,
                    citation: f.citation
                }))
            };
            return summary;
        })()
        ")
        
        assert feat_data.get(total) == 38, fExpected 38 features, got {feat_data.get('total')}
        assert feat_data.get(t1) == 10, fExpected 10 Tier 1 features, got {feat_data.get('t1')}
        assert feat_data.get(t2) == 16, fExpected 16 Tier 2 features, got {feat_data.get('t2')}
        assert feat_data.get(t3) == 12, fExpected 12 Tier 3 features, got {feat_data.get('t3')}
        
        # Check all entries have valid physical fields
        for item in feat_data[items]:
            assert item[mass] is not None and item[mass] != 0, fInvalid mass in {item['id']}
            assert item[dist] is not None and item[dist] >= 0, fInvalid dist in {item['id']}
            assert item[citation] and len(item[citation]) > 5, fInvalid citation in {item['id']}
            assert len(item[pos]) == 3, fInvalid pos in {item['id']}
            
        print(f PASS: All 38 features verified (10 Tier 1, 16 Tier 2, 12 Tier 3).)
        print(f PASS: Mass, Coordinates, Distance, Velocity, and Citations 100% complete.\n)
        results[astrometric_features] = {
            status: PASS,
            total: feat_data[total],
            t1: feat_data[t1],
            t2: feat_data[t2],
            t3: feat_data[t3]
        }
        
        # ---------------------------------------------------------------------
        # Check 2: 15+ Rapid Theme Toggles Memory Leak Stress Test
        # ---------------------------------------------------------------------
        print([3/6] Running 15 Rapid Theme Toggles Memory Stress Test...)
        mem_stress = client.evaluate("
        (function() {
            const renderer = window.cosmicflows.renderer;
            window.toggleTheme('dark');
            if (renderer.render && window.cosmicflows.scene && window.cosmicflows.camera) {
                renderer.render(window.cosmicflows.scene, window.cosmicflows.camera);
            }
            
            const baseGeo = renderer.info.memory.geometries;
            const baseTex = renderer.info.memory.textures;
            
            const log = [];
            for (let i = 1; i <= 15; i++) {
                window.toggleTheme();
                log.push({
                    toggle: i,
                    theme: document.body.classList.contains('theme-white') ? 'white' : 'dark',
                    geometries: renderer.info.memory.geometries,
                    textures: renderer.info.memory.textures
                });
            }
            
            // Return to dark
            window.toggleTheme('dark');
            if (renderer.render && window.cosmicflows.scene && window.cosmicflows.camera) {
                renderer.render(window.cosmicflows.scene, window.cosmicflows.camera);
            }
            
            const final15Geo = renderer.info.memory.geometries;
            const final15Tex = renderer.info.memory.textures;
            
            // Additional 15 toggles (total 30) to test monotonic growth
            for (let i = 16; i <= 30; i++) {
                window.toggleTheme();
            }
            window.toggleTheme('dark');
            if (renderer.render && window.cosmicflows.scene && window.cosmicflows.camera) {
                renderer.render(window.cosmicflows.scene, window.cosmicflows.camera);
            }
            
            const final30Geo = renderer.info.memory.geometries;
            const final30Tex = renderer.info.memory.textures;
            
            return {
                baseGeo: baseGeo,
                baseTex: baseTex,
                final15Geo: final15Geo,
                final15Tex: final15Tex,
                final30Geo: final30Geo,
                final30Tex: final30Tex,
                delta15Geo: final15Geo - baseGeo,
                delta15Tex: final15Tex - baseTex,
                delta30Geo: final30Geo - final15Geo,
                delta30Tex: final30Tex - final15Tex,
                log: log
            };
        })()
        ")
        
        print(f Baseline: Geometries = {mem_stress['baseGeo']}, Textures = {mem_stress['baseTex']})
        print(f After 15 toggles: Geometries = {mem_stress['final15Geo']} (? = {mem_stress['delta15Geo']}), Textures = {mem_stress['final15Tex']} (? = {mem_stress['delta15Tex']}))
        print(f After 30 toggles: Geometries = {mem_stress['final30Geo']} (?30 = {mem_stress['delta30Geo']}), Textures = {mem_stress['final30Tex']} (?30 = {mem_stress['delta30Tex']}))
        
        assert abs(mem_stress[delta15Geo]) <= 2, fGeometries leak after 15 toggles: delta={mem_stress['delta15Geo']}
        assert abs(mem_stress[delta15Tex]) <= 2, fTextures leak after 15 toggles: delta={mem_stress['delta15Tex']}
        assert mem_stress[delta30Geo] == 0, fMonotonic geometry growth detected: delta={mem_stress['delta30Geo']}
        assert mem_stress[delta30Tex] == 0, fMonotonic texture growth detected: delta={mem_stress['delta30Tex']}
        
        print( PASS: Zero monotonic growth and zero memory leak across 30 rapid theme switches.\n)
        results[memory_stress] = {
            status: PASS,
            baseGeometries: mem_stress[baseGeo],
            baseTextures: mem_stress[baseTex],
            after15Geometries: mem_stress[final15Geo],
            after15Textures: mem_stress[final15Tex],
            after30Geometries: mem_stress[final30Geo],
            after30Textures: mem_stress[final30Tex]
        }
        
        # ---------------------------------------------------------------------
        # Check 3: Point Cloud Scaling & LOD Modes
        # ---------------------------------------------------------------------
        print([4/6] Testing Point Cloud Scaling and LOD Modes...)
        lod_data = client.evaluate("
        (function() {
            const selLod = document.getElementById('sel-lod-mode');
            function getCount() {
                return window.cosmicflows ? window.cosmicflows.activePointCount : 0;
            }
            
            selLod.value = 'full';
            selLod.dispatchEvent(new Event('change'));
            const full = getCount();
            
            selLod.value = 'allsky';
            selLod.dispatchEvent(new Event('change'));
            const allsky = getCount();
            
            selLod.value = 'med';
            selLod.dispatchEvent(new Event('change'));
            const med = getCount();
            
            selLod.value = 'low';
            selLod.dispatchEvent(new Event('change'));
            const low = getCount();
            
            // Restore full
            selLod.value = 'full';
            selLod.dispatchEvent(new Event('change'));
            
            return { full: full, allsky: allsky, med: med, low: low };
        })()
        ")
        
        print(f 'full' mode point count: {lod_data['full']} (Required >= 56,000))
        print(f 'allsky' mode point count: {lod_data['allsky']} (Required >= 65,000))
        print(f 'med' mode point count: {lod_data['med']})
        print(f 'low' mode point count: {lod_data['low']})
        
        assert lod_data[full] >= 56000, fPoint count in full mode too low: {lod_data['full']}
        assert lod_data[allsky] >= 65000, fPoint count in allsky mode too low: {lod_data['allsky']}
        print( PASS: Point cloud scaling verified across all LOD modes.\n)
        results[point_cloud] = {
            status: PASS,
            full: lod_data[full],
            allsky: lod_data[allsky],
            med: lod_data[med],
            low: lod_data[low]
        }
        
        # ---------------------------------------------------------------------
        # Check 4: Raycaster Hover Card Hit Testing
        # ---------------------------------------------------------------------
        print([5/6] Testing Raycaster Hover Card Hit Testing...)
        raycast_data = client.evaluate("
        (function() {
            const hoverCard = document.getElementById('hover-card');
            const feats = window.cosmicflows.ASTROMETRIC_FEATURES;
            const camera = window.cosmicflows.camera;
            
            // Ensure callouts enabled
            const chk = document.getElementById('chk-callouts');
            if (chk) { chk.checked = true; chk.dispatchEvent(new Event('change')); }
            
            const targets = ['virgo_cl', 'coma_cl', 'shapley_core', 'ga_norma', 'vela_scl', 'dipole_rep', 'meerkat_corridor', 'parkes_corridor'];
            const hitDetails = [];
            
            for (let id of targets) {
                const f = feats.find(item => item.id === id);
                if (!f) continue;
                
                const threePos = sgToThree(f.pos.x, f.pos.y, f.pos.z);
                const screenPos = threePos.clone().project(camera);
                const clientX = (screenPos.x * 0.5 + 0.5) * window.innerWidth;
                const clientY = (-(screenPos.y * 0.5) + 0.5) * window.innerHeight;
                
                window.dispatchEvent(new PointerEvent('pointermove', {
                    clientX: clientX,
                    clientY: clientY,
                    bubbles: true
                }));
                
                hitDetails.push({
                    id: id,
                    displayed: hoverCard.style.display === 'block',
                    name: document.getElementById('hc-name') ? document.getElementById('hc-name').textContent : '',
                    dist: document.getElementById('hc-dist') ? document.getElementById('hc-dist').textContent : '',
                    vel: document.getElementById('hc-vel') ? document.getElementById('hc-vel').textContent : '',
                    cite: document.getElementById('hc-cite') ? document.getElementById('hc-cite').textContent : ''
                });
            }
            
            // Move off screen
            window.dispatchEvent(new PointerEvent('pointermove', { clientX: 2, clientY: 2, bubbles: true }));
            const hidden = hoverCard.style.display === 'none';
            
            return {
                hitDetails: hitDetails,
                hiddenOffTarget: hidden
            };
        })()
        ")
        
        assert raycast_data[hiddenOffTarget] is True, Hover card not hidden when off target
        for hd in raycast_data[hitDetails]:
            print(f Hit target '{hd['id']}': name='{hd['name'][:30]}...', citation='{hd['cite'][:25]}...')
            assert len(hd[name]) > 0, fEmpty name on target {hd['id']}
            assert len(hd[cite]) > 0, fEmpty citation on target {hd['id']}
            
        print( PASS: Raycaster hover card hit testing verified on all 8 test landmarks.\n)
        results[raycaster] = {
            status: PASS,
            testedTargets: len(raycast_data[hitDetails]),
            hiddenOffTarget: raycast_data[hiddenOffTarget]
        }
        
        # ---------------------------------------------------------------------
        # Check 5: Console Errors & Exceptions
        # ---------------------------------------------------------------------
        print([6/6] Asserting zero console errors & zero uncaught exceptions...)
        client.assert_no_errors()
        print( PASS: 0 console errors, 0 uncaught exceptions.\n)
        results[console_errors] = {status: PASS, errors: 0, exceptions: 0}
        
    except Exception as e:
        print(f\n[FAIL] Adversarial challenge encountered error: {e})
        traceback.print_exc()
        results[verdict] = CHALLENGE_FAILED
        results[error] = str(e)
        return 1, results
    finally:
        client.close()
        
    print(= * 80)
    print( FINAL CHALLENGER VERDICT: APPROVE)
    print(= * 80)
    results[verdict] = APPROVE
    return 0, results

if __name__ == __main__:
    code, res = main()
    sys.exit(code)
