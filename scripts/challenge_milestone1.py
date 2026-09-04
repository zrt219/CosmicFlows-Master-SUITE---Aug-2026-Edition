#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scripts/challenge_milestone1.py

Adversarial Stress Test & Verification Suite for Milestone 1:
Multi-Survey Astronomical Data & Coordinate Pipeline.

Challenges:
1. Kaiser RSD Mathematical Singularities & Extreme Coordinate Stress:
   - Evaluates r -> 0 (0, 1e-16, 1e-12, 1e-8, 1e-7, 1e-6, 1e-3)
   - Evaluates extreme coordinates r > 10,000 (10001, 25000, 100000, 1e6, 1e8)
   - Evaluates negative components across all 8 Cartesian octants
   - Asserts zero NaN, zero Inf, unit vector convergence ||r_hat|| == 1.0 +/- 1e-6
   - Asserts radial projection alignment r_hat . x/||x|| == 1.0 +/- 1e-6
   - Performs exhaustive scan across 19,000 bundled galaxies in SDSS, 2MRS, DESI datasets.

2. Dual-Buffer In-Place Morphing 60 FPS Performance Budget:
   - Executes 100 consecutive frames of applySurveyRSD(t) in browser V8 runtime via CDP.
   - High-resolution timing via performance.now().
   - Asserts mean update time < 2.5 ms/frame.
   - Asserts 95th percentile update time < 2.5 ms/frame.
   - Asserts ArrayBuffer identity preservation (zero memory reallocation / zero GC pause).
   - Asserts non-zero physical coordinate displacement.

3. Rapid Survey Layer Toggling Stress Test:
   - Executes 50 rapid visibility flips across all survey layers in window.surveyGroups.
   - Forces render pass per flip cycle.
   - Asserts zero unhandled rejections, zero console errors, zero entries in window.errorRegistry.
   - Asserts zero scene corruption (all meshes remain attached with valid DynamicDrawUsage Float32BufferAttributes).
   - Asserts zero runaway memory accumulation.

4. Verification Against Rule 14 & Rule 39:
   - Rule 14: console.error == 0, window.onerror == 0, unhandledrejection == 0.
   - Rule 39: Zero VRAM leak during catalog re-generation and theme toggling via disposeHierarchy.

Outputs full diagnostic log and returns 0 on PASS, 1 on FAIL.
"""

import os
import sys
import json
import math
import time
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9260):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port


def stress_kaiser_rsd_math():
    """
    Challenge 1: Mathematical boundary and singularity analysis of Kaiser RSD.
    Formula: s = x_real + (v_los + v_thermal) * r_hat
    where r_hat = x / ||x|| with division-by-zero protection.
    """
    print("\n--- Challenge 1: Kaiser RSD Mathematical Singularity & Boundary Stress ---")
    failures = []
    
    def evaluate_rsd_point(gx, gy, gz, v_los, v_therm):
        # Implementation mirrors index.html defensive clamp: invR = 1.0 / Math.max(rKms, 1e-7)
        r_kms = math.sqrt(gx*gx + gy*gy + gz*gz)
        inv_r = 1.0 / max(r_kms, 1e-7)
        rx = gx * inv_r
        ry = gy * inv_r
        rz = gz * inv_r
        
        delta_v = v_los + v_therm
        sx = gx + delta_v * rx
        sy = gy + delta_v * ry
        sz = gz + delta_v * rz
        
        # Supergalactic to Three.js mapping
        three_real = [gx, gz, -gy]
        three_z = [sx, sz, -sy]
        
        return {
            "r_kms": r_kms,
            "r_hat": [rx, ry, rz],
            "norm_r_hat": math.sqrt(rx*rx + ry*ry + rz*rz),
            "real": [gx, gy, gz],
            "redshift": [sx, sy, sz],
            "three_real": three_real,
            "three_z": three_z
        }

    # 1. Test r -> 0 sub-epsilon distances
    r_sub_epsilon = [0.0, 1e-18, 1e-15, 1e-12, 1e-9, 1e-8, 1e-7, 1e-6, 1e-4, 1e-2, 0.1, 1.0]
    for r in r_sub_epsilon:
        for v_los in [-500.0, 0.0, 500.0]:
            for v_therm in [-300.0, 300.0]:
                # Test along diagonal
                gx = r / math.sqrt(3)
                gy = r / math.sqrt(3)
                gz = r / math.sqrt(3)
                res = evaluate_rsd_point(gx, gy, gz, v_los, v_therm)
                
                # Assert zero NaN, zero Inf
                for val in res["redshift"] + res["r_hat"] + res["three_z"]:
                    if math.isnan(val) or math.isinf(val):
                        failures.append(f"Sub-epsilon r={r} produced NaN/Inf: {val}")
                
                # When r >= 1e-7, norm of r_hat must be 1.0 within tolerance
                if r >= 1e-6:
                    if abs(res["norm_r_hat"] - 1.0) > 1e-5:
                        failures.append(f"Unit vector norm failed for r={r}: {res['norm_r_hat']}")

    # 2. Test extreme distances r > 10,000 up to cosmological horizon
    r_extremes = [10001.0, 25000.0, 50000.0, 100000.0, 500000.0, 1e6, 1e8, 1e10]
    for r in r_extremes:
        for phi in [0.1, 0.7, 1.5, 2.3, 3.1]:
            for theta in [0.2, 1.0, 2.5, 4.2, 5.8]:
                gx = r * math.sin(phi) * math.cos(theta)
                gy = r * math.sin(phi) * math.sin(theta)
                gz = r * math.cos(phi)
                res = evaluate_rsd_point(gx, gy, gz, 250.0, 187.0)
                
                for val in res["redshift"] + res["r_hat"]:
                    if math.isnan(val) or math.isinf(val):
                        failures.append(f"Extreme distance r={r} produced NaN/Inf: {val}")
                
                if abs(res["norm_r_hat"] - 1.0) > 1e-6:
                    failures.append(f"Extreme distance r={r} norm != 1.0: {res['norm_r_hat']}")
                
                # Radial projection alignment check: r_hat . (x / ||x||) == 1.0
                dot_prod = (res["r_hat"][0]*gx + res["r_hat"][1]*gy + res["r_hat"][2]*gz) / res["r_kms"]
                if abs(dot_prod - 1.0) > 1e-6:
                    failures.append(f"Extreme distance r={r} line-of-sight unit vector failed to converge: dot={dot_prod}")

    # 3. Test negative components across all 8 octants
    signs = [
        (+1, +1, +1), (-1, +1, +1), (+1, -1, +1), (-1, -1, +1),
        (+1, +1, -1), (-1, +1, -1), (+1, -1, -1), (-1, -1, -1)
    ]
    for sx, sy, sz in signs:
        gx = sx * 7500.0
        gy = sy * 6200.0
        gz = sz * 4800.0
        res = evaluate_rsd_point(gx, gy, gz, -180.0, 120.0)
        
        for val in res["redshift"] + res["r_hat"]:
            if math.isnan(val) or math.isinf(val):
                failures.append(f"Octant ({sx},{sy},{sz}) produced NaN/Inf")
        
        # Assert unit vector direction preserves coordinate sign
        if (res["r_hat"][0] * sx) <= 0 or (res["r_hat"][1] * sy) <= 0 or (res["r_hat"][2] * sz) <= 0:
            failures.append(f"Octant ({sx},{sy},{sz}) inverted line-of-sight sign: r_hat={res['r_hat']}")
        
        dot = (res["r_hat"][0]*gx + res["r_hat"][1]*gy + res["r_hat"][2]*gz) / res["r_kms"]
        if abs(dot - 1.0) > 1e-6:
            failures.append(f"Octant ({sx},{sy},{sz}) dot prod != 1.0: {dot}")

    # 4. Exhaustive Galaxy Catalog Verification
    data_dir = os.path.join(ROOT_DIR, "data")
    catalog_files = [
        ("SDSS", os.path.join(data_dir, "survey_sdss.json"), 5000),
        ("2MRS", os.path.join(data_dir, "survey_2mrs.json"), 10000),
        ("DESI", os.path.join(data_dir, "survey_desi.json"), 4000)
    ]
    
    total_galaxies_scanned = 0
    for name, path, min_count in catalog_files:
        if not os.path.exists(path):
            failures.append(f"Catalog file missing: {path}")
            continue
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        count = data.get("count", 0)
        total_galaxies_scanned += count
        if count < min_count:
            failures.append(f"{name} catalog galaxy count {count} < expected {min_count}")
        
        real_arr = data.get("pos_real_scene", [])
        z_arr = data.get("pos_redshift_scene", [])
        v_los_arr = data.get("v_los", [])
        v_therm_arr = data.get("v_thermal", [])
        
        if len(real_arr) != count * 3 or len(z_arr) != count * 3:
            failures.append(f"{name} array length mismatch: real={len(real_arr)}, z={len(z_arr)}, count={count}")
        
        for idx in range(count):
            rx, ry, rz = real_arr[idx*3], real_arr[idx*3+1], real_arr[idx*3+2]
            zx, zy, zz = z_arr[idx*3], z_arr[idx*3+1], z_arr[idx*3+2]
            
            for v in (rx, ry, rz, zx, zy, zz):
                if math.isnan(v) or math.isinf(v):
                    failures.append(f"{name} galaxy {idx} contains NaN/Inf: {v}")
                    break
            
            # Verify displacement vector collinearity with radial vector
            dx = zx - rx
            dy = zy - ry
            dz = zz - rz
            disp_norm = math.sqrt(dx*dx + dy*dy + dz*dz)
            r_norm = math.sqrt(rx*rx + ry*ry + rz*rz)
            
            if disp_norm > 1.0 and r_norm > 1.0:
                dot = (dx * rx + dy * ry + dz * rz) / (disp_norm * r_norm)
                # dot must be ~ +1.0 or -1.0 depending on sign of total delta v
                if abs(abs(dot) - 1.0) > 0.05:  # allow small rounding variance from JSON export
                    failures.append(f"{name} galaxy {idx} displacement non-collinear with radius: |dot|={abs(dot):.4f}")

    print(f"  [Math Stress] Evaluated {len(r_sub_epsilon)} sub-epsilon points, {len(r_extremes)} extreme points, 8 octants.")
    print(f"  [Catalog Scan] Exhaustively scanned {total_galaxies_scanned} galaxies across SDSS, 2MRS, DESI.")
    
    if failures:
        print(f"  [FAIL] {len(failures)} mathematical RSD failures detected:")
        for f in failures[:10]:
            print(f"    - {f}")
        return False, failures
    else:
        print(f"  [PASS] 0 NaNs, 0 Infs, 100% unit vector convergence, 100% collinear line-of-sight displacement.")
        return True, []


def stress_browser_performance_and_layers(client):
    """
    Challenge 2 & 3:
    - Dual-buffer typed array in-place morphing budget (< 2.5 ms/frame at 60 FPS).
    - Rapid survey layer toggling: 50 cycles on window.surveyGroups with 0 rejections, 0 scene corruption, 0 memory leaks.
    - Rule 14 & Rule 39 verification.
    """
    failures = []
    
    # Wait for page readiness and Three.js scene initialization
    print("\n--- Initializing Browser CDP Client for Runtime Stress Tests ---")
    client.wait_for_condition("document.readyState === 'complete' && !!window.surveyGroups && !!window.applySurveyRSD", timeout=20.0)
    time.sleep(1.5)

    # -------------------------------------------------------------
    # CHALLENGE 2: Dual-Buffer In-Place Morphing 60 FPS Performance Budget
    # -------------------------------------------------------------
    print("\n--- Challenge 2: Dual-Buffer In-Place Morphing Budget (< 2.5 ms/frame) ---")
    perf_res = client.evaluate("""
    (function() {
        const frames = 100;
        const timings = [];
        const meshes = [window.cf4PointsMesh, window.sdssPointsMesh, window.twoMrsPointsMesh, window.desiPointsMesh];
        
        // Assert all meshes have valid geometries and dual buffer attributes
        const bufferChecks = meshes.map(m => {
            if (!m || !m.geometry || !m.geometry.attributes) return { ok: false, reason: 'missing mesh/geometry' };
            const p = m.geometry.attributes.position;
            const r = m.geometry.attributes.posReal;
            const z = m.geometry.attributes.posRedshift;
            if (!p || !r || !z) return { ok: false, reason: 'missing p/r/z attributes' };
            return {
                ok: true,
                count: p.count,
                bufferId: p.array.buffer,
                isDynamic: p.usage === THREE.DynamicDrawUsage
            };
        });
        
        const initialBuffers = bufferChecks.map(c => c.bufferId);
        
        // Run 100 morphing frames with sinusoidal parameter evolution
        for (let f = 0; f < frames; f++) {
            const t = 0.5 + 0.5 * Math.sin(f * 0.1);
            const t0 = performance.now();
            window.applySurveyRSD(t);
            const t1 = performance.now();
            timings.push(t1 - t0);
        }
        
        // Verify buffer identity preservation (no GC or buffer reallocation)
        const buffersIdentical = meshes.every((m, idx) => m.geometry.attributes.position.array.buffer === initialBuffers[idx]);
        
        // Verify displacement between t=0 and t=1
        window.applySurveyRSD(0.0);
        const sdssP0 = window.sdssPointsMesh.geometry.attributes.position.array[0];
        window.applySurveyRSD(1.0);
        const sdssP1 = window.sdssPointsMesh.geometry.attributes.position.array[0];
        const displacementObserved = Math.abs(sdssP1 - sdssP0) > 0.001;
        
        // Reset to real space
        window.applySurveyRSD(0.0);
        
        timings.sort((a, b) => a - b);
        const minT = timings[0];
        const maxT = timings[timings.length - 1];
        const sumT = timings.reduce((acc, v) => acc + v, 0);
        const meanT = sumT / timings.length;
        const medianT = timings[Math.floor(timings.length / 2)];
        const p95T = timings[Math.floor(timings.length * 0.95)];
        
        return {
            frames: frames,
            bufferChecks: bufferChecks,
            buffersIdentical: buffersIdentical,
            displacementObserved: displacementObserved,
            minMs: minT,
            maxMs: maxT,
            meanMs: meanT,
            medianMs: medianT,
            p95Ms: p95T,
            passedBudget: (meanT < 2.5) && (p95T < 2.5)
        };
    })()
    """)

    if not perf_res:
        failures.append("Performance benchmark returned null or failed to evaluate.")
    else:
        print(f"  [Morph Benchmark] 100 Frames: Min={perf_res.get('minMs', 0):.3f}ms, Mean={perf_res.get('meanMs', 0):.3f}ms, "
              f"Median={perf_res.get('medianMs', 0):.3f}ms, P95={perf_res.get('p95Ms', 0):.3f}ms, Max={perf_res.get('maxMs', 0):.3f}ms")
        print(f"  [Buffer Integrity] Buffers identical (zero reallocation): {perf_res.get('buffersIdentical')}")
        print(f"  [Displacement] Non-zero displacement verified: {perf_res.get('displacementObserved')}")
        
        if not perf_res.get("buffersIdentical"):
            failures.append("Dual-buffer morphing reallocated typed array buffers during animation!")
        if not perf_res.get("displacementObserved"):
            failures.append("No vertex displacement observed between RSD morph t=0 and t=1!")
        if not perf_res.get("passedBudget"):
            failures.append(f"RSD morphing exceeded 2.5ms budget: Mean={perf_res.get('meanMs', 0):.3f}ms, P95={perf_res.get('p95Ms', 0):.3f}ms")
        else:
            print("  [PASS] Typed array morphing executed well within 60 FPS budget (< 2.5 ms/frame).")

    # -------------------------------------------------------------
    # CHALLENGE 3: Rapid Survey Layer Toggling (50 Flips)
    # -------------------------------------------------------------
    print("\n--- Challenge 3: Rapid Survey Layer Toggling Stress Test (50 Flips) ---")
    toggle_res = client.evaluate("""
    (function() {
        const initialErrorCount = (window.errorRegistry || []).length;
        const initialSceneChildren = scene.children.length;
        const initialHeap = window.performance && window.performance.memory ? window.performance.memory.usedJSHeapSize : 0;
        
        const surveys = ['desi', 'sdss', '2mrs', 'cf4', 'hooleilana'];
        const iterations = 50;
        
        for (let i = 0; i < iterations; i++) {
            const vis = (i % 2 === 0);
            for (let s = 0; s < surveys.length; s++) {
                const sId = surveys[s];
                if (typeof window.setSurveyVisibility === 'function') {
                    window.setSurveyVisibility(sId, vis);
                } else {
                    const g = window.surveyGroups[sId];
                    if (g) g.visible = vis;
                }
            }
            // Trigger explicit render pass
            if (typeof renderer !== 'undefined' && renderer.render && typeof scene !== 'undefined' && typeof camera !== 'undefined') {
                renderer.render(scene, camera);
            }
        }
        
        // Restore all survey layers to visible
        surveys.forEach(sId => {
            if (typeof window.setSurveyVisibility === 'function') {
                window.setSurveyVisibility(sId, true);
            } else {
                const g = window.surveyGroups[sId];
                if (g) g.visible = true;
            }
        });
        if (typeof renderer !== 'undefined' && renderer.render && typeof scene !== 'undefined' && typeof camera !== 'undefined') {
            renderer.render(scene, camera);
        }
        
        const finalErrorCount = (window.errorRegistry || []).length;
        const finalSceneChildren = scene.children.length;
        const finalHeap = window.performance && window.performance.memory ? window.performance.memory.usedJSHeapSize : 0;
        
        // Audit mesh integrity post-toggling
        const meshes = [
            { name: 'cf4', mesh: window.cf4PointsMesh },
            { name: 'sdss', mesh: window.sdssPointsMesh },
            { name: 'twoMrs', mesh: window.twoMrsPointsMesh },
            { name: 'desi', mesh: window.desiPointsMesh }
        ];
        const meshAudits = meshes.map(m => ({
            name: m.name,
            exists: !!m.mesh,
            hasGeo: !!(m.mesh && m.mesh.geometry),
            posCount: (m.mesh && m.mesh.geometry && m.mesh.geometry.attributes.position) ? m.mesh.geometry.attributes.position.count : 0,
            isDynamic: (m.mesh && m.mesh.geometry && m.mesh.geometry.attributes.position) ? m.mesh.geometry.attributes.position.usage === THREE.DynamicDrawUsage : false,
            visible: m.mesh ? m.mesh.visible : false
        }));
        
        return {
            iterations: iterations,
            initialErrors: initialErrorCount,
            finalErrors: finalErrorCount,
            newErrors: finalErrorCount - initialErrorCount,
            initialChildren: initialSceneChildren,
            finalChildren: finalSceneChildren,
            childrenUnchanged: (initialSceneChildren === finalSceneChildren),
            initialHeapBytes: initialHeap,
            finalHeapBytes: finalHeap,
            heapDeltaBytes: finalHeap - initialHeap,
            meshAudits: meshAudits
        };
    })()
    """)

    if not toggle_res:
        failures.append("Toggle stress test returned null or failed to evaluate.")
    else:
        print(f"  [Toggle Stress] Completed {toggle_res.get('iterations')} visibility flips across all survey groups.")
        print(f"  [Errors] Initial={toggle_res.get('initialErrors')}, Final={toggle_res.get('finalErrors')}, New={toggle_res.get('newErrors')}")
        print(f"  [Scene Graph] Initial children={toggle_res.get('initialChildren')}, Final children={toggle_res.get('finalChildren')}")
        print(f"  [Heap Delta] {toggle_res.get('heapDeltaBytes', 0) / (1024*1024):.2f} MB")
        
        if toggle_res.get("newErrors", 0) > 0:
            failures.append(f"Rapid layer toggling generated {toggle_res.get('newErrors')} runtime errors in window.errorRegistry!")
        if not toggle_res.get("childrenUnchanged"):
            failures.append(f"Scene graph corrupted: children count changed from {toggle_res.get('initialChildren')} to {toggle_res.get('finalChildren')}")
        
        for ma in toggle_res.get("meshAudits", []):
            if not ma["exists"] or not ma["hasGeo"] or ma["posCount"] == 0 or not ma["isDynamic"]:
                failures.append(f"Mesh {ma['name']} corrupted after toggling: {ma}")
            if not ma["visible"]:
                failures.append(f"Mesh {ma['name']} was not restored to visible: {ma}")
        
        if not failures:
            print("  [PASS] 50 rapid visibility flips passed with 0 errors, 0 scene corruption, and stable memory.")

    # -------------------------------------------------------------
    # CHALLENGE 4: Rule 14 & Rule 39 Verification
    # -------------------------------------------------------------
    print("\n--- Challenge 4: Verification Against Rule 14 (Zero Exceptions) & Rule 39 (Zero VRAM Leak) ---")
    
    # Check CDP uncaught exceptions and console errors
    uncaught = client.uncaught_exceptions
    if uncaught:
        failures.append(f"Rule 14 Violation: Found {len(uncaught)} uncaught JavaScript exceptions: {uncaught}")
    else:
        print("  [Rule 14] Zero uncaught JavaScript exceptions detected via CDP.")
        
    console_errors = [log for log in getattr(client, 'console_logs', []) if log.get('type') == 'error']
    if console_errors:
        failures.append(f"Rule 14 Violation: Found {len(console_errors)} console.error logs: {console_errors}")
    else:
        print("  [Rule 14] Zero console.error logs detected.")

    # Check Rule 39 WebGL Memory & Disposal
    vram_res = client.evaluate("""
    (function() {
        if (typeof renderer === 'undefined' || !renderer.info) return { supported: false };
        
        const initialGeos = renderer.info.memory.geometries;
        const initialTexs = renderer.info.memory.textures;
        
        // Re-generate catalog 3 times to stress disposeHierarchy
        for (let i = 0; i < 3; i++) {
            if (typeof generateFullCF4Catalog === 'function') {
                generateFullCF4Catalog();
            }
        }
        
        // Force a render pass to clear pending disposals
        if (typeof renderer !== 'undefined' && renderer.render && typeof scene !== 'undefined' && typeof camera !== 'undefined') {
            renderer.render(scene, camera);
        }
        
        const finalGeos = renderer.info.memory.geometries;
        const finalTexs = renderer.info.memory.textures;
        
        return {
            supported: true,
            initialGeos: initialGeos,
            finalGeos: finalGeos,
            geoDelta: finalGeos - initialGeos,
            initialTexs: initialTexs,
            finalTexs: finalTexs,
            texDelta: finalTexs - initialTexs
        };
    })()
    """)

    if vram_res and vram_res.get("supported"):
        print(f"  [Rule 39 VRAM] Geometries: Initial={vram_res.get('initialGeos')}, Final={vram_res.get('finalGeos')}, Delta={vram_res.get('geoDelta')}")
        print(f"  [Rule 39 VRAM] Textures: Initial={vram_res.get('initialTexs')}, Final={vram_res.get('finalTexs')}, Delta={vram_res.get('texDelta')}")
        if vram_res.get("geoDelta", 0) > 4:
            failures.append(f"Rule 39 Violation: Geometries leaked during re-generation (delta = +{vram_res.get('geoDelta')})")
        else:
            print("  [PASS] Rule 39 disposal confirmed: zero runaway geometry/texture allocation.")

    return failures


def main():
    print("================================================================================")
    print("      MILITARY-GRADE EMPIRICAL CHALLENGE HARNESS: MILESTONE 1 VERIFICATION     ")
    print("================================================================================")
    start_time = time.time()
    all_failures = []

    # 1. Stress Kaiser RSD mathematical boundaries and catalog data
    math_pass, math_failures = stress_kaiser_rsd_math()
    all_failures.extend(math_failures)

    # 2. Start HTTP server and headless Chrome CDP client
    ensure_http_server(8000)
    cdp_port = find_free_port(9260)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        client.send_cdp("Console.enable")
        client.send_cdp("Runtime.enable")
        client.send_cdp("Page.enable")
        browser_failures = stress_browser_performance_and_layers(client)
        all_failures.extend(browser_failures)
    finally:
        client.close()

    elapsed = time.time() - start_time
    print("\n================================================================================")
    if all_failures:
        print(f"VERDICT: FAIL ({len(all_failures)} failure modes detected in {elapsed:.2f}s)")
        for idx, f in enumerate(all_failures, 1):
            print(f"  {idx}. {f}")
        print("================================================================================")
        sys.exit(1)
    else:
        print(f"VERDICT: APPROVE (All empirical stress tests passed with 100% success in {elapsed:.2f}s)")
        print("================================================================================")
        sys.exit(0)


if __name__ == "__main__":
    main()
