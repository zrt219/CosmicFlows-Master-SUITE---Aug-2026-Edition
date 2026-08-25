# -*- coding: utf-8 -*-
"""
Milestone M1 Adversarial Stress Test Harness
Challenger 2: Engine Kinematics & Stress Verifier

Tests:
1. 7 Cosmographic Engines Kinematic Differentiability Matrix across 30 spatial sample points.
2. 100 Ultra-rapid engine switches + theme toggle interleaving with zero console errors and zero WebGL context losses.
3. Dynamic parameter sweeps and randomized fuzzing (gamma, soften, turb, dt, linesCount, colormaps).
4. Cash-Karp RK45 singularity & extreme boundary box stress testing.
5. Physical divergence and vorticity vector calculus verification across 8 critical hubs.
"""

import subprocess
import time
import json
import urllib.request
import websocket
import sys
import os
import http.server
import socketserver
import threading
import base64
import math
import random

PORT = 8035
DEBUG_PORT = 9249

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

def start_server():
    os.chdir(r"c:\Users\Zhane\Documents\antigravity\gallant-newton")
    httpd = socketserver.TCPServer(("127.0.0.1", PORT), QuietHandler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    return httpd

def run_stress_suite():
    print("=" * 80, flush=True)
    print("  CHALLENGER 2: M1 ADVERSARIAL ENGINE KINEMATICS & STRESS TEST HARNESS", flush=True)
    print("=" * 80, flush=True)
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}", flush=True)
    print(f"Serving HTTP on port {PORT}, Chrome Remote Debugging on port {DEBUG_PORT}", flush=True)

    httpd = start_server()
    chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    user_data = f"C:\\Users\\Zhane\\AppData\\Local\\Temp\\chrome_m1_chal2_{int(time.time())}"
    cmd = [
        chrome_path,
        "--headless=new",
        f"--remote-debugging-port={DEBUG_PORT}",
        "--remote-debugging-address=0.0.0.0",
        "--remote-allow-origins=*",
        f"--user-data-dir={user_data}",
        "--window-size=1920,1080",
        "--enable-webgl",
        "--use-gl=angle",
        "--no-first-run",
        "--no-default-browser-check",
        "--disable-background-networking",
        f"http://localhost:{PORT}/index.html"
    ]

    proc = subprocess.Popen(cmd)
    time.sleep(2.5)

    console_errors = []
    console_warnings = []
    exceptions = []

    try:
        req = urllib.request.urlopen(f"http://localhost:{DEBUG_PORT}/json/list")
        pages = json.loads(req.read().decode())
        target = next((p for p in pages if f"{PORT}" in p.get("url", "")), None)
        if not target:
            raise RuntimeError(f"Target page on port {PORT} not found in Chrome list: {pages}")

        ws = websocket.create_connection(target["webSocketDebuggerUrl"], timeout=30.0)
        msg_id = 0

        def send_cmd(method, params=None):
            nonlocal msg_id
            msg_id += 1
            payload = {"id": msg_id, "method": method}
            if params:
                payload["params"] = params
            ws.send(json.dumps(payload))
            while True:
                resp = json.loads(ws.recv())
                if resp.get("method") == "Runtime.consoleAPICalled":
                    args = resp.get("params", {}).get("args", [])
                    t = resp.get("params", {}).get("type", "")
                    text = " ".join(str(a.get("value", a.get("description", ""))) for a in args)
                    if t == "error":
                        console_errors.append(text)
                    elif t == "warning":
                        console_warnings.append(text)
                elif resp.get("method") == "Runtime.exceptionThrown":
                    exc_details = resp.get("params", {}).get("exceptionDetails", {})
                    exc_msg = exc_details.get("text", "")
                    if "exception" in exc_details:
                        exc_msg += " " + str(exc_details["exception"].get("description", ""))
                    exceptions.append(exc_msg)
                if resp.get("id") == msg_id:
                    if "error" in resp:
                        raise RuntimeError(f"CDP error in {method}: {resp['error']}")
                    return resp.get("result", {})

        def eval_js(expr):
            res = send_cmd("Runtime.evaluate", {"expression": expr, "returnByValue": True, "awaitPromise": True})
            if "exceptionDetails" in res:
                exc = res["exceptionDetails"]
                exc_str = f"JS Exception: {exc.get('text', '')} {exc.get('exception', {}).get('description', '')}"
                exceptions.append(exc_str)
                raise RuntimeError(f"Evaluation failed in '{expr[:100]}...': {exc_str}")
            return res.get("result", {}).get("value")

        send_cmd("Page.enable")
        send_cmd("Runtime.enable")
        send_cmd("DOM.enable")
        time.sleep(1.0)

        results = {
            "differentiability": False,
            "rapid_switching_100x": False,
            "parameter_sweeps": False,
            "rk45_singularity_bounds": False,
            "vector_calculus_hubs": False,
            "webgl_context_stable": False,
            "console_exceptions_zero": False
        }

        # -------------------------------------------------------------
        # 1. API Bridge Integrity Check
        # -------------------------------------------------------------
        print("\n[TEST 1/6] Verifying API Bridge Interface...", flush=True)
        api_check = eval_js("""
            (() => {
                return {
                    defined: typeof window.cosmicflows !== 'undefined',
                    version: window.cosmicflows ? window.cosmicflows.version : null,
                    hasVel: typeof window.cosmicflows.getVelocitySg === 'function',
                    hasDiv: typeof window.cosmicflows.getDivergenceSg === 'function',
                    hasVort: typeof window.cosmicflows.getVorticitySg === 'function',
                    hasRK4: typeof window.cosmicflows.integrateRK4 === 'function',
                    hasSwitch: typeof window.cosmicflows.switchScienceEngine === 'function',
                    hasTheme: typeof window.cosmicflows.toggleTheme === 'function',
                    hasSimState: typeof window.cosmicflows.simState === 'object',
                    engines: Object.keys(window.cosmicflows.ENGINES || {})
                };
            })()
        """)
        assert api_check["defined"], "window.cosmicflows is missing!"
        assert api_check["version"] == "2026.1", f"Unexpected version: {api_check['version']}"
        assert api_check["hasVel"] and api_check["hasDiv"] and api_check["hasVort"] and api_check["hasRK4"]
        print(f"  -> API bridge verified: version={api_check['version']}, engines={api_check['engines']}", flush=True)

        # -------------------------------------------------------------
        # 2. 7 Cosmographic Engines Differentiability Matrix
        # -------------------------------------------------------------
        print("\n[TEST 2/6] Evaluating 7 Cosmographic Engines Differentiability Matrix...", flush=True)
        engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026']

        diff_data = eval_js("""
            (() => {
                const engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026'];
                const testPoints = [
                    { name: 'Local Group Origin', p: new THREE.Vector3(0, 0, 0) },
                    { name: 'Shapley Core', p: new THREE.Vector3(7200, -8600, -2400) },
                    { name: 'Laniakea Core GA', p: new THREE.Vector3(15000, -14000, 6800) },
                    { name: 'Centaurus Cluster', p: new THREE.Vector3(-4200, 1200, 3100) },
                    { name: 'Coma Supercluster', p: new THREE.Vector3(500, 7000, 1500) },
                    { name: 'Perseus-Pisces Spine', p: new THREE.Vector3(4500, -3000, 0) },
                    { name: 'Southern Infall Sink', p: new THREE.Vector3(-4500, -5000, -14500) },
                    { name: 'Vela ZoA Core', p: new THREE.Vector3(-8500, -12000, -3200) },
                    { name: 'Dipole Repeller Void', p: new THREE.Vector3(-10000, 10000, 12000) },
                    { name: 'Sloan Great Wall Basin', p: new THREE.Vector3(12000, 14000, 11000) },
                    { name: 'ZoA Corridor Midpoint', p: new THREE.Vector3(-6350, -5400, -50) },
                    { name: 'ZoA Filament Bridge', p: new THREE.Vector3(-7000, -8000, -1500) },
                    { name: 'Perseus Axis Point', p: new THREE.Vector3(8500, -4500, 2500) },
                    { name: 'Void Deep Interior', p: new THREE.Vector3(-12000, 12000, 13000) },
                    { name: 'Radial Shell 4 Point', p: new THREE.Vector3(-6200, -7400, 2600) },
                    { name: 'Radial Shell 2 Point', p: new THREE.Vector3(2000, 2000, 2000) },
                    { name: 'Corner Positive', p: new THREE.Vector3(14000, 14000, 14000) },
                    { name: 'Corner Negative', p: new THREE.Vector3(-14000, -14000, -14000) },
                    { name: 'Midplane X', p: new THREE.Vector3(5000, 0, 0) },
                    { name: 'Midplane Y', p: new THREE.Vector3(0, -6000, 0) },
                    { name: 'Midplane Z', p: new THREE.Vector3(0, 0, 8000) },
                    { name: 'Virial Shapley Halo', p: new THREE.Vector3(7200 + 1500, -8600, -2400) },
                    { name: 'Virial Coma Halo', p: new THREE.Vector3(500, 7000 + 1200, 1500) },
                    { name: 'Virial Vela Halo', p: new THREE.Vector3(-8500 + 1000, -12000, -3200) },
                    { name: 'Arbitrary Field 1', p: new THREE.Vector3(-3500, 4500, -6500) },
                    { name: 'Arbitrary Field 2', p: new THREE.Vector3(8500, -3500, 9500) },
                    { name: 'Arbitrary Field 3', p: new THREE.Vector3(-11000, -2000, 4000) },
                    { name: 'Arbitrary Field 4', p: new THREE.Vector3(6000, 9000, -8000) },
                    { name: 'Arbitrary Field 5', p: new THREE.Vector3(-1500, -11500, 7500) },
                    { name: 'Arbitrary Field 6', p: new THREE.Vector3(10500, -7500, -10500) }
                ];

                const engineVelocities = {};
                for (let eng of engines) {
                    window.cosmicflows.switchScienceEngine(eng);
                    engineVelocities[eng] = [];
                    for (let pt of testPoints) {
                        const v = new THREE.Vector3();
                        window.cosmicflows.getVelocitySg(pt.p, v);
                        engineVelocities[eng].push({
                            name: pt.name,
                            x: v.x, y: v.y, z: v.z,
                            spd: v.length(),
                            isFinite: isFinite(v.x) && isFinite(v.y) && isFinite(v.z)
                        });
                    }
                }
                return { pointsCount: testPoints.length, velocities: engineVelocities };
            })()
        """)

        num_points = diff_data["pointsCount"]
        vels = diff_data["velocities"]

        # Check finiteness across all points and engines
        for eng in engines:
            for pt in vels[eng]:
                assert pt["isFinite"], f"Engine {eng} produced non-finite velocity at {pt['name']}: ({pt['x']}, {pt['y']}, {pt['z']})"

        # Compute Pairwise Difference Matrix across all 21 pairs
        pair_distances = {}
        min_pairwise_distance = 1e9
        print(f"\n  Pairwise Mean Velocity Differences across {num_points} 3D locations (km/s):", flush=True)
        print(f"  {'Engine Pair':<35} | {'Mean Delta V':<12} | {'Max Delta V':<12} | {'Min Delta V':<12}", flush=True)
        print("  " + "-" * 75, flush=True)

        for i in range(len(engines)):
            for j in range(i + 1, len(engines)):
                e1, e2 = engines[i], engines[j]
                vlist1, vlist2 = vels[e1], vels[e2]
                deltas = []
                for k in range(num_points):
                    p1, p2 = vlist1[k], vlist2[k]
                    d = math.sqrt((p1['x'] - p2['x'])**2 + (p1['y'] - p2['y'])**2 + (p1['z'] - p2['z'])**2)
                    deltas.append(d)
                mean_d = sum(deltas) / len(deltas)
                max_d = max(deltas)
                min_d = min(deltas)
                pair_name = f"{e1} vs {e2}"
                pair_distances[pair_name] = mean_d
                min_pairwise_distance = min(min_pairwise_distance, mean_d)
                print(f"  {pair_name:<35} | {mean_d:<12.2f} | {max_d:<12.2f} | {min_d:<12.2f}", flush=True)
                assert mean_d >= 10.0, f"Engines {e1} and {e2} are mathematically indistinguishable! Mean Delta = {mean_d:.2f}"

        print(f"\n  -> Minimum pairwise engine distinction across all pairs: {min_pairwise_distance:.2f} km/s (Threshold: >= 10.0 km/s)", flush=True)
        results["differentiability"] = True

        # -------------------------------------------------------------
        # 3. 100 Ultra-Rapid Engine Switches Stress Test
        # -------------------------------------------------------------
        print("\n[TEST 3/6] Running 100 Ultra-Rapid Engine Switches with Theme Interleaving...", flush=True)
        t0_switch = time.time()

        switch_res = eval_js("""
            (() => {
                const engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026'];
                const colormaps = ['ink', 'doppler', 'divergence', 'vorticity', 'velocity', 'vweb'];
                const themes = ['white', 'dark'];
                const history = [];

                // Use 250 streamlines during 100 rapid switches
                window.cosmicflows.simState.linesCount = 250;

                // 1. Sequential Round-Robin Switches (50 switches)
                for (let i = 0; i < 50; i++) {
                    const eng = engines[i % engines.length];
                    window.cosmicflows.switchScienceEngine(eng);
                    if (i % 5 === 0) {
                        window.cosmicflows.toggleTheme(themes[i % 2]);
                    }
                    if (i % 7 === 0) {
                        window.cosmicflows.simState.colormap = colormaps[i % colormaps.length];
                        const selCol = document.getElementById('sel-colormap');
                        if (selCol) selCol.value = colormaps[i % colormaps.length];
                    }
                    history.push(eng);
                }

                // 2. Randomized Fuzz Switches (50 switches)
                for (let i = 0; i < 50; i++) {
                    const randEng = engines[Math.floor(Math.random() * engines.length)];
                    window.cosmicflows.switchScienceEngine(randEng);
                    if (Math.random() < 0.2) {
                        window.cosmicflows.toggleTheme();
                    }
                    history.push(randEng);
                }

                // Restore full density standard state
                window.cosmicflows.simState.linesCount = 2600;
                window.cosmicflows.switchScienceEngine('cf4-wf');
                window.cosmicflows.toggleTheme('white');

                // Check WebGL Context Status
                const gl = window.cosmicflows.renderer.getContext();
                const isLost = gl ? gl.isContextLost() : false;

                return {
                    totalSwitches: history.length,
                    isContextLost: isLost,
                    currentTurb: window.cosmicflows.simState.turb,
                    currentSoften: window.cosmicflows.simState.soften,
                    currentGamma: window.cosmicflows.simState.gamma
                };
            })()
        """)

        switch_time = time.time() - t0_switch
        print(f"  -> Completed {switch_res['totalSwitches']} engine switches in {switch_time:.2f}s ({switch_res['totalSwitches'] / switch_time:.1f} switches/sec)", flush=True)
        print(f"  -> WebGL context lost: {switch_res['isContextLost']}", flush=True)
        print(f"  -> Active parameters synchronized: turb={switch_res['currentTurb']}, soften={switch_res['currentSoften']}, gamma={switch_res['currentGamma']}", flush=True)
        assert not switch_res["isContextLost"], "WebGL context was lost during rapid engine switching!"
        assert switch_res["totalSwitches"] == 100, f"Expected 100 switches, got {switch_res['totalSwitches']}"
        results["rapid_switching_100x"] = True
        results["webgl_context_stable"] = not switch_res["isContextLost"]

        # -------------------------------------------------------------
        # 4. Dynamic Parameter Sweeps & Randomized Fuzzing
        # -------------------------------------------------------------
        print("\n[TEST 4/6] Running Dynamic Parameter Sweeps & Fuzzing...", flush=True)
        param_fuzz_res = eval_js("""
            (() => {
                const engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026'];
                const soften_vals = [10, 50, 200, 500, 1000, 2000, 5000, 10000, 50000];
                const turb_vals = [0, 1, 5, 15, 30, 60, 120];
                const gamma_vals = [0.1, 0.5, 1.0, 1.55, 2.5, 5.0];

                let allFinite = true;
                let evalCount = 0;

                // Systematic sweep across parameter extremes on field evaluation methods
                for (let eng of engines) {
                    window.cosmicflows.switchScienceEngine(eng);
                    for (let s of soften_vals) {
                        window.cosmicflows.simState.soften = s;
                        for (let g of gamma_vals) {
                            window.cosmicflows.simState.gamma = g;
                            for (let t of turb_vals) {
                                window.cosmicflows.simState.turb = t;
                                const v = new THREE.Vector3();
                                window.cosmicflows.getVelocitySg(new THREE.Vector3(1000, -2000, 3000), v);
                                const div = window.cosmicflows.getDivergenceSg(new THREE.Vector3(1000, -2000, 3000));
                                const vort = new THREE.Vector3();
                                window.cosmicflows.getVorticitySg(new THREE.Vector3(1000, -2000, 3000), vort);

                                evalCount += 3;
                                if (!isFinite(v.x) || !isFinite(v.y) || !isFinite(v.z) || !isFinite(div) || !isFinite(vort.length())) {
                                    allFinite = false;
                                }
                            }
                        }
                    }
                }

                // 80 iterations of randomized fuzzing with RK45 integration
                let rk4AllBounded = true;
                for (let i = 0; i < 80; i++) {
                    const eng = engines[i % engines.length];
                    window.cosmicflows.switchScienceEngine(eng);
                    window.cosmicflows.simState.soften = 50 + Math.random() * 5000;
                    window.cosmicflows.simState.turb = Math.random() * 80;
                    window.cosmicflows.simState.gamma = 0.2 + Math.random() * 3.0;
                    const dt = 5 + Math.random() * 110;

                    const rx = (Math.random() * 2 - 1) * 14000;
                    const ry = (Math.random() * 2 - 1) * 14000;
                    const rz = (Math.random() * 2 - 1) * 14000;
                    const pts = window.cosmicflows.integrateRK4(new THREE.Vector3(rx, ry, rz), 1, 50, dt);

                    for (let pt of pts) {
                        if (Math.abs(pt.x) > 15000.1 || Math.abs(pt.y) > 15000.1 || Math.abs(pt.z) > 15000.1) {
                            rk4AllBounded = false;
                        }
                        if (!isFinite(pt.x) || !isFinite(pt.y) || !isFinite(pt.z)) {
                            allFinite = false;
                        }
                    }
                }

                // Reset standard parameters
                window.cosmicflows.switchScienceEngine('cf4-wf');

                return {
                    evalCount: evalCount,
                    allFinite: allFinite,
                    rk4AllBounded: rk4AllBounded
                };
            })()
        """)

        print(f"  -> Parameter sweep executed {param_fuzz_res['evalCount']} field evaluations across all 7 engines.", flush=True)
        print(f"  -> All evaluation outputs finite: {param_fuzz_res['allFinite']}", flush=True)
        print(f"  -> All randomized RK45 trajectories bounded in [-15000, 15000]: {param_fuzz_res['rk4AllBounded']}", flush=True)
        assert param_fuzz_res["allFinite"], "Non-finite values encountered during parameter sweep!"
        assert param_fuzz_res["rk4AllBounded"], "RK45 trajectory exceeded bounding box during fuzzing!"
        results["parameter_sweeps"] = True

        # -------------------------------------------------------------
        # 5. Cash-Karp RK45 Numerical Singularity & Boundary Stress
        # -------------------------------------------------------------
        print("\n[TEST 5/6] Testing Cash-Karp RK45 Singularity & Boundary Extremes...", flush=True)
        rk4_stress = eval_js("""
            (() => {
                const singularities = [
                    { name: 'Shapley Center', p: new THREE.Vector3(7200, -8600, -2400) },
                    { name: 'Great Attractor Center', p: new THREE.Vector3(15000, -14000, 6800) },
                    { name: 'Virgo Origin Center', p: new THREE.Vector3(0, 0, 0) },
                    { name: 'Coma Center', p: new THREE.Vector3(500, 7000, 1500) },
                    { name: 'Dipole Repeller Center', p: new THREE.Vector3(-10000, 10000, 12000) }
                ];

                const boundaries = [
                    { name: 'Box Face +X', p: new THREE.Vector3(15000, 0, 0) },
                    { name: 'Box Face -X', p: new THREE.Vector3(-15000, 0, 0) },
                    { name: 'Box Face +Y', p: new THREE.Vector3(0, 15000, 0) },
                    { name: 'Box Face -Y', p: new THREE.Vector3(0, -15000, 0) },
                    { name: 'Box Face +Z', p: new THREE.Vector3(0, 0, 15000) },
                    { name: 'Box Face -Z', p: new THREE.Vector3(0, 0, -15000) },
                    { name: 'Far Outside Domain', p: new THREE.Vector3(35000, -40000, 50000) }
                ];

                const singResults = [];
                for (let s of singularities) {
                    const fwd = window.cosmicflows.integrateRK4(s.p, 1, 100, 46);
                    const bwd = window.cosmicflows.integrateRK4(s.p, -1, 100, 46);
                    singResults.push({
                        name: s.name,
                        fwdLen: fwd.length,
                        bwdLen: bwd.length,
                        valid: fwd.every(pt => isFinite(pt.x) && isFinite(pt.y) && isFinite(pt.z)) &&
                               bwd.every(pt => isFinite(pt.x) && isFinite(pt.y) && isFinite(pt.z))
                    });
                }

                const bndResults = [];
                for (let b of boundaries) {
                    const fwd = window.cosmicflows.integrateRK4(b.p, 1, 100, 46);
                    bndResults.push({
                        name: b.name,
                        len: fwd.length,
                        allWithinBox: fwd.every(pt => Math.abs(pt.x) <= 15000.1 && Math.abs(pt.y) <= 15000.1 && Math.abs(pt.z) <= 15000.1)
                    });
                }

                return { singResults: singResults, bndResults: bndResults };
            })()
        """)

        for s in rk4_stress["singResults"]:
            print(f"  -> Singularity at {s['name']}: fwd_pts={s['fwdLen']}, bwd_pts={s['bwdLen']}, finite={s['valid']}", flush=True)
            assert s["valid"], f"Non-finite RK45 integration at singularity {s['name']}"

        for b in rk4_stress["bndResults"]:
            print(f"  -> Boundary at {b['name']}: pts={b['len']}, within_box={b['allWithinBox']}", flush=True)
            assert b["allWithinBox"], f"Boundary clipping violated at {b['name']}"

        results["rk45_singularity_bounds"] = True

        # -------------------------------------------------------------
        # 6. Physical Divergence & Vorticity Across All Attractors & Voids
        # -------------------------------------------------------------
        print("\n[TEST 6/6] Testing Physical Vector Calculus across Landmarks & Engines...", flush=True)
        vec_calc = eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const landmarks = [
                    { id: 'shapley', name: 'Shapley Supercluster', pos: new THREE.Vector3(7200, -8600, -2400), expectedDiv: 'negative' },
                    { id: 'laniakea', name: 'Laniakea Core / GA', pos: new THREE.Vector3(15000, -14000, 6800), expectedDiv: 'negative' },
                    { id: 'centaurus', name: 'Centaurus Cluster', pos: new THREE.Vector3(-4200, 1200, 3100), expectedDiv: 'negative' },
                    { id: 'coma', name: 'Coma Supercluster', pos: new THREE.Vector3(500, 7000, 1500), expectedDiv: 'negative' },
                    { id: 'southern', name: 'Southern Infall Sink', pos: new THREE.Vector3(-4500, -5000, -14500), expectedDiv: 'negative' },
                    { id: 'vela', name: 'Vela Supercluster', pos: new THREE.Vector3(-8500, -12000, -3200), expectedDiv: 'negative' },
                    { id: 'sloan', name: 'Sloan Great Wall Basin', pos: new THREE.Vector3(12000, 14000, 11000), expectedDiv: 'negative' },
                    { id: 'dipole', name: 'Dipole Repeller Great Void', pos: new THREE.Vector3(-10000, 10000, 12000), expectedDiv: 'positive' }
                ];

                const divResults = [];
                for (let lm of landmarks) {
                    const div = window.cosmicflows.getDivergenceSg(lm.pos);
                    divResults.push({
                        name: lm.name,
                        div: div,
                        expected: lm.expectedDiv,
                        passed: lm.expectedDiv === 'positive' ? div > 0 : div < 0
                    });
                }

                // Vorticity check: cf4-wf (pure potential) vs gadget4 (N-body halo spin)
                const pt = new THREE.Vector3(7200 + 1200, -8600, -2400);
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const vWF = new THREE.Vector3();
                window.cosmicflows.getVorticitySg(pt, vWF);

                window.cosmicflows.switchScienceEngine('gadget4');
                const vG4 = new THREE.Vector3();
                window.cosmicflows.getVorticitySg(pt, vG4);

                return {
                    divResults: divResults,
                    vortWF: vWF.length(),
                    vortG4: vG4.length()
                };
            })()
        """)

        print(f"\n  Divergence Signatures at Critical Cosmographic Nodes:", flush=True)
        for r in vec_calc["divResults"]:
            print(f"  - {r['name']:<30}: div = {r['div']:<10.6f} (Expected: {r['expected'].upper()}) -> {'PASS' if r['passed'] else 'FAIL'}", flush=True)
            assert r["passed"], f"Divergence sign mismatch for {r['name']}: expected {r['expected']}, got {r['div']}"

        print(f"\n  Vorticity Comparison at Shapley Halo (r = 1200 km/s):", flush=True)
        print(f"  - CF4 WF/CR Potential Flow Vorticity: {vec_calc['vortWF']:.8f} (Expected: ~0, irrotational)", flush=True)
        print(f"  - Gadget-4 Nonlinear Halo Vorticity:   {vec_calc['vortG4']:.6f} (Expected: > 0.01, physical spin)", flush=True)
        assert vec_calc["vortWF"] < 1e-4, f"Potential flow must have near-zero vorticity, got {vec_calc['vortWF']}"
        assert vec_calc["vortG4"] > 0.01, f"Gadget-4 halo must have positive vorticity, got {vec_calc['vortG4']}"
        results["vector_calculus_hubs"] = True

        # -------------------------------------------------------------
        # 7. Final Console Errors & Exceptions Audit
        # -------------------------------------------------------------
        print("\n[AUDIT] Checking Console Logs & JS Exceptions...", flush=True)
        print(f"  Total Console Errors: {len(console_errors)}", flush=True)
        for err in console_errors:
            print(f"    ERR: {err}", flush=True)
        print(f"  Total Uncaught Exceptions: {len(exceptions)}", flush=True)
        for exc in exceptions:
            print(f"    EXC: {exc}", flush=True)

        assert len(console_errors) == 0, f"Detected {len(console_errors)} console errors!"
        assert len(exceptions) == 0, f"Detected {len(exceptions)} uncaught exceptions!"
        results["console_exceptions_zero"] = True

        print("\n" + "=" * 80, flush=True)
        print("  SUMMARY OF EMPIRICAL ADVERSARIAL CHALLENGE FOR M1", flush=True)
        print("=" * 80, flush=True)
        for k, v in results.items():
            status = "PASS" if v else "FAIL"
            print(f"  - {k:<30}: {status}", flush=True)
        print("=" * 80, flush=True)
        print("  FINAL VERDICT: APPROVE (100% Tests Passed without Flaws)", flush=True)
        print("=" * 80 + "\n", flush=True)
        return 0

    finally:
        proc.terminate()
        try:
            httpd.shutdown()
        except Exception:
            pass

if __name__ == "__main__":
    sys.exit(run_stress_suite())
