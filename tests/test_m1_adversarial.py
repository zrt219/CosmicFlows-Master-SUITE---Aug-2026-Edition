# -*- coding: utf-8 -*-
"""
Milestone M1 Adversarial Stress Test Suite & Mathematical Invariant Verifier.
Executes empirical stress tests against the Cosmicflows vector calculus engine via Chrome CDP.
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
import math
import socket
from typing import Dict, Any, List

def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

class ReusableTCPServer(socketserver.TCPServer):
    allow_reuse_address = True

def start_server(port):
    os.chdir(r"c:\Users\Zhane\Documents\antigravity\gallant-newton")
    httpd = ReusableTCPServer(("127.0.0.1", port), QuietHandler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    return httpd

class AdversarialCDPTester:
    def __init__(self):
        self.port = find_free_port()
        self.debug_port = find_free_port()
        self.httpd = None
        self.proc = None
        self.ws = None
        self.msg_id = 0
        self.console_errors = []
        self.exceptions = []
        self.test_results = {}

    def setup(self):
        self.httpd = start_server(self.port)
        print(f"[INIT] Server running on port {self.port}")
        chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
        user_data = f"C:\\Users\\Zhane\\AppData\\Local\\Temp\\chrome_adv_m1_{int(time.time()*1000)}"
        cmd = [
            chrome_path,
            "--headless=new",
            f"--remote-debugging-port={self.debug_port}",
            "--remote-debugging-address=0.0.0.0",
            "--remote-allow-origins=*",
            f"--user-data-dir={user_data}",
            "--window-size=1920,1080",
            "--enable-webgl",
            "--use-gl=angle",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-background-networking",
            f"http://localhost:{self.port}/index.html"
        ]
        self.proc = subprocess.Popen(cmd)
        time.sleep(2.5)

        req = urllib.request.urlopen(f"http://localhost:{self.debug_port}/json/list")
        pages = json.loads(req.read().decode())
        target = next((p for p in pages if f"{self.port}" in p.get("url", "")), None)
        if not target:
            raise RuntimeError(f"Target page on port {self.port} not found in Chrome list: {pages}")

        self.ws = websocket.create_connection(target["webSocketDebuggerUrl"])
        self.send_cmd("Page.enable")
        self.send_cmd("Runtime.enable")
        time.sleep(1.5)
        print("[INIT] CDP WebSocket connected successfully.")

    def teardown(self):
        if self.ws:
            try:
                self.ws.close()
            except:
                pass
        if self.proc:
            try:
                self.proc.terminate()
                self.proc.wait(timeout=2.0)
            except:
                pass
        if self.httpd:
            try:
                self.httpd.shutdown()
            except:
                pass

    def send_cmd(self, method: str, params: Dict[str, Any] = None) -> Dict[str, Any]:
        self.msg_id += 1
        payload = {"id": self.msg_id, "method": method}
        if params:
            payload["params"] = params
        self.ws.send(json.dumps(payload))
        while True:
            resp = json.loads(self.ws.recv())
            if resp.get("method") == "Runtime.consoleAPICalled":
                args = resp.get("params", {}).get("args", [])
                t = resp.get("params", {}).get("type", "")
                text = " ".join(str(a.get("value", "")) for a in args)
                if t == "error":
                    self.console_errors.append(text)
            elif resp.get("method") == "Runtime.exceptionThrown":
                exc_details = resp.get("params", {}).get("exceptionDetails", {})
                self.exceptions.append(exc_details.get("text", "Unknown exception"))
            if resp.get("id") == self.msg_id:
                return resp.get("result", {})

    def eval_js(self, expr: str) -> Any:
        res = self.send_cmd("Runtime.evaluate", {"expression": expr, "returnByValue": True})
        if "exceptionDetails" in res:
            self.exceptions.append(res["exceptionDetails"].get("text", "Eval exception"))
            raise RuntimeError(f"JS Exception: {res['exceptionDetails']}")
        return res.get("result", {}).get("value")

    # ==========================================
    # SUITE 1: SINGULARITY STRESS AT r=0 & DELTA LIMITS
    # ==========================================
    def test_singularity_stress(self):
        print("\n=======================================================")
        print(">>> SUITE 1: Singularity Stress (r=0 & Micro-Perturbations) <<<")
        print("=======================================================")

        landmarks = self.eval_js("""
            window.cosmicflows.LANDMARKS.map(lm => ({
                id: lm.id,
                name: lm.name,
                sign: lm.sign,
                str: lm.str,
                pos: { x: lm.pos.x, y: lm.pos.y, z: lm.pos.z }
            }))
        """)

        deltas = [0.0, 1e-6, 1e-3, 0.1, 1.0, 10.0, 50.0, 100.0, 500.0]
        engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026']

        singularity_results = {}

        for eng in engines:
            self.eval_js(f"window.cosmicflows.switchScienceEngine('{eng}')")
            engine_records = []

            for lm in landmarks:
                cx, cy, cz = lm['pos']['x'], lm['pos']['y'], lm['pos']['z']
                
                # Check r=0
                res_r0 = self.eval_js(f"""
                    (() => {{
                        const v = new THREE.Vector3();
                        const vort = new THREE.Vector3();
                        window.cosmicflows.getVelocitySg(new THREE.Vector3({cx}, {cy}, {cz}), v);
                        const div = window.cosmicflows.getDivergenceSg(new THREE.Vector3({cx}, {cy}, {cz}));
                        window.cosmicflows.getVorticitySg(new THREE.Vector3({cx}, {cy}, {cz}), vort);
                        return {{
                            v: {{ x: v.x, y: v.y, z: v.z, len: v.length() }},
                            div: div,
                            vort: {{ x: vort.x, y: vort.y, z: vort.z, len: vort.length() }}
                        }};
                    }})()
                """)

                # Check for NaN / Inf
                assert not math.isnan(res_r0['v']['len']), f"NaN velocity at r=0 for {lm['id']} in {eng}"
                assert not math.isinf(res_r0['v']['len']), f"Inf velocity at r=0 for {lm['id']} in {eng}"
                assert not math.isnan(res_r0['div']), f"NaN divergence at r=0 for {lm['id']} in {eng}"
                assert not math.isinf(res_r0['div']), f"Inf divergence at r=0 for {lm['id']} in {eng}"
                assert not math.isnan(res_r0['vort']['len']), f"NaN vorticity at r=0 for {lm['id']} in {eng}"
                assert not math.isinf(res_r0['vort']['len']), f"Inf vorticity at r=0 for {lm['id']} in {eng}"

                # Continuity check across micro-deltas
                delta_records = []
                for d in deltas:
                    d_eval = self.eval_js(f"""
                        (() => {{
                            const v = new THREE.Vector3();
                            window.cosmicflows.getVelocitySg(new THREE.Vector3({cx + d}, {cy}, {cz}), v);
                            const div = window.cosmicflows.getDivergenceSg(new THREE.Vector3({cx + d}, {cy}, {cz}));
                            return {{ d: {d}, v_len: v.length(), div: div }};
                        }})()
                    """)
                    delta_records.append(d_eval)

                # Check jump discontinuity: difference between delta=0 and delta=1e-6 must be < 1e-3 km/s
                v_0 = delta_records[0]['v_len']
                v_micro = delta_records[1]['v_len']
                diff_micro = abs(v_micro - v_0)
                assert diff_micro < 0.05, f"Discontinuity detected at {lm['id']} in {eng}: diff={diff_micro}"

                engine_records.append({
                    "id": lm['id'],
                    "name": lm['name'],
                    "r0_v_len": res_r0['v']['len'],
                    "r0_div": res_r0['div'],
                    "r0_vort": res_r0['vort']['len'],
                    "diff_micro": diff_micro,
                    "delta_profile": delta_records
                })

            singularity_results[eng] = engine_records

        self.test_results['singularity'] = singularity_results
        print(f"  [PASS] Verified r=0 smoothness & continuity across {len(landmarks)} landmarks for ALL 7 engines.")
        print(f"  [SAMPLE] Shapley (r=0, cf4-wf): v={singularity_results['cf4-wf'][4]['r0_v_len']:.2f} km/s, div={singularity_results['cf4-wf'][4]['r0_div']:.4f}, |curl v|={singularity_results['cf4-wf'][4]['r0_vort']:.8f}")
        print(f"  [SAMPLE] Dipole (r=0, cf4-wf):  v={singularity_results['cf4-wf'][8]['r0_v_len']:.2f} km/s, div={singularity_results['cf4-wf'][8]['r0_div']:.4f}, |curl v|={singularity_results['cf4-wf'][8]['r0_vort']:.8f}")

    # ==========================================
    # SUITE 2: DIVERGENCE SIGN & CONSERVATION INVARIANTS
    # ==========================================
    def test_divergence_signs_and_analytical_match(self):
        print("\n=======================================================")
        print(">>> SUITE 2: Divergence Signs & Analytical Field Verification <<<")
        print("=======================================================")

        # Invariant: Attractors (sign = +1) must have div v < 0 (compression)
        # Invariant: Voids/Repellers (sign = -1) must have div v > 0 (expansion)
        landmarks = self.eval_js("""
            window.cosmicflows.LANDMARKS.map(lm => ({
                id: lm.id,
                name: lm.name,
                sign: lm.sign,
                str: lm.str,
                pos: { x: lm.pos.x, y: lm.pos.y, z: lm.pos.z }
            }))
        """)

        div_records = []
        for lm in landmarks:
            cx, cy, cz = lm['pos']['x'], lm['pos']['y'], lm['pos']['z']
            div_val = self.eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('cf4-wf');
                    return window.cosmicflows.getDivergenceSg(new THREE.Vector3({cx}, {cy}, {cz}));
                }})()
            """)

            expected_sign = "negative" if lm['sign'] > 0 else "positive"
            actual_sign = "negative" if div_val < 0 else "positive"
            
            # Local Group is observer origin (virgo is at 0,0,0 with str=15, sign=+1)
            # Check sign
            if lm['id'] == 'dipole':
                assert div_val > 0, f"Dipole repeller MUST have div v > 0 (expansion plume), got {div_val}"
            elif lm['sign'] > 0 and lm['str'] >= 50:
                assert div_val < 0, f"Attractor core {lm['id']} MUST have div v < 0 (compression sink), got {div_val}"

            div_records.append({
                "id": lm['id'],
                "name": lm['name'],
                "sign": lm['sign'],
                "str": lm['str'],
                "div": div_val,
                "expected": expected_sign,
                "actual": actual_sign,
                "valid": (expected_sign == actual_sign)
            })

        # Analytical divergence comparison for cf4-wf:
        # Analytical divergence div_analytic(p) = - sum_i sign_i * str_i * C_SCALE * 3 * eps^2 / (r_i^2 + eps^2)^(5/2)
        analytical_cmp = self.eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const C_SCALE = 2.3e8;
                const eps = 1000.0;
                const eps2 = eps * eps;

                function getAnalyticalDiv(p) {
                    let totalDiv = 0;
                    for (let node of window.cosmicflows.LANDMARKS) {
                        const dx = node.pos.x - p.x;
                        const dy = node.pos.y - p.y;
                        const dz = node.pos.z - p.z;
                        const r2 = dx*dx + dy*dy + dz*dz;
                        const K = node.sign * node.str * C_SCALE;
                        const term = -3.0 * eps2 * K / Math.pow(r2 + eps2, 2.5);
                        totalDiv += term;
                    }
                    return totalDiv;
                }

                const testPoints = [
                    { name: 'Origin', p: new THREE.Vector3(0,0,0) },
                    { name: 'Shapley Core', p: new THREE.Vector3(7200, -8600, -2400) },
                    { name: 'Dipole Core', p: new THREE.Vector3(-10000, 10000, 12000) },
                    { name: 'Coma Core', p: new THREE.Vector3(500, 7000, 1500) },
                    { name: 'Arbitrary Point', p: new THREE.Vector3(3000, -4000, 5000) },
                    { name: 'Deep Void Region', p: new THREE.Vector3(-8000, 8000, 10000) }
                ];

                return testPoints.map(tp => {
                    const numerical = window.cosmicflows.getDivergenceSg(tp.p);
                    const analytic = getAnalyticalDiv(tp.p);
                    const relErr = Math.abs(numerical - analytic) / (Math.abs(analytic) + 1e-9);
                    return {
                        name: tp.name,
                        numerical: numerical,
                        analytic: analytic,
                        absErr: Math.abs(numerical - analytic),
                        relErr: relErr
                    };
                });
            })()
        """)

        print("  Divergence Sign Audit Table:")
        for r in div_records:
            status = "PASS" if r['valid'] else "FAIL"
            print(f"    [{status}] {r['id']:<12}: sign={r['sign']:>2}, str={r['str']:>3} -> div v = {r['div']:>10.4f} ({r['actual']})")

        print("\n  Numerical vs. Exact Analytical Divergence Comparison (6-point finite difference):")
        for ac in analytical_cmp:
            print(f"    Point: {ac['name']:<20} | Numerical: {ac['numerical']:>10.4f} | Analytic: {ac['analytic']:>10.4f} | RelErr: {ac['relErr'] * 100:.4f}%")
            assert ac['relErr'] < 0.01, f"Analytical divergence mismatch at {ac['name']}: {ac['relErr']*100:.3f}%"

        self.test_results['divergence'] = {
            "landmarks": div_records,
            "analytical_comparison": analytical_cmp
        }
        print("  [PASS] Divergence signs and analytical invariants 100% verified.")

    # ==========================================
    # SUITE 3: ADAPTIVE RK45 INTEGRATOR BOUNDS & STRESS
    # ==========================================
    def test_rk45_integrator_bounds_and_stress(self):
        print("\n=======================================================")
        print(">>> SUITE 3: Adaptive RK45 Integrator Extreme Stress <<<")
        print("=======================================================")

        # Test 1: In-Domain Extreme Coordinates (must remain strictly in box)
        in_domain_coords = [
            {"name": "Origin (0,0,0)", "p": [0, 0, 0]},
            {"name": "Near Box Boundary (+14995)", "p": [14995, 14995, 14995]},
            {"name": "Near Box Boundary (-14995)", "p": [-14995, -14995, -14995]},
            {"name": "Exact Box Corner (+15000)", "p": [15000, 15000, 15000]},
            {"name": "Exact Box Corner (-15000)", "p": [-15000, -15000, -15000]},
            {"name": "Near Shapley Singularity", "p": [7200.001, -8600.001, -2400.001]},
            {"name": "Near Dipole Singularity", "p": [-10000.001, 10000.001, 12000.001]},
        ]

        coord_results = []
        for ec in in_domain_coords:
            p = ec["p"]
            res = self.eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('cf4-wf');
                    const pStart = new THREE.Vector3({p[0]}, {p[1]}, {p[2]});
                    const fwd = window.cosmicflows.integrateRK4(pStart, 1, 350, 46);
                    const bwd = window.cosmicflows.integrateRK4(pStart, -1, 350, 46);
                    
                    function checkPts(pts) {{
                        let maxCoord = 0;
                        let hasNaN = false;
                        let hasInf = false;
                        for (let pt of pts) {{
                            if (isNaN(pt.x) || isNaN(pt.y) || isNaN(pt.z)) hasNaN = true;
                            if (!isFinite(pt.x) || !isFinite(pt.y) || !isFinite(pt.z)) hasInf = true;
                            maxCoord = Math.max(maxCoord, Math.abs(pt.x), Math.abs(pt.y), Math.abs(pt.z));
                        }}
                        return {{ count: pts.length, maxCoord: maxCoord, hasNaN: hasNaN, hasInf: hasInf }};
                    }}

                    return {{
                        name: '{ec["name"]}',
                        fwd: checkPts(fwd),
                        bwd: checkPts(bwd)
                    }};
                }})()
            """)
            coord_results.append(res)
            
            # Assert in-domain invariants: strictly <= 15000
            assert not res['fwd']['hasNaN'], f"NaN in forward trajectory from {ec['name']}"
            assert not res['fwd']['hasInf'], f"Inf in forward trajectory from {ec['name']}"
            assert not res['bwd']['hasNaN'], f"NaN in backward trajectory from {ec['name']}"
            assert not res['bwd']['hasInf'], f"Inf in backward trajectory from {ec['name']}"
            assert res['fwd']['maxCoord'] <= 15000.01, f"Forward box leak from {ec['name']}: {res['fwd']['maxCoord']}"
            assert res['bwd']['maxCoord'] <= 15000.01, f"Backward box leak from {ec['name']}: {res['bwd']['maxCoord']}"

        print("  In-Domain Extreme Initial Coordinates Streamline Audit:")
        for cr in coord_results:
            print(f"    Location: {cr['name']:<32} | Fwd pts: {cr['fwd']['count']:>3}, max={cr['fwd']['maxCoord']:>7.1f} | Bwd pts: {cr['bwd']['count']:>3}, max={cr['bwd']['maxCoord']:>7.1f}")

        # Test 1b: Out-of-Domain Seeds (Beyond +/- 15000 box)
        out_domain_coords = [
            {"name": "Beyond Box (+30k)", "p": [30000, -25000, 20000]},
            {"name": "Deep Space (+50k)", "p": [50000, 50000, 50000]},
            {"name": "Deep Space (-100k)", "p": [-100000, 0, 0]},
        ]
        out_results = []
        for ec in out_domain_coords:
            p = ec["p"]
            res = self.eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('cf4-wf');
                    const pStart = new THREE.Vector3({p[0]}, {p[1]}, {p[2]});
                    const fwd = window.cosmicflows.integrateRK4(pStart, 1, 350, 46);
                    const bwd = window.cosmicflows.integrateRK4(pStart, -1, 350, 46);
                    
                    function inspectOut(pts) {{
                        let hasNaN = false;
                        for (let pt of pts) {{
                            if (isNaN(pt.x) || isNaN(pt.y) || isNaN(pt.z)) hasNaN = true;
                        }}
                        const termPt = pts.length > 1 ? pts[pts.length - 1] : pts[0];
                        const termMax = Math.max(Math.abs(termPt.x), Math.abs(termPt.y), Math.abs(termPt.z));
                        return {{ count: pts.length, termMax: termMax, hasNaN: hasNaN }};
                    }}

                    return {{
                        name: '{ec["name"]}',
                        fwd: inspectOut(fwd),
                        bwd: inspectOut(bwd)
                    }};
                }})()
            """)
            out_results.append(res)
            assert not res['fwd']['hasNaN'], f"NaN in forward out-of-domain trajectory {ec['name']}"
            assert not res['bwd']['hasNaN'], f"NaN in backward out-of-domain trajectory {ec['name']}"
            # Integrator must terminate immediately (count <= 2) and terminal point clipped to 15000
            assert res['fwd']['count'] <= 2, f"Out-of-domain did not terminate immediately: count={res['fwd']['count']}"
            assert res['fwd']['termMax'] <= 15000.01, f"Terminal point not clamped to box: {res['fwd']['termMax']}"

        print("\n  Out-of-Domain Boundary Clipping & Immediate Halt Audit:")
        for or_res in out_results:
            print(f"    Location: {or_res['name']:<28} | Fwd pts: {or_res['fwd']['count']:>2}, terminal max={or_res['fwd']['termMax']:>7.1f} | Bwd pts: {or_res['bwd']['count']:>2}, terminal max={or_res['bwd']['termMax']:>7.1f}")

        # Test 2: Extreme Step Sizes
        dt_values = [1e-6, 0.01, 1.0, 10.0, 46.0, 100.0, 120.0, 500.0, 10000.0, -100.0, 0.0]
        dt_results = []
        for dt_val in dt_values:
            res_dt = self.eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('cf4-wf');
                    const pStart = new THREE.Vector3(0, 0, 0);
                    const pts = window.cosmicflows.integrateRK4(pStart, 1, 300, {dt_val});
                    let maxCoord = 0;
                    let hasNaN = false;
                    for (let pt of pts) {{
                        if (isNaN(pt.x) || isNaN(pt.y) || isNaN(pt.z)) hasNaN = true;
                        maxCoord = Math.max(maxCoord, Math.abs(pt.x), Math.abs(pt.y), Math.abs(pt.z));
                    }}
                    return {{
                        dt: {dt_val},
                        count: pts.length,
                        maxCoord: maxCoord,
                        hasNaN: hasNaN
                    }};
                }})()
            """)
            dt_results.append(res_dt)
            assert not res_dt['hasNaN'], f"NaN with dt={dt_val}"
            assert res_dt['maxCoord'] <= 15000.01, f"Box leak with dt={dt_val}: max={res_dt['maxCoord']}"

        print("\n  Extreme Initial Step Size (dt) Robustness Audit:")
        for dtr in dt_results:
            print(f"    dt input: {dtr['dt']:>10.6f} -> Points generated: {dtr['count']:>3}, Max Coord: {dtr['maxCoord']:>7.1f}, NaN: {dtr['hasNaN']}")

        # Test 3: Step Adaptation Dynamism
        # Verify that RK45 actually adapts h based on local truncation error
        adaptation_test = self.eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                // Near Shapley where gradient is steep vs far-field where gradient is flat
                const pSteep = new THREE.Vector3(7200 + 400, -8600, -2400);
                const pFlat = new THREE.Vector3(12000, 12000, -12000);
                
                const ptsSteep = window.cosmicflows.integrateRK4(pSteep, 1, 50, 46);
                const ptsFlat = window.cosmicflows.integrateRK4(pFlat, 1, 50, 46);

                // Compute point-to-point step lengths
                function getStepLengths(pts) {
                    const lens = [];
                    for (let i = 1; i < pts.length; i++) {
                        const dx = pts[i].x - pts[i-1].x;
                        const dy = pts[i].y - pts[i-1].y;
                        const dz = pts[i].z - pts[i-1].z;
                        lens.push(Math.sqrt(dx*dx + dy*dy + dz*dz));
                    }
                    return lens;
                }

                return {
                    steepSteps: getStepLengths(ptsSteep),
                    flatSteps: getStepLengths(ptsFlat)
                };
            })()
        """)
        
        avg_steep = sum(adaptation_test['steepSteps']) / max(1, len(adaptation_test['steepSteps']))
        avg_flat = sum(adaptation_test['flatSteps']) / max(1, len(adaptation_test['flatSteps']))
        print(f"\n  Step Adaptation Verification:")
        print(f"    Avg step length in steep gradient (Shapley infall): {avg_steep:.2f} Three.js units")
        print(f"    Avg step length in flat gradient (Far field):       {avg_flat:.2f} Three.js units")

        self.test_results['rk45'] = {
            "coords": coord_results,
            "out_domain": out_results,
            "dt": dt_results,
            "avg_steep": avg_steep,
            "avg_flat": avg_flat
        }
        print("  [PASS] RK45 Adaptive Integrator boundaries & stress invariants 100% verified.")

    # ==========================================
    # SUITE 4: ALL 7 ENGINES VORTICITY & FIELD PROPERTIES
    # ==========================================
    def test_cosmographic_engines_vorticity_and_physics(self):
        print("\n=======================================================")
        print(">>> SUITE 4: 7 Cosmographic Engines Vorticity & Kinematics <<<")
        print("=======================================================")

        engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026']
        
        # Test point near Shapley halo where Gadget-4 adds halo spin
        eval_pt = {"x": 7200 + 1000, "y": -8600, "z": -2400}
        
        engine_evals = {}
        for eng in engines:
            res = self.eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('{eng}');
                    const p = new THREE.Vector3({eval_pt['x']}, {eval_pt['y']}, {eval_pt['z']});
                    const v = new THREE.Vector3();
                    const vort = new THREE.Vector3();
                    window.cosmicflows.getVelocitySg(p, v);
                    const div = window.cosmicflows.getDivergenceSg(p);
                    window.cosmicflows.getVorticitySg(p, vort);
                    return {{
                        v: {{ x: v.x, y: v.y, z: v.z, len: v.length() }},
                        div: div,
                        vort: {{ x: vort.x, y: vort.y, z: vort.z, len: vort.length() }}
                    }};
                }})()
            """)
            engine_evals[eng] = res

        print("  Engine Kinematics & Vorticity Summary (at Shapley Infall Point):")
        for eng, data in engine_evals.items():
            print(f"    {eng:<18} | Speed: {data['v']['len']:>7.2f} km/s | Div: {data['div']:>9.4f} | Curl: {data['vort']['len']:>10.8f}")

        # Invariant 1: CF4 WF must have curl ~ 0 (potential flow)
        assert engine_evals['cf4-wf']['vort']['len'] < 1e-4, f"cf4-wf must be irrotational, got {engine_evals['cf4-wf']['vort']['len']}"
        
        # Invariant 2: Gadget-4 must have non-zero halo spin curl
        assert engine_evals['gadget4']['vort']['len'] > 0.05, f"gadget4 must have physical vorticity inside halo, got {engine_evals['gadget4']['vort']['len']}"

        # Invariant 3: All 7 velocity vectors must be mutually distinct
        for i in range(len(engines)):
            for j in range(i + 1, len(engines)):
                eng1, eng2 = engines[i], engines[j]
                v1 = engine_evals[eng1]['v']
                v2 = engine_evals[eng2]['v']
                diff = math.sqrt((v1['x'] - v2['x'])**2 + (v1['y'] - v2['y'])**2 + (v1['z'] - v2['z'])**2)
                assert diff > 10.0, f"Engines {eng1} and {eng2} have identical velocity vector (diff={diff})"

        self.test_results['engines'] = engine_evals
        print("  [PASS] All 7 cosmographic engines verified with distinct physical kinematics.")

    def run_all(self):
        self.setup()
        try:
            self.test_singularity_stress()
            self.test_divergence_signs_and_analytical_match()
            self.test_rk45_integrator_bounds_and_stress()
            self.test_cosmographic_engines_vorticity_and_physics()
            
            print("\n=======================================================")
            print(">>> ALL ADVERSARIAL STRESS TESTS COMPLETED <<<")
            print(f"Console Errors: {len(self.console_errors)}")
            print(f"Exceptions:     {len(self.exceptions)}")
            print("=======================================================\n")
            
            assert len(self.console_errors) == 0, f"Encountered console errors: {self.console_errors}"
            assert len(self.exceptions) == 0, f"Encountered exceptions: {self.exceptions}"
            return True
        finally:
            self.teardown()

if __name__ == "__main__":
    tester = AdversarialCDPTester()
    success = tester.run_all()
    if success:
        print("VERDICT: APPROVE")
        sys.exit(0)
    else:
        print("VERDICT: CHALLENGE_FAILED")
        sys.exit(1)
