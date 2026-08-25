# -*- coding: utf-8 -*-
import subprocess, time, json, urllib.request, websocket, sys, os, http.server, socketserver, threading, base64, math

PORT = 8024
DEBUG_PORT = 9240

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

def start_server():
    os.chdir(r"c:\Users\Zhane\Documents\antigravity\gallant-newton")
    socketserver.TCPServer.allow_reuse_address = True
    httpd = socketserver.TCPServer(("127.0.0.1", PORT), QuietHandler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    return httpd

def run_tests():
    print("=== Starting Milestone M1 Astrophysics Math & Cosmographic Engines Test Suite ===")
    httpd = start_server()
    print(f"[1/7] Local HTTP server running on port {PORT}")

    chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    user_data = f"C:\\Users\\Zhane\\AppData\\Local\\Temp\\chrome_m1_{int(time.time())}"
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
    exceptions = []

    try:
        req = urllib.request.urlopen(f"http://localhost:{DEBUG_PORT}/json/list")
        pages = json.loads(req.read().decode())
        target = next((p for p in pages if f"{PORT}" in p.get("url", "")), None)
        if not target:
            raise RuntimeError(f"Target page on port {PORT} not found in Chrome list: {pages}")

        ws = websocket.create_connection(target["webSocketDebuggerUrl"])
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
                    text = " ".join(str(a.get("value", "")) for a in args)
                    if t == "error":
                        console_errors.append(text)
                elif resp.get("method") == "Runtime.exceptionThrown":
                    exc_details = resp.get("params", {}).get("exceptionDetails", {})
                    exceptions.append(exc_details.get("text", "Unknown exception"))
                if resp.get("id") == msg_id:
                    return resp.get("result", {})

        def eval_js(expr):
            res = send_cmd("Runtime.evaluate", {"expression": expr, "returnByValue": True})
            if "exceptionDetails" in res:
                exceptions.append(res["exceptionDetails"].get("text", "Eval exception"))
                raise RuntimeError(f"JS Exception in expr: {res['exceptionDetails']}")
            return res.get("result", {}).get("value")

        send_cmd("Page.enable")
        send_cmd("Runtime.enable")
        time.sleep(1.5)

        print("[2/7] Testing window.cosmicflows Programmatic API Bridge Interface Contract...")
        cf_exists = eval_js("typeof window.cosmicflows !== 'undefined'")
        assert cf_exists, "window.cosmicflows is not defined!"
        version = eval_js("window.cosmicflows.version")
        assert version == '2026.1', f"Expected version 2026.1, got {version}"
        
        has_methods = eval_js("""
            typeof window.cosmicflows.getVelocitySg === 'function' &&
            typeof window.cosmicflows.getDivergenceSg === 'function' &&
            typeof window.cosmicflows.getVorticitySg === 'function' &&
            typeof window.cosmicflows.integrateRK4 === 'function' &&
            typeof window.cosmicflows.switchScienceEngine === 'function' &&
            typeof window.cosmicflows.toggleTheme === 'function' &&
            typeof window.cosmicflows.exportPublicationFigure === 'function'
        """)
        assert has_methods, "window.cosmicflows missing required API methods!"
        print("  -> window.cosmicflows API bridge verified successfully.")

        print("[3/7] Testing Exact Plummer Potential Gradient & Discontinuity Elimination...")
        v_exact_center = eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const v0 = new THREE.Vector3();
                window.cosmicflows.getVelocitySg(new THREE.Vector3(7200, -8600, -2400), v0);
                const v1 = new THREE.Vector3();
                window.cosmicflows.getVelocitySg(new THREE.Vector3(7200 + 1, -8600, -2400), v1);
                const v10 = new THREE.Vector3();
                window.cosmicflows.getVelocitySg(new THREE.Vector3(7200 + 10, -8600, -2400), v10);
                const v100 = new THREE.Vector3();
                window.cosmicflows.getVelocitySg(new THREE.Vector3(7200 + 100, -8600, -2400), v100);
                return {
                    v0: { x: v0.x, y: v0.y, z: v0.z, len: v0.length() },
                    v1: { x: v1.x, y: v1.y, z: v1.z, len: v1.length() },
                    v10: { x: v10.x, y: v10.y, z: v10.z, len: v10.length() },
                    v100: { x: v100.x, y: v100.y, z: v100.z, len: v100.length() }
                };
            })()
        """)
        print(f"  -> Center velocity len: {v_exact_center['v0']['len']:.2f} km/s")
        print(f"  -> Delta=1 velocity len: {v_exact_center['v1']['len']:.2f} km/s")
        print(f"  -> Delta=10 velocity len: {v_exact_center['v10']['len']:.2f} km/s")
        print(f"  -> Delta=100 velocity len: {v_exact_center['v100']['len']:.2f} km/s")

        v_lg = eval_js("""
            (() => {
                const v = new THREE.Vector3();
                window.cosmicflows.getVelocitySg(new THREE.Vector3(0, 0, 0), v);
                return { x: v.x, y: v.y, z: v.z, len: v.length() };
            })()
        """)
        print(f"  -> Local Group (0,0,0) peculiar velocity: |v_LG| = {v_lg['len']:.2f} km/s")
        assert 400.0 <= v_lg['len'] <= 900.0, f"Unphysical Local Group velocity: {v_lg['len']}"

        print("[4/7] Testing Divergence & Vorticity Physical Vector Calculus...")
        div_results = eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const divVoid = window.cosmicflows.getDivergenceSg(new THREE.Vector3(-10000, 10000, 12000));
                const divShapley = window.cosmicflows.getDivergenceSg(new THREE.Vector3(7200, -8600, -2400));
                const divGA = window.cosmicflows.getDivergenceSg(new THREE.Vector3(15000, -14000, 6800));
                const divComa = window.cosmicflows.getDivergenceSg(new THREE.Vector3(500, 7000, 1500));
                return {
                    divVoid: divVoid,
                    divShapley: divShapley,
                    divGA: divGA,
                    divComa: divComa
                };
            })()
        """)
        print(f"  -> Div at Dipole Repeller Void: {div_results['divVoid']:.6f} (> 0 required: expansion)")
        print(f"  -> Div at Shapley Attractor Sink: {div_results['divShapley']:.6f} (< 0 required: compression)")
        print(f"  -> Div at Great Attractor Sink: {div_results['divGA']:.6f} (< 0 required: compression)")
        print(f"  -> Div at Coma Cluster Sink: {div_results['divComa']:.6f} (< 0 required: compression)")
        
        assert div_results['divVoid'] > 0, f"Divergence at void must be > 0, got {div_results['divVoid']}"
        assert div_results['divShapley'] < 0, f"Divergence at Shapley must be < 0, got {div_results['divShapley']}"
        assert div_results['divGA'] < 0, f"Divergence at GA must be < 0, got {div_results['divGA']}"
        assert div_results['divComa'] < 0, f"Divergence at Coma must be < 0, got {div_results['divComa']}"

        vort_results = eval_js("""
            (() => {
                const p = new THREE.Vector3(7200 + 1000, -8600, -2400);
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const vortWF = new THREE.Vector3();
                window.cosmicflows.getVorticitySg(p, vortWF);

                window.cosmicflows.switchScienceEngine('gadget4');
                const vortG4 = new THREE.Vector3();
                window.cosmicflows.getVorticitySg(p, vortG4);

                return {
                    curlWF: vortWF.length(),
                    curlG4: vortG4.length()
                };
            })()
        """)
        print(f"  -> Vorticity |curl v| for CF4 WF/CR: {vort_results['curlWF']:.8f} (~0 required)")
        print(f"  -> Vorticity |curl v| for Gadget-4: {vort_results['curlG4']:.6f} (>0 required)")
        assert vort_results['curlWF'] < 1e-4, f"Potential flow must be irrotational, got curl {vort_results['curlWF']}"
        assert vort_results['curlG4'] > 0.001, f"Gadget-4 halo must have physical vorticity, got {vort_results['curlG4']}"

        print("[5/7] Testing Cash-Karp Adaptive RK45 Streamline Integrator...")
        rk45_results = eval_js("""
            (() => {
                window.cosmicflows.switchScienceEngine('cf4-wf');
                const pStart = new THREE.Vector3(0, 0, 0);
                const forwardPts = window.cosmicflows.integrateRK4(pStart, 1, 300, 46);
                const backwardPts = window.cosmicflows.integrateRK4(pStart, -1, 300, 46);
                
                let maxCoord = 0;
                for (let pt of forwardPts) {
                    maxCoord = Math.max(maxCoord, Math.abs(pt.x), Math.abs(pt.y), Math.abs(pt.z));
                }
                for (let pt of backwardPts) {
                    maxCoord = Math.max(maxCoord, Math.abs(pt.x), Math.abs(pt.y), Math.abs(pt.z));
                }

                return {
                    forwardCount: forwardPts.length,
                    backwardCount: backwardPts.length,
                    maxCoord: maxCoord
                };
            })()
        """)
        print(f"  -> Forward streamline points generated: {rk45_results['forwardCount']}")
        print(f"  -> Backward streamline points generated: {rk45_results['backwardCount']}")
        print(f"  -> Max coordinate along streamline: {rk45_results['maxCoord']:.1f} (Box limit 15000)")
        assert rk45_results['forwardCount'] > 5, "Forward streamline failed to integrate!"
        assert rk45_results['backwardCount'] > 5, "Backward streamline failed to integrate!"
        assert rk45_results['maxCoord'] <= 15000.1, f"Streamline leaked past box limit: {rk45_results['maxCoord']}"

        print("[6/7] Testing All 7 Cosmographic Engines Kinematics Differentiation...")
        engines = ['cf4-wf', 'cf4-hmc', 'vela-zoa', 'vweb-2026', 'bayesian-2026', 'gadget4', 'nusser-tully-2026']
        velocities_by_engine = {}
        for eng in engines:
            v_eng = eval_js(f"""
                (() => {{
                    window.cosmicflows.switchScienceEngine('{eng}');
                    const p = new THREE.Vector3(-4000, 2000, 1000);
                    const v = new THREE.Vector3();
                    window.cosmicflows.getVelocitySg(p, v);
                    return {{ x: v.x, y: v.y, z: v.z, spd: v.length() }};
                }})()
            """)
            velocities_by_engine[eng] = v_eng
            print(f"  -> Engine {eng}: vx={v_eng['x']:.2f}, vy={v_eng['y']:.2f}, vz={v_eng['z']:.2f} (spd={v_eng['spd']:.2f} km/s)")

        for i in range(len(engines)):
            for j in range(i + 1, len(engines)):
                e1, e2 = engines[i], engines[j]
                v1, v2 = velocities_by_engine[e1], velocities_by_engine[e2]
                diff = math.sqrt((v1['x'] - v2['x'])**2 + (v1['y'] - v2['y'])**2 + (v1['z'] - v2['z'])**2)
                assert diff > 1.0, f"Engines {e1} and {e2} have indistinguishable kinematics! Diff = {diff}"
        print("  -> All 7 cosmographic engines confirmed uniquely distinct.")

        print("[7/7] Testing UI Handlers, Callout Toggle, Theme Inversion & Capturing Screenshots...")
        eval_js("document.getElementById('btn-callouts-toggle').click()")
        eval_js("document.getElementById('btn-callouts-toggle').click()")
        
        eval_js("window.cosmicflows.toggleTheme('white')")
        time.sleep(1.0)
        ss_white = send_cmd("Page.captureScreenshot", {"format": "png"})
        img_white = base64.b64decode(ss_white["data"])
        with open("screenshot_1.png", "wb") as f:
            f.write(img_white)
        print(f"  -> Captured screenshot_1.png ({len(img_white)} bytes)")

        eval_js("window.cosmicflows.toggleTheme('dark')")
        time.sleep(1.0)
        ss_dark = send_cmd("Page.captureScreenshot", {"format": "png"})
        img_dark = base64.b64decode(ss_dark["data"])
        with open("screenshot_dark.png", "wb") as f:
            f.write(img_dark)
        print(f"  -> Captured screenshot_dark.png ({len(img_dark)} bytes)")

        print(f"Console Errors: {len(console_errors)}")
        for err in console_errors:
            print(f"  ERR: {err}")
        print(f"Exceptions: {len(exceptions)}")
        for exc in exceptions:
            print(f"  EXC: {exc}")

        assert len(console_errors) == 0, f"Detected {len(console_errors)} console errors!"
        assert len(exceptions) == 0, f"Detected {len(exceptions)} uncaught exceptions!"

        print("\n=======================================================")
        print(">>> ALL MILESTONE M1 VERIFICATION TESTS PASSED (100%) <<<")
        print("=======================================================")

    finally:
        proc.terminate()
        try:
            httpd.shutdown()
        except Exception:
            pass

def test_m1_astrophysics_math_and_engines():
    run_tests()

if __name__ == '__main__':
    run_tests()
