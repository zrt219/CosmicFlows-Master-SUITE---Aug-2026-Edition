import subprocess, time, json, urllib.request, websocket, sys, base64

chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
user_data = f"C:\\Users\\Zhane\\AppData\\Local\\Temp\\chrome_ss_{int(time.time())}"
cmd = [
    chrome_path,
    "--headless=new",
    "--remote-debugging-port=9275",
    "--remote-debugging-address=0.0.0.0",
    "--remote-allow-origins=*",
    f"--user-data-dir={user_data}",
    "--window-size=1920,1080",
    "--enable-webgl",
    "--use-gl=angle",
    "--no-first-run",
    "--no-default-browser-check",
    "--disable-background-networking",
    "http://localhost:8000"
]

proc = subprocess.Popen(cmd)
time.sleep(2.5)

try:
    req = urllib.request.urlopen("http://localhost:9275/json/list")
    pages = json.loads(req.read().decode())
    target = next((p for p in pages if "8000" in p.get("url", "")), None)
    if not target:
        print("Page not found:", pages)
        sys.exit(1)

    ws = websocket.create_connection(target["webSocketDebuggerUrl"])
    ws.send(json.dumps({"id": 1, "method": "Page.enable"}))
    ws.recv()
    ws.send(json.dumps({"id": 2, "method": "Runtime.enable"}))
    ws.recv()

    # Diagnostics and click
    diag_script = """
    (() => {
      // 1. Hide loader
      const loader = document.getElementById('loading');
      if (loader) loader.style.display = 'none';

      // 2. Select compute tab
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.pane').forEach(p => p.classList.remove('active'));
      
      const computeBtn = document.querySelector('[data-tab="tab-compute"]');
      if (computeBtn) computeBtn.classList.add('active');
      const computePane = document.getElementById('tab-compute');
      if (computePane) computePane.classList.add('active');

      return {
        computeBtnFound: !!computeBtn,
        computePaneFound: !!computePane,
        computePaneDisplay: computePane ? window.getComputedStyle(computePane).display : null,
        provenanceCardFound: !!document.getElementById('compute-provenance')
      };
    })()
    """
    ws.send(json.dumps({
        "id": 10,
        "method": "Runtime.evaluate",
        "params": {"expression": diag_script, "returnByValue": True}
    }))
    while True:
        r = json.loads(ws.recv())
        if r.get("id") == 10:
            print("Diagnostics result:", r)
            break

    time.sleep(1.0)

    ws.send(json.dumps({"id": 3, "method": "Page.captureScreenshot", "params": {"format": "png"}}))
    while True:
        resp = json.loads(ws.recv())
        if resp.get("id") == 3:
            img_data = base64.b64decode(resp["result"]["data"])
            with open("screenshot_compute_tab.png", "wb") as f:
                f.write(img_data)
            print(f"SUCCESS! Saved screenshot_compute_tab.png ({len(img_data)} bytes)")
            break
finally:
    proc.terminate()
