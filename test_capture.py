import subprocess, time, json, urllib.request, websocket, sys

chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
user_data = f"C:\\Users\\Zhane\\AppData\\Local\\Temp\\chrome_ss_{int(time.time())}"
cmd = [
    chrome_path,
    "--headless=new",
    "--remote-debugging-port=9222",
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
time.sleep(2.0)

try:
    req = urllib.request.urlopen("http://localhost:9222/json/list")
    pages = json.loads(req.read().decode())
    target = next((p for p in pages if "8000" in p.get("url", "")), None)
    if not target:
        print("Page not found in list:", pages)
        sys.exit(1)

    ws = websocket.create_connection(target["webSocketDebuggerUrl"])
    ws.send(json.dumps({"id": 1, "method": "Page.enable"}))
    ws.recv()
    ws.send(json.dumps({"id": 2, "method": "Runtime.enable"}))
    ws.recv()

    # Wait for rendering to complete
    time.sleep(1.5)

    ws.send(json.dumps({"id": 3, "method": "Page.captureScreenshot", "params": {"format": "png"}}))
    while True:
        resp = json.loads(ws.recv())
        if resp.get("method") == "Runtime.consoleAPICalled" or resp.get("method") == "Runtime.exceptionThrown":
            print("BROWSER CONSOLE/EXCEPTION:", resp)
        if resp.get("id") == 3:
            import base64
            img_data = base64.b64decode(resp["result"]["data"])
            with open("screenshot_1.png", "wb") as f:
                f.write(img_data)
            print(f"SUCCESS! Saved screenshot_1.png ({len(img_data)} bytes)")
            break
finally:
    proc.terminate()
