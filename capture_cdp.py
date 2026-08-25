import subprocess
import time
import json
import base64
import requests
import websocket
import os
import sys

def take_screenshot(output_path="screenshot.png", wait_scene=3.5):
    user_data_dir = os.path.join(os.environ.get("TEMP", "C:\\Temp"), f"chrome_ss_{int(time.time())}")
    chrome_cmd = [
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        "--headless=new",
        "--remote-debugging-port=9222",
        "--remote-allow-origins=*",
        f"--user-data-dir={user_data_dir}",
        "--window-size=1920,1080",
        "--enable-webgl",
        "--use-gl=angle",
        "--no-first-run",
        "--no-default-browser-check",
        "http://localhost:8000"
    ]
    print("Spawning Chrome process...")
    proc = subprocess.Popen(chrome_cmd)
    
    ws_url = None
    for i in range(25):
        time.sleep(0.3)
        try:
            r = requests.get("http://127.0.0.1:9222/json/list", timeout=1)
            targets = r.json()
            for t in targets:
                if t.get("type") == "page" and "8000" in t.get("url", ""):
                    ws_url = t.get("webSocketDebuggerUrl")
                    break
            if ws_url:
                break
        except Exception:
            pass
            
    if not ws_url:
        print("Failed to find target page ws_url")
        proc.kill()
        sys.exit(1)
        
    print("Connecting to WebSocket:", ws_url)
    ws = websocket.create_connection(ws_url, timeout=10)
    
    # Enable Page
    ws.send(json.dumps({"id": 1, "method": "Page.enable"}))
    time.sleep(wait_scene)
    
    # Capture Screenshot
    ws.send(json.dumps({
        "id": 2,
        "method": "Page.captureScreenshot",
        "params": {
            "format": "png",
            "captureBeyondViewport": False
        }
    }))
    
    while True:
        resp = json.loads(ws.recv())
        if resp.get("id") == 2:
            data = resp["result"]["data"]
            img_bytes = base64.b64decode(data)
            with open(output_path, "wb") as f:
                f.write(img_bytes)
            print(f"Successfully captured and saved {output_path} ({len(img_bytes)} bytes)")
            break
            
    ws.close()
    proc.kill()

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "screenshot_1.png"
    take_screenshot(out)
