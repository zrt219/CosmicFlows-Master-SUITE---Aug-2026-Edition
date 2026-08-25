import os
import sys
import time
import json
import base64
import subprocess
import urllib.request
import urllib.error
import websocket
from typing import Any, Dict, List, Optional, Tuple

class ChromeCDPClient:
    """
    Chrome DevTools Protocol (CDP) WebSocket client for automated testing
    of the Cosmicflows / Wiener Filter Master Suite.
    """
    def __init__(self, port: int = 9222, url: str = "http://localhost:8000", spawn_headless: bool = True):
        self.port = port
        self.target_url = url
        self.spawn_headless = spawn_headless
        self.proc: Optional[subprocess.Popen] = None
        self.ws: Optional[websocket.WebSocket] = None
        self.msg_id = 0
        self.console_logs: List[Dict[str, Any]] = []
        self.console_errors: List[str] = []
        self.console_warnings: List[str] = []
        self.uncaught_exceptions: List[str] = []

    def clear_logs(self):
        self.console_logs.clear()
        self.console_errors.clear()
        self.console_warnings.clear()
        self.uncaught_exceptions.clear()

    def find_chrome_path(self) -> str:
        candidates = [
            r"C:\Program Files\Google\Chrome\Application\chrome.exe",
            r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
            os.path.expanduser(r"~\AppData\Local\Google\Chrome\Application\chrome.exe"),
            r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
            r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
            r"C:\Program Files\BraveSoftware\Brave-Browser\Application\brave.exe",
            "/usr/bin/google-chrome",
            "/usr/bin/chromium-browser",
            "/usr/bin/chromium",
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        ]
        for c in candidates:
            if os.path.exists(c):
                return c
        return "chrome"

    def is_port_open(self) -> bool:
        try:
            req = urllib.request.urlopen(f"http://localhost:{self.port}/json/version", timeout=1.0)
            return req.status == 200
        except Exception:
            return False

    def start(self):
        if not self.is_port_open():
            if self.spawn_headless:
                chrome_exe = self.find_chrome_path()
                user_data = os.path.join(
                    os.environ.get("TEMP", r"C:\Windows\Temp"),
                    f"chrome_cdp_test_{int(time.time() * 1000)}"
                )
                cmd = [
                    chrome_exe,
                    "--headless=new",
                    f"--remote-debugging-port={self.port}",
                    "--remote-debugging-address=0.0.0.0",
                    "--remote-allow-origins=*",
                    f"--user-data-dir={user_data}",
                    "--window-size=1920,1080",
                    "--enable-webgl",
                    "--use-gl=angle",
                    "--no-first-run",
                    "--no-default-browser-check",
                    "--disable-background-networking",
                    self.target_url
                ]
                self.proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                time.sleep(2.0)
            else:
                raise RuntimeError(f"CDP Port {self.port} is not accessible and spawn_headless=False")

        # Connect to page target
        for attempt in range(10):
            try:
                req = urllib.request.urlopen(f"http://localhost:{self.port}/json/list", timeout=2.0)
                pages = json.loads(req.read().decode())
                target = next((p for p in pages if "8000" in p.get("url", "") or p.get("type") == "page"), None)
                if target and "webSocketDebuggerUrl" in target:
                    self.ws = websocket.create_connection(target["webSocketDebuggerUrl"], timeout=30.0)
                    break
            except Exception as e:
                time.sleep(0.5)
        
        if not self.ws:
            raise RuntimeError(f"Failed to connect to Chrome WebSocket on port {self.port}")

        # Enable domains
        self.send_cdp("Page.enable")
        self.send_cdp("Runtime.enable")
        self.send_cdp("DOM.enable")
        time.sleep(1.0)

    def send_cdp(self, method: str, params: Optional[Dict[str, Any]] = None) -> Any:
        self.msg_id += 1
        req_id = self.msg_id
        payload = {"id": req_id, "method": method}
        if params is not None:
            payload["params"] = params
        self.ws.send(json.dumps(payload))

        while True:
            raw = self.ws.recv()
            if not raw:
                continue
            data = json.loads(raw)
            if data.get("method") == "Runtime.consoleAPICalled":
                params_data = data.get("params", {})
                msg_type = params_data.get("type", "log")
                args = params_data.get("args", [])
                text = " ".join(str(a.get("value", a.get("description", ""))) for a in args)
                log_entry = {"type": msg_type, "text": text, "timestamp": time.time()}
                self.console_logs.append(log_entry)
                if msg_type == "error":
                    self.console_errors.append(text)
                elif msg_type == "warning":
                    self.console_warnings.append(text)
            elif data.get("method") == "Runtime.exceptionThrown":
                exc_details = data.get("params", {}).get("exceptionDetails", {})
                exc_msg = exc_details.get("text", "")
                if "exception" in exc_details:
                    exc_msg += " " + str(exc_details["exception"].get("description", ""))
                self.uncaught_exceptions.append(exc_msg)

            if data.get("id") == req_id:
                if "error" in data:
                    raise RuntimeError(f"CDP Error in {method}: {data['error']}")
                return data.get("result", {})

    def evaluate(self, expression: str, return_by_value: bool = True, await_promise: bool = False) -> Any:
        res = self.send_cdp("Runtime.evaluate", {
            "expression": expression,
            "returnByValue": return_by_value,
            "awaitPromise": await_promise
        })
        if "exceptionDetails" in res:
            exc = res["exceptionDetails"]
            raise RuntimeError(f"JS Evaluation Error in '{expression}': {exc.get('text', '')} {exc.get('exception', {}).get('description', '')}")
        result_obj = res.get("result", {})
        if result_obj.get("subtype") == "error":
            raise RuntimeError(f"JS Evaluation Exception: {result_obj.get('description', 'Unknown error')}")
        if "value" in result_obj:
            return result_obj["value"]
        if result_obj.get("type") == "undefined":
            return None
        return result_obj

    def wait_for_condition(self, js_predicate: str, timeout: float = 10.0, poll_interval: float = 0.1) -> bool:
        start_time = time.time()
        while time.time() - start_time < timeout:
            try:
                res = self.evaluate(f"Boolean({js_predicate})")
                if res is True:
                    return True
            except Exception:
                pass
            time.sleep(poll_interval)
        return False

    def capture_screenshot(self, filepath: Optional[str] = None, format_type: str = "png") -> bytes:
        res = self.send_cdp("Page.captureScreenshot", {"format": format_type})
        data_base64 = res.get("data", "")
        img_bytes = base64.b64decode(data_base64)
        if filepath:
            with open(filepath, "wb") as f:
                f.write(img_bytes)
        return img_bytes

    def assert_no_errors(self):
        if self.uncaught_exceptions:
            raise AssertionError(f"Uncaught JavaScript exceptions detected: {self.uncaught_exceptions}")
        if self.console_errors:
            # Filter out non-fatal or favicon errors if any
            critical_errors = [e for e in self.console_errors if "favicon" not in e.lower()]
            if critical_errors:
                raise AssertionError(f"Console errors detected: {critical_errors}")

    def close(self):
        if self.ws:
            try:
                self.ws.close()
            except Exception:
                pass
            self.ws = None
        if self.proc:
            try:
                self.proc.terminate()
                self.proc.wait(timeout=2.0)
            except Exception:
                try:
                    self.proc.kill()
                except Exception:
                    pass
            self.proc = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
