#!/usr/bin/env python3
"""
Gate 4: Dual-Theme PIL Non-Blank Visual Variance Check
Renders canvas in White Mode and Dark Mode, asserting standard deviation sigma > 12.0 and size > 35 KB.
"""
import os
import sys
import io
import time
import base64
import socket
from PIL import Image
import numpy as np

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9350):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Testing Dual-Theme Visual PIL Non-Blank Variance...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9350)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        ready = client.wait_for_condition(
            "document.readyState === 'complete' && !!window.camera && !!window.controls && typeof (window.toggleTheme || (window.cosmicflows && window.cosmicflows.toggleTheme)) === 'function'",
            timeout=15.0
        )
        if not ready:
            raise RuntimeError("Timed out waiting for camera, controls, and toggleTheme initialization")
            
        client.evaluate("if (window.splashController) window.splashController.dismiss(true);")
        time.sleep(0.5)
        
        # 1. White Mode Check
        client.evaluate("(window.toggleTheme || (window.cosmicflows && window.cosmicflows.toggleTheme))('white');")
        time.sleep(0.5)
        res_w = client.send_cdp("Page.captureScreenshot", {"format": "png"})
        img_w = Image.open(io.BytesIO(base64.b64decode(res_w["data"]))).convert("RGB")
        arr_w = np.array(img_w)
        std_w = float(np.std(arr_w))
        sz_w = len(base64.b64decode(res_w["data"])) / 1024.0
        
        # 2. Dark Mode Check
        client.evaluate("(window.toggleTheme || (window.cosmicflows && window.cosmicflows.toggleTheme))('dark');")
        time.sleep(0.5)
        res_d = client.send_cdp("Page.captureScreenshot", {"format": "png"})
        img_d = Image.open(io.BytesIO(base64.b64decode(res_d["data"]))).convert("RGB")
        arr_d = np.array(img_d)
        std_d = float(np.std(arr_d))
        sz_d = len(base64.b64decode(res_d["data"])) / 1024.0
    finally:
        client.close()
    
    print(f"  White Mode: size={sz_w:.1f} KB, std={std_w:.1f}")
    print(f"  Dark Mode:  size={sz_d:.1f} KB, std={std_d:.1f}")
    
    if std_w > 12.0 and sz_w > 35.0 and std_d > 12.0 and sz_d > 35.0:
        print("[PASS] Both White and Dark modes render high-variance non-blank visual content.")
        sys.exit(0)
    else:
        print("[FAIL] Visual variance check failed: blank or flat frame detected!")
        sys.exit(1)

if __name__ == "__main__":
    main()
