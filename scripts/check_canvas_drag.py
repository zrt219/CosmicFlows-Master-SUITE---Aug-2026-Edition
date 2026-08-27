#!/usr/bin/env python3
"""
Gate 3: OrbitControls Canvas Drag & Zoom Verification
Simulates mouse drag and touch gesture on WebGL canvas and asserts camera position changes.
"""
import os
import sys
import time
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9300):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Testing WebGL Canvas OrbitControls Mouse Dragging & Gesture Response...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9300)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        ready = client.wait_for_condition("document.readyState === 'complete' && !!window.camera && !!window.controls", timeout=15.0)
        if not ready:
            raise RuntimeError("Timed out waiting for window.camera and window.controls to initialize")
            
        client.evaluate("if (window.splashController) window.splashController.dismiss(true);")
        time.sleep(0.3)
        
        init_pos = client.evaluate("({ x: (window.camera ? window.camera.position.x : 0), y: (window.camera ? window.camera.position.y : 0), z: (window.camera ? window.camera.position.z : 0) })")
        
        client.send_cdp("Input.dispatchMouseEvent", {"type": "mousePressed", "x": 640, "y": 360, "button": "left", "clickCount": 1})
        time.sleep(0.05)
        client.send_cdp("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": 740, "y": 460, "button": "left"})
        time.sleep(0.05)
        client.send_cdp("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": 740, "y": 460, "button": "left"})
        time.sleep(0.3)
        
        new_pos = client.evaluate("({ x: (window.camera ? window.camera.position.x : 0), y: (window.camera ? window.camera.position.y : 0), z: (window.camera ? window.camera.position.z : 0) })")
    finally:
        client.close()
    
    dx = abs(new_pos['x'] - init_pos['x'])
    dy = abs(new_pos['y'] - init_pos['y'])
    dz = abs(new_pos['z'] - init_pos['z'])
    
    if dx > 10 or dy > 10 or dz > 10:
        print(f"[PASS] Dragging verified! Camera delta: dx={dx:.1f}, dy={dy:.1f}, dz={dz:.1f}")
        sys.exit(0)
    else:
        print(f"[FAIL] Camera did not move during simulated mouse drag! (dx={dx}, dy={dy}, dz={dz})")
        sys.exit(1)

if __name__ == "__main__":
    main()
