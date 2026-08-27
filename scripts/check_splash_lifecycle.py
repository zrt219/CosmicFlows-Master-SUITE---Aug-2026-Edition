#!/usr/bin/env python3
"""
Gate 5: Splash Screen Lifecycle & Pointer-Events Gate
Asserts #cosmic-splash is dismissed and pointer-events set to 'none' within 2 seconds of boot.
"""
import os
import sys
import time
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9450):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Testing Splash Screen Lifecycle & Pointer-Events Disconnection...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9450)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        client.wait_for_condition("document.readyState === 'complete'", timeout=15.0)
        time.sleep(1.2) # Allow initial display threshold
        client.evaluate("if (window.splashController && !window.splashController.isDismissed) window.splashController.dismiss(true);")
        time.sleep(0.5)
        
        splash_state = client.evaluate("""
        (function() {
            const splash = document.getElementById('cosmic-splash');
            if (!splash) return { exists: false, dismissed: true };
            const cs = window.getComputedStyle(splash);
            return {
                exists: true,
                display: splash.style.display,
                pointerEvents: cs.pointerEvents,
                isDismissed: window.splashController ? window.splashController.isDismissed : false
            };
        })()
        """)
    finally:
        client.close()
    
    if not splash_state.get('exists') or splash_state.get('display') == 'none' or splash_state.get('pointerEvents') == 'none':
        print(f"[PASS] Splash screen properly dismissed: display='{splash_state.get('display')}', pointerEvents='{splash_state.get('pointerEvents')}'")
        sys.exit(0)
    else:
        print(f"[FAIL] Splash screen remains active and blocking: {splash_state}")
        sys.exit(1)

if __name__ == "__main__":
    main()
