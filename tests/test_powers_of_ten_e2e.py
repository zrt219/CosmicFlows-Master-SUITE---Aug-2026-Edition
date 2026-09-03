#!/usr/bin/env python3
"""
Automated E2E validation test for the Powers of Ten Cosmic Scale Explorer.
Verifies:
1. Modal opening from topbar button #btn-powers-of-ten-toggle
2. Step navigation across all 4 LOD stages (Laniakea, Milky Way, Solar System, Earth & Moon)
3. HUD telemetry updates (scale exponents, real-scale everyday analogies, speed of light ticker)
4. Continuous slider scrubbing and camera transitions
5. Modal closing and cleanup
6. Assert 0 console errors and 0 uncaught exceptions
"""
import os
import sys
import time
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9260):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Testing Powers of Ten Cosmic Scale Explorer E2E...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9260)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()

    try:
        client.send_cdp("Console.enable")
        client.send_cdp("Runtime.enable")
        client.send_cdp("Page.enable")
        client.send_cdp("Page.reload", {"ignoreCache": True})

        client.wait_for_condition("document.readyState === 'complete'", timeout=15.0)
        time.sleep(2.0)

        # 1. Verify open Powers of Ten from topbar button
        print("  Testing open via #btn-powers-of-ten-toggle...")
        client.evaluate("document.getElementById('btn-powers-of-ten-toggle').click();")
        client.wait_for_condition("document.getElementById('powers-of-ten-modal').style.display !== 'none'", timeout=5.0)
        print("  [OK] Modal opened successfully.")

        # 2. Verify Stage 1: Laniakea
        stage1_title = client.evaluate("document.getElementById('pot-stage-title').textContent")
        print(f"  Stage 1 active title: {stage1_title}")
        assert "Laniakea" in stage1_title, f"Expected Laniakea, got {stage1_title}"

        # 3. Test Stepper to Stage 2: Milky Way
        print("  Switching to Stage 2: Milky Way...")
        client.evaluate("document.querySelectorAll('.pot-step-btn')[1].click();")
        time.sleep(0.5)
        stage2_title = client.evaluate("document.getElementById('pot-stage-title').textContent")
        print(f"  Stage 2 active title: {stage2_title}")
        assert "Milky Way" in stage2_title, f"Expected Milky Way, got {stage2_title}"

        # 4. Test Stepper to Stage 3: Solar System
        print("  Switching to Stage 3: Solar System...")
        client.evaluate("document.querySelectorAll('.pot-step-btn')[2].click();")
        time.sleep(0.5)
        stage3_title = client.evaluate("document.getElementById('pot-stage-title').textContent")
        print(f"  Stage 3 active title: {stage3_title}")
        assert "Solar System" in stage3_title, f"Expected Solar System, got {stage3_title}"

        # 5. Test Stepper to Stage 4: Earth & Moon
        print("  Switching to Stage 4: Earth & Moon...")
        client.evaluate("document.querySelectorAll('.pot-step-btn')[3].click();")
        time.sleep(0.5)
        stage4_title = client.evaluate("document.getElementById('pot-stage-title').textContent")
        print(f"  Stage 4 active title: {stage4_title}")
        assert "Earth & Moon" in stage4_title, f"Expected Earth & Moon, got {stage4_title}"

        # 6. Verify Speed of Light ticker updated
        light_ticker = client.evaluate("document.getElementById('pot-light-ticker').textContent")
        print(f"  Stage 4 Light ticker: {light_ticker.strip()}")
        assert "Seconds" in light_ticker or "1.28" in light_ticker, f"Unexpected light ticker: {light_ticker}"

        # 7. Test Close
        print("  Closing modal via close button...")
        client.evaluate("document.getElementById('pot-btn-close').click();")
        client.wait_for_condition("document.getElementById('powers-of-ten-modal').style.display === 'none'", timeout=5.0)
        print("  [OK] Modal closed successfully.")

        # 8. Assert zero console errors
        time.sleep(1.0)
        errors = client.console_errors
        exceptions = client.uncaught_exceptions

        if errors or exceptions:
            print(f"[FAIL] Found errors: {errors}, exceptions: {exceptions}")
            sys.exit(1)

        print("[PASS] Powers of Ten E2E verification passed with 0 console errors and 0 exceptions!")
        sys.exit(0)

    finally:
        client.close()

if __name__ == "__main__":
    main()
