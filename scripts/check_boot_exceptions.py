#!/usr/bin/env python3
"""
Gate 2: Headless CDP Boot-Time Zero-Exception Gate
Starts index.html on port 8000, listens for Runtime exceptions and console errors for 3s.
"""
import os
import sys
import time
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9250):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Running Headless Chrome CDP Boot-Time Zero-Exception Gate...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9250)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        client.send_cdp("Console.enable")
        client.send_cdp("Runtime.enable")
        client.send_cdp("Page.enable")
        client.send_cdp("Page.reload", {"ignoreCache": True})
        
        client.wait_for_condition("document.readyState === 'complete'", timeout=15.0)
        time.sleep(2.0)
        exceptions = client.uncaught_exceptions
    finally:
        client.close()
    
    if exceptions:
        print(f"[FAIL] Found {len(exceptions)} uncaught boot exceptions:")
        for ex in exceptions:
            print("  -", ex)
        sys.exit(1)
    else:
        print("[PASS] Boot sequence completed with 0 uncaught exceptions.")
        sys.exit(0)

if __name__ == "__main__":
    main()
