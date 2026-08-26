#!/usr/bin/env python3
"""
ZRT CosmicFlows-4 Master Pre-Commit Quality Gate
Executes 6 automated verification gates before any commit or push:
1. JavaScript AST / Static TDZ & Variable Scope Check
2. Headless Chrome CDP Boot-Time Zero-Exception Gate
3. OrbitControls Canvas Drag & Zoom Verification
4. Dual-Theme PIL Non-Blank Visual Variance Check
5. Splash Screen Lifecycle & Pointer-Events Gate
6. KaTeX In-Browser Mathematical Formula Linter

Usage:
  python scripts/pre_commit_gate.py
"""

import sys
import os
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPTS_DIR)

GATES = [
    {
        "name": "Gate 1: JavaScript AST / TDZ Static Analysis",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_js_tdz.py")],
        "timeout": 15
    },
    {
        "name": "Gate 2: Headless CDP Boot-Time Zero-Exception Gate",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_boot_exceptions.py")],
        "timeout": 30
    },
    {
        "name": "Gate 3: OrbitControls Canvas Drag & Zoom Verification",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_canvas_drag.py")],
        "timeout": 30
    },
    {
        "name": "Gate 4: Dual-Theme PIL Non-Blank Visual Variance Check",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_image_variance.py")],
        "timeout": 30
    },
    {
        "name": "Gate 5: Splash Screen Lifecycle & Pointer-Events Gate",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_splash_lifecycle.py")],
        "timeout": 30
    },
    {
        "name": "Gate 6: KaTeX In-Browser Mathematical Formula Linter",
        "cmd": [sys.executable, os.path.join(SCRIPTS_DIR, "check_katex_invariants.py")],
        "timeout": 30
    }
]

def run_gate(gate):
    start = time.time()
    try:
        proc = subprocess.run(
            gate["cmd"],
            cwd=ROOT_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=gate["timeout"]
        )
        elapsed = time.time() - start
        passed = (proc.returncode == 0)
        return {
            "name": gate["name"],
            "passed": passed,
            "elapsed": elapsed,
            "output": proc.stdout.strip()
        }
    except Exception as e:
        elapsed = time.time() - start
        return {
            "name": gate["name"],
            "passed": False,
            "elapsed": elapsed,
            "output": f"ERROR: {str(e)}"
        }

def main():
    print("=" * 80)
    print("   ZRT COSMICFLOWS-4 AUTOMATED PRE-COMMIT QUALITY GATE")
    print("=" * 80)
    print(f"Running {len(GATES)} automated gates...\n")

    start_all = time.time()
    results = []
    
    # Run with 2 workers to avoid Windows browser port/process resource contention
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(run_gate, g) for g in GATES]
        for f in futures:
            results.append(f.result())

    total_elapsed = time.time() - start_all
    all_passed = all(r["passed"] for r in results)

    print("-" * 80)
    for r in results:
        status_tag = "[PASS] PASSED" if r["passed"] else "[FAIL] FAILED"
        print(f"{r['name']:<58} {r['elapsed']:>6.2f}s | {status_tag}")
        if not r["passed"]:
            print(f"\n--- Output from {r['name']} ---\n{r['output']}\n{'-'*80}")

    print("=" * 80)
    if all_passed:
        print(f"ALL {len(GATES)} GATES PASSED CLEANLY in {total_elapsed:.2f}s! READY FOR COMMIT & PUSH.")
        print("=" * 80)
        sys.exit(0)
    else:
        print(f"QUALITY GATE FAILED in {total_elapsed:.2f}s. ABORTING COMMIT.")
        print("=" * 80)
        sys.exit(1)

if __name__ == "__main__":
    main()
