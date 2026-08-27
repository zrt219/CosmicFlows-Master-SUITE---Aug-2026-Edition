#!/usr/bin/env python3
"""
Gate 6: KaTeX In-Browser Mathematical Formula Linter
Validates all LaTeX formulas in README.md and markdown files for KaTeX 0-error parse rate.
"""
import os
import sys
import json
import re
import socket

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

from tests.cdp_client import ChromeCDPClient
from tests.conftest import ensure_http_server

def find_free_port(start_port=9500):
    for p in range(start_port, start_port + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    return start_port

def main():
    print("Testing KaTeX Mathematical Syntax across README.md...")
    ensure_http_server(8000)
    cdp_port = find_free_port(9500)
    client = ChromeCDPClient(port=cdp_port, url="http://localhost:8000", spawn_headless=True)
    client.start()
    
    try:
        client.evaluate("""
        (function() {
            if (typeof katex === 'undefined') {
                const s = document.createElement('script');
                s.src = 'https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.js';
                document.head.appendChild(s);
            }
        })()
        """)
        ready = client.wait_for_condition("typeof katex !== 'undefined'", timeout=15.0)
        if not ready:
            raise RuntimeError("Timed out waiting for KaTeX library to load in browser")
            
        readme_path = os.path.join(ROOT_DIR, "README.md")
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()
            
        display_eqs = re.findall(r'\$\$([\s\S]*?)\$\$', content)
        display_eqs_clean = [eq.strip() for eq in display_eqs if eq.strip()]
        
        # Also extract inline math (single dollar)
        content_no_display = re.sub(r'\$\$[\s\S]*?\$\$', '', content)
        inline_eqs = re.findall(r'(?<!\$)\$(?!\$)([^\$\n]+?)(?<!\$)\$(?!\$)', content_no_display)
        inline_eqs_clean = [eq.strip() for eq in inline_eqs if eq.strip()]

        eqs_json = json.dumps(display_eqs_clean)
        inline_json = json.dumps(inline_eqs_clean)
        
        js_eval = f"""
        (function() {{
            if (typeof katex === 'undefined') return {{ available: false }};
            const displayEqs = {eqs_json};
            const inlineEqs = {inline_json};
            let errors = [];
            displayEqs.forEach((eq, idx) => {{
                try {{
                    katex.renderToString(eq, {{ displayMode: true, throwOnError: true }});
                }} catch (e) {{
                    errors.push({{ type: 'display', index: idx, eq: eq, error: e.message }});
                }}
            }});
            inlineEqs.forEach((eq, idx) => {{
                try {{
                    katex.renderToString(eq, {{ displayMode: false, throwOnError: true }});
                }} catch (e) {{
                    errors.push({{ type: 'inline', index: idx, eq: eq, error: e.message }});
                }}
            }});
            return {{ available: true, displayTotal: displayEqs.length, inlineTotal: inlineEqs.length, errors: errors }};
        }})()
        """
        res = client.evaluate(js_eval)
    finally:
        client.close()
    
    if not res:
        print("[FAIL] Failed to evaluate KaTeX render in browser.")
        sys.exit(1)
        
    if res.get('available') and len(res.get('errors', [])) > 0:
        print(f"[FAIL] Found {len(res['errors'])} KaTeX parse errors:")
        for err in res['errors']:
            print("  -", err)
        sys.exit(1)
    else:
        print(f"[PASS] {len(display_eqs_clean)} display + {len(inline_eqs_clean)} inline equations parsed with 0 KaTeX errors.")
        sys.exit(0)

if __name__ == "__main__":
    main()
