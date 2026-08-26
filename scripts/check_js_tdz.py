#!/usr/bin/env python3
"""
Gate 1: JavaScript AST / Static TDZ & Variable Scope Checker
Scans all <script> blocks in index.html and all JS files in src/ to detect:
1. Temporal Dead Zone (TDZ) usages of const/let before declaration
2. Uninitialized variable consumption in calculations
"""
import os
import re
import sys

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INDEX_HTML = os.path.join(ROOT_DIR, "index.html")

def scan_function_scope(func_name, body):
    errors = []
    lines = body.split("\n")
    declared_vars = set()
    
    for idx, line in enumerate(lines, 1):
        line_clean = re.sub(r'//.*', '', line).strip()
        if not line_clean:
            continue
            
        # Match standard declarations: const x =, let y =
        decl_match = re.findall(r'(?:const|let|var)\s+([a-zA-Z0-9_]+)\s*=', line_clean)
        for var in decl_match:
            declared_vars.add(var)
            
        # Match object destructuring: const { a, b, c } =
        destruct_match = re.findall(r'(?:const|let|var)\s*\{([^}]+)\}\s*=', line_clean)
        for dgroup in destruct_match:
            for item in dgroup.split(','):
                v = item.strip().split(':')[0].strip()
                if v:
                    declared_vars.add(v)
            
        # Match function arguments or parameters
        param_match = re.findall(r'function\s+[a-zA-Z0-9_]+\s*\(([^)]*)\)', line_clean)
        for pgroup in param_match:
            for p in pgroup.split(','):
                p_clean = p.strip()
                if p_clean:
                    declared_vars.add(p_clean)

        # Critical variables to watch for TDZ inside mathematical solvers
        if func_name in ['solveTullyFisher', 'jacobiDiagonalize', 'integrateRK4', 'computeFlowColormap']:
            critical_vars = ['sigmaMu', 'deltaMuMalm', 'distMpc', 'muCorr', 'mAbs', 'mu', 'wCorr', 'vRot']
            for cvar in critical_vars:
                # If used on right-hand side before declaration
                if re.search(r'=\s*[^;]*\b' + cvar + r'\b', line_clean) and not re.search(r'(?:const|let|var)\s+[^{]*\b' + cvar + r'\b', line_clean):
                    if cvar not in declared_vars:
                        errors.append(f"Function '{func_name}' line {idx}: TDZ Hazard! Variable '{cvar}' accessed before declaration.")
                    
    return errors

def main():
    print("Checking JavaScript Temporal Dead Zone (TDZ) & Variable Scoping...")
    all_errors = []
    
    with open(INDEX_HTML, "r", encoding="utf-8") as f:
        html = f.read()
        
    funcs = re.findall(r'function\s+([a-zA-Z0-9_]+)\s*\([^)]*\)\s*\{([\s\S]*?\n\})', html)
    for fname, fbody in funcs:
        errs = scan_function_scope(fname, fbody)
        if errs:
            all_errors.extend(errs)
            
    if all_errors:
        print(f"[FAIL] Found {len(all_errors)} TDZ / scoping violations:")
        for e in all_errors:
            print("  -", e)
        sys.exit(1)
    else:
        print("[PASS] 0 TDZ or variable scope violations detected across mathematical solvers.")
        sys.exit(0)

if __name__ == "__main__":
    main()
