#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys

if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ALLOWED_EXTS = ('.js', '.mjs', '.ts', '.tsx', '.py')
EXCLUDE_DIRS = {
    '.git', '__pycache__', 'node_modules', '.gemini', '.agents', '.pytest_cache',
    'legacy-prototype', 'dist', 'build', '.idea', '.vscode'
}

def count_meaningful_file_loc(filepath):
    count = 0
    in_block_comment = False
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue
                if in_block_comment:
                    if '*/' in line:
                        in_block_comment = False
                        remaining = line.split('*/', 1)[1].strip()
                        if remaining and not remaining.startswith('//'):
                            count += 1
                    continue
                if line.startswith('/*'):
                    if '*/' in line:
                        remaining = line.split('*/', 1)[1].strip()
                        if remaining and not remaining.startswith('//'):
                            count += 1
                    else:
                        in_block_comment = True
                    continue
                if line.startswith('//') or line.startswith('#'):
                    continue
                count += 1
    except Exception as e:
        pass
    return count

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subsystem_counts = {}
    total_loc = 0
    file_count = 0

    for dp, dn, fn in os.walk(root_dir):
        rel_dir = os.path.relpath(dp, root_dir)
        parts = rel_dir.split(os.sep)
        # Exclude the root data/ folder (which holds raw FITS/BIN data), but DO NOT exclude src/data or tests/data
        if parts[0] == 'data':
            continue
        if any(p in EXCLUDE_DIRS for p in parts if p != '.'):
            continue

        for f in fn:
            if f.endswith(ALLOWED_EXTS) and f != 'count-meaningful-loc.py':
                filepath = os.path.join(dp, f)
                subsystem = parts[0] if parts[0] != '.' else 'root'
                if len(parts) > 1 and parts[0] in ('src', 'tests'):
                    subsystem = f"{parts[0]}/{parts[1]}"

                loc = count_meaningful_file_loc(filepath)
                if loc > 0:
                    subsystem_counts[subsystem] = subsystem_counts.get(subsystem, 0) + loc
                    total_loc += loc
                    file_count += 1

    print('=' * 68)
    print('ZRT COSMICFLOWS RESEARCH WORKBENCH - MEANINGFUL NON-HTML LOC AUDIT')
    print('=' * 68)
    print(f"{'Subsystem / Module Path':<42} | {'Meaningful LOC':>15}")
    print('-' * 68)
    for sub, count in sorted(subsystem_counts.items(), key=lambda x: x[0]):
        print(f"{sub:<42} | {count:>15,d}")
    print('-' * 68)
    print(f"{'TOTAL MEANINGFUL NON-HTML LOC':<42} | {total_loc:>15,d}")
    print(f"{'TOTAL NUMBER OF CODE FILES':<42} | {file_count:>15,d}")
    print('=' * 68)

    target_floor = 80000
    if total_loc >= target_floor:
        print(f'[PASS] LOC REQUIREMENT MET: {total_loc:,d} >= {target_floor:,d}')
        return 0
    else:
        print(f'[FAIL] LOC REQUIREMENT NOT YET MET: {total_loc:,d} < {target_floor:,d} (Remaining: {target_floor - total_loc:,d} LOC)')
        return 1

if __name__ == '__main__':
    sys.exit(main())
