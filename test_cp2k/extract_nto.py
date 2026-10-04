#!/usr/bin/env python3
"""
Extract a single MO block from a molden file.

Usage: extract_nto.py source.molden output.molden <selector>

Selector is either a Sym tag (e.g. "2001a") or a 1-based index.

Block delimiters:
  - If the file has Sym= lines, each Sym= starts a new block.
  - Otherwise each Ene= line starts a new block.
"""
import sys


def find_blocks(lines):
    """Return list of (start, end) index pairs for each MO block."""
    mo_start = None
    for i, l in enumerate(lines):
        if l.strip() == "[MO]":
            mo_start = i
            break
    if mo_start is None:
        return []

    sym_lines = [i for i in range(mo_start, len(lines))
                 if lines[i].strip().startswith("Sym=")]
    ene_lines = [i for i in range(mo_start, len(lines))
                 if lines[i].strip().startswith("Ene=")]

    delims = sym_lines if sym_lines else ene_lines
    if not delims:
        return []

    blocks = []
    for k, s in enumerate(delims):
        e = delims[k + 1] if k + 1 < len(delims) else len(lines)
        blocks.append((s, e))
    return blocks


def get_sym(lines, start, end):
    for j in range(start, min(start + 5, end)):
        if lines[j].strip().startswith("Sym="):
            return lines[j].split()[1]
    return None


def extract(src, dst, selector):
    lines = open(src).readlines()

    # Header: everything up to and including the first [MO]
    header = []
    for i, l in enumerate(lines):
        header.append(l)
        if l.strip() == "[MO]":
            break

    blocks = find_blocks(lines)
    print(f"Found {len(blocks)} MO blocks in {src}")

    chosen = None
    for s, e in blocks:
        sym = get_sym(lines, s, e)
        if sym == selector:
            chosen = (s, e)
            print(f"Matched by Sym={selector}")
            break

    if chosen is None:
        try:
            idx = int(selector)
            if 1 <= idx <= len(blocks):
                chosen = blocks[idx - 1]
                print(f"Matched by index {selector}")
        except ValueError:
            pass

    if chosen is None:
        raise SystemExit(f"Selector {selector!r} not found in {src}")

    s, e = chosen
    with open(dst, "w") as f:
        f.writelines(header)
        f.writelines(lines[s:e])
    print(f"Wrote {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: extract_nto.py source.molden output.molden <sym|index>")
        sys.exit(1)
    extract(sys.argv[1], sys.argv[2], sys.argv[3])
