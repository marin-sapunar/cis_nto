#!/usr/bin/env python3
"""
Extract a single MO block from a molden file.

Usage: extract_nto.py source.molden output.molden <selector>

Selector is either a 1-based block index or a Sym tag (e.g. "2001a").

If the selector parses as an integer in [1, nblocks], that block is used.
Otherwise a Sym-tag lookup is attempted; if exactly one block matches, it
is used; if more than one matches, the program aborts (ambiguous Sym tags
appear in NTO molden files where the same label is emitted once per
excited state).
"""
import sys


def find_blocks(lines):
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

    header = []
    for i, l in enumerate(lines):
        header.append(l)
        if l.strip() == "[MO]":
            break

    blocks = find_blocks(lines)
    print(f"Found {len(blocks)} MO blocks in {src}")

    chosen = None

    try:
        idx = int(selector)
        if 1 <= idx <= len(blocks):
            chosen = blocks[idx - 1]
            print(f"Matched by index {selector}")
    except ValueError:
        pass

    if chosen is None:
        matches = [b for b in blocks if get_sym(lines, b[0], b[1]) == selector]
        if len(matches) == 1:
            chosen = matches[0]
            print(f"Matched by Sym={selector}")
        elif len(matches) > 1:
            raise SystemExit(
                f"Selector {selector!r} matches {len(matches)} blocks in {src}. "
                f"Pass a 1-based block index instead."
            )

    if chosen is None:
        raise SystemExit(f"Selector {selector!r} not found in {src}")

    s, e = chosen
    with open(dst, "w") as f:
        f.writelines(header)
        f.writelines(lines[s:e])
    print(f"Wrote {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: extract_nto.py source.molden output.molden <index|sym>")
        sys.exit(1)
    extract(sys.argv[1], sys.argv[2], sys.argv[3])
