#!/usr/bin/env python3
"""Combine VMD-rendered TGAs into a 4x4 grid for the adamantane series."""
from PIL import Image
import os

mols = ['adamantane', '1methyladamantane', '13dimethyladamantane', '135trimethyladamantane']
orbs = ['homo', 'lumo', 's1_occ', 's1_virt']

tile = (700, 700)
canvas_w = tile[0] * len(orbs)
canvas_h = tile[1] * len(mols)

canvas = Image.new('RGB', (canvas_w, canvas_h), 'white')

for i, mol in enumerate(mols):
    for j, orb in enumerate(orbs):
        fname = f'img_{orb}_{mol}.tga'
        if not os.path.exists(fname):
            print(f"missing: {fname}")
            continue
        im = Image.open(fname).convert('RGB').resize(tile)
        canvas.paste(im, (j * tile[0], i * tile[1]))

canvas.save('adamantane_nto_grid.png')
print('Wrote adamantane_nto_grid.png')

# Also save individual PNGs
for mol in mols:
    for orb in orbs:
        fname = f'img_{orb}_{mol}.tga'
        if os.path.exists(fname):
            Image.open(fname).convert('RGB').save(f'img_{orb}_{mol}.png')

print('Converted individual TGAs to PNG')