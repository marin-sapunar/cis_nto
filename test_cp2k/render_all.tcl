# ============================================================================
# Render HOMO/LUMO/S1-NTO isosurfaces for the adamantane series.
#
# Output: one TGA per (molecule, orbital) pair, plus a combined 4x4 grid.
#
# Usage: vmd -dispdev text -e render_all.tcl
# ============================================================================

display projection Orthographic
axes location Off
color Display Background white
menu files off
menu graphics off
menu render off
menu color off

set iso 0.030
set color_pos 24   ;# positive lobe: orange-red
set color_neg 30   ;# negative lobe: steel blue

# Molecule tags (must match the ones in run_all.sh)
set mols {adamantane 1methyladamantane 13dimethyladamantane 135trimethyladamantane}

# Orbital tags
set orbs {homo lumo s1_occ s1_virt}

foreach mol $mols {
    foreach orb $orbs {
        set inp "cubes/${orb}_${mol}.cube"
        set out "img_${orb}_${mol}.tga"

        if {![file exists $inp]} {
            puts "SKIP: $inp not found"
            continue
        }

        puts "Rendering $inp -> $out"
        mol delete all
        mol new $inp type cube first 0 last -1 step 1 waitfor 1

        display resetview
        display resetview

        # CPK atom representation
        mol delrep 0 top
        mol addrep top
        mol modstyle 0 top CPK 1.0 0.3 12.0 12.0
        mol modcolor 0 top Name
        mol modmaterial 0 top Opaque

        # Positive isosurface
        mol addrep top
        mol modstyle 1 top Isosurface $iso 0 0 0 1 1
        mol modcolor 1 top ColorID $color_pos
        mol modmaterial 1 top Transparent

        # Negative isosurface
        mol addrep top
        mol modstyle 2 top Isosurface [expr {-$iso}] 0 0 0 1 1
        mol modcolor 2 top ColorID $color_neg
        mol modmaterial 2 top Transparent

        scale by 1.6

        material change opacity Transparent 0.55
        material change shininess Transparent 0.4

        render TachyonInternal $out -aasamples 24 -format TGA -res 1600 1600
        mol delete top
    }
}

puts "All renders complete."

# ---- Assemble the 4x4 grid (molecules as rows, orbitals as columns) ----
# Produces one PNG per molecule and one summary PNG if PIL is available
exec python3 make_grid.py

quit