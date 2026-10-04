cat > /home/hamid/cis_nto/test_cp2k/render_nto.tcl << 'TCLEOF'
# Render adamantane HOMO/LUMO/NTOs as publication-quality isosurfaces.
# Usage: vmd -dispdev text -e render_nto.tcl

display projection Orthographic
axes location Off
color Display Background white

set iso 0.030
set color_pos 24
set color_neg 30

foreach pair {
    {homo.cube     homo.tga     "S0 HOMO"}
    {lumo.cube     lumo.tga     "S0 LUMO"}
    {s1_occ.cube   s1_occ.tga   "S1 occupied NTO"}
    {s1_virt.cube  s1_virt.tga  "S1 virtual NTO"}
} {
    lassign $pair inp out label

    if {![file exists $inp]} {
        puts "WARNING: $inp not found, skipping."
        continue
    }

    puts "Processing $inp -> $out"

    mol delete all
    mol new $inp type cube first 0 last -1 step 1 waitfor 1

    display resetview
    display resetview

    mol delrep 0 top
    mol addrep top
    mol modstyle 0 top CPK 1.0 0.3 12.0 12.0
    mol modcolor 0 top Name
    mol modmaterial 0 top Opaque

    mol addrep top
    mol modstyle 1 top Isosurface $iso 0 0 0 1 1
    mol modcolor 1 top ColorID $color_pos
    mol modmaterial 1 top Transparent

    mol addrep top
    mol modstyle 2 top Isosurface [expr {-$iso}] 0 0 0 1 1
    mol modcolor 2 top ColorID $color_neg
    mol modmaterial 2 top Transparent

    scale by 1.6

    material change opacity Transparent 0.55
    material change shininess Transparent 0.4

    render TachyonInternal $out -aasamples 24 -format TGA -res 2048 2048

    mol delete top
}

puts "All images rendered."
quit
TCLEOF