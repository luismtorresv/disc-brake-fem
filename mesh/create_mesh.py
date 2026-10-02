import gmsh
import meshio

# Example values
ri, ro, e = (
    0.07,  # Inner radius
    0.15,  # Outer radius
    0.01,  # Half-thickness [m]
)
ra, rb = 0.08, 0.14  # Start and end of the friction track [m]
h_fino, h_grueso = 0.0004, 0.004  # Element size near and far from the track [m]

gmsh.initialize()
gmsh.model.add("disc")
occ = gmsh.model.geo

p = [
    occ.addPoint(x, z, 0)
    for x, z in [(ri, 0), (ro, 0), (ro, e), (rb, e), (ra, e), (ri, e)]
]
l = [occ.addLine(p[i], p[(i + 1) % 6]) for i in range(6)]
# l[0]: z=0 (Gamma_N) | l[1]: r=ro | l[2]: upper right face | l[3]: friction track | l[4]: upper left face | l[5]: r=ri
s = occ.addPlaneSurface([occ.addCurveLoop(l)])
occ.synchronize()

gmsh.model.addPhysicalGroup(1, [l[3]], 1, "brake")
gmsh.model.addPhysicalGroup(1, [l[1], l[2], l[4], l[5]], 2, "conv")
gmsh.model.addPhysicalGroup(1, [l[0]], 3, "neumann")
gmsh.model.addPhysicalGroup(2, [s], 4, "disc")

# Element size: fine near the friction track, coarsening toward the interior
f = gmsh.model.mesh.field
f.add("Distance", 1)
f.setNumbers(1, "CurvesList", [l[3]])
f.setNumber(1, "Sampling", 200)
f.add("Threshold", 2)
f.setNumber(2, "InField", 1)
f.setNumber(2, "SizeMin", h_fino)
f.setNumber(2, "SizeMax", h_grueso)
f.setNumber(2, "DistMin", 0.001)
f.setNumber(2, "DistMax", e)
f.setAsBackgroundMesh(2)

gmsh.model.mesh.generate(2)
gmsh.write("disc.msh")
gmsh.finalize()
