module WSI3DMesh1

using PreconditionedVLFS
using DrWatson
using GridapGmsh: gmsh
using PartitionedArrays

function create_mesh(ranks, height::Float64)

    path = mkpath(datadir("wsi_3d", "model"))

    mesh_file = "$path/mesh_wsi_3d_1.msh"

    map_main(ranks) do rank

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("3D_domain")

    # Outter box dimensions
    length_x = 9 * pi * height
    length_y = pi * height
    length_z = height

    # Membrane dimensions
    length_mem = 2 * pi * height
    width = 0.5 * pi * height

    # Mesh partitioning per unit length
    float_x_mesh = 3
    float_y_mesh = 3
    top_x_mesh = 3
    top_y_mesh = 3
    bottom_x_mesh = 3
    bottom_y_mesh = 3
    sides_mesh = 3

    # Membrane origin coordinates
    membrane_begin = length_x / 2 - length_mem / 2
    membrane_end = membrane_begin + length_mem
    membrane_left = length_y / 2 - width / 2
    membrane_right = membrane_left + width

    # Points: gmsh.model.geo.addPoint(x, y, z, meshSize)
    p1 = gmsh.model.geo.addPoint(membrane_begin, membrane_left, length_z, 1.0)
    p2 = gmsh.model.geo.addPoint(membrane_end,   membrane_left, length_z, 1.0)
    p3 = gmsh.model.geo.addPoint(membrane_end,   membrane_right, length_z, 1.0)
    p4 = gmsh.model.geo.addPoint(membrane_begin, membrane_right, length_z, 1.0)
    p5  = gmsh.model.geo.addPoint(0.0,      0.0,      length_z, 1.0)
    p6  = gmsh.model.geo.addPoint(length_x, 0.0,      length_z, 1.0)
    p7  = gmsh.model.geo.addPoint(length_x, length_y, length_z, 1.0)
    p8  = gmsh.model.geo.addPoint(0.0,      length_y, length_z, 1.0)
    p9  = gmsh.model.geo.addPoint(0.0,      length_y, 0.0,      1.0)
    p10 = gmsh.model.geo.addPoint(0.0,      0.0,      0.0,      1.0)
    p11 = gmsh.model.geo.addPoint(length_x, 0.0,      0.0,      1.0)
    p12 = gmsh.model.geo.addPoint(length_x, length_y, 0.0,      1.0)

    # Lines: gmsh.model.geo.addLine(start_point_tag, end_point_tag)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5  = gmsh.model.geo.addLine(p5, p6)
    l6  = gmsh.model.geo.addLine(p6, p7)
    l7  = gmsh.model.geo.addLine(p7, p8)
    l8  = gmsh.model.geo.addLine(p8, p5)
    l9  = gmsh.model.geo.addLine(p5, p10)
    l10 = gmsh.model.geo.addLine(p10, p11)
    l11 = gmsh.model.geo.addLine(p11, p12)
    l12 = gmsh.model.geo.addLine(p12, p9)
    l13 = gmsh.model.geo.addLine(p9, p10)
    l14 = gmsh.model.geo.addLine(p6, p11)
    l15 = gmsh.model.geo.addLine(p7, p12)
    l16 = gmsh.model.geo.addLine(p8, p9)

    # Surfaces (Curve Loops & Plane Surfaces)
    loop1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf1 = gmsh.model.geo.addPlaneSurface([loop1])

    loop2 = gmsh.model.geo.addCurveLoop([l8, l5, l6, l7])
    surf2 = gmsh.model.geo.addPlaneSurface([loop2, loop1])

    loop3 = gmsh.model.geo.addCurveLoop([l13, l10, l11, l12])
    surf3 = gmsh.model.geo.addPlaneSurface([loop3])

    loop4 = gmsh.model.geo.addCurveLoop([l14, l11, -l15, -l6])
    surf4 = gmsh.model.geo.addPlaneSurface([loop4])

    loop5 = gmsh.model.geo.addCurveLoop([l12, -l16, -l7, l15])
    surf5 = gmsh.model.geo.addPlaneSurface([loop5])

    loop6 = gmsh.model.geo.addCurveLoop([l13, -l9, -l8, l16])
    surf6 = gmsh.model.geo.addPlaneSurface([loop6])

    loop7 = gmsh.model.geo.addCurveLoop([l10, -l14, -l5, l9])
    surf7 = gmsh.model.geo.addPlaneSurface([loop7])

    # Volume (Surface Loops & volumes)
    shell1 = gmsh.model.geo.addSurfaceLoop([surf7, surf3, surf6, surf2, surf4, surf5, surf1])
    vol1 = gmsh.model.geo.addVolume([shell1])

    # ---------------------------------------------------------
    # Transfinite curve definitions
    # ---------------------------------------------------------

    # Vertical partitions (Passing positive l5 since progression is 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, max(2, Int(round(float_x_mesh * length_mem))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, max(2, Int(round(float_x_mesh * length_mem))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, max(2, Int(round(float_y_mesh * width))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, max(2, Int(round(float_y_mesh * width))), "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteCurve(l5, max(2, Int(round(length_x * top_x_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l7, max(2, Int(round(length_x * top_x_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l8, max(2, Int(round(length_y * top_y_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l6, max(2, Int(round(length_y * top_y_mesh))), "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteCurve(l9,  max(2, Int(round(length_z * sides_mesh))), "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(l14, max(2, Int(round(length_z * sides_mesh))), "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(l15, max(2, Int(round(length_z * sides_mesh))), "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(l16, max(2, Int(round(length_z * sides_mesh))), "Progression", 0.9)

    gmsh.model.geo.mesh.setTransfiniteCurve(l10, max(2, Int(round(length_x * bottom_x_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l12, max(2, Int(round(length_x * bottom_x_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l13, max(2, Int(round(length_y * bottom_y_mesh))), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l11, max(2, Int(round(length_y * bottom_y_mesh))), "Progression", 1.0)

    # Transfinite surface definitions

    # -----------------------------
    # Physical group tags and names
    # -----------------------------
    pg1 = gmsh.model.addPhysicalGroup(2, [surf1])
    gmsh.model.setPhysicalName(2, pg1, "FloatingSolid")

    pg2 = gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4])
    gmsh.model.setPhysicalName(1, pg2, "FloatingSolid")

    pg3 = gmsh.model.addPhysicalGroup(0, [p1, p2, p3, p4])
    gmsh.model.setPhysicalName(0, pg3, "FloatingSolid")

    pg4 = gmsh.model.addPhysicalGroup(2, [surf2])
    gmsh.model.setPhysicalName(2, pg4, "FreeSurface")

    pg5 = gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l7])
    gmsh.model.setPhysicalName(1, pg5, "FreeSurface")

    pg6 = gmsh.model.addPhysicalGroup(0, [p1, p2, p3, p4])
    gmsh.model.setPhysicalName(0, pg6, "FreeSurface")

    pg7 = gmsh.model.addPhysicalGroup(2, [surf4])
    gmsh.model.setPhysicalName(2, pg7, "Outlet")

    pg8 = gmsh.model.addPhysicalGroup(2, [surf6])
    gmsh.model.setPhysicalName(2, pg8, "Inlet")

    pg9 = gmsh.model.addPhysicalGroup(2, [surf7, surf5])
    gmsh.model.setPhysicalName(2, pg9, "SideWalls")

    pg10 = gmsh.model.addPhysicalGroup(1, [l9, l14, l15, l16, l5, l7])
    gmsh.model.setPhysicalName(1, pg10, "SideWalls")

    pg11 = gmsh.model.addPhysicalGroup(2, [surf3])
    gmsh.model.setPhysicalName(2, pg11, "Bed")

    pg12 = gmsh.model.addPhysicalGroup(1, [l13, l11])
    gmsh.model.setPhysicalName(1, pg12, "Bed")

    pg13 = gmsh.model.addPhysicalGroup(1, [l12, l10])
    gmsh.model.setPhysicalName(1, pg13, "BedInt")

    pg14 = gmsh.model.addPhysicalGroup(0, [p9, p11, p10, p12])
    gmsh.model.setPhysicalName(0, pg14, "BedInt")

    pg15 = gmsh.model.addPhysicalGroup(1, [l8])
    gmsh.model.setPhysicalName(1, pg15, "LeftPoint")

    pg16 = gmsh.model.addPhysicalGroup(0, [p5, p8])
    gmsh.model.setPhysicalName(0, pg16, "LeftPoint")

    pg17 = gmsh.model.addPhysicalGroup(1, [l6])
    gmsh.model.setPhysicalName(1, pg17, "RightPoint")

    pg18 = gmsh.model.addPhysicalGroup(0, [p6, p7])
    gmsh.model.setPhysicalName(0, pg18, "RightPoint")

    pg19 = gmsh.model.addPhysicalGroup(3, [vol1])
    gmsh.model.setPhysicalName(3, pg19, "Domain")

    gmsh.model.geo.synchronize()

    # Generate 3D mesh
    gmsh.model.mesh.generate(3)
    # Finalize
    gmsh.write(mesh_file)
    gmsh.finalize()
    end
    sum(ranks) # Adding this line to create a MPI Barrier using PartitinedArrays
    return mesh_file
end
end
