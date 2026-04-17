module WSI2DMesh_warmup

using PreconditionedVLFS
using DrWatson
using GridapGmsh: gmsh
using PartitionedArrays

function create_mesh(ranks)

    path = mkpath(datadir("wsi_2d", "model"))

    mesh_file = "$path/mesh_wsi_2d_warmup.msh"
    
    map_main(ranks) do rank
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    gmsh.model.add("2D_domain")
    
    height = 1.0;
    mem_length = 1;
    Outter = 4 ;
    membrane_begin = 2 ; 
    membrane_end = membrane_begin + mem_length;

    meshPartFree = 1.0;
    meshPartFloat = 1.0;
    z_partitions = 1;

    # Points: gmsh.model.geo.addPoint(x, y, z, meshSize)
    p1 = gmsh.model.geo.addPoint(0.0,            0.0,    0.0, 1.0)
    p2 = gmsh.model.geo.addPoint(0.0,            height, 0.0, 1.0)
    p3 = gmsh.model.geo.addPoint(Outter,         0.0,    0.0, 1.0)
    p4 = gmsh.model.geo.addPoint(Outter,         height, 0.0, 1.0)
    p5 = gmsh.model.geo.addPoint(membrane_begin, height, 0.0, 1.0)
    p6 = gmsh.model.geo.addPoint(membrane_begin, 0.0,    0.0, 1.0)
    p7 = gmsh.model.geo.addPoint(membrane_end,   0.0,    0.0, 1.0)
    p8 = gmsh.model.geo.addPoint(membrane_end,   height, 0.0, 1.0)

    # Lines: gmsh.model.geo.addLine(start_point_tag, end_point_tag)
    l1  = gmsh.model.geo.addLine(p1, p2)
    l2  = gmsh.model.geo.addLine(p2, p5)
    l3  = gmsh.model.geo.addLine(p5, p8)
    l4  = gmsh.model.geo.addLine(p8, p4)
    l5  = gmsh.model.geo.addLine(p4, p3)
    l6  = gmsh.model.geo.addLine(p3, p7)
    l7  = gmsh.model.geo.addLine(p7, p6)
    l8  = gmsh.model.geo.addLine(p6, p1)
    l9  = gmsh.model.geo.addLine(p6, p5)
    l10 = gmsh.model.geo.addLine(p7, p8)

    # Surfaces (Curve Loops & Plane Surfaces)
    loop1 = gmsh.model.geo.addCurveLoop([l8, l1, l2, -l9])
    surf1 = gmsh.model.geo.addPlaneSurface([loop1])

    loop2 = gmsh.model.geo.addCurveLoop([l7, l9, l3, -l10])
    surf2 = gmsh.model.geo.addPlaneSurface([loop2])

    loop3 = gmsh.model.geo.addCurveLoop([l6, l10, l4, l5])
    surf3 = gmsh.model.geo.addPlaneSurface([loop3])

    # ---------------------------------------------------------
    # Transfinite curve definitions
    # ---------------------------------------------------------

    # Vertical partitions (Passing positive l5 since progression is 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l1,  z_partitions, "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(l9,  z_partitions, "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(l10, z_partitions, "Progression", 0.9)
    gmsh.model.geo.mesh.setTransfiniteCurve(-l5,  z_partitions, "Progression", 0.9)

    # Horizontal partitions
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, Int(round(membrane_begin * meshPartFree)), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l8, Int(round(membrane_begin * meshPartFree)), "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteCurve(l4, Int(round((Outter - membrane_end) * meshPartFree)), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l6, Int(round((Outter - membrane_end) * meshPartFree)), "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteCurve(l3, Int(round(mem_length * meshPartFloat)), "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(l7, Int(round(mem_length * meshPartFloat)), "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteSurface(surf1, [p1, p6, p5, p2])
    gmsh.model.geo.mesh.setTransfiniteSurface(surf2, [p6, p7, p8, p5])
    gmsh.model.geo.mesh.setTransfiniteSurface(surf3, [p7, p3, p4, p8])

    gmsh.model.geo.mesh.setRecombine(2, surf1)
    gmsh.model.geo.mesh.setRecombine(2, surf2)
    gmsh.model.geo.mesh.setRecombine(2, surf3)

    # Transfinite surface definitions

    # -----------------------------
    # Physical group tags and names
    # -----------------------------
    pg1 = gmsh.model.addPhysicalGroup(0, [p1, p6, p7, p3])
    gmsh.model.setPhysicalName(0, pg1, "Bed")

    pg2 = gmsh.model.addPhysicalGroup(1, [l8, l7, l6])
    gmsh.model.setPhysicalName(1, pg2, "Bed")

    pg3 = gmsh.model.addPhysicalGroup(1, [l3])
    gmsh.model.setPhysicalName(1, pg3, "FloatingSolid")

    pg4 = gmsh.model.addPhysicalGroup(1, [l2, l4])
    gmsh.model.setPhysicalName(1, pg4, "FreeSurface")

    pg5 = gmsh.model.addPhysicalGroup(0, [p2])
    gmsh.model.setPhysicalName(0, pg5, "LeftPoint")

    pg6 = gmsh.model.addPhysicalGroup(0, [p4])
    gmsh.model.setPhysicalName(0, pg6, "RightPoint")

    pg7 = gmsh.model.addPhysicalGroup(1, [l1])
    gmsh.model.setPhysicalName(1, pg7, "Inlet")

    pg8 = gmsh.model.addPhysicalGroup(1, [l5])
    gmsh.model.setPhysicalName(1, pg8, "Outlet")

    pg9 = gmsh.model.addPhysicalGroup(2, [surf1, surf3])
    gmsh.model.setPhysicalName(2, pg9, "FreeSurfaceBulk")

    pg10 = gmsh.model.addPhysicalGroup(2, [surf2])
    gmsh.model.setPhysicalName(2, pg10, "FloatingSolidBulk")

    pg11 = gmsh.model.addPhysicalGroup(1, [l9, l10])
    gmsh.model.setPhysicalName(1, pg11, "Domain") 

    pg12 = gmsh.model.addPhysicalGroup(2, [surf1, surf2, surf3])
    gmsh.model.setPhysicalName(2, pg12, "Domain")

    gmsh.model.geo.synchronize()

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Finalize
    gmsh.write(mesh_file)
    gmsh.finalize()
    end

    sum(ranks) # Adding this line to create a MPI Barrier using PartitinedArrays
    return mesh_file
end

end
