// Gmsh project created on Mon Sep 28 10:10:31 2026
SetFactory("OpenCASCADE");
// Inputs Box length
Lref = 35;
length_x = 1000/Lref;
length_y = 400/Lref;
length_z = 50/Lref;
// Inputs (membrane)
radius = 35/Lref;
space = 100/Lref; // diameter plus spacing
//+
// Inputs - Partition per unit length
mesh_sides       = 0.30*Lref;
mesh_fs_x_top    = 1.00*Lref;
mesh_fs_x_bottom = 0.5*Lref;
mesh_fs_y_top    = 1.00*Lref;
mesh_fs_y_bottom = 0.5*Lref;
mesh_mem         = 1.00*Lref;
//+
// Outter box points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, length_y, 0, 1.0};
Point(3) = {0, length_y, length_z, 1.0};
Point(4) = {length_x, length_y, length_z, 1.0};
Point(5) = {length_x, 0, 0, 1.0};
Point(6) = {0, 0, length_z, 1.0};
Point(7) = {length_x, 0, length_z, 1.0};
Point(8) = {length_x, length_y, 0, 1.0};
//+
// Outter box lengths
Line(1) = {1, 2};
Line(2) = {2, 8};
Line(3) = {8, 5};
Line(4) = {5, 1};
Line(5) = {1, 6};
Line(6) = {6, 3};
Line(7) = {3, 2};
Line(8) = {3, 4};
Line(9) = {4, 8};
Line(10) = {4, 7};
Line(11) = {7, 5};
Line(12) = {7, 6};
//+
Circle(13) = {length_x/2, length_y/2 - space/2, length_z, radius, 0, 2*Pi};
Circle(14) = {length_x/2, length_y/2 + space/2, length_z, radius, 0, 2*Pi};
Circle(15) = {length_x/2 + space, length_y/2 - space/2, length_z, radius, 0, 2*Pi};
Circle(16) = {length_x/2 + space, length_y/2 + space/2, length_z, radius, 0, 2*Pi};
Circle(17) = {length_x/2 - space, length_y/2 - space/2, length_z, radius, 0, 2*Pi};
Circle(18) = {length_x/2 - space, length_y/2 + space/2, length_z, radius, 0, 2*Pi};
Circle(19) = {length_x/2 + 2*space, length_y/2 - space/2, length_z, radius, 0, 2*Pi};
Circle(20) = {length_x/2 + 2*space, length_y/2 + space/2, length_z, radius, 0, 2*Pi};
Circle(21) = {length_x/2 - 2*space, length_y/2 - space/2, length_z, radius, 0, 2*Pi};
Circle(22) = {length_x/2 - 2*space, length_y/2 + space/2, length_z, radius, 0, 2*Pi};
// Surfaces
//+
Curve Loop(1) = {21};
Plane Surface(1) = {1};
Curve Loop(2) = {17};
Plane Surface(2) = {2};
Curve Loop(3) = {13};
Plane Surface(3) = {3};
Curve Loop(4) = {15};
Plane Surface(4) = {4};
Curve Loop(5) = {19};
Plane Surface(5) = {5};
Curve Loop(6) = {20};
Plane Surface(6) = {6};
Curve Loop(7) = {16};
Plane Surface(7) = {7};
Curve Loop(8) = {14};
Plane Surface(8) = {8};
Curve Loop(9) = {18};
Plane Surface(9) = {9};
Curve Loop(10) = {22};
Plane Surface(10) = {10};
Curve Loop(11) = {12, 6, 8, 10};
Curve Loop(12) = {21};
Curve Loop(13) = {17};
Curve Loop(14) = {13};
Curve Loop(15) = {15};
Curve Loop(16) = {19};
Curve Loop(17) = {20};
Curve Loop(18) = {16};
Curve Loop(19) = {14};
Curve Loop(20) = {18};
Curve Loop(21) = {22};
//+
Plane Surface(11) = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
Curve Loop(22) = {4, 5, -12, 11};
Plane Surface(12) = {22};
Curve Loop(23) = {1, -7, -6, -5};
Plane Surface(13) = {23};
Curve Loop(24) = {2, -9, -8, 7};
Plane Surface(14) = {24};
Curve Loop(25) = {3, -11, -10, 9};
Plane Surface(15) = {25};
Curve Loop(26) = {4, 1, 2, 3};
Plane Surface(16) = {26};
//+
Surface Loop(1) = {12, 16, 13, 14, 15, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Volume(1) = {1};
//+
Transfinite Curve {5, 11, 9, 7} = mesh_sides*length_z Using Progression 1.2; // Sides
Transfinite Curve {6, 10} = mesh_fs_y_top*length_y Using Progression 1; // freesurface_y_top
Transfinite Curve {1, 3} = mesh_fs_y_bottom*length_y Using Progression 1; // freesurface_y_bottom
Transfinite Curve {12, 8} = mesh_fs_x_top*length_x Using Progression 1; // freesurface_x_top
Transfinite Curve {4, 2} = mesh_fs_x_bottom*length_x Using Progression 1; // freesurface_x_bottom
Transfinite Curve {21, 22, 18, 17, 13, 14, 16, 15, 19, 20} = mesh_mem*2*Pi*radius Using Progression 1; // floating_membrane
//+
Physical Surface("FloatingSolid", 37) = {6, 7, 8, 9, 10, 1, 2, 3, 4, 5};
Physical Curve("FloatingSolid", 38) = {20, 16, 14, 18, 22, 21, 17, 13, 15, 19};
//+
Physical Surface("FreeSurface", 39) = {17};
Physical Curve("FreeSurface", 39) = {19, 20, 16, 15, 13, 14, 18, 17, 21, 22, 12, 8};
//+
Physical Surface("Inlet", 40) = {13};
Physical Surface("Outlet", 41) = {15};
//+
Physical Surface("SideWalls", 42) = {12, 14};
Physical Curve("SideWalls", 43) = {5, 12, 11, 7, 8, 9};
//+
Physical Surface("Bed", 44) = {16};
Physical Curve("Bed", 45) = {1,3};
Physical Curve("BedInt", 46) = {4, 2};
Physical Point("BedInt", 47) = {1, 2, 8, 5};
//+
Physical Point("LeftPoint", 48) = {6, 3};
Physical Curve("LeftPoint", 49) = {6};
//+
Physical Point("RightPoint", 50) = {7, 4};
Physical Curve("RightPoint", 51) = {10};
//+
Physical Volume("Domain", 52) = {1};
