Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {1, 1, 0, 0.1};
Point(4) = {0, 1, 0, 0.1};

Line(1) = {1, 2};
Transfinite Line {1} = 3 Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = 3 Using Progression 1;
Line(3) = {3, 4};
Transfinite Line {3} = 3 Using Progression 1;
Line(4) = {4, 1};
Transfinite Line {4} = 3 Using Progression 1;

Line Loop(6) = {1, 2, 3, 4};

Plane Surface(6) = {6};


Recombine Surface {6};
Mesh.Algorithm = 8; // Delaunay for quads
Mesh.MshFileVersion = 2.2; // Mesh Format version 2.2
Mesh 2;
Physical Surface("FLUID") = {6};
Physical Line("WALL") = {1, 3}; // Solid Inviscid Wall (Euler) [TYPE 0]:
Physical Line("INLET") = {4}; // Supersonic Inlet BC or Far Field (Dirichlet BC) 
Physical Line("OUTLET") = {2}; // Supersonic Outlet BC (Neumann BC)
Coherence Mesh;
