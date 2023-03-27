// Gen mesh with:
// gmsh rect.geo -parse_and_exit

lc = 0.05;

Point (1) = {-0.1, 0, 0, lc};
Point (2) = {0, 0, 0, lc};
Point (3) = {1, 0, 0, lc};
Point (4) = {1, 0.2, 0, lc};
Point (5) = {0, 0.2, 0, lc};
Point (6) = {-0.1, 0.2, 0, lc};

Line (1) = {1, 2};
Line (2) = {2, 5};
Line (3) = {5, 6};
Line (4) = {6, 1};
Line (5) = {2, 3};
Line (6) = {3, 4};
Line (7) = {4, 5};

Curve Loop (1) = {1, 2, 3, 4};
Curve Loop (2) = {5, 6, 7, -2};

Plane Surface (1) = {1};
Plane Surface (2) = {2};

Physical Curve("top") = {3, 7};
Physical Curve("bot0") = {1};
Physical Curve("bot1") = {5};
Physical Curve("left") = {4};
Physical Curve("right") = {6};
Physical Surface("internal") = {1, 2};

Transfinite Curve {5, -7} = 60 Using Progression 1.05;   // 90
Transfinite Curve {-1, 3} = 20 Using Progression 1.1;   // 30
Transfinite Curve {2, -4, 6} = 90 Using Progression 1.08;  // 110

Transfinite Surface {1};
Transfinite Surface {2};

Recombine Surface (1);
Recombine Surface (2);

// Mesh the surface
Mesh 2;

Save "flat_plate.msh";