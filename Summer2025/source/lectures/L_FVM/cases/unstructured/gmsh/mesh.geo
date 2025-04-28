lc = 0.1;
zmin=0;
zmax=0.05;
width=0.25;
x0=0;
y0=0;

Point(1) = {x0, y0, zmin, lc};
Point(2) = {x0-width, y0, zmin, lc};
Point(3) = {x0-width, y0+width, zmin, lc};
Point(4) = {x0+width, y0+width, zmin, lc};
Point(5) = {x0+width, y0-width, zmin, lc};
Point(6) = {x0, y0-width, zmin, lc};
Line(1) = {6, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 6};
//+
Curve Loop(1) = {4, 5, 6, 1, 2, 3};
//+
Plane Surface(1) = {1};

// Transfinite Surface {1};
// Recombine Surface {1};

Extrude {0, 0, zmax} {
Surface{1};
Layers{1};
Recombine;
}

Physical Volume("main") = {1};
Physical Surface("frontAndBack") = {1, 38};
Physical Surface("left") = {37, 33, 29};
Physical Surface("right") = {21};
Physical Surface("top") = {17};
Physical Surface("bottom") = {25};

Color Gray{Surface{1, 38};}
Color Green{Surface{37, 33, 29};}
Color Blue{Surface{21};}
Color Purple{Surface{17};}
Color Yellow{Surface{25};}
