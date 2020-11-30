SetFactory("OpenCASCADE");
lc = 100;
z=0;
// zmax=1;
xmin=-1700;
xmax=1700;
ymax=0;
ymin=-1000;
len_x = xmax-xmin;
len_y = ymax-ymin;
s_main=news;
Rectangle(s_main) = {xmin,ymin,z, len_x, len_y};


// //points
// Point(1) = {xmin,ymin,z, lc};
// Point(2) = {xmax,ymin,z, lc} ;
// Point(3) = {xmax,ymax,z, lc} ;
// Point(4) = {xmin,ymax,z, lc} ;
// //new points 
// Point(5) = {xmin, y_pip_top, z, lc};
// Point(6) = {x_pip_left, y_pip_top, z, lc};
// Point(7) = {x_pip_left, ymin, z, lc};
// Point(8) = {x_pip_right, ymin, z, lc};
// Point(9) = {x_pip_right, y_pip_top, z, lc};
// Point(10) = {xmax, y_pip_top, z, lc};
// Point(11) = {x_pip_left, ymax, z, lc};
// Point(12) = {x_pip_right, ymax, z, lc};


// Line(1) = {7, 8};

// Line(2) = {8, 9};

// Line(3) = {9, 10};

// Line(4) = {10, 3};

// Line(5) = {3, 12};

// Line(6) = {12, 11};

// Line(7) = {11, 4};

// Line(8) = {4, 5};

// Line(9) = {5, 6};

// Line(10) = {6, 7};

// Line(11) = {11, 6};

// Line(12) = {12, 9};

// Line(13) = {6, 9};

// Line Loop(1) = {7, 8, 9, -11};

// Plane Surface(1) = {1};

// Line Loop(2) = {11, 13, -12, 6};

// Plane Surface(2) = {2};

// Line Loop(3) = {12, 3, 4, 5};

// Plane Surface(3) = {3};

// Line Loop(4) = {13, -2, -1, -10};

// Plane Surface(4) = {4};

//refine
Point(200) = {(xmin+xmax)/2, ymin, z, lc};
Point(201) = {(xmin+xmax)/2, ymax, z, lc};
Line(200)={200,201};
Field[1] = Attractor;
Field[1].NNodesByEdge = 400;
Field[1].EdgesList = {200};

// Field[2] = Threshold;
// Field[2].IField = 1;
// Field[2].LcMin = 2; // minimum cell size is 2m
// Field[2].LcMax = lc;
// Field[2].DistMin = (x_pip_right-x_pip_left)*1;
// Field[2].DistMax = (x_pip_right-x_pip_left)*3;

// Field[7] = Min;
// Field[7].FieldsList = {2};
// Background Field = 7;

l() = Unique(Abs(Boundary{ Surface{s_main}; }));
p() = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc;

// Transfinite Surface {4,2};
// Recombine Surface {4,2};

// Extrude {0, 0, zmax} {
// Surface{1,2,3,4};
// Layers{1};
// Recombine;
// }

// Physical Volume("layer2A") = {1,2,3};
// Physical Volume("layer2B") = {4};
// Physical Surface("frontAndBack") = {222,244,266,288,  1,2,4,3};
// Physical Surface("wall") = {217,257};
// Physical Surface("bottom") = {283};
// Physical Surface("left") = {213,287};
// Physical Surface("top") = {209,243,265};
// Physical Surface("right") = {261,279};