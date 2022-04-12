SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 4, 3, 0};
// refine
l() = Unique(Abs(Boundary{ Surface{1}; }));
p() = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = 1;
// make structured mesh
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("main") = {1};
Physical Line("bottom") = {1};
Physical Line("right") = {2};
Physical Line("top") = {3};
Physical Line("left") = {4};