SetFactory("OpenCASCADE");
lc = 50;
z=0;
zmax=1;
xmin=-1700;
xmax=1700;
ymax=0;
ymin=-1000;
len_x = xmax-xmin;
len_y = ymax-ymin;
width_refine = 450;
s_main=news;
Rectangle(s_main) = {xmin,ymin,z, len_x, len_y};

//refine
p1=newp; Point(p1) = {(xmin+xmax)/2, ymin, z, lc};
p2=newp; Point(p2) = {(xmin+xmax)/2, ymax, z, lc};
l1=newl; Line(l1)={p1,p2};
Field[1] = Attractor;
Field[1].NNodesByEdge = 400;
Field[1].EdgesList = {l1};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 5;
Field[2].LcMax = lc;
Field[2].DistMin = width_refine;
Field[2].DistMax = width_refine*1.2;
Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

l() = Unique(Abs(Boundary{ Surface{s_main}; }));
p() = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc;
// extrude to 3d
Extrude {0, 0, zmax} {
Surface{s_main};
Layers{1};
Recombine;
}

Physical Volume("crust") = {1};
Physical Surface("frontAndBack") = {s_main, 6};
Physical Surface("bottom") = {2};
Physical Surface("left") = {5};
Physical Surface("right") = {3};
Physical Surface("top") = {4};


Color Red{Surface{2};}
Color Gray{Surface{s_main, 6};}
Color Blue{Surface{4};}
Color Green{Surface{3, 5};}