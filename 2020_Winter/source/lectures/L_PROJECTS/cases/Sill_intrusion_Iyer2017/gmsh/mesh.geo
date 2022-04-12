SetFactory("OpenCASCADE");
xmin=0;
xmax=50000;
ymin=-9000;
ymax=0;
zmin=0;
zmax=1000;
lc=200;
// s_main=news;
// Rectangle(s_main) = {xmin, ymin,zmin, xmax-xmin,ymax-ymin,zmin};
Point(1) = {xmin,ymin,zmin,lc};
Point(2) = {xmax,ymin,zmin,lc};
Point(3) = {xmax,ymax,zmin,lc};
Point(4) = {xmin,ymax,zmin,lc};
Line(1)={1,2};
Line(2)={3,4};
// 1. Tertiary layer
Include "Tertiary.geo";
l_tmp1=newl; Line(l_tmp1)={4, pts_Tertiary[0]};
l_tmp2=newl; Line(l_tmp2)={3, pts_Tertiary[#pts_Tertiary[]-1]};
s_Tertiary=news; Curve Loop(s_Tertiary) = {l_tmp1, lines_Tertiary[], -l_tmp2, 2};
Plane Surface(s_Tertiary) = {s_Tertiary};

// 2. Upper Cretaceous layer
Include "UpperCretaceous.geo";
l_tmp1=newl; Line(l_tmp1)={pts_Tertiary[0], pts_UpperCretaceous[0]};
l_tmp2=newl; Line(l_tmp2)={pts_Tertiary[#pts_Tertiary[]-1], pts_UpperCretaceous[#pts_UpperCretaceous[]-1]};
s_UpperCretaceous=news; Curve Loop(s_UpperCretaceous) = {l_tmp1, lines_UpperCretaceous[], -l_tmp2, -lines_Tertiary[]};
Plane Surface(s_UpperCretaceous) = {s_UpperCretaceous};

// 3. lower Cretaceous layer
l_tmp1=newl; Line(l_tmp1)={pts_UpperCretaceous[0], 1};
l_tmp2=newl; Line(l_tmp2)={pts_UpperCretaceous[#pts_UpperCretaceous[]-1], 2};
s_LowerCretaceous=news; Curve Loop(s_LowerCretaceous) = {l_tmp1, 1, -l_tmp2, -lines_UpperCretaceous[]};
Plane Surface(s_LowerCretaceous) = {s_LowerCretaceous};

// 4. sills
Include "sill1.geo";
Include "sill2.geo";
Include "sill3.geo";
Include "sill4.geo";
Include "sill5.geo";

// refinement of sills
Field[1] = Attractor;
// Field[1].NNodesByEdge = 10;
Field[1].EdgesList = {lines_sill1[], lines_sill2[], lines_sill3[], lines_sill4[], lines_sill5[]};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 20; // minimum cell size is 2m
Field[2].LcMax = lc;
Field[2].DistMin = 100;
Field[2].DistMax = 300;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

Extrude {0, 0, zmax} {
Surface{s_Tertiary, s_UpperCretaceous, s_LowerCretaceous};
Layers{1};
Recombine;
}
front={847, 1129, 1125};
back={s_Tertiary, s_UpperCretaceous, s_LowerCretaceous};
sidewalls={566, 1126, 848, 845, 1128, 1124};
bottom={1127};
seafloor={846};
Physical Volume("layer1") = {1};
Physical Volume("layer2") = {2};
Physical Volume("layer3") = {3};
Physical Surface("front") = {front[]};
Physical Surface("back") = {back[]};
Physical Surface("sidewalls") = {sidewalls[]};
Physical Surface("bottom") = {bottom[]};
Physical Surface("seafloor") = {seafloor[]};

Color Grey{Surface{front[],back[]};}
Color Pink{Surface{bottom[]};}
Color Green{Surface{sidewalls[]};}
Color Blue{Surface{seafloor[]};}

Color Yellow{Volume{1};}
Color Purple{Volume{2};}
Color Green{Volume{3};}