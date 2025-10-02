SetFactory("OpenCASCADE");
xmin=0;
xmax=50000;
ymin=-9000;
ymax=0;
zmin=-1000;
zmax=2000;
lc=20;
width_sill=100;
// 
Point(1) = {xmin,ymin,zmin,lc};
Point(2) = {xmax,ymin,zmin,lc};
Point(3) = {xmax,ymax,zmin,lc};
Point(4) = {xmin,ymax,zmin,lc};
Line(1)={1,2};
Line(2)={3,4};

// 4. sills
Include "sill1.geo";
Include "sill2.geo";
Include "sill3.geo";
Include "sill4.geo";
Include "sill5.geo";
curve_sill1=newl; Spline(curve_sill1) = {pts_sill1[]};
curve_sill2=newl; Spline(curve_sill2) = {pts_sill2[]};
curve_sill3=newl; Spline(curve_sill3) = {pts_sill3[]};
curve_sill4=newl; Spline(curve_sill4) = {pts_sill4[]};
curve_sill5=newl; Spline(curve_sill5) = {pts_sill5[]};
Delete { Line{lines_sill1[]}; Line{lines_sill2[]}; Line{lines_sill3[]};Line{lines_sill4[]};Line{lines_sill5[]}; }
// make surface
curve_sills={curve_sill1, curve_sill2, curve_sill3, curve_sill4,curve_sill5};
Translate {0, -width_sill/2, 0} {Curve{curve_sills[]}; }
Extrude {0, width_sill, 0} {Curve{curve_sills[]}; }

Extrude {0, 0, zmax-zmin} {
Surface{1:5};
Layers{1};
Recombine;
}
// refinement of sills
s()=Unique(Abs(Boundary{ Volume{1:5}; }));
l() = Unique(Abs(Boundary{ Surface{s()}; }));
p[] = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc;

// front={847, 1129, 1125};
// back={s_Tertiary, s_UpperCretaceous, s_LowerCretaceous};
// sidewalls={566, 1126, 848, 845, 1128, 1124};
// bottom={1127};
// seafloor={846};
Physical Volume("sill1") = {1};
Physical Volume("sill2") = {2};
Physical Volume("sill3") = {3};
Physical Volume("sill4") = {4};
Physical Volume("sill5") = {5};
// Physical Surface("front") = {front[]};
// Physical Surface("back") = {back[]};
// Physical Surface("sidewalls") = {sidewalls[]};
// Physical Surface("bottom") = {bottom[]};
// Physical Surface("seafloor") = {seafloor[]};

// Color Grey{Surface{front[],back[]};}
// Color Pink{Surface{bottom[]};}
// Color Green{Surface{sidewalls[]};}
// Color Blue{Surface{seafloor[]};}

// Color Yellow{Volume{1};}
// Color Purple{Volume{2};}
// Color Green{Volume{3};}//+

//+

