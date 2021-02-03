SetFactory("OpenCASCADE");
lc=200;                     // mesh size control region
lc_df=25;                   // mesh size control sill
z1=-250;                    // coordinates controlling pseudo-thickness of 2D slice 
z2=250;
xmin=0;                     // x left
xmax=5000;                  // x right
ymin=-5000;                 // y bottom
ymax=0;                     // y top

s_thick=200;
s_length=4000;
s_depth=4000;
s_tilt=10*Pi/180;

// modeling domain
s_box = news;   // next availabe surface tag
Rectangle(s_box) = {xmin, ymin, z1, (xmax-xmin), (ymax-ymin), 0};

// detachment fault
s_df = news;   // next availabe surface tag
Rectangle(s_df) = {xmin-200, ymax-s_depth, z1, s_length, s_thick, 0};
Rotate {{0,0,1}, {xmin, ymax-s_depth, z1}, s_tilt} { Surface{s_df}; }

// intersection of surfaces
surfaces_df_merge[]=BooleanIntersection{ Surface{s_box};}{ Surface{s_df}; Delete;};
Printf("DF: %g", surfaces_df_merge[0]);
s_df = surfaces_df_merge[0];
s_crust=s_box;
// make intersection unique
s() = BooleanFragments{ Surface{s_df}; Delete; }{ Surface{s_box};Delete; };//+
s_df = s[0];
s_box = s[1];

// 2d to 3D
ss[]=Extrude {0, 0, z2-z1} {
    Surface{s_box,s_df};
    Layers{1}; //set layer number to 1 for 2D model
    Recombine;
};

v_box=1;
v_df=2;

Physical Volume("crust") = {v_box};
Physical Volume("intrusion") = {v_df};
s_back={s_df,s_box};
s_front={15,17};
s_left={11,16,7};
s_right=13;
s_bottom=14;
s_seafloor={12};

s()=Unique(Abs(Boundary{ Volume{v_box}; }));
l() = Unique(Abs(Boundary{ Surface{s()}; }));
p[] = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc;

s()=Unique(Abs(Boundary{ Volume{v_df}; }));
l() = Unique(Abs(Boundary{ Surface{s()}; }));
p[] = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc_df;


Physical Surface("front") = {s_front[]};
Physical Surface("back") = {s_back[]};
Physical Surface("sidewalls") = {s_left[], s_right};
Physical Surface("seafloor") = {s_seafloor[]};
Physical Surface("bottom") = {s_bottom[]};

Color Grey{Surface{s_back[],s_front[]};}
Color Pink{Surface{s_bottom[]};}
Color Green{Surface{s_left[],s_right};}
Color Blue{Surface{s_seafloor[]};}

// Color Yellow{Volume{v_box[]};}
// Color Purple{Volume{v_df};}//+
// Hide "*";
// //+
// Show {
// Point{5,6,7,8,13,14,15,16};
// Curve{5,7,17,19,20,21,23,24};
// Surface{11,13};
// }

