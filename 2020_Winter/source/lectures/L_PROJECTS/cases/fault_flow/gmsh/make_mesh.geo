SetFactory("OpenCASCADE");
lc=200;                     // mesh size control region
lc_df=25;                   // mesh size control detachment
z1=-250;                    // coordinates controlling pseudo-thickness of 2D slice 
z2=250;
xmin=0;                     // x left
xmax=6000;                  // x right
ymin=-7000;                 // y bottom
ymax=-3000;                     // y top
width_df = 75;              // fault width
ymin_df  = -6000;
xmin_df  = 2134;            // (xmax-min)/tan(60)/2, fault is in middle of box
theta    = 60*Pi/180;
height_df = 8000;

// modeling domain
s_box = news;   // next availabe surface tag
Rectangle(s_box) = {xmin, ymin, z1, (xmax-xmin), (ymax-ymin), 0};

// detachment fault
s_df = news;   // next availabe surface tag
Rectangle(s_df) = {xmin_df-width_df/2.0, ymin_df, z1, width_df, height_df+1000, 0};
Rotate {{0,0,1}, {xmin_df, ymin_df, z1}, -Pi/2.0+theta} { Surface{s_df}; }

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
Physical Volume("detachment") = {v_df};
s_back={s_df,s_box};
s_front={15, 17};
s_left=7;
s_right=13;
s_bottom=14;
s_seafloor={8, 16, 12};

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
Physical Surface("sidewalls") = {s_left, s_right};
Physical Surface("seafloor") = {s_seafloor[]};
Physical Surface("bottom") = {s_bottom[]};

Color Grey{Surface{s_back[],s_front[]};}
Color Pink{Surface{s_bottom[]};}
Color Green{Surface{s_left,s_right};}
Color Blue{Surface{s_seafloor[]};}

Color Yellow{Volume{v_box[]};}
Color Purple{Volume{v_df};}