SetFactory("OpenCASCADE");
lc=50;                     // mesh size control region
lc_df=25;                   // mesh size control intrusion
z1=-250;                    // coordinates controlling pseudo-thickness of 2D slice 
z2=250;
xmin=0;                     // x left
xmax=4000;                  // x right
ymin=-6000;                 // y bottom
ymax=-3000;                     // y top
e_center_x=(xmax-xmin)/2;
e_center_y=ymax-2000;
e_center_z=z1;
e_len_x=750;
e_len_y=150;
// modeling domain
s_box = news;   // next availabe surface tag
Rectangle(s_box) = {xmin, ymin, z1, (xmax-xmin), (ymax-ymin), 0};

// intrusion
s_in = news;   // next availabe surface tag
Disk(s_in) = {e_center_x, e_center_y, e_center_z, e_len_x, e_len_y};

// intersection of surfaces
surfaces_merge[]=BooleanIntersection{ Surface{s_box};}{ Surface{s_in}; Delete;};
s_in = surfaces_merge[0];

// make intersection unique
s() = BooleanFragments{ Surface{s_in}; Delete; }{ Surface{s_box};Delete; };//+
s_in = s[0];
s_box = s[1];

// 2d to 3D
ss[]=Extrude {0, 0, z2-z1} {
    Surface{s_box,s_in};
    Layers{1}; //set layer number to 1 for 2D model
    Recombine;
};

v_box=1;
v_in=2;

Physical Volume("crust") = {v_box};
Physical Volume("intrusion") = {v_in};
s_back={s_in,s_box};
s_front={13, 12};
s_left=8;
s_right=9;
s_bottom=7;
s_seafloor={10};

s()=Unique(Abs(Boundary{ Volume{v_box}; }));
l() = Unique(Abs(Boundary{ Surface{s()}; }));
p[] = Unique(Abs(Boundary{ Line{l()}; }));
Characteristic Length{p()} = lc;

s()=Unique(Abs(Boundary{ Volume{v_in}; }));
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
Color Purple{Volume{v_in};}