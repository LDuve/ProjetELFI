Lx = 1.0;
Ly = 1.0;
LxBar = 0.5;
LyBar = 2.0;   
h  = LyBar + Ly;

lc = 0.1;

Point(1) = {-LxBar/2, Ly, 0., lc};
Point(2) = {-LxBar/2,  h, 0., lc};
Point(3) = { LxBar/2,  h, 0., lc};
Point(4) = { LxBar/2, Ly, 0., lc};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Curve Loop(1) = {1:5, -6, 7:8};
Plane Surface(1) = {1};