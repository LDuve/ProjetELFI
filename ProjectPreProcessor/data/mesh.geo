Lx = 1.0;
Ly = 1.0;
LxBar = 0.5;
LyBar = 2.0;   


h  =  Ly;

lc = 0.1;

Point(1) = { -2.5*Lx, 0., 0., lc};
Point(2) = { -2.5*Lx,  Ly, 0., lc};
Point(3) = {-LxBar/2, Ly, 0., lc};
Point(4) = {-LxBar/2,  Ly+LyBar, 0., lc};
Point(5) = { LxBar/2,  Ly+LyBar  , 0., lc};
Point(6) = { LxBar/2,  Ly  , 0., lc};
Point(7) = {  2.5*Lx, Ly  , 0., lc};
Point(8) = {  2.5*Lx, 0.  , 0., lc};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};


Plane Surface(1) = {1};



