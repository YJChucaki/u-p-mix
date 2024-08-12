
c = 0.1;
Point(1) = {5.0,0.0,0.0,c};
Point(2) = {11.0,0.0,0.0,c};
Point(3) = {14.0,3.0,0.0,c};
Point(4) = {3/10*6^0.5+46/5,9/10*6^0.5+18/5,0.0,c};
Point(5) = {0.0,5.0,0.0,c};
Point(6) = {5.0,5.0,0.0,c};
Point(7) = {11.0,3.0,0.0,c};

Line(1) = {1,2};
Circle(2) = {2,7,3};
Circle(3) = {3,7,4};
Circle(4) = {4,6,5};
Circle(5) = {5,6,1};

Curve Loop(6) = {1,2,3,4,5};

Plane Surface(1) = {6};
Physical Curve("Γ") = {1,2,3,4,5};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;