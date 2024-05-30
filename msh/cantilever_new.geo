
a = 48.0;
b = 12.0;
n = 8;

Point(1) = {0.0, -b/2, 0.0};
Point(2) = {  a, -b/2, 0.0};
Point(3) = {  a,  b/2, 0.0};
Point(4) = {0.0,  b/2, 0.0};
Point(5) = {8.0,  -b/2, 0.0};
Point(6) = {8.0,  b/2, 0.0};


Line(1) = {5,2};
Line(2) = {2,3};
Line(3) = {3,6};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {6,4};
Line(7) = {5,6};

Curve Loop(1) = {1,2,3,-7};
Curve Loop(2) = {5,7,6,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Curve{1,3} = 4*n+1;
Transfinite Curve{2} = n+1;
Transfinite Curve{4,5,6,7} = n/2+1;



Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Curve("Γ") = {1,3,5,6,7};
Physical Surface("Ω") = {2};

//Transfinite Surface{2};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//Mesh.SecondOrderIncomplete = 1;
//SetOrder 2;
//RecombineMesh;
