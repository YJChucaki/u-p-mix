
a = 1.0;
b = 5.0;
n = 8;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  b, 0.0, 0.0};
Point(4) = {  b,   b, 0.0};
Point(5) = {0.0,   b, 0.0};
Point(6) = {0.0,   a, 0.0};


Circle(1) = {6,1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};


Curve Loop(1) = {1,2,3,4,5};

Plane Surface(1) = {16};
Plane Surface(2) = {17};
Plane Surface(3) = {18};
Plane Surface(4) = {19};
Plane Surface(5) = {20};

//Plane Surface(1) = {1};
Transfinite Curve{3,4} = 4*n+1;
Transfinite Curve{2,5} = 3*n+1;
Transfinite Curve{1} = n+1;


Physical Curve("Γᵗ₁") = {4};
Physical Curve("Γᵗ₂") = {1};
Physical Curve("Γᵍ₁") = {2};
Physical Curve("Γᵍ₂") = {3};
Physical Curve("Γᵍ₃") = {5};
//Physical Surface("Ω") = {1,2,3,4,5};

//Transfinite Surface{1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
