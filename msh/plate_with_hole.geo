
a = 1.0;
b = 5.0;
// n = 8;
// c = 0.099970;
c = 0.45;

Point(1) = {0.0, 0.0, 0.0, c};
Point(2) = {  a, 0.0, 0.0, c};
Point(3) = {  b, 0.0, 0.0, c};
Point(4) = {  b,   b, 0.0, c};
Point(5) = {0.0,   b, 0.0, c};
Point(6) = {0.0,   a, 0.0, c};
Point(7) = {2*a, 0.0, 0.0, c};
Point(8) = {  b, 1.6*a, 0.0, c};
Point(9) = {1.6*a,   b, 0.0, c};
Point(10) = {0.0, 2*a, 0.0, c};
Point(11) = {a*2^0.5,a*2^0.5, 0.0, c};
Point(12) = {a*2^0.5/2,a*2^0.5/2, 0.0, c};

Circle(1) = {12,1,2};
Line(2) = {2,7};
Line(3) = {7,3};
Line(4) = {3,8};
Line(5) = {8,4};
Line(6) = {4,9};
Line(7) = {9,5};
Line(8) = {5,10};
Line(9) = {10,6};
Circle(10) = {6,1,12};
Circle(11) = {10,1,11};
Circle(12) = {11,1,7};
Line(13) = {8,11};
Line(14) = {9,11};
Line(15) = {11,12};

Curve Loop(16) = {15,1,2,-12};
Curve Loop(17) = {12,3,4,13};
Curve Loop(18) = {5,6,14,-13};
Curve Loop(19) = {8,11,-14,7};
Curve Loop(20) = {10,-15,-11,9};
Curve Loop(21) = {1,2,3,4,5};

Plane Surface(1) = {16};
Plane Surface(2) = {17};
Plane Surface(3) = {18};
Plane Surface(4) = {19};
Plane Surface(5) = {20};


// Transfinite Curve{1,2,4,7,9,10,11,12,15} = n+1;
// Transfinite Curve{3,5,6,8,13,14} = 2*n+1;


Physical Curve("Γᵍ₃") = {4,5};
Physical Curve("Γᵗ₁") = {6,7};
Physical Curve("Γᵗ₂") = {1,10};
Physical Curve("Γᵍ₁") = {2,3};
Physical Curve("Γᵍ₂") = {8,9};
//Physical Curve("Γ") = {11,12,13,14,15};
Physical Curve("Γ") = {1,2,3,4,5};
Physical Surface("Ω") = {1,2,3,4,5};

// Transfinite Surface{1,2,3,4,5};

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
