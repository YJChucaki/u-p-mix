
a = 48.0;
b = 44.0;
c = 16.0;
n =33;


Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, b, 0.0};
Point(3) = {  a,  b+c, 0.0};
Point(4) = {0.0,  b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};

Transfinite Curve{-1,3} = 4*n+1 Using Progression 1.00;
Transfinite Curve{4} = 2*n+1;
Transfinite Curve{2} = n+1;

Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Curve("Γ") = {1,3};
//Physical Surface("Ω") = {1};

//Transfinite Surface{1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
Mesh.SecondOrderIncomplete = 1;
SetOrder 2;
//RecombineMesh;
