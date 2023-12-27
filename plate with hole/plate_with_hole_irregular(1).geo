
a = 1.0;
b = 5.0;
c = 0.2;
d = 0.8;

// Mesh.CharacteristicLengthMin = c;
// Mesh.CharacteristicLengthMax = d;
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

Curve Loop(6) = {1,2,3,4,5};

Plane Surface(1) = {6};

Physical Curve("Γᵗ₁") = {3};
Physical Curve("Γᵗ₂") = {4};
Physical Curve("Γᵗ₃") = {1};
Physical Curve("Γᵍ₁") = {2};
Physical Curve("Γᵍ₂") = {5};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 7;
Mesh.MshFileVersion = 2;
Mesh 2;
