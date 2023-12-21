
R = 1.0;
n = 28;

Point(1) = {0, 0, 0};
Point(2) = {R, R, 0};
Point(3) = {R, 0, 0};

Circle(1) = {1,3,2};

Transfinite Curve{1} = n+1;
Physical Curve("Ω") = {1};
Physical Point("Γᵍ") = {1};
Physical Point("Γᵗ") = {2};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;
