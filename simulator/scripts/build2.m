clear()
id = 1;

fileName = sprintf('../../output/1102/information1102.csv');
M = readmatrix(fileName);
[row,col] = size(M);
est_coe = M(id, 66:71);

A = est_coe(1);
B = est_coe(2);
C = est_coe(3);
D = est_coe(4);
F = est_coe(5);
G = est_coe(6);

y0 = 247;
x0 = 329.5;

G = A*x0^2 + B*x0*y0 + C*y0^2 + D*x0 + F*y0 +G;
F = B*x0 + 2*C*y0 + F;
D = 2*A*x0 + B*y0 + D;

