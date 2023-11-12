clear()
date = "1102";
ID = 3;
fileName = "../../output/"+date + "/information.csv";
M = readmatrix(fileName);

row = M(ID,:);
correct_ellipse = row(45:50);
correct_terminator = row(51:56);
estimate_ellipse = row(66:71);

A = correct_ellipse(1);
B = correct_ellipse(2);
C = correct_ellipse(3);
D = correct_ellipse(4);
F = correct_ellipse(5);
G = correct_ellipse(6);

A_term = correct_terminator(1);
B_term = correct_terminator(2);
C_term = correct_terminator(3);
D_term = correct_terminator(4);
F_term = correct_terminator(5);
G_term = correct_terminator(6);

A_estimate = estimate_ellipse(1);
B_estimate = estimate_ellipse(2);
C_estimate = estimate_ellipse(3);
D_estimate = estimate_ellipse(4);
F_estimate = estimate_ellipse(5);
G_estimate = estimate_ellipse(6);

fnc = @(X,Y) A * (X).^2 + B * (X) .* (Y) + C * (Y).^2 + D * X + F * Y + G;
fnc_term = @(X,Y) A_term * (X).^2 + B_term * (X) .* (Y) + C_term * (Y).^2 + D_term * X + F_term * Y + G_term;
fnc_est = @(X,Y) A_estimate * (X).^2 + B_estimate * (X) .* (Y) + C_estimate * (Y).^2 + D_estimate * X + F_estimate * Y + G_estimate;

figure
fimplicit(fnc)
hold on
fimplicit(fnc_term)
hold on
fimplicit(fnc_est)
xlim([0, 659])
ylim([0, 494])
grid on
daspect([1 1 1])

ax = gca;
% ax.XDir = 'reverse';
ax.YDir = 'reverse';

% ellipseFile = sprintf("../../output//ellipseImages/%d.png",id);
% saveas(figure,ellipseFile)
% print(ellipseFile,'-dpng')