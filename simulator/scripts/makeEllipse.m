id = 65439;
[shou, data_id] = quorem(sym(id),sym(10000));
csv_id = shou + 1;

csv_name = sprintf("./informations/output/output_%d.csv",csv_id);



M = readmatrix(csv_name);
A = M(data_id, 45);
B = M(data_id, 46);
C = M(data_id, 47);
D = M(data_id, 48);
F = M(data_id, 49);
G = M(data_id, 50);

disp(M(data_id, 1));

fnc = @(X,Y) A * (X/50).^2 + B * (X/50) .* (Y/50) + C * (Y/50).^2 + D * X/50 + F * Y/50 + G;

figure
fimplicit(fnc)
xlim([-2.44, 2.44])
ylim([-1.83, 1.83])
grid on
daspect([1 1 1])
% ellipseFile = sprintf("./ellipseImages/%d.png",id);
% % saveas(figure,ellipseFile)
%  print(ellipseFile,'-dpng')
