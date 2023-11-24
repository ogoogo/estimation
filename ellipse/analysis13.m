% 0.034に収まってる図を作る
clear()
date = "1121";
celestial = "moon";



f = 50;
phi = atan(0.494*7.4/2/f);
k = 1000000*tan(phi);

fileName = "../output/"+date + "/information.csv";
M = readmatrix(fileName);

nonZeroRows = all(M ~= 0, 2);
M = M(nonZeroRows,:);


[row,col] = size(M);

data = createData2(M);
data = filterData(data);
% 
% figure(1)
% scatter(data(13,:),data(14,:),50,data(16,:),"filled")
% colorbar;
% hold on 
% scatter(0,0,200,"magenta","filled")
% scatter(384400,0,200,"magenta","filled")
% daspect([1 1 1])

figure(1)
scatter(data(13,:),data(17,:),10,"filled")
title("angular error - x")
yline(0.017)
xlabel("x in earth-moon fixed frame")
ylabel("angular error")

figure(2)
scatter(data(14,:),data(17,:),10,"filled")
title("angular error - y")
yline(0.017)
xlabel("y in earth-moon fixed frame")
ylabel("angular error")

saveas(1, '../output/graph2/0.034_earth/limit_x.png')
saveas(2, '../output/graph2/0.034_earth/limit_y.png')


