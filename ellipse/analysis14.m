% 川端先生の図を作る
clear()
date = "1121";
celestial = "moon";
limit = 0.2;
limit2 = "normal";


f = 50;
phi = atan(0.494*7.4/2/f);
k = 1000000*tan(phi);

fileName = "../output/"+date + "/information.csv";
M = readmatrix(fileName);

nonZeroRows = all(M ~= 0, 2);
M = M(nonZeroRows,:);


[row,col] = size(M);

data = createData2(M);
data = filterData2(data,limit);

figure(1)
scatter(data(13,:),data(14,:),30,data(16,:),"filled")
colorbar;
hold on 
scatter(0,0,200,"magenta","filled")
scatter(384400,0,200,"magenta","filled")
daspect([1 1 1])
xlabel("x(km)")
ylabel("y(km)")
title("angular error in earth-moon fixed frame")
saveas(1, '../output/graph2/kawabata_deg/'+limit2+'.png')


