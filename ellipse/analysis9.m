% サンセンサのずれ比較
clear()

date = "1115";

fileName = "../output/"+date + "/information.csv";
M = readmatrix(fileName);

M1 = M(1:41,:);
M2 = M(42:82,:);
M3 = M(83:123,:);

m1 = createData(M1);
m2 = createData(M2);
m3 = createData(M3);

figure(1)
plot(m1(13,:),m1(2,:))
hold on
plot(m2(13,:),m2(2,:))
hold on
plot(m3(13,:),m3(2,:))
legend("45°","90°","135°")