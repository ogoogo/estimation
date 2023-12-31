% 45度ずつに分けて視線ベクトルx,y,z方向誤差を出す
clear()
date = "1107";
fileName = "../output/"+date + "/information01.csv";
M = readmatrix(fileName);
% [row,col] = size(M);



M1_ = M(1:100,:);
M2_ = M(101:200,:);
M3_ = M(201:300,:);
M4_ = M(301:400,:);


rows1 = M1_(:, 72) ~= 0;
rows2 = M2_(:, 72) ~= 0;
rows3 = M3_(:, 72) ~= 0;
rows4 = M4_(:, 72) ~= 0;


M1 = M1_(rows1, :);
M2 = M2_(rows2, :);
M3 = M3_(rows3, :);
M4 = M4_(rows4, :);

[row1,col1] = size(M1);
[row2,col2] = size(M2);
[row3,col3] = size(M3);
[row4,col4] = size(M4);

data1 = zeros(4,row1);
data2 = zeros(4,row2);
data3 = zeros(4,row3);
data4 = zeros(4,row4);

data1(1,:) = M1(:,44)';
data1(2,:) = M1(:,72)'- M1(:,42)';
data1(3,:) = M1(:,73)'- M1(:,43)';
data1(4,:) = M1(:,74)'- M1(:,44)';

data2(1,:) = M2(:,44)';
data2(2,:) = M2(:,72)'- M2(:,42)';
data2(3,:) = M2(:,73)'- M2(:,43)';
data2(4,:) = M2(:,74)'- M2(:,44)';

data3(1,:) = M3(:,44)';
data3(2,:) = M3(:,72)'- M3(:,42)';
data3(3,:) = M3(:,73)'- M3(:,43)';
data3(4,:) = M3(:,74)'- M3(:,44)';

data4(1,:) = M4(:,44)';
data4(2,:) = M4(:,72)'- M4(:,42)';
data4(3,:) = M4(:,73)'- M4(:,43)';
data4(4,:) = M4(:,74)'- M4(:,44)';


fileName = "../output/"+date + "/information.csv";
m = readmatrix(fileName);
% [row,col] = size(M);



m1_ = m(1:100,:);
m2_ = m(101:200,:);
m3_ = m(201:300,:);
m4_ = m(301:400,:);


rows1 = m1_(:, 72) ~= 0;
rows2 = m2_(:, 72) ~= 0;
rows3 = m3_(:, 72) ~= 0;
rows4 = m4_(:, 72) ~= 0;


m1 = m1_(rows1, :);
m2 = m2_(rows2, :);
m3 = m3_(rows3, :);
m4 = m4_(rows4, :);

[row1,col1] = size(m1);
[row2,col2] = size(m2);
[row3,col3] = size(m3);
[row4,col4] = size(m4);

Data1 = zeros(4,row1);
Data2 = zeros(4,row2);
Data3 = zeros(4,row3);
Data4 = zeros(4,row4);

Data1(1,:) = m1(:,44)';
Data1(2,:) = m1(:,72)'- m1(:,42)';
Data1(3,:) = m1(:,73)'- m1(:,43)';
Data1(4,:) = m1(:,74)'- m1(:,44)';

Data2(1,:) = m2(:,44)';
Data2(2,:) = m2(:,72)'- m2(:,42)';
Data2(3,:) = m2(:,73)'- m2(:,43)';
Data2(4,:) = m2(:,74)'- m2(:,44)';

Data3(1,:) = m3(:,44)';
Data3(2,:) = m3(:,72)'- m3(:,42)';
Data3(3,:) = m3(:,73)'- m3(:,43)';
Data3(4,:) = m3(:,74)'- m3(:,44)';

Data4(1,:) = m4(:,44)';
Data4(2,:) = m4(:,72)'- m4(:,42)';
Data4(3,:) = m4(:,73)'- m4(:,43)';
Data4(4,:) = m4(:,74)'- m4(:,44)';

figure(1)
scatter(data1(1,:), data1(2,:), 25, 'blue','filled')
hold on
scatter(Data1(1,:), Data1(2,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data1(2,:))), max(abs(data1(2,:)))])
title("x Error (Angle with the Sun 0-45 Deg)")
xlabel("true distance (km)")
ylabel("x error (km)")

figure(2)
scatter(data1(1,:), data1(3,:), 25, 'blue','filled')
hold on
scatter(Data1(1,:), Data1(3,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data1(3,:))), max(abs(data1(3,:)))])
title("y Error (Angle with the Sun 0-45 Deg)")
xlabel("true distance (km)")
ylabel("y error (km)")

figure(3)
scatter(data1(1,:), data1(4,:), 25, 'blue','filled')
hold on
scatter(Data1(1,:), Data1(4,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data1(4,:))), max(abs(data1(4,:)))])
title("z Error (Angle with the Sun 0-45 Deg)")
xlabel("true distance (km)")
ylabel("z error (km)")


figure(4)
scatter(data2(1,:), data2(2,:), 25, 'blue','filled')
hold on
scatter(Data2(1,:), Data2(2,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data2(2,:))), max(abs(data2(2,:)))])
title("x Error (Angle with the Sun 45-90 Deg)")
xlabel("true distance (km)")
ylabel("x error (km)")

figure(5)
scatter(data2(1,:), data2(3,:), 25, 'blue','filled')
hold on
scatter(Data2(1,:), Data2(3,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data2(3,:))), max(abs(data2(3,:)))])
title("y Error (Angle with the Sun 45-90 Deg)")
xlabel("true distance (km)")
ylabel("y error (km)")

figure(6)
scatter(data2(1,:), data2(4,:), 25, 'blue','filled')
hold on
scatter(Data2(1,:), Data2(4,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data2(4,:))), max(abs(data2(4,:)))])
title("z Error (Angle with the Sun 45-90 Deg)")
xlabel("true distance (km)")
ylabel("z error (km)")


figure(7)
scatter(data3(1,:), data3(2,:), 25, 'blue','filled')
hold on
scatter(Data3(1,:), Data3(2,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data3(2,:))), max(abs(data3(2,:)))])
title("x Error (Angle with the Sun 90-135 Deg)")
xlabel("true distance (km)")
ylabel("x error (km)")

figure(8)
scatter(data3(1,:), data3(3,:), 25, 'blue','filled')
hold on
scatter(Data3(1,:), Data3(3,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data3(3,:))), max(abs(data3(3,:)))])
title("y Error (Angle with the Sun 90-135 Deg)")
xlabel("true distance (km)")
ylabel("y error (km)")

figure(9)
scatter(data3(1,:), data3(4,:), 25, 'blue','filled')
hold on
scatter(Data3(1,:), Data3(4,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data3(4,:))), max(abs(data3(4,:)))])
title("z Error (Angle with the Sun 90-135 Deg)")
xlabel("true distance (km)")
ylabel("z error (km)")


figure(10)
scatter(data4(1,:), data4(2,:), 25, 'blue','filled')
hold on
scatter(Data4(1,:), Data4(2,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data4(2,:))), max(abs(data4(2,:)))])
title("x Error (Angle with the Sun 135-180 Deg)")
xlabel("true distance (km)")
ylabel("x error (km)")

figure(11)
scatter(data4(1,:), data4(3,:), 25, 'blue','filled')
hold on
scatter(Data4(1,:), Data4(3,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data4(3,:))), max(abs(data4(3,:)))])
title("y Error (Angle with the Sun 135-180 Deg)")
xlabel("true distance (km)")
ylabel("y error (km)")

figure(12)
scatter(data4(1,:), data4(4,:), 25, 'blue','filled')
hold on
scatter(Data4(1,:), Data4(4,:), 25, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
ylim([-max(abs(data4(4,:))), max(abs(data4(4,:)))])
title("z Error (Angle with the Sun 135-180 Deg)")
xlabel("true distance (km)")
ylabel("z error (km)")

saveas(1, '../output/1107/graph/normal/x_0_45.png')
saveas(2, '../output/1107/graph/normal/y_0_45.png')
saveas(3, '../output/1107/graph/normal/z_0_45.png')

saveas(4, '../output/1107/graph/normal/x_45_90.png')
saveas(5, '../output/1107/graph/normal/y_45_90.png')
saveas(6, '../output/1107/graph/normal/z_45_90.png')

saveas(7, '../output/1107/graph/normal/x_90_135.png')
saveas(8, '../output/1107/graph/normal/y_90_135.png')
saveas(9, '../output/1107/graph/normal/z_90_135.png')

saveas(10, '../output/1107/graph/normal/x_135_180.png')
saveas(11, '../output/1107/graph/normal/y_135_180.png')
saveas(12, '../output/1107/graph/normal/z_135_180.png')



