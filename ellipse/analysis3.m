% 画角の比較
clear()

fileName = "../output/50/information_e0.csv";
M = readmatrix(fileName);
data50_e0 = createData(M);

fileName = "../output/50/information.csv";
M = readmatrix(fileName);
data50_e = createData(M);

fileName = "../output/28/information_e0.csv";
M = readmatrix(fileName);
data28_e0 = createData(M);

fileName = "../output/28/information.csv";
M = readmatrix(fileName);
data28_e = createData(M);

fileName = "../output/10/information_e0.csv";
M = readmatrix(fileName);
data10_e0 = createData(M);

fileName = "../output/10/information.csv";
M = readmatrix(fileName);
data10_e = createData(M);

fileName = "../output/5/information_e0.csv";
M = readmatrix(fileName);
data5_e0 = createData(M);

fileName = "../output/5/information.csv";
M = readmatrix(fileName);
data5_e = createData(M);

figure(1)
scatter(data50_e0(7,:), data50_e0(2,:), 10, 'red','filled')
hold on
scatter(data50_e(7,:), data50_e(2,:), 10, 'green','filled')
legend('e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
% xlim([0,1000000])
ylim([0,max([max(data50_e0(2,:)),max(data50_e(2,:))])])
% ylim([0,20])
title("projected radius - error (f=50mm)")
xlabel("projected radius(px)")
ylabel("error(deg)")


figure(2)
scatter(data28_e0(7,:), data28_e0(2,:), 10, 'red','filled')
hold on
scatter(data28_e(7,:), data28_e(2,:), 10, 'green','filled')
legend('e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
% xlim([0,1000000])
ylim([0,max([max(data28_e0(2,:)),max(data28_e(2,:))])])
% ylim([0,20])
title("projected radius - error (f=28mm)")
xlabel("projected radius(px)")
ylabel("error(deg)")


figure(3)
scatter(data10_e0(7,:), data10_e0(2,:), 10, 'red','filled')
hold on
scatter(data10_e(7,:), data10_e(2,:), 10, 'green','filled')
legend('e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
% xlim([0,1000000])
ylim([0,max([max(data10_e0(2,:)),max(data10_e(2,:))])])
% ylim([0,20])
title("projected radius - error (f=10mm)")
xlabel("projected radius(px)")
ylabel("error(deg)")


figure(4)
scatter(data5_e0(7,:), data5_e0(2,:), 10, 'red','filled')
hold on
scatter(data5_e(7,:), data5_e(2,:), 10, 'green','filled')
legend('e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
% xlim([0,1000000])
ylim([0,max([max(data5_e0(2,:)),max(data5_e(2,:))])])
% ylim([0,20])
title("projected radius - error (f=5mm)")
xlabel("projected radius(px)")
ylabel("error(deg)")

saveas(1, '../output/graph/focus_compare/f=50.png')
saveas(2, '../output/graph/focus_compare/f=28.png')
saveas(3, '../output/graph/focus_compare/f=10.png')
saveas(4, '../output/graph/focus_compare/f=5.png')
