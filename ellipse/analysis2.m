% エラーはデータがないやつのみ
clear()
date = "1114";
celestial = "moon";

fileName = "../output/"+date + "/information_fitz.csv";
M1 = readmatrix(fileName);

fileName = "../output/"+date + "/information_3.csv";
M2 = readmatrix(fileName);

fileName = "../output/"+date + "/information_4.csv";
M3 = readmatrix(fileName);

[row,col] = size(M1);

nondata1 = createNonData(M1);
nondata2 = createNonData(M2);
nondata3 = createNonData(M3);

rows = M1(:, 72).*M2(:, 72).*M3(:, 72) ~= 0;
M1 = M1(rows,:);
M2 = M2(rows,:);
M3 = M3(rows,:);

data1 = createData(M1);
data2 = createData(M2);
data3 = createData(M3);





figure(1)
scatter3(data1(1,:), data1(3,:),data1(2,:), 10, 'blue','filled')
hold on
scatter3(data2(1,:), data2(3,:), data2(2,:), 10, 'red','filled')
hold on
scatter3(data3(1,:),data3(3,:), data3(2,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
xlim([0,180])
ylim([0,1000000])
zlim([0,max([max(data1(2,:)),max(data2(2,:)),max(data3(2,:))])])
% zlim([0,20])
title("sun-camera angular - true distance - error")
xlabel("sun-camera angular(deg)")
ylabel("true distance")
zlabel("error(deg)")

figure(2)
scatter(data1(1,:),data1(2,:), 10, 'blue','filled')
hold on
scatter(data2(1,:), data2(2,:), 10, 'red','filled')
hold on
scatter(data3(1,:), data3(2,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
xlim([0,180])
% ylim([0,1000000])
ylim([0,max([max(data1(2,:)),max(data2(2,:)),max(data3(2,:))])])
% ylim([0,20])
title("sun-camera angular - error")
xlabel("sun-camera angular(deg)")
ylabel("error(deg)")

figure(3)
scatter(data1(3,:),data1(2,:), 10, 'blue','filled')
hold on
scatter(data2(3,:), data2(2,:), 10, 'red','filled')
hold on
scatter(data3(3,:), data3(2,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
% yline(0)
xlim([0,1000000])
ylim([0,max([max(data1(2,:)),max(data2(2,:)),max(data3(2,:))])])
% ylim([0,20])
title("true-distance - error")
xlabel("true distance(km)")
ylabel("error(deg)")

figure(4)
scatter(data1(3,:),data1(4,:), 10, 'blue','filled')
hold on
scatter(data2(3,:), data2(4,:), 10, 'red','filled')
hold on
scatter(data3(3,:), data3(4,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(4,:))),max(abs(data2(4,:))),max(abs(data3(4,:)))]);
ylim([-limit, limit])
title("x Error")
xlabel("true distance (km)")
ylabel("x error (km)")

figure(5)
scatter(data1(3,:),data1(5,:), 10, 'blue','filled')
hold on
scatter(data2(3,:), data2(5,:), 10, 'red','filled')
hold on
scatter(data3(3,:), data3(5,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(5,:))),max(abs(data2(5,:))),max(abs(data3(5,:)))]);
ylim([-limit, limit])
title("y Error")
xlabel("true distance (km)")
ylabel("y error (km)")

figure(6)
scatter(data1(3,:),data1(6,:), 10, 'blue','filled')
hold on
scatter(data2(3,:), data2(6,:), 10, 'red','filled')
hold on
scatter(data3(3,:), data3(6,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(6,:))),max(abs(data2(6,:))),max(abs(data3(6,:)))]);
ylim([-limit, limit])
title("z Error")
xlabel("true distance (km)")
ylabel("z error (km)")



saveas(1, '../output/graph/normal_all_existed/sun-cam_true-z_error.png')
saveas(2, '../output/graph/normal_all_existed/sun-cam_error.png')
saveas(3, '../output/graph/normal_all_existed/true-z_error.png')
saveas(4, '../output/graph/normal_all_existed/true-z_x.png')
saveas(5, '../output/graph/normal_all_existed/true-z_y.png')
saveas(6, '../output/graph/normal_all_existed/true-z_z.png')




% 
% function data =  makeData(date,filename,celestial)
% fileName = "../output/"+date + filename;
% M = readmatrix(fileName);
% [row,col] = size(M);
% rows = M(:, 72) ~= 0;
% 
% m = M(rows,:);
% [rows2,col] = size(m);
% 
% 
% data = zeros(2,rows2);
% for i = 1:rows2
%     r_correct = m(i,42:44);
%     r_estimate = m(i,72:74);
%     radError = acos(dot(r_correct,r_estimate)/(norm(r_correct)*norm(r_estimate)));
% %     disp(radError)
%     data(2,i) = rad2deg(radError);
% 
%     if celestial == "moon"
%         l_cele = m(i,6:8);
%     else
%         l_cele = zeros(1,3);
%     end
%     
%     l_sun = m(i,9:11);
%     l_equ = m(i,3:5);
% 
%     rad_sun_cele = acos(dot(l_equ-l_cele, l_sun-l_cele)/(norm(l_equ-l_cele)*norm(l_sun-l_cele)));
%     data(1,i) = rad2deg(rad_sun_cele);
%     data(3,i) = m(i,44);
%     data(4,i) = m(i,72)-m(i,42);
%     data(5,i) = m(i,73)-m(i,43);
%     data(6,i) = m(i,74)-m(i,44);
% 
% 
% 
% end
% 
% 
% end