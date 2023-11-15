clear()
date = "1114";
celestial = "moon";

fileName = "../output/"+date + "/information_fitz.csv";
M1 = readmatrix(fileName);

fileName = "../output/"+date + "/information_3.csv";
M2 = readmatrix(fileName);

fileName = "../output/"+date + "/information_4.csv";
M3 = readmatrix(fileName);

rows = M1(:, 72).*M2(:, 72).*M3(:, 72) ~= 0;
M1 = M1(rows,:);
M2 = M2(rows,:);
M3 = M3(rows,:);

data1 = createData(M1);
data2 = createData(M2);
data3 = createData(M3);

data1 = filterData(data1);
data2 = filterData(data2);
data3 = filterData(data3);






% figure(1)
% scatter3(data1(1,:), data1(8,:),data1(12,:), 10, 'blue','filled')
% hold on
% scatter3(data2(1,:), data2(8,:), data2(12,:), 10, 'red','filled')
% hold on
% scatter3(data3(1,:),data3(8,:), data3(12,:), 10, 'green','filled')
% legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
% grid on;
% % yline(0)
% xlim([0,180])
% ylim([0,1000000])
% zlim([0,max([max(data1(12,:)),max(data2(12,:)),max(data3(12,:))])])
% % zlim([0,20])
% title("sun-camera angular - true distance - error")
% xlabel("sun-camera angular(deg)")
% ylabel("true distance")
% zlabel("error(deg)")
figure(1)
scatter(data1(3,:),data1(12,:), 10, 'blue','filled')
hold on
scatter(data2(3,:), data2(12,:), 10, 'red','filled')
hold on
scatter(data3(3,:), data3(12,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0.034)
xlim([0,1000000])
ylim([0,max([max(data1(12,:)),max(data2(12,:)),max(data3(12,:))])])
ylim([0,0.2])
title("true-distance(camera to moon) - error")
xlabel("true distance(camera to moon)(km)")
ylabel("error(deg)")

figure(2)
scatter(data1(1,:),data1(12,:), 10, 'blue','filled')
hold on
scatter(data2(1,:), data2(12,:), 10, 'red','filled')
hold on
scatter(data3(1,:), data3(12,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0.034)
xlim([0,180])
% ylim([0,1000000])
ylim([0,max([max(data1(12,:)),max(data2(12,:)),max(data3(12,:))])])
ylim([0,0.2])
title("sun-camera angular - error")
xlabel("sun-camera angular(deg)")
ylabel("error(deg)")

figure(3)
scatter(data1(8,:),data1(12,:), 10, 'blue','filled')
hold on
scatter(data2(8,:), data2(12,:), 10, 'red','filled')
hold on
scatter(data3(8,:), data3(12,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0.034)
xlim([0,1000000])
ylim([0,max([max(data1(12,:)),max(data2(12,:)),max(data3(12,:))])])
ylim([0,0.2])
title("true-distance(earth to camera) - error")
xlabel("true distance(earth to camera)(km)")
ylabel("error(deg)")

figure(4)
scatter(data1(8,:),data1(9,:), 10, 'blue','filled')
hold on
scatter(data2(8,:), data2(9,:), 10, 'red','filled')
hold on
scatter(data3(8,:), data3(9,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(9,:))),max(abs(data2(9,:))),max(abs(data3(9,:)))]);
ylim([-limit, limit])
title("x Error")
xlabel("true distance(earth to camera) (km)")
ylabel("x error (km)")

figure(5)
scatter(data1(8,:),data1(10,:), 10, 'blue','filled')
hold on
scatter(data2(8,:), data2(10,:), 10, 'red','filled')
hold on
scatter(data3(8,:), data3(10,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(10,:))),max(abs(data2(10,:))),max(abs(data3(10,:)))]);
ylim([-limit, limit])
title("y Error")
xlabel("true distance(earth to camera) (km)")
ylabel("y error (km)")

figure(6)
scatter(data1(8,:),data1(11,:), 10, 'blue','filled')
hold on
scatter(data2(8,:), data2(11,:), 10, 'red','filled')
hold on
scatter(data3(8,:), data3(11,:), 10, 'green','filled')
legend('ellipse constraint', 'e=0 constraint','e constraint','AutoUpdate','off')
grid on;
yline(0)
limit = max([max(abs(data1(11,:))),max(abs(data2(11,:))),max(abs(data3(11,:)))]);
ylim([-limit, limit])
title("z Error")
xlabel("true distance(earth to camera) (km)")
ylabel("z error (km)")


% 
% saveas(1, '../output/graph/j2000_all_existed/true-z(cam to moon)_error.png')
% saveas(2, '../output/graph/j2000_all_existed/sun-cam_error.png')
% saveas(3, '../output/graph/j2000_all_existed/true-z_error.png')
% saveas(4, '../output/graph/j2000_all_existed/true-z_x.png')
% saveas(5, '../output/graph/j2000_all_existed/true-z_y.png')
% saveas(6, '../output/graph/j2000_all_existed/true-z_z.png')

% 
% saveas(1, '../output/graph/j2000_each_existed/true-z(cam to moon)_error.png')
% saveas(2, '../output/graph/j2000_each_existed/sun-cam_error.png')
% saveas(3, '../output/graph/j2000_each_existed/true-z_error.png')
% saveas(4, '../output/graph/j2000_each_existed/true-z_x.png')
% saveas(5, '../output/graph/j2000_each_existed/true-z_y.png')
% saveas(6, '../output/graph/j2000_each_existed/true-z_z.png')

% saveas(1, '../output/graph/j2000_filtered/true-z(cam to moon)_error.png')
% saveas(2, '../output/graph/j2000_filtered/sun-cam_error.png')
% saveas(3, '../output/graph/j2000_filtered/true-z_error.png')
% saveas(4, '../output/graph/j2000_filtered/true-z_x.png')
% saveas(5, '../output/graph/j2000_filtered/true-z_y.png')
% saveas(6, '../output/graph/j2000_filtered/true-z_z.png')

saveas(1, '../output/graph/j2000_filtered_limit/true-z(cam to moon)_error.png')
saveas(2, '../output/graph/j2000_filtered_limit/sun-cam_error.png')
saveas(3, '../output/graph/j2000_filtered_limit/true-z_error.png')
saveas(4, '../output/graph/j2000_filtered_limit/true-z_x.png')
saveas(5, '../output/graph/j2000_filtered_limit/true-z_y.png')
saveas(6, '../output/graph/j2000_filtered_limit/true-z_z.png')
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