clear()
date = "1107";
celestial = "moon";
fileName = "../output/"+date + "/information01.csv";
M = readmatrix(fileName);
% [row,col] = size(M);
rows = M(:, 72) ~= 0;

m = M(rows,:);
[rows2,col] = size(m);

degErrors = zeros(1,rows2);
degSuns = zeros(1,rows2);

for i = 1:rows2
    r_correct = m(i,42:44);
    r_estimate = m(i,72:74);
    radError = acos(dot(r_correct,r_estimate)/(norm(r_correct)*norm(r_estimate)));
%     disp(radError)
    degErrors(i) = rad2deg(radError);

    if celestial == "moon"
        l_cele = m(i,6:8);
    else
        l_cele = zeros(1,3);
    end
    
    l_sun = m(i,9:11);
    l_equ = m(i,3:5);

    rad_sun_cele = acos(dot(l_equ-l_cele, l_sun-l_cele)/(norm(l_equ-l_cele)*norm(l_sun-l_cele)));
    degSuns(i) = rad2deg(rad_sun_cele);

end

fileName = "../output/"+date + "/information.csv";
M2 = readmatrix(fileName);
% [row,col] = size(M);
rows = M2(:, 72) ~= 0;

m2 = M2(rows,:);
[rows2,col] = size(m2);

degErrors2 = zeros(1,rows2);
degSuns2 = zeros(1,rows2);

for i = 1:rows2
    r_correct = m2(i,42:44);
    r_estimate = m2(i, 72:74);
    radError2 = acos(dot(r_correct,r_estimate)/(norm(r_correct)*norm(r_estimate)));
    degErrors2(i) = rad2deg(radError2);

    if celestial == "moon"
        l_cele = m2(i,6:8);
    else
        l_cele = zeros(1,3);
    end
    
    l_sun = m2(i,9:11);
    l_equ = m2(i,3:5);

    rad_sun_cele2 = acos(dot(l_equ-l_cele, l_sun-l_cele)/(norm(l_equ-l_cele)*norm(l_sun-l_cele)));
    degSuns2(i) = rad2deg(rad_sun_cele2);

end



figure()
scatter(degSuns, degErrors, 10, 'blue','filled')
hold on
scatter(degSuns2, degErrors2, 10, 'red','filled')
legend('ellipse constraint', 'e constraint','AutoUpdate','off')
grid on;
yline(0)
% xlim([0,180])
% ylim([-90, 90])
title("sun-camera angular - error")
xlabel("sun-camera angular(deg)")
ylabel("error(deg)")


saveas(1, '../output/1107/graph/normal/sun-cam_error.png')


