function data =  createData(M)

% [row,col] = size(M);

rows = M(:, 72) ~= 0;

m = M(rows,:);
[rows2,col] = size(m);

m = M;
[rows2,col] = size(m);


data = zeros(2,rows2);
for i = 1:rows2
    r_correct = m(i,42:44);
    r_estimate = m(i,96:98);
    radError = acos(dot(r_correct,r_estimate)/(norm(r_correct)*norm(r_estimate)));
%     disp(radError)
% 角度誤差
    data(2,i) = rad2deg(radError);

%     if celestial == "moon"
%         
%     else
%         l_cele = zeros(1,3);
%     end
    l_cele = m(i,6:8);
    
    l_sun = m(i,9:11);
    l_equ = m(i,3:5);

    rad_sun_cele = acos(dot(l_equ-l_cele, l_sun-l_cele)/(norm(l_equ-l_cele)*norm(l_sun-l_cele)));
%     太陽との角度
    data(1,i) = rad2deg(rad_sun_cele);
%     true z
    data(3,i) = m(i,44);
%     xの誤差(視線)
    data(4,i) = m(i,72)-m(i,42);
%     yの誤差(視線)
    data(5,i) = m(i,73)-m(i,43);
%     zの誤差(視線)
    data(6,i) = m(i,74)-m(i,44);

    a = m(i,45);
    b = m(i,46);
    c = m(i,47);
    d = m(i,48);
    f = m(i,49);
    g = m(i,50);
%     半径
    data(7,i) = sqrt(2*(a*f^2 + c*d^2 - b*d*f + g*(b^2 - 4*a*c))/((b^2 - 4*a*c)*(sqrt((a-c)^2 * b^2) - a - c)));

%     衛星の位置
    dcm = [m(i,57:59);m(i,60:62);m(i,63:65)];
    r_w = r_estimate*dcm;
    l_equ_e = l_cele - r_w;
    l_equ_c = m(i,3:5);
    error_w = l_equ_c - l_equ_e;
%     disp(l_equ_e)
%     disp(l_equ_c)

%     地球から衛星の距離(true)
    data(8,i) = norm(l_equ_c);
%     xyzの誤差(eclipj2000)
    data(9:11,i) = error_w';

%     角度誤差(eclipj2000)
    radError2 = acos(dot(l_equ_c,l_equ_e)/(norm(l_equ_c)*norm(l_equ_e)));
    data(12,i) = rad2deg(radError2);

%     川端先生図の座標(x,y)
    coordinate = transform_coordinates([0,0,0],l_cele,l_equ_c);
    data(13:14,i) = coordinate';
%     絶対値誤差
    l_equ_e_wearth = m(i,99:101);
    error = l_equ_e_wearth - l_equ_c;
    abserror = norm(error);
    data(15,i) = abserror;
%     角度誤差新(全体)
    radError3 = acos(dot(l_equ_c,l_equ_e_wearth)/(norm(l_equ_c)*norm(l_equ_e_wearth)));
    data(16,i) = rad2deg(radError3);


    r_correct_e = m(i,87:89);
    r_estimate_e = m(i,102:104);
    radError_e = acos(dot(r_correct_e,r_estimate_e)/(norm(r_correct_e)*norm(r_estimate_e)));
%     disp(radError)
% 角度誤差(地球)
    data(17,i) = rad2deg(radError_e);

end


end