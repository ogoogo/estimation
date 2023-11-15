function data =  createNonData(M)

% [row,col] = size(M);
m = M(M(:, 72) == 0, :);
% rows = M(:, 72) == 0.0;
% disp(rows)

% m = M(rows,:);
% disp(m)
[rows2,col] = size(m);


data = zeros(2,rows2);
for i = 1:rows2

    l_cele = m(i,6:8);
    
    l_sun = m(i,9:11);
    l_equ = m(i,3:5);

    rad_sun_cele = acos(dot(l_equ-l_cele, l_sun-l_cele)/(norm(l_equ-l_cele)*norm(l_sun-l_cele)));
%     太陽との角度
    data(1,i) = rad2deg(rad_sun_cele);
%     true z
    data(2,i) = m(i,44);
%     半径
    a = m(i,45);
    b = m(i,46);
    c = m(i,47);
    d = m(i,48);
    f = m(i,49);
    g = m(i,50);

    data(3,i) = sqrt(2*(a*f^2 + c*d^2 - b*d*f + g*(b^2 - 4*a*c))/((b^2 - 4*a*c)*(sqrt((a-c)^2 * b^2) - a - c)));




end


end