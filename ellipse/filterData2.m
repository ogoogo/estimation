function data =  filterData2(Data,limit)
f = 50;
phi = atan(0.494*7.4/2/f);
phi_deg = rad2deg(phi);
% [row,col] = size(M);
index = Data(16, :)<= phi_deg & Data(17, :)<= phi_deg;
% index = Data(16, :)<= 0.2;
% index = Data(2, :)<= phi_deg;


% index = Data(16, :) <= limit;

% 条件を満たす列だけを取り出す
data = Data(:, index);

end