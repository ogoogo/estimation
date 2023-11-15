function data =  filterData(Data)
f = 50;
phi = atan(0.494*7.4/2/f);
phi_deg = rad2deg(phi);
% [row,col] = size(M);
index = Data(2, :) <= phi_deg;

% 条件を満たす列だけを取り出す
data = Data(:, index);

end