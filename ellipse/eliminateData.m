function data =  eliminateData(Data)
f = 50;
phi = atan(0.494*7.4/2/f);
phi_deg = rad2deg(phi);
% [row,col] = size(M);
index = Data(2, :) > phi_deg |  Data(3, :) + Data(6,:) == 0;

% 条件を満たす列だけを取り出す
data(1,:) = Data(1, index);
data(2,:) = Data(3,index);
data(3,:) = Data(7,index);

end