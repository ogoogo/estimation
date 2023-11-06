function calculate_data(equ_source_name,output_name,celestial, id)


rng('shuffle');

equ_source = readmatrix(equ_source_name);
[row2,col2] = size(equ_source);

source_i = randi(row2);
r_equ = equ_source(source_i,2:4);

et = equ_source(source_i,1);

moon = cspice_spkezr('MOON', et, 'ECLIPJ2000','NONE','EARTH');
sun = cspice_spkezr('SUN', et, 'ECLIPJ2000','NONE','EARTH');

l_moon = moon(1:3).';
l_sun = sun(1:3).';

if celestial == "moon"
    l_cele = l_moon;
    R = 1737.4;
else
    l_cele = [0,0,0];
    R = 6378.1;
end

cele_vector = l_cele - r_equ;



% グラムシュミット法
a1 = [1,0,0];
a2 = [0,1,0];
dcm1_y = cele_vector/norm(cele_vector);
b1 = a1 - dot(a1,dcm1_y)*dcm1_y;
dcm1_x = b1/norm(b1);
dcm1_z = cross(dcm1_x,dcm1_y);

dcm1 = [dcm1_x; dcm1_y; dcm1_z];


% ランダムにずらす
l = norm(cele_vector);

% 回転させる
dcm2 = cspice_rotmat(dcm1, 2*pi*randi(100)/100, 2);

%縦方向にずらす
psi = deg2rad(2.09);
dpsi_max = atan((l*tan(psi) - R)/l);
dpsi = -dpsi_max + dpsi_max * randi(200) / 100;

dcm3 = cspice_rotmat(dcm2, dpsi, 1);

% 横方向にずらす
phi = deg2rad(2.79);
dphi_max = atan((l*tan(phi) - R)/l);
dphi = -dphi_max + dphi_max * randi(200) / 100;

dcm4 = cspice_rotmat(dcm3, dphi, 3);

% dlpから見た姿勢
dlp_dcm = cspice_rotmat(dcm4, pi, 2);
dlp_dcm2 = cspice_rotmat(dlp_dcm,-pi/2,1);

% dlpから見た太陽方向
sun_dlp = l_sun*dlp_dcm2';
% disp(sun_dlp)

% 月の座標変換
moon_i = l_moon/norm(l_moon);
z = acos(norm([l_moon(1),l_moon(2)])/norm(l_moon));
if l_moon(3)>0
    moon_dk = [0, 0, norm(l_moon)/sin(z)] - l_moon;
else
    moon_dk = -([0, 0, -norm(l_moon)/sin(z)] - l_moon);
end
moon_k = moon_dk/norm(moon_dk);
moon_j = cross(moon_k,moon_i);

dcm_moon1 = [moon_i;moon_j;moon_k];
% dcm_moon2 = [1,0,0;0,0,1;0,-1,0]*[moon_i;moon_j;moon_k];
% 真の月のdcm
dcm_moon2 = cspice_rotmat(dcm_moon1,pi/2,1);

% 月面図調整
dcm_moon0 = eye(3);
dcm_moon3 = cspice_rotmat(dcm_moon0, -23.4/180*pi, 1);

dcm_moon6 = dcm_moon3*dcm_moon2;

dcm_moon = [dcm_moon6(1,:),dcm_moon6(2,:),dcm_moon6(3,:)];
writematrix(dcm_moon,".moondcm.txt", 'Delimiter',',')

if celestial == "moon"
    cele_dlp = (l_moon - r_equ)*dlp_dcm2';
else
    cele_dlp = (- r_equ)*dlp_dcm2';
end
% disp(cele_dlp')

% 楕円係数を求める
K = diag([50*494/3.66, 50*494/3.66, 1]);

Tcp = dcm_moon1/dlp_dcm2;
Tcp = cspice_rotmat(Tcp, pi/2, 1);

% A_mat = Tcp'*(1/(R^2))*eye(3)*Tcp;
A_mat = (1/(R^2))*eye(3);
disp("A")
disp(A_mat)
r_vec = cele_dlp';
M_mat =K\((A_mat*r_vec*(r_vec')*A_mat - (r_vec'*A_mat*r_vec - 1)*A_mat))/K;
disp(M_mat);

% ターミネータ楕円係数を求める

% t_c = cele_dlp';

% グラムシュミット法

n = ((l_sun - l_cele)/norm(l_sun - l_cele))';
c1 = a1' - dot(a1',n)*n;
p_tilda = c1/norm(c1);
s_tilda = cross(n,p_tilda);



Q4 = inv(R^2*eye(3));


omega = (1/2) * atan((2 * p_tilda'* Q4 * s_tilda)/(p_tilda' * Q4 * p_tilda - s_tilda'* Q4 * s_tilda));
p = cos(omega)*p_tilda + sin(omega)*s_tilda;
s = -sin(omega)*p_tilda + cos(omega)*s_tilda;
% disp("confi")
% disp(cross(p,s))
% disp(n)



M_term = sqrt(1/(p'* Q4 *p));
m_term = sqrt(1/(s'* Q4 *s));

R_tw = [p,s,n];
% disp(R_tw')
% disp(inv(R_tw))

T = diag([1/M_term^2, 1/m_term^2, -1]);

% dlp_dcm2 = [1 0 0; 0 0 -1; 0 1 0]*dlp_dcm;

% dlp_dcm3 = cspice_rotmat(dlp_dcm2, pi, 2);

R_ct = dlp_dcm2*R_tw;
t_c = dlp_dcm2*cele_vector';


% t_t = R_tw * t_w;
% R_wt = inv(R_tw);
H = K * [R_ct(:,1), R_ct(:,2), t_c];
T_dash = inv(H')* T * inv(H);
% T_tilda_star = R_wt * H_tilda * T_star * H_tilda' * R_wt';
% D_star = R * T_tilda_star * R';


y0 = 247;
x0 = 329.5;

% A(x-x0)^2 + B(x-x0)(y-y0) + C(y-y0)^2 + D(x-x0) + F(y-y0) + G = 0
A = M_mat(1,1);
B = 2*M_mat(1,2);
C = M_mat(2,2);
D = 2*M_mat(1,3);
F = 2*M_mat(2,3);
G = M_mat(3,3);

G = A*x0^2 + B*x0*y0 + C*y0^2 -D*x0 -F*y0 +G;
F = -B*x0 - 2*C*y0 + F;
D = -2*A*x0 -B*y0 + D;

coefficient = [A,B,C,D,F,G];

A_term = T_dash(1,1);
B_term = 2*T_dash(1,2);
C_term = T_dash(2,2);
D_term = 2*T_dash(1,3);
F_term = 2*T_dash(2,3);
G_term = T_dash(3,3);

G_term = A_term*x0^2 + B_term*x0*y0 + C_term*y0^2 -D_term*x0 -F_term*y0 +G_term;
F_term = -B_term*x0 - 2*C_term*y0 + F_term;
D_term= -2*A_term*x0 -B_term*y0 + D_term;
coefficient_term = [A_term,B_term,C_term,D_term,F_term,G_term];


inforow = [id, et, r_equ, l_moon, l_sun, dcm4(1,:), dcm4(2,:), dcm4(3,:), dlp_dcm(1,:), dlp_dcm(2,:), dlp_dcm(3,:), sun_dlp, dcm_moon, cele_dlp, coefficient, coefficient_term, dlp_dcm2(1,:), dlp_dcm2(2,:), dlp_dcm2(3,:),];

writematrix(inforow, output_name, 'WriteMode', 'append')



end







