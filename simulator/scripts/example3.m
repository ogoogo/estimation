clear()

% 編集引数
min_distance = 70000;
max_distance = 100000;
min_deg = 45;
max_deg = 90;

celestial = "moon";
% celestial = "earth";

f = 50;


equ_source_name = sprintf("./orbit_equ_output/orbit_%d_%d_%d_%d.csv", min_distance/10000, max_distance/10000, min_deg, max_deg);


addpath('~/issl/research/simulator/mice/src/mice/')
addpath('~/issl/research/simulator/mice/lib/')

% addpath('~/issl/research/simulator/mice/src/mice/')
% addpath('~/issl/research/simulator/mice/lib/')

if cspice_ktotal( 'ALL' ) >= 1
else
    disp('kernels settings')
    cspice_furnsh('../kernel/naif0012.tls');
    cspice_furnsh('../kernel/de440.bsp');
end


if exist('./information.csv') == 0
    id = 1;
else
    info = readmatrix('./information.csv');
    [row1,col1] = size(info);
    id = info(row1,1) + 1;
end





rng('shuffle');
% et_min = 7.2247686918289995e+08;
% et_max = 7.5962046918519998e+08;
%
% et = et_min + randi(round(et_max-et_min));
%
% date = cspice_et2utc(et,'C',6);
% disp(date)



% writematrix(l_moon,'l_moon.txt','Delimiter',',')
% writematrix(l_sun,'l_sun.txt','Delimiter',',')






if exist(equ_source_name,"file") == 0
    for j = 1:85
        disp(j)
        fileName = sprintf('./orbit_equ/orbit_equ%d.dat',j);
        M = readmatrix(fileName);
        [row,col] = size(M);
        for i = 1:row

            r_equ = M(i,2:4);
            et = M(i,1);
            cspice_sun = cspice_spkezr('SUN', et, 'ECLIPJ2000','NONE','EARTH');
            r_sun = cspice_sun(1:3)';
            if celestial == "moon"
                cspice_cele = cspice_spkezr('MOON', et, 'ECLIPJ2000','NONE','EARTH');
                r_cele = cspice_cele(1:3)';
            else
                r_cele = [0,0,0];
            end
            norm_cele = norm(r_cele - r_equ);
            rad_cele = acos(dot(r_equ-r_cele,r_sun-r_cele)/(norm(r_equ-r_cele)*norm(r_sun-r_cele)));
            deg_cele = rad2deg(rad_cele);
            if min_distance < norm_cele && norm_cele < max_distance
                if min_deg < deg_cele && deg_cele < max_deg
                    writematrix(M(i,:), equ_source_name, 'WriteMode', 'append')
                end

            end
        end
    end

end

equ_source = readmatrix(equ_source_name);
[row2,col2] = size(equ_source);

source_i = randi(row2);
r_equ = equ_source(source_i,2:4);

et = equ_source(source_i,1);


moon = cspice_spkezr('MOON', et, 'ECLIPJ2000','NONE','EARTH');
sun = cspice_spkezr('SUN', et, 'ECLIPJ2000','NONE','EARTH');

% moon = cspice_spkezr('MOON', et, 'J2000','NONE','EARTH');
% sun = cspice_spkezr('SUN', et, 'J2000','NONE','EARTH');

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


% while 1
%     for i = 1:row
%         if M(i,1) < et && et < M(i+1,1)
%             r_equ = M(i,2:4);
%             break
%         end
%     end
%     cele_vector = l_cele - r_equ;
%     norm_cele = norm(cele_vector);
%     if 0 < norm_cele && norm_cele < 700000
%
%         break
%     end
%
%
% end



% cele_vector = l_cele - r_equ;



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
% psi = deg2rad(2.09);
psi = atan(0.494*7.4/2/f);
dpsi_max = atan((l*tan(psi) - R)/l);
dpsi = -dpsi_max + dpsi_max * randi(200) / 100;

dcm3 = cspice_rotmat(dcm2, dpsi, 1);

% 横方向にずらす
% phi = deg2rad(2.79);
phi = atan(0.659*7.4/2/f);
dphi_max = atan((l*tan(phi) - R)/l);
dphi = -dphi_max + dphi_max * randi(200) / 100;

dcm4 = cspice_rotmat(dcm3, dphi, 3);


% writematrix(dcm4(1,:),'dcm1.txt','Delimiter',',')
% writematrix(dcm4(2,:),'dcm2.txt','Delimiter',',')
% writematrix(dcm4(3,:),'dcm3.txt','Delimiter',',')

% dlpから見た姿勢
dlp_dcm = cspice_rotmat(dcm4, pi, 2);
dlp_dcm2 = cspice_rotmat(dlp_dcm,-pi/2,1);

% writematrix(dlp_dcm(1,:),'dcm1.txt','Delimiter',',')
% writematrix(dlp_dcm(2,:),'dcm2.txt','Delimiter',',')
% writematrix(dlp_dcm(3,:),'dcm3.txt','Delimiter',',')

% dlpから見た太陽方向
sun_dlp = l_sun*dlp_dcm2';
sun_dlp = sun_dlp/norm(sun_dlp);
disp(sun_dlp)

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
%     t_w = ((l_moon - r_equ)*dlp_dcm')';
%     cele_w = l_moon';
else
    cele_dlp = (- r_equ)*dlp_dcm2';
%     t_w = ((- r_equ)*dlp_dcm')';
%     cele_w = [0;0;0];
end
disp(cele_dlp')

% 楕円係数を求める
% K = diag([f*1000/7.4, f*494/3.66, 1]);
K = diag([f*1000/7.4, f*1000/7.4, 1]);


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

sun_cele = dlp_dcm2*((l_sun - l_cele)/norm(l_sun - l_cele))';
% disp(sun_cele(1:2))
% disp(sun_dlp(1:2)/norm(sun_dlp(1:2)))
n = dlp_dcm2*((l_sun - l_cele)/norm(l_sun - l_cele))';
% n = (sun_dlp/norm(sun_dlp))';
% disp(n1)
% disp(n)
c1 = a1' - dot(a1',n)*n;
p = c1/norm(c1);
s = cross(n,p);



Q4 = inv(R^2*eye(3));


% omega = (1/2) * atan((2 * p_tilda'* Q4 * s_tilda)/(p_tilda' * Q4 * p_tilda - s_tilda'* Q4 * s_tilda));
% p = cos(omega)*p_tilda + sin(omega)*s_tilda;
% s = -sin(omega)*p_tilda + cos(omega)*s_tilda;
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

R_ct = R_tw;
% R_ct = dlp_dcm2*R_tw;

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


% disp([A,B,C,D,F,G])



fnc = @(X,Y) A * (X).^2 + B * (X) .* (Y) + C * (Y).^2 + D * X + F * Y + G;
fnc_term = @(X,Y) A_term * (X).^2 + B_term * (X) .* (Y) + C_term * (Y).^2 + D_term * X + F_term * Y + G_term;

figure
fimplicit(fnc)
hold on
fimplicit(fnc_term)
xlim([0, 659])
ylim([0, 494])
% xlim([-329.5, 329.5])
% ylim([-247, 247])
grid on
daspect([1 1 1])

ax = gca;
% ax.XDir = 'reverse';
ax.YDir = 'reverse';

ellipseFile = sprintf("./ellipseImages/%d.png",id);
% saveas(figure,ellipseFile)
print(ellipseFile,'-dpng')


inforow = [id, et, r_equ, l_moon, l_sun,  dcm4(1,:), dcm4(2,:), dcm4(3,:), dlp_dcm(1,:), dlp_dcm(2,:), dlp_dcm(3,:), sun_dlp, dcm_moon, cele_dlp, coefficient];

writematrix(inforow, "information.csv", 'WriteMode', 'append')
writematrix(inforow,"inforow.txt", 'Delimiter',',')








