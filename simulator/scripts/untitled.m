clear()
dcm = eye(3);

dcm2 = cspice_rotmat(dcm, pi/6, 3);
dcm1 = cspice_rotmat(dcm2,pi/2,1);
dcm3 = dcm1*dcm2;
