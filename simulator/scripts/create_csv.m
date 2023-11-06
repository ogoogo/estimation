function equ_source_name = create_csv(min_distance, max_distance, min_deg, max_deg, celestial)

equ_source_name = sprintf("./orbit_equ_output/orbit_%d_%d_%d_%d.csv", min_distance, max_distance/10000, min_deg, max_deg);


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
end