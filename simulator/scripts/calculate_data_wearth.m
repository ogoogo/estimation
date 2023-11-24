function calculate_data_wearth(equ_source_name,output_name,celestial, id,f)


rng('shuffle');

equ_source = readmatrix(equ_source_name);
[row2,col2] = size(equ_source);

source_i = randi(row2);
[et, r_equ, l_moon, l_sun, dcm4_m, dlp_dcm_pov_m, sun_dlp_m, dcm_moon, cele_dlp_m, coefficient_m, coefficient_term_m, dlp_dcm_m] = calculate_data_sub(equ_source(source_i,:),"moon",f);


[et, r_equ, l_moon, l_sun, dcm4_e, dlp_dcm_pov_e, sun_dlp_e, dcm_moon_e, cele_dlp_e, coefficient_e, coefficient_term_e, dlp_dcm_e] = calculate_data_sub(equ_source(source_i,:),"earth",f);

inforow = [id, et, r_equ, l_moon, l_sun, dcm4_m(1,:), dcm4_m(2,:), dcm4_m(3,:), dlp_dcm_pov_m(1,:), dlp_dcm_pov_m(2,:), dlp_dcm_pov_m(3,:), sun_dlp_m, dcm_moon, cele_dlp_m, coefficient_m, coefficient_term_m, dlp_dcm_m(1,:), dlp_dcm_m(2,:), dlp_dcm_m(3,:), dlp_dcm_pov_e(1,:),dlp_dcm_pov_e(2,:),dlp_dcm_pov_e(3,:),dlp_dcm_e(1,:),dlp_dcm_e(2,:),dlp_dcm_e(3,:),sun_dlp_e,cele_dlp_e,coefficient_e];

writematrix(inforow, output_name, 'WriteMode', 'append')



end







