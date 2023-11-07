function build1(min_distance,max_distance, min_deg, max_deg, n)

% n = 5;
% min_distance = 70000;
% max_distance = 200000;
% min_deg = 0;
% max_deg = 90;

situation = "mac";
celestial = "moon";

output_name = "../../output/1107/information.csv";

path_setting(situation)

equ_source_name = create_csv(min_distance,max_distance,min_deg,max_deg, celestial);


if exist(output_name) == 0
    id = 1;
else
    info = readmatrix(output_name);
    [row1,col1] = size(info);
    id = info(row1,1) + 1;
end



for i = id:id+n-1
    calculate_data(equ_source_name,output_name, celestial, i)
    disp(i+"回")
end