clear()
min_distance = 70000;
max_distance = 100000;
min_deg = 0;
max_deg = 45;
n=10;

for i = 1:4
    for j = 1:10
        build1(min_distance,max_distance, min_deg, max_deg,n)
        min_distance = max_distance;
        max_distance = max_distance + 100000;
        disp(min_distance)
        disp(max_distance)
        disp(min_deg)
        disp(max_deg)

    end
    min_deg = max_deg;
    max_deg = max_deg + 45;
    min_distance = 70000;
    max_distance = 100000;
end