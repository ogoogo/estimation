clear()
distances = 10000*[15,16];
min_deg = 45;
max_deg = 90;
deg_start = [44.5, 89.5,134.5];
deg_end = [45.5, 90.5, 135.5];
n=1;
f = 50;

for i = 1:1
    for j = 1:1
        build1(distances(j),distances(j+1), deg_start(i), deg_end(i),n,f)
%         min_distance = max_distance;
%         max_distance = max_distance + 100000;
        disp(distances(j))
        disp(distances(j+1))
        disp(deg_start(1))
        disp(deg_end)

    end
%     min_deg = max_deg;
%     max_deg = max_deg + 45;
%     min_distance = 70000;
%     max_distance = 100000;
end