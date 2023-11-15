clear()
distances = 10000*[0.8,2,4,7,10];
min_deg = 45;
max_deg = 90;
n=20;
f = 5;

for i = 1:1
    for j = 1:4
        build1(distances(j),distances(j+1), min_deg, max_deg,n,f)
%         min_distance = max_distance;
%         max_distance = max_distance + 100000;
        disp(distances(j))
        disp(distances(j+1))
        disp(min_deg)
        disp(max_deg)

    end
    min_deg = max_deg;
    max_deg = max_deg + 45;
    min_distance = 70000;
    max_distance = 100000;
end