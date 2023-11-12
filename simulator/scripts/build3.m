clear()
distances = 10000*[7,10,20,30,40,50,60,70,80,90,100];
min_deg = 0;
max_deg = 45;
n=50;

for i = 1:4
    for j = 1:10
        build1(distances(j),distances(j+1), min_deg, max_deg,n)
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