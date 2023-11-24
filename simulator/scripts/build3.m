clear()
distances = 10000*[15,16];
min_deg = 45;
max_deg = 90;
deg_start = [0, 45,90,135];
deg_end = [45, 90, 135,180];
n=10;
f = 50;

for i = 1:4

    for j = 7:99
        build1(j*10000,(j+1)*10000, deg_start(i), deg_end(i),n,f)
        disp(j*10000)
    end

end