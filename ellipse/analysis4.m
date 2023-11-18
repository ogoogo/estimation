% データないのがエラーのヒストグラム
clear()

date = "1114";

fileName = "../output/"+date + "/information_fitz.csv";
M1 = readmatrix(fileName);

fileName = "../output/"+date + "/information_3.csv";
M2 = readmatrix(fileName);

fileName = "../output/"+date + "/information_4.csv";
M3 = readmatrix(fileName);

[row,col] = size(M1);

nondata1 = createNonData(M1);
nondata2 = createNonData(M2);
nondata3 = createNonData(M3);


figure(1)
h1 = histogram(nondata1(1,:));
h1.BinWidth = 10;
xlim([0,180])
title("sun angular - fitz")
xlabel("sun angular")
ylabel("number of error data(fitz)")

figure(2)
h2 = histogram(nondata2(1,:));
h2.BinWidth = 10;
xlim([0,180])
title("sun angular - f0")
xlabel("sun angular")
ylabel("number of error data(f0)")

figure(3)
h3 = histogram(nondata3(1,:));
h3.BinWidth = 10;
xlim([0,180])
title("sun angular - ef")
xlabel("sun angular")
ylabel("number of error data(ef)")


figure(4)
h1 = histogram(nondata1(2,:));
h1.BinWidth = 50000;
xlim([0,1000000]);
ylim([0,20])
title("true z - fitz")
xlabel("true z")
ylabel("number of error data(fitz)")

figure(5)
h1 = histogram(nondata2(2,:));
h1.BinWidth = 50000;
xlim([0,1000000]);
ylim([0,20])
title("true z - f0")
xlabel("true z")
ylabel("number of error data(f0)")

figure(6)
h1 = histogram(nondata2(2,:));
h1.BinWidth = 50000;
xlim([0,1000000]);
ylim([0,20])
title("true z - ef")
xlabel("true z")
ylabel("number of error data(ef)")

saveas(1, '../output/graph/nonDataHist/sun-angular_fitz.png')
saveas(6, '../output/graph/nonDataHist/sun-angular_f0.pngg')
saveas(6, '../output/graph/nonDataHist/sun-angular_ef.png')
saveas(6, '../output/graph/nonDataHist/true-z_fitz.png')
saveas(6, '../output/graph/nonDataHist/true-z_f0.png')
saveas(6, '../output/graph/nonDataHist/true-z_ef.png')