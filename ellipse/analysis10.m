% 計算量比較
clear()

date1 = "1118_mac";
date2 = "1118_ras";

fileName = "../output/"+date1 + "/information.csv";
M1 = readmatrix(fileName);

m1 = createData(M1);

fileName = "../output/"+date2 + "/information.csv";
M2 = readmatrix(fileName);

m2 = createData(M2);

% 項目の作成
items = {'fitz','a','b'};

% カテゴリカルデータに変換
categories = categorical(items,items);

figure(1)
% 棒グラフの描画
bar(categories, m1(14,:));

% グラフにタイトルと軸ラベルを追加
title('Computational Performance Comparison on Mac');
xlabel('method');
ylabel('time(s)');
% ylim([0,8])

figure(2)
% 棒グラフの描画
bar(categories, m2(14,:));

% グラフにタイトルと軸ラベルを追加
title('Computational Performance Comparison on Raspberry Pi 4B');
xlabel('method');
ylabel('time(s)');
ylim([0,8])

saveas(1, '../output/graph/calculation/mac2.png')
saveas(2, '../output/graph/calculation/ras.png')