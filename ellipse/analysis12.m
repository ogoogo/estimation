% 計算量比較(zernike)
clear()

date1 = "1118_mac";


fileName = "../output/"+date1 + "/information.csv";
M1 = readmatrix(fileName);

m1 = createData(M1);



% 項目の作成
items = {'b','b with zernike'};

% カテゴリカルデータに変換
categories = categorical(items,items);

figure(1)
% 棒グラフの描画
bar(categories, [0.5126640796661377,28.53331995010376]);

% グラフにタイトルと軸ラベルを追加
title('Computational Performance Comparison');
xlabel('method');
ylabel('time(s)');
% ylim([0,8])


saveas(1, '../output/graph2/zernike/zernike.png')
