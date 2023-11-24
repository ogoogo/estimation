function coordinate = transform_coordinates(p1, p2, p3)
    % p1: 原点の座標
    % p2: x軸方向の座標
    % p3: 平面上の座標を求める対象の点
    
    % x軸の方向ベクトルを計算
    v1 = (p2 - p1) / norm(p2 - p1);
    v2 = p3 - p1;
    v3 = dot(v2,v1)*v1;
    
    coordinate = [dot(v2,v1),norm(v2-v3)];
%     disp(coordinate)

end





