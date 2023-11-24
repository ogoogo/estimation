import numpy as np
import sympy

def find_point_with_min_distance(line1, line2):
    p,v = line1
    q,w = line2

    #変数s, tを用意
    s = sympy.Symbol('s')
    t = sympy.Symbol('t')
    #点Pと点Qの距離の２乗
    PQ2 = ( (q[0]+t*w[0]) - (p[0]+s*v[0]) )**2\
        +( (q[1]+t*w[1]) - (p[1]+s*v[1]) )**2\
        +( (q[2]+t*w[2]) - (p[2]+s*v[2]) )**2
    #これを偏微分
    dPQ2_ds = sympy.diff(PQ2, s)
    dPQ2_dt = sympy.diff(PQ2, t)
    ans = sympy.solve([dPQ2_ds, dPQ2_dt])
    # print('ans = {}'.format(ans))
    s, t = ans[s], ans[t] #ここでs, tは通常の変数に戻る
    # print('(s, t) = ({}, {})'.format(s, t))
    #それは直線１上のどこなのか
    x1, y1, z1 = p[0]+s*v[0], p[1]+s*v[1], p[2]+s*v[2]
    # print('On line 1 : (x, y, z) = ({}, {}, {})'.format(x1, y1, z1))
    #それは直線２上のどこなのか
    x2, y2, z2 = q[0]+t*w[0], q[1]+t*w[1], q[2]+t*w[2]
    # print('On line 2 : (x, y, z) = ({}, {}, {})'.format(x2, y2, z2))

    return np.array([(x1+x2)/2,(y1+y2)/2,(z1+z2)/2])



# 例として、直線1が原点を通り、方向ベクトルが[1, 1, 1]で、直線2が原点を通り、方向ベクトルが[1, -1, 0]の場合
line1 = (np.array([0, 0, 0]), np.array([1, 1, 1]))
line2 = (np.array([1, 1, 0]), np.array([-1, -1, 0]))

min_distance_point = find_point_with_min_distance(line1, line2)

print("Point with minimum distance:", min_distance_point)