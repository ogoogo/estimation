import numpy as np
from scipy.optimize import minimize

def slsqp(x,y,e):
    num = np.size(x)
    if num != np.size(y):
        print("Function sizes incorrect")
        exit()

    # build design and constraint matrix
    D = np.zeros((num,6))
    for i in range(0,num):
        D[i,0] = x[i]**2
        D[i,1] = x[i]*y[i]
        D[i,2] = y[i]**2
        D[i,3] = x[i]
        D[i,4] = y[i]
        D[i,5] = 1
    C = np.array([[0, 0, 2,0,0,0], [0, -1, 0,0,0,0],[2, 0, 0, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]],dtype=np.float64)
    e_calc = (2-e**2)**2-e**4
    B = np.array([[e_calc, 0, -e_calc,0,0,0], [0, (2-e**2)**2, 0,0,0,0],[-e_calc, 0, e_calc, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]],dtype=np.float64)


    # 目的関数の定義
    def objective_function(a):
        sum = 0
        for i in range(0,num):
            sum += a[0]*x[i]**2 + a[1]*x[i]*y[i] + a[2]*y[i]**2 + a[3]*x[i] + a[4]*y[i] + a[5]
        
        return sum
        # return np.dot(a.T, np.dot(D.T, np.dot(D, a))) 

    # 制約条件
    def constraint1(a):
        # return np.dot(a.T, np.dot(C, a))
        return 4*a[0]*a[2] - a[1]**2 -1
    

    def constraint2(a):

        # build design and constraint matrix

        # return np.dot(a.T, np.dot(D.T, np.dot(D, a)) + l1*(1 - np.dot(a.T, np.dot(C, a)) ) + l2*np.dot(a.T, np.dot(B, a)))
        return np.dot(a.T, np.dot(B, a))



    # 初期値と制約条件を設定
    a0 = np.array([1000,50,1000,1,1,-1000], dtype=np.float64)
    # constraints = ({'type': 'eq', 'fun': constraint1},
    #             {'type': 'eq', 'fun': constraint2})

    constraints = ({'type': 'eq', 'fun': constraint1})
    # ラグランジュの未定乗数法を使って最適化問題を解く
    result = minimize(objective_function, a0, constraints=constraints, method='SLSQP')

    # 結果の表示
    print("最適解:", result.x)
    print("最小値:", result.fun)
    a = result.x[0]
    b = result.x[1]
    c = result.x[2]
    d = result.x[3]
    e = result.x[4]
    f = result.x[5]
    return np.array([a, b, c, d, e, f])