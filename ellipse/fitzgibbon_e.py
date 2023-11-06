import numpy as np
def bfgs_fit(x, y, e):
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
    B = np.array([[(2-e**2)**2-e**4, 0, -(2-e**2)**2-e**4,0,0,0], [0, (2-e**2)**2, 0,0,0,0],[-(2-e**2)**2-e**4, 0, (2-e**2)**2-e**4, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]],dtype=np.float64)

    def obj_func(al):
        """評価関数"""
        a = al[0:6]
        l1 = al[6]
        l2 = al[7]
        return np.dot(a.T, np.dot(D.T, np.dot(D, a)) + l1*(1 - np.dot(a.T, np.dot(C, a)) ) + l2*np.dot(a.T, np.dot(B, a)))
        # return np.dot(a.T, np.dot(D.T, np.dot(D, a))) + l1*(1 - np.dot(a.T, np.dot(C, a)) ) 


    def obj_grad(al):
        """評価関数の勾配"""
        a = al[0:6]
        l1 = al[6]
        l2 = al[7]
        # tmp = 2*np.dot(D.T, np.dot(D,a)) - 2*l1*np.dot(C,a) + 2*l2*np.dot(B,a)
        tmp = 2*np.dot(D.T, np.dot(D,a)) - 2*l1*np.dot(C,a)

        return np.array([tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], 1-np.dot(a.T, np.dot(C,a)), np.dot(a.T, np.dot(B,a))])
        # return np.array([tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], 1-np.dot(a.T, np.dot(C,a))])


    al = np.array([1,0.3,1,0.3,0.3,-1,10,10], dtype=np.float64)
    # al = np.array([1,0.3,1,0.3,0.3,-1,10], dtype=np.float64)

    B_ = np.eye(len(al)) # ヘッセ行列の近似行列
    tau = 0.9  # 直線探索のパラメータ
    xi = 0.3  # 直線探索のパラメータ

    grad = obj_grad(al)
    # print(grad)

    for i in range(10):
        f = obj_func(al)
        d = -np.dot(np.linalg.inv(B_), grad)  # 探索方向
        # print(d)

        # Armijo条件を用いた直線探索
        alpha = 1
        temp = xi * np.dot(grad, d)
        # print(al + alpha * d)
        # print(obj_func(x,y,e,al + alpha * d))
        # print(f + alpha * temp)
        while obj_func(al + alpha * d) > f + alpha * temp:
            alpha *= tau

        al += alpha * d  # 解の更新

        old_grad = grad
        grad = obj_grad(al)

        s = (alpha * d).reshape([-1, 1])  # 解の変化量
        b_ = (grad - old_grad).reshape([-1, 1])  # 勾配の変化量
        Bs = np.dot(B_, s)

        B_ -= (Bs * Bs.T) / np.dot(s.T, Bs) - (b_ * b_.T) / np.dot(s.T, b_)
        # ヘッセ行列の近似行列を更新
        
        print()
        print(al)
    # C = np.array([[0, 0, 2,0,0,0], [0, -1, 0,0,0,0],[2, 0, 0, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]])
    print(al[0:6]@C@al[0:6])

    # determine ellipse parameters
    a = al[0]
    b = al[1]
    c = al[2]
    d = al[3]
    e = al[4]
    f = al[5]
    x_c = (2 * c * d - b * e) / (b * b - 4 * a * c)
    y_c = (2 * a * e - b * d) / (b * b - 4 * a * c)
    major = np.sqrt(
        2
        * (a * e * e + c * d * d + f * b * b - b * d * e - a * c * f)
        / (b - 4 * a * c)
        / (np.sqrt((a - c) ** 2 + b**2) - (a + c))
    )
    minor = np.sqrt(
        2
        * (a * e * e + c * d * d + f * b * b - b * d * e - a * c * f)
        / (4 * a * c - b)
        / (np.sqrt((a - c) ** 2 + b**2) + a + c)
    )
    if b == 0 and a < c:
        theta = 0
    elif b == 0 and c < a:
        theta = np.pi / 2
    elif a < c:
        theta = 0.5 * np.arctan(2 * b / (a - c))
    else:
        theta = np.pi + 0.5 * np.arctan(2 * b / (a - c))

    return np.array([a, b, c, d, e, f])
