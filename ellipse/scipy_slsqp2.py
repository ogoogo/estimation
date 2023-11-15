import numpy as np
from scipy.optimize import minimize
import estimate_r
import matplotlib.pyplot

def slsqp2(x,y,e,f,sun_dlp,x_t,y_t,width,height,x_3_e):
    R = 1737.4
    num = np.size(x)
    num_t = np.size(x_t)
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
    D_t = np.zeros((num_t,6))
    for i in range(0,num_t):
        D_t[i,0] = x_t[i]**2
        D_t[i,1] = x_t[i]*y_t[i]
        D_t[i,2] = y_t[i]**2
        D_t[i,3] = x_t[i]
        D_t[i,4] = y_t[i]
        D_t[i,5] = 1
    C = np.array([[0, 0, 2,0,0,0], [0, -1, 0,0,0,0],[2, 0, 0, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]],dtype=np.float64)
    e_calc = (2-e**2)**2-e**4
    B = np.array([[(2-e**2)**2-e**4, 0, -(2-e**2)**2-e**4,0,0,0], [0, (2-e**2)**2, 0,0,0,0],[-(2-e**2)**2-e**4, 0, (2-e**2)**2-e**4, 0, 0, 0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]],dtype=np.float64)


    # 目的関数の定義
    def objective_function(a):
        # sum = 0
        # for i in range(0,num):
        #     sum += a[0]*x[i]**2 + a[1]*x[i]*y[i] + a[2]*y[i]**2 + a[3]*x[i] + a[4]*y[i] + a[5]
        
        # return sum
        return np.dot(a, np.dot(D.T, np.dot(D, a))) 
    def objective_grad(a):
        # sum = 0
        # for i in range(0,num):
        #     sum += a[0]*x[i]**2 + a[1]*x[i]*y[i] + a[2]*y[i]**2 + a[3]*x[i] + a[4]*y[i] + a[5]
        
        # return sum
        return 2*np.dot(D.T, np.dot(D, a))
    # 制約条件

    def constraint1(a):
        # return np.dot(a.T, np.dot(C, a))
        return np.dot(a, np.dot(B, a))
    

    def constraint2(a):
        r = estimate_r.estimate_r(a)
        print(r)
        n = sun_dlp/np.linalg.norm(sun_dlp)
        a1 = np.array([1,0,0])
        c1 = a1 - np.dot(a1,n)*n
        p = c1/np.linalg.norm(c1)
        s = np.cross(n,p)

        Rct = np.array([[p[0],s[0],n[0]],
                        [p[1],s[1],n[1]],
                        [p[2],s[2],n[2]]],dtype=np.float64) 

        
        Q = np.array([[1/R**2,0,0],
                        [0,1/R**2,0],
                        [0,0,1/R**2]],dtype=np.float64) 
        
        M = (1/np.dot(p,np.dot(Q,p)))**0.5
        m = (1/np.dot(s,np.dot(Q,s)))**0.5

        T = np.array([[1/M**2,0,0],
                        [0,1/m**2,0],
                        [0,0,-1]],dtype=np.float64) 
        
        H1 = np.array([[Rct[0][0], Rct[0][1], r[0]],
                      [Rct[1][0], Rct[1][1], r[1]],
                      [Rct[2][0], Rct[2][1], r[2]]],dtype=np.float64)
        
        # K = np.diag([f*1000/7.4,f*1000/7.4,1])
        focus = 50
        K = np.array([[focus*1000/7.4,0,0],
                      [0,focus*1000/7.4,0],
                      [0,0,1]],dtype=np.float64)
        H = K@H1
        T_dash = np.linalg.inv(H.T)@T@np.linalg.inv(H)

        y0 = 247
        x0 = 329.5
        a_t = T_dash[0][0]
        b_t = 2*T_dash[0][1]
        c_t = T_dash[1][1]
        d_t = 2*T_dash[0][2]
        f_t = 2*T_dash[1][2]
        g_t = T_dash[2][2]
        g_t = a_t*x0**2 + b_t*x0*y0 + c_t*y0**2 -d_t*x0 -f_t*y0 +g_t
        f_t = -b_t*x0 - 2*c_t*y0 + f_t
        d_t= -2*a_t*x0 -b_t*y0 + d_t

        b = np.array([a_t,b_t,c_t,d_t,f_t,g_t]) 
        x = np.arange(0,width,1)
        y = np.arange(0,height,1)
        x,y = np.meshgrid(x,y)
        z = b[0]*x**2 + b[1]*x*y + b[2]*y**2 + b[3]*x + b[4]*y + b[5]
        plt = matplotlib.pyplot.contour(z,[0])
        x = plt.collections[0].get_paths()[0].vertices[:,0]
        y = plt.collections[0].get_paths()[0].vertices[:,1]
        edge_t = np.array([x,y]).T
        edge_r = np.zeros((np.size(x),2))
        mat_r = np.array([[sun_dlp[1],-sun_dlp[0]],[sun_dlp[0],sun_dlp[1]]])
        xmax = -659
        xmin = 659
        for i in range(np.size(x)):
            # print(mat_r)
            edge_r[i] = np.dot(edge_t[i],mat_r.T)
            # print(edge_r[i])
            if xmax < edge_r[i,1]:
                xmax = edge_r[i,1]
                # print(xmax)
                end = i
            if xmin > edge_r[i,1]:
                xmin = edge_r[i,1]
                # print(xmin)
                start = i
        
        x_3 = np.array([edge_t[start][0], edge_t[start][1]],dtype=np.float64)
        print("x_3")
        print(x_3)
        print(x_3_e)
        return np.linalg.norm(x_3 - x_3_e)
        # return np.dot(b, np.dot(D_t.T, np.dot(D_t, b))) 


        



    # 初期値と制約条件を設定
    a0 = np.array(f, dtype=np.float64)
    # constraints = ({'type': 'eq', 'fun': constraint1},
    #             {'type': 'eq', 'fun': constraint2})

    constraints = ({'type': 'eq', 'fun': constraint2})
    # ラグランジュの未定乗数法を使って最適化問題を解く
    result = minimize(objective_function, a0, jac = objective_grad, constraints=constraints, method='SLSQP')

    # 結果の表示
    print("最適解:", result.x)
    print("最小値:", result.fun)
    a = result.x[0]
    b = result.x[1]
    c = result.x[2]
    d = result.x[3]
    e = result.x[4]
    f = result.x[5]
    # print(np.array([a, b, c, d, e, f]))
    return np.array([a, b, c, d, e, f])