import numpy as np
from scipy.optimize import least_squares

def estimate_r(coe):
    R = 1737.4
    y0 = 247
    x0 = 329.5
    (a,b,c,d,f,g) = coe
    g = a*x0**2 + b*x0*y0 + c*y0**2 + d*x0 + f*y0 +g
    f = b*x0 + 2*c*y0 + f
    d = 2*a*x0 +b*y0 + d

    # 行列MとAを定義
    M1 = np.array([[a,b/2,d/2], [b/2, c, f/2], [d/2, f/2, g]])
    K = np.array([[50*494/3.66,0,0],[0,50*494/3.66,0],[0,0,1]])  # 3x3の行列
    C = K@M1@K
    A = np.array([[1/R**2,0,0],[0,1/R**2,0],[0,0,1/R**2]]) # 3x3の行列

    w,V = np.linalg.eig(np.linalg.inv(A)@C)
    if (w[0] < 0 and w[1] > 0 and w[2] > 0) or (w[0] > 0 and w[1] < 0 and w[2] < 0):
        v = V[:,0]/np.linalg.norm(V[:,0])
        l = w[0]
    elif (w[1] < 0 and w[0] > 0 and w[2] > 0) or (w[1] > 0 and w[0] < 0 and w[2] < 0):
        v = V[:,1]/np.linalg.norm(V[:,1])
        l = w[1]
    elif(w[2] < 0 and w[1] > 0 and w[0] > 0) or (w[2] > 0 and w[1] < 0 and w[0] < 0):
        v = V[:,2]/np.linalg.norm(V[:,2])
        l = w[2]
    else:
        print("esitimate-r failed")
        exit()
    


    p = ((np.trace(C)-l*np.trace(A))/(l*(np.dot(v, np.dot(A, np.dot(A, v)))) - l*np.dot(v, np.dot(A,v))*np.trace(A)))**0.5
    r = p*v
    # print(r)
    return r
    # # 初期推定値としてrを設定
    # r0 = np.random.rand(3)  # ランダムな初期値を使用する例
    # # print(r0)
    # # print(r0.reshape([3,1]))
    # # 最小二乗法の目的関数を定義
    # def fun(r):
    #     r_ = np.array([[r[0],r[1],r[2]]])
    #     return np.linalg.norm( - A@r_.reshape([3,1])@r_@A + (r_@A@r_.reshape([3,1])-1)@A, 'fro')**2

    # # 最小二乗法を実行
    # result = least_squares(fun, r0, method='trf')

    # # 結果を表示
    # r_solution = result.x
    # print('解(r):')
    # print(r_solution)

if __name__ == "__main__":
    estimate_r([0.692516954302689,0.001525634136022375,0.7213999864456683,-414.35377901181386,-448.84188056892737,114170.30887439639])