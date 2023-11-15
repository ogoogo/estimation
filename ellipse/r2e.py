import numpy as np


def r2e(r):
    R = 1737.4
    A = np.array([[1/R**2, 0, 0],[0, 1/R**2, 0], [0, 0, 1/R**2]])
    r_vec = np.array([[r[0]],[r[1]],[r[2]]])
    M = A @ r_vec @ r_vec.T @ A - (r_vec.T @ A @ r_vec - 1) * A
    # M = np.dot(A,np.dot(r, np.dot(r,A))) - (np.dot(r, np.dot(A,r)) - 1)*A
    a = M[0,0]
    b = 2*M[0,1]
    c = M[1,1]

    print([a,b,c])
    # print(((a - c)**2 - b**2))
    print("a")
    print((a - c)**2 - b**2)
    e = abs((2*abs((a - c)**2 - b**2)**0.5)/(abs((a - c)**2 - b**2)**0.5 + a + c))**0.5
    print("e")
    print(e)
    return e
                                             

                                             
