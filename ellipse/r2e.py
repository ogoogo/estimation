import numpy as np


def r2e(r):
    R = 1737.4
    A = np.array([[1/R**2, 0, 0],[0, 1/R**2, 0], [0, 0, 1/R**2]])
    M = np.dot(A,np.dot(r, np.dot(r,A))) - (np.dot(r, np.dot(A,r)) - 1)*A
    a = M[0,0]
    b = 2*M[0,1]
    c = M[1,1]

    print([a,b,c])

    e = ((2*((a - c)**2 - b**2)**0.5)/(((a - c)**2 - b**2)**0.5 + a + c))**0.5
    print(e)
                                             

                                             
