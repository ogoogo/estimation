import numpy as np
def fitzgibbon_fit(x,y,e):
    
    num = np.size(x)
    if num != np.size(y):
        print("Function sizes incorrect")
        exit()
    
    # build design and constraint matrix
    # C = np.array([[0,0,2],[0,-1,0],[2,0,0]])
    # C = np.array([[0,0,2],[0,-1,0],[2,0,0]])
    C = np.array([[(2-e**2)**2-e**4, 0, -(2-e**2)**2-e**4], [0, (2-e**2)**2, 0],[-(2-e**2)**2-e**4, 0, (2-e**2)**2-e**4]],dtype=np.float64)
    D1 = np.zeros((num,3))
    D2 = np.ones((num,3))
    for i in range(0,num):
        D1[i,0] = x[i]**2
        D1[i,1] = x[i]*y[i]
        D1[i,2] = y[i]**2
        D2[i,0] = x[i]
        D2[i,1] = y[i]
        
    # calculate scatter matrix
    S1 = np.matmul(np.transpose(D1),D1)
    S2 = np.matmul(np.transpose(D1),D2)
    S3 = np.matmul(np.transpose(D2),D2)
    
    # calculate eigenvalue M matrix
    T = -np.matmul(np.linalg.inv(S3),np.transpose(S2))
    M = np.matmul(S2,T)
    M = np.add(S1,M)
    M = np.matmul(np.linalg.inv(C),M)
    w,V = np.linalg.eig(M)
    
    # calculate constraint and identify ellipse coefficients
    if w[0] < 0 and (4*V[0,0]*V[2,0] - V[1,0]**2) > 0 and \
        (w[1] > 0 or w[0] > w[1]) and (w[2] > 1 or w[0] > w[2]):
        A1 = V[:,0]
    elif w[1] < 0 and (4*V[0,1]*V[2,1] - V[1,1]**2) > 0 and \
        (w[2] > 1 or w[1] > w[2]):
        A1 = V[:,1]
    elif w[2] < 0 and (4*V[0,2]*V[2,2] - V[1,2]**2) > 0:
        A1 = V[:,2]
    else:
        print("Fitzgibbon failed")
        exit()
    # A1 = V[:,0]
    
    A2 = np.matmul(T,A1)
    
    # determine ellipse parameters
    a = A1[0]
    b = A1[1]
    c = A1[2]
    d = A2[0]
    e = A2[1]
    f = A2[2]

    print((2*((a-c)**2 + b**2)**0.5) / (((a-c)**2 + b**2)**0.5 + a + c))

    print(a,b,c)
    x_c = (2*c*d - b*e)/(b*b - 4*a*c) 
    y_c = (2*a*e - b*d)/(b*b - 4*a*c)
    major = np.sqrt(2*(a*e*e + c*d*d + f*b*b - b*d*e - a*c*f) / \
                (b - 4*a*c)/(np.sqrt((a-c)**2 + b**2) - (a + c))) 
    minor = np.sqrt(2*(a*e*e + c*d*d + f*b*b - b*d*e - a*c*f) / \
                (4*a*c - b)/(np.sqrt((a-c)**2 + b**2) + a + c)) 
    if b == 0 and a < c:
        theta = 0
    elif b == 0 and c < a:
        theta = np.pi/2
    elif a < c:
        theta = 0.5*np.arctan(2*b/(a-c))
    else:
        theta = np.pi + 0.5*np.arctan(2*b/(a-c))
    
    return np.array([a,b,c,d,e,f])
