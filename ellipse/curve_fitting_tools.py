import numpy as np

def fitzgibbon_hyp_fit(x,y):
    
    num = np.size(x)
    if num != np.size(y):
        print("Function sizes incorrect")
        exit()
    
    # build design and constraint matrix
    C = np.array([[0,0,2],[0,-1,0],[2,0,0]])
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
    M = -np.matmul(np.linalg.inv(C),M)
    w,V = np.linalg.eig(M)
    
    # calculate constraint and identify hyperbola coefficients
    if w[0] > 0 and (4*V[0,0]*V[2,0] - V[1,0]**2) < 0 and \
        (w[1] < 0 or w[0] < w[1]) and (w[2] < 1 or w[0] < w[2]):
        A1 = V[:,0]
    elif w[1] > 0 and (4*V[0,1]*V[2,1] - V[1,1]**2) < 0 and \
        (w[2] < 1 or w[1] < w[2]):
        A1 = V[:,1]
    elif w[2] > 0 and (4*V[0,2]*V[2,2] - V[1,2]**2) < 0:
        A1 = V[:,2]
    else:
        print("Fitzgibbon hyperbola failed")
        raise ValueError("error!")
        # exit()
    A2 = np.matmul(T,A1)
    
    return np.concatenate((A1,A2))

def fitzgibbon_fit(x,y):
    
    num = np.size(x)
    if num != np.size(y):
        print("Function sizes incorrect")
        exit()
    
    # build design and constraint matrix
    # C = np.array([[0,0,2],[0,-1,0],[2,0,0]])
    C = np.array([[0,0,2],[0,-1,0],[2,0,0]])
    # C = np.array([[(2-e**2)**2-e**4, 0, -(2-e**2)**2-e**4], [0, (2-e**2)**2, 0],[-(2-e**2)**2-e**4, 0, (2-e**2)**2-e**4]],dtype=np.float64)
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
    if w[0] > 0 and (4*V[0,0]*V[2,0] - V[1,0]**2) > 0 and \
        (w[1] < 0 or w[0] < w[1]) and (w[2] < 1 or w[0] < w[2]):
        A1 = V[:,0]
    elif w[1] > 0 and (4*V[0,1]*V[2,1] - V[1,1]**2) > 0 and \
        (w[2] < 1 or w[1] < w[2]):
        A1 = V[:,1]
    elif w[2] > 0 and (4*V[0,2]*V[2,2] - V[1,2]**2) > 0:
        A1 = V[:,2]
    else:
        print("Fitzgibbon failed")
        raise ValueError("error!")
        # exit()
    A2 = np.matmul(T,A1)
    
    # determine ellipse parameters
    a = A1[0]
    b = A1[1]
    c = A1[2]
    d = A2[0]
    e = A2[1]
    f = A2[2]
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



    
def taubin_fit(x, y):
    
    num = np.size(x)
    if num != np.size(y):
        print("Function sizes incorrect")
        exit()
    
    bar_x = np.mean(x)
    bar_y = np.mean(y)
    z = (x-bar_x)**2 + (y - bar_y)**2
    bar_z = np.mean(z)
    Z0 = np.zeros((num,3))
    for i in range(0,num):
        Z0[i,1] = x[i] - bar_x
        Z0[i,2] = y[i] - bar_y
        Z0[i,0] = (z[i] - bar_z)/2/np.sqrt(bar_z)
    u,s,v = np.linalg.svd(Z0) # TODO: investigate characteristic equation speed and accuracy
    A = np.transpose(v)[:,2]
    a = A[0]/2/np.sqrt(bar_z)
    b = A[1]
    c = A[2]
    d = -a*bar_z
    
    # calculate circle dimensions (in px)
    x_c = -b/2/a + bar_x
    y_c = -c/2/a + bar_y
    R = np.sqrt(b**2 + c**2 - 4*a*d)/2/a
    return x_c, y_c, R

def sigmoid_circle(x, y, g, x_c, y_c, R, g_max=255, g_min=1, k=0.5):
    
    # set NLLS max boundaries
    tol = 0.03/50
    max_iter = 50
    
    # initialise sigmoid fit
    num = np.size(x)
    if num != np.size(y) or num != np.size(g):
        print("Sigmoid function sizes incorrect")
        exit()
    delta = np.ones((6,1))
    count = 0
    J = np.zeros((num,6))
    deltaf = np.zeros((num,1))
    
    # run NLLS
    while count < max_iter and np.linalg.norm(delta[:,0]) > tol:
        for i in range(0,num):
            beta = np.sqrt((x_c-x[i])**2 + (y_c-y[i])**2)
            alpha = np.exp(k*(R-beta))
            J[i,0] = (g_min - g_max)*alpha/((1+alpha)**2)*k*(x_c-x[i])/beta
            J[i,1] = (g_min - g_max)*alpha/((1+alpha)**2)*k*(y_c-y[i])/beta
            J[i,2] = -(g_min - g_max)*alpha/((1+alpha)**2)*k
            J[i,3] = -(g_min - g_max)*alpha/((1+alpha)**2)*(R-beta)
            J[i,4] = alpha/(1+alpha)
            J[i,5] = 1/(1+alpha)
            f = g_max + (g_min - g_max)/(1 + alpha)
            deltaf[i,0] = g[i] - f
    
        delta = np.matmul(np.matmul(np.linalg.pinv(np.matmul(J.T,J)),J.T),deltaf)
        x_c = x_c + delta[0,0]
        y_c = y_c + delta[1,0]
        R = R + delta[2,0]
        k = k + delta[3,0]
        g_max = g_max + delta[4,0]
        g_min = g_min + delta[5,0]
        
    return x_c, y_c, R

def sigmoid_conical(x, y, g, x_c, y_c, R, foc, g_max=255, g_min=1, k=0.5):
    
    # check size
    num = np.size(x)
    if num != np.size(y) or num != np.size(g):
        print("Sigmoid function sizes incorrect")
        exit()
    
    # set NLLS max boundaries
    tol = 0.3/50
    max_iter = 50
    
    # calculate unit vector set
    b = np.zeros((num,3))
    for i in range(0,num):
        b[i,:] = 1/np.sqrt(x[i]**2 + y[i]**2 +foc**2)*np.array([-x[i],-y[i],foc])
    
    # calculate mean direction angles and angle radius
    xbar = np.mean(x)
    ybar = np.mean(y)
    theta = np.arccos(foc/np.sqrt(xbar**2 + ybar**2 + foc**2))
    phi = np.arctan(ybar/xbar)
    theta_r = np.arctan(np.sqrt(x_c**2 + y_c**2)/foc)
    
    # initialise sigmoid fit
    delta = np.ones((6,1))
    count = 0
    J = np.zeros((num,6))
    deltaf = np.zeros((num,1))
    
    # run NLLS
    while count < max_iter and np.linalg.norm(delta[:,0]) > tol:
        for i in range(0,num):
            
            beta = np.cos(theta_r) - b[i,0]*np.sin(theta)*np.cos(phi) - \
                    b[i,1]*np.sin(theta)*np.sin(phi) - b[i,2]*np.cos(theta)
            alpha = np.exp(k*beta)
            J[i,0] = -(g_min - g_max)*alpha/((1+alpha)**2)*k* \
                        (b[i,0]*np.sin(phi) - b[i,1]*np.cos(phi))*np.sin(theta)
            J[i,1] = -(g_min - g_max)*alpha/((1+alpha)**2)*k* \
                        (-b[i,0]*np.cos(theta)*np.cos(phi) - \
                          b[i,1]*np.cos(theta)*np.sin(phi) + b[i,2]*np.sin(theta))
            J[i,2] = (g_min - g_max)*alpha/((1+alpha)**2)*k*np.sin(theta_r)
            J[i,3] = (g_min - g_max)*alpha/((1+alpha)**2)*beta
            J[i,4] = alpha/(1+alpha)
            J[i,5] = 1/(1+alpha)
            f = g_max + (g_min - g_max)/(1 + alpha)
            deltaf[i,0] = g[i] - f
    
        delta = np.matmul(np.matmul(np.linalg.pinv(np.matmul(J.T,J)),J.T),deltaf)
        phi = phi + delta[0,0]
        theta = theta + delta[1,0]
        theta_r = theta_r + delta[2,0]
        k = k + delta[3,0]
        g_max = g_max + delta[4,0]
        g_min = g_min + delta[5,0]
        
    x_c = foc*np.tan(theta)*np.sin(phi)
    y_c = foc*np.tan(theta)*np.cos(phi)
    R = np.sqrt(x_c**2 + y_c**2)*np.tan(theta_r)

    
    return x_c, y_c, R
    