# Applies horizon estimation model according to:
# Christian, John A. "Accurate planetary limb localization 
# for image-based spacecraft navigation." 
# Journal of Spacecraft and Rockets 54.3 (2017): 708-730.
import numpy as np
import pandas as pd
from PIL import Image
import scipy.integrate
import math
import warnings
import curve_fitting_tools
import fitzgibbon_e_2
import scipy_slsqp
import random
import matplotlib.pyplot
import estimate_r
import christian_robinson
import scipy_slsqp2
import scipy_powell

def main_christian(ID,date,option,e,focus,F):

    # date = "1102"
    # ID = 2
    # e = 0.1
    # constants
    R_EARTH = 6371E3

    # image settings
    IMAGE_FILE_NAME = '../output/' +str(date)+'/images/raw/' + str(ID) +'.png'
    D_TRUE = 408+R_EARTH # true distance to Earth [m]

    # set sun earth direction (TODO: ephemeris-based approach)
    u = np.array([1,0])

    # number of illumination scan lines
    m = 100 # number of illumination scan lines across image

    # processing settings
    step_size = 2 # image scan step size
    window_radius = 4 # Sobel window size
    zernike_radius = 4 # Zernike moment size (NOTE: cannot be larger than window_size) 
    sigma_psf = 0.2 # pixel point spread function
    d_min = 100 # usually 5*radius_px_body/100
    N_T = 500 # maximum RANSAC tests
    shape = 1 # hyperbola = 0, ellipse = 1
    line_thickness = 1 # pixel line fit thinkness
    bright_thresh = 20 # bright pixel threshold

    # algorithm parameter
    gradient_ratio = 0.6
    direction_ratio = 0.3



    # CSVファイルのパスを指定
    csv_file_path = "../output/" + str(date) + "/information.csv"

    # CSVファイルを読み込む
    df = pd.read_csv(csv_file_path, header=None)

    # 6行目のデータを取得
    row = df.iloc[ID-1]  # インデックスは0から始まるため

    # 10番目と11番目の値を取得し、配列に格納
    value_10 = -row[29]  # インデックスは0から始まるため
    value_11 = -row[30]
    values = np.array([value_10, value_11])
    dis = np.linalg.norm(values)
    u = values/dis
    u2 = -u

    sun_dlp = np.array([row[29],row[30],row[31]])
    sun_dlp = sun_dlp/np.linalg.norm(sun_dlp)
    f_answer = np.array(row[44:50])
    # print(f_answer)

    # 結果を表示
    print("取得した値:", u)

    
        

    # read image from file and convert to black and white
    original_image = Image.open(IMAGE_FILE_NAME).convert('RGB')
    image = np.array(Image.open(IMAGE_FILE_NAME).convert('L'))
    height, width = image.shape
    
    def insert_pixel(image,x,y,color,size):
        
        height,width = np.array(image.convert('L')).shape
        
        for i in range(-size,size+1):
            for j in range(-size,size+1):
                
                if int(y)+j >= height-1 or int(y)+j <= 0 or \
                    int(x)+i >= width-1 or int(x)+i <= 0:
                    continue
                image.putpixel((int(x)+i,int(y)+j), color)
    
        return image
    
    # perform Gaussian blur on image
    # image = scipy.ndimage.filters.gaussian_filter(image, sigma=2.5, order=0,
    #                                                                     output=np.uint8,mode="constant", cval=0.0,
    #                                                                     truncate=2.0)
    # image = image.astype(npbright_thresh32)
    
    # determine bright pixel threshold
    bright_thresh = 50
    
    def find_edge(u,bright_thresh):
        # generate horizon scan lines
        b = np.array([[1,0],[-1,0],[0,1],[0,-1]])
        edge = np.zeros((2*m,2))
        pos = 0
        for i in range(0,4,1):
            
            # check if border is appropriate
            if np.dot(b[i],u) > 0:
                
                # determine each corner start position
                if b[i,0] == 0:
                    yy = np.ones(m)*(0 if b[i,1] > 0 else height-1)
                    xx = np.arange(0, width, width/m)
                else:
                    yy = np.arange(0, height, height/m)
                    xx = np.ones(m)*(0 if b[i,0] > 0 else width-1)
                    
        
                # move along each strip and determine if encounter limb
                for j in range(0, m, 1):
                    
                    y = int(yy[j])
                    x = int(xx[j])
                    
                    # move along the strip and break if suspected limb
                    count = 0
                    
                    while y < height and y >= 0 and x < width and x >= 0:
                        
                        # might be a limb, increase counter
                        if image[y,x] > bright_thresh:
                            count = count + 1
                        
                        # this could be a limb, hence include in measured position
                        if count > window_radius:
                            y = y - int(window_radius*step_size*u[1])
                            x = x - int(window_radius*step_size*u[0])
                            edge[pos,0] = np.float64(y)
                            edge[pos,1] = np.float64(x)
                            pos = pos + 1

                            break

                        # increase along strip
                        y = y + int(step_size*u[1])
                        x = x + int(step_size*u[0])
 
        
        # reduce suspected positions to pos
        edge = edge[:pos,:]
        
        # apply Sobel-based pixel edge localisation
        sobel = np.array([[-1,0,1],[-2,0,2],[-1,0,1]])
        for i in range(0, pos):
            y = int(edge[i,0])
            x = int(edge[i,1])
            
            # consider across each window mass
            g_max = 0
            for j in range(-window_radius,window_radius):
                for k in range(-window_radius,window_radius):
                    
                    # check if image corner encountered
                    if y+j+2 >= height or y+j-1 < 0 or x+k+2 >= width or x+k-1 < 0:
                        continue
                    
                    # determine gradient
                    g_x = np.sum(np.multiply(sobel,image[y+j-1:y+j+2,x+k-1:x+k+2]))
                    g_y = np.sum(np.multiply(np.transpose(sobel),image[y+j-1:y+j+2,x+k-1:x+k+2]))
                    g = np.sqrt(g_x**2 + g_y**2)
                    
                    
                    # determine max
                    if g > g_max:
                        edge[i,0] = np.float64(y + j)
                        edge[i,1]= np.float64(x + k)
                        g_max = g 
        return edge,pos
                    
    edge,pos = find_edge(u,50)
    edge_t,pos_t = find_edge(u2,30)
    edge_r = np.zeros((pos_t,2))
    xmax = -659
    xmin = 659
    mat_r = np.array([[u2[1],-u2[0]],u2])
    for i in range(pos_t):
        
        # print(mat_r)
        edge_r[i] = np.dot(edge_t[i],mat_r)
        # print(edge_r[i])
        if xmax < edge_r[i,1]:
            xmax = edge_r[i,1]
            # print(xmax)
            end = i
        if xmin > edge_r[i,1]:
            xmin = edge_r[i,1]
            # print(xmin)
            start = i
    
    # mean = int((xmin + xmax)/2)
    # print(mean)
    # x_3_r = edge_r[edge_r[:, 1] == mean]
    # absolute_differences = np.abs(edge_r[:, 1] - mean)
    # closest_row_index = np.argmin(absolute_differences)
    # x_3_r = edge_r[closest_row_index]
    # x_3 = np.dot(x_3_r,mat_r.T)
    # x_3 = np.array([x_3[1],x_3[0]],dtype=np.float64)
    
    # original_image = insert_pixel(original_image,x_3[0],x_3[1],(0,255,255),3)
    
    x_1 = np.array([edge_t[start,1],edge_t[start,0]])
    original_image = insert_pixel(original_image,x_1[0],x_1[1],(0,255,255),3)
    # print(x_1)
    x_2 = np.array([edge_t[end,1],edge_t[end,0]])
    original_image = insert_pixel(original_image,x_2[0],x_2[1],(0,255,255),3)
    # print(x_2)
    x_c = np.array([(x_1[0]+x_2[0])/2, (x_1[1]+x_2[1])/2])
    # print(x_c)
    radius = np.linalg.norm(x_1-x_2)/2
    # print(radius)
    a_i = 1.1
    b_i = 0
    c_i = 1
    d_i = -2*x_c[0]
    f_i = -2*x_c[1]
    g_i = x_c[0]**2 + x_c[1]**2 -radius**2
    F_ini = np.array([a_i,b_i,c_i,d_i,f_i,g_i])
    print(F_ini)


    
    # for sub-pixel edge localisation, apply Zernike moments
    def zernike(n,m,window_radius,j,k):
        # determine mask parameters
        N = window_radius*2 + 1
        u_i = 2*k/N
        v_i = 2*j/N

        # define integral
        def integ(v,u,n,m,mode):
            
            # calculate radius from centre and check if outside circle
            r = np.sqrt(u**2 + v**2)
            if r > 1:
                return 0.0

            # calculate radial polynomial
            theta = np.arctan2(v,u)
            R_nm = 0
            temp = np.arange(0,(n-np.abs(m))/2 + 1,1)
            for s in temp:
                R_nm = R_nm + (-1)**s*math.gamma(n-s+1)*r**(n-2*s)/ \
                        math.gamma(s+1)*math.gamma((n + abs(m))/2 - s+1)* \
                        math.gamma((n - abs(m))/2 - s+1)
            
            # calculate polynomial
            if mode == 'real':
                T_nm = R_nm*np.cos(m*theta)
            else:
                T_nm = R_nm*np.sin(m*theta)

            return T_nm
        
        options={'limit':5}
        warnings.filterwarnings('ignore')
        M_nm_real = scipy.integrate.nquad(integ, [[u_i-1/N, u_i+1/N], [v_i-1/N, v_i+1/N]], \
                                        args=(n,m,'real'),opts=[options,options])
        M_nm_imag = scipy.integrate.nquad(integ, [[u_i-1/N, u_i+1/N], [v_i-1/N, v_i+1/N]], \
                                        args=(n,m,'complex'),opts=[
                                            options,options])
        
        return M_nm_real[0] + 1j*M_nm_imag[0]

    # construct Zernike moments
    # M_11 = np.zeros((zernike_radius*2+1,zernike_radius*2+1),dtype=complex)
    # M_20 = np.zeros((zernike_radius*2+1,zernike_radius*2+1),dtype=complex)    
    # for j in range(-zernike_radius,zernike_radius+1):
    #     for k in range(-zernike_radius,zernike_radius+1):
    #         M_11[zernike_radius+j,zernike_radius+k] = zernike(1,1,zernike_radius,k,j)
    #         M_20[zernike_radius+j,zernike_radius+k] = zernike(2,0,zernike_radius,k,j)

    # perform sub-pixel determination by Zernike
    for i in range(0, pos):
        # y = int(edge[i,0])
        # x = int(edge[i,1])
        # A_11 = 0
        # A_20 = 0
        # for j in range(-zernike_radius,zernike_radius+1):
        #     for k in range(-zernike_radius,zernike_radius+1):
        #         A_11 = A_11 + image[y+j,x+k]*M_11[j+zernike_radius,k+zernike_radius]
        #         A_20 = A_20 + image[y+j,x+k]*M_20[j+zernike_radius,k+zernike_radius]
                
        # # determine orientation and length from estimated pixel position
        # A_20 = np.real(A_20)
        # psi = np.arctan2(np.imag(A_11),np.real(A_11))
        # Aprime_11 = np.real(A_11)*np.cos(psi) + np.imag(A_11)*np.sin(psi)
        # w = 1.66*sigma_psf
        # l = (1 - w**2 - np.sqrt((w**2-1)**2 - 2*w**2*A_20/Aprime_11))/w**2
        
        # # adjust estimated pixel position
        # edge[i,0] = edge[i,0] + (zernike_radius*2+1)*l*np.cos(psi)
        # edge[i,1] = edge[i,1] + (zernike_radius*2+1)*l*np.sin(psi)
        try:
            original_image.putpixel((int(edge[i,1])+1,int(edge[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge[i,1])+1,int(edge[i,0])), (255,0,0))
            original_image.putpixel((int(edge[i,1])+1,int(edge[i,0])+1), (255,0,0))
            original_image.putpixel((int(edge[i,1]),int(edge[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge[i,1]),int(edge[i,0])), (255,0,0))
            original_image.putpixel((int(edge[i,1]),int(edge[i,0])+1), (255,0,0))
            original_image.putpixel((int(edge[i,1])-1,int(edge[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge[i,1])-1,int(edge[i,0])), (255,0,0))
            original_image.putpixel((int(edge[i,1])-1,int(edge[i,0])+1), (255,0,0))
        except Exception as e:
            print(e)
    for i in range(0, pos_t):
        try:
            original_image.putpixel((int(edge_t[i,1])+1,int(edge_t[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge_t[i,1])+1,int(edge_t[i,0])), (255,0,0))
            original_image.putpixel((int(edge_t[i,1])+1,int(edge_t[i,0])+1), (255,0,0))
            original_image.putpixel((int(edge_t[i,1]),int(edge_t[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge_t[i,1]),int(edge_t[i,0])), (255,0,0))
            original_image.putpixel((int(edge_t[i,1]),int(edge_t[i,0])+1), (255,0,0))
            original_image.putpixel((int(edge_t[i,1])-1,int(edge_t[i,0])-1), (255,0,0))
            original_image.putpixel((int(edge_t[i,1])-1,int(edge_t[i,0])), (255,0,0))
            original_image.putpixel((int(edge_t[i,1])-1,int(edge_t[i,0])+1), (255,0,0))
        except Exception as e:
            print(e)
    # perform fitzgibbon fit using RANSAC test
    N_test = 0
    N_min = pos
    while N_test < N_T:
        
        # increase the test count
        N_test = N_test + 1
        
        # random sample of points
        # m = random.sample(range(pos),6)
        
        # determine fit
        # if shape == 0:
        #     f = curve_fitting_tools.fitzgibbon_hyp_fit(edge[m,1], edge[m,0])
        # else:
        #     f = curve_fitting_tools.fitzgibbon_fit(edge[:,1], edge[:,0])
        #     # f = fitzgibbon_e_2.fitzgibbon_fit(edge[:,1], edge[:,0], e)
        #     # f = scipy_slsqp.slsqp(edge[:,1], edge[:,0], e)

        if option == "fitz":
            f = curve_fitting_tools.fitzgibbon_fit(edge[:,1], edge[:,0])
        elif option == "e":
            f = fitzgibbon_e_2.fitzgibbon_fit(edge[:,1], edge[:,0], e)
        elif option == "slsqp":
            f = scipy_slsqp.slsqp(edge[:,1], edge[:,0], e, F_ini)
        elif option == "slsqp2":
            f = scipy_slsqp.slsqp(edge[:,1], edge[:,0], e, F)
        elif option == "term":
            f = scipy_slsqp2.slsqp2(edge[:,1], edge[:,0], e, F_ini, sun_dlp, edge_t[:,1], edge_t[:,0],width,height,x_3)
        elif option == "powell":
            f = scipy_powell.powell(edge[:,1], edge[:,0], e, F_ini, sun_dlp, edge_t[:,1], edge_t[:,0],width,height,x_3)

        else:
            raise ValueError("option is wrong")
        
        # check distance with all points
        N_k = 0
        elem = []
        for i in range(0,pos):
            y = edge[i,0]
            x = edge[i,1]
            d = (f[0]*x**2 + f[1]*x*y + f[2]*y**2 + \
                f[3]*x + f[4]*y + f[5]) / \
                np.sqrt((2*f[0]*x + f[1]*y + f[3])**2 + \
                        (2*f[2]*y + f[1]*x + f[4])**2)
            if d_min > d:
                elem.append(i)
                N_k = N_k + 1
        
        # check if number of valid points exceed 
        if N_k < N_min:
            N_min = N_k
        else: 
            break
    
    # perform best fit using valid image set (using either hyperbola or ellipsoid method)
    # TODO: Check if image estimated altitude implies hyperbola or ellipse
    print("number of identified horizon points : " + str(N_k))
    edge = edge[elem,:]
    # if shape == 0:
    #     f = curve_fitting_tools.fitzgibbon_hyp_fit(edge[:,1], edge[:,0])
    if option == "fitz":
        f = curve_fitting_tools.fitzgibbon_fit(edge[:,1], edge[:,0])
    elif option == "e":
        f = fitzgibbon_e_2.fitzgibbon_fit(edge[:,1], edge[:,0], e)
    elif option == "slsqp":
        f = scipy_slsqp.slsqp(edge[:,1], edge[:,0], e,F_ini)
    elif option == "slsqp2":
        f = scipy_slsqp.slsqp(edge[:,1], edge[:,0], e,F)
    elif option == "term":
        f = scipy_slsqp2.slsqp2(edge[:,1], edge[:,0], e, F_ini, sun_dlp, edge_t[:,1], edge_t[:,0],width,height,x_3)
    elif option == "powell":
        f = scipy_powell.powell(edge[:,1], edge[:,0], e, F_ini, sun_dlp, edge_t[:,1], edge_t[:,0],width,height,x_3)

    else:
        raise ValueError("option is wrong")
    # print(f.tolist())
    # print(row.tolist())
    # print(df)
    if df.shape[1] == 65:
        df[65] = 0
        df[66] = 0
        df[67] = 0
        df[68] = 0
        df[69] = 0
        df[70] = 0
    df.iloc[ID-1,65:71] = f
    
    df.to_csv(csv_file_path, index=False, header=False)  # index=Falseを指定すると、インデックス列が保存されません

    
    
    # plot implicit function
    x = np.arange(0,width,1)
    y = np.arange(0,height,1)
    x,y = np.meshgrid(x,y)
    z = f[0]*x**2 + f[1]*x*y + f[2]*y**2 + f[3]*x + f[4]*y + f[5]
    plt = matplotlib.pyplot.contour(z,[0])
    x = plt.collections[0].get_paths()[0].vertices[:,0]
    y = plt.collections[0].get_paths()[0].vertices[:,1]

    x2 = np.array(x, dtype=np.float64)
    y2 = np.array(y, dtype=np.float64)

    x_ans = np.arange(0,width,1)
    y_ans = np.arange(0,height,1)
    x_ans,y_ans = np.meshgrid(x_ans,y_ans)

    z_ans = f_answer[0]*(x_ans)**2 + f_answer[1]*(x_ans)*(y_ans) + f_answer[2]*(y_ans)**2 + f_answer[3]*(x_ans) + f_answer[4]*(y_ans) + f_answer[5]
    plt_ans = matplotlib.pyplot.contour(z_ans,[0])
    x_ans = plt_ans.collections[0].get_paths()[0].vertices[:,0] 
    y_ans = plt_ans.collections[0].get_paths()[0].vertices[:,1] 
    # print(x)


    for i in range(0,len(x_ans)):
        if int(y_ans[i]) >= height-1 or int(y_ans[i]) <= 0 or \
            int(x_ans[i]) >= width-1 or int(x_ans[i]) <= 0:
            continue
        original_image = insert_pixel(original_image,x_ans[i],y_ans[i],(0,255,0),line_thickness)  

    for i in range(0,len(x)):
        if int(y[i]) >= height-1 or int(y[i]) <= 0 or \
            int(x[i]) >= width-1 or int(x[i]) <= 0:
            continue
        original_image = insert_pixel(original_image,x[i],y[i],(255,0,0),line_thickness)

    # calculate error
    d_err = np.zeros(N_k)
    for i in range(0,N_k):
        y = edge[i,0]
        x = edge[i,1]
        d_err[i] = (f[0]*x**2 + f[1]*x*y + f[2]*y**2 + \
                f[3]*x + f[4]*y + f[5]) / \
                np.sqrt((2*f[0]*x + f[1]*y + f[3])**2 + \
                        (2*f[2]*y + f[1]*x + f[4])**2)
    d_std = np.std(d_err)
    print("line fit error : " + str(d_std) + " px")
    
    # draw images
    # original_image.show()
    original_image.save("../output/" + str(date) + "/images/fitted/" + str(ID) + ".png")


    y0 = 247
    x0 = 329.5


    # f[5] = f[0]*x0**2 + f[1]*x0*y0 + f[2]*y0**2 +f[3]*x0 +f[4]*y0 +f[5]
    # f[4] = f[1]*x0 + 2*f[2]*y0 + f[4]
    # f[3] = 2*f[0]*x0 + f[1]*y0 + f[3]

    # x2 = np.arange(-width/2,width/2,1)
    # y2 = np.arange(-height/2,height/2,1)
    # x2,y2 = np.meshgrid(x2,y2)
    # z2 = f[0]*(x2)**2 + f[1]*x2*y2 + f[2]*y2**2 + f[3]*x2 + f[4]*y2 + f[5]
    # plt2 = matplotlib.pyplot.contour(z2,[0])
    # x2 = plt2.collections[0].get_paths()[0].vertices[:,0]
    # y2 = plt2.collections[0].get_paths()[0].vertices[:,1]
    # print(x2)
    # print(y2)
    x2 = x2- x0
    y2 = y2 - y0
    # print(x2)
    m2 = x2.shape
    # print(m2)
    f_column = np.full((m2[0], 1),focus*494/3.66)
    # print(f_column)
    x2 = np.array(x2).reshape(-1, 1)
    y2 = np.array(y2).reshape(-1, 1)
    f_column = f_column.reshape(-1, 1)
    s = np.hstack((x2,y2, f_column))
    # s = s/(50*494/3.66)
    # print(s)

    r = christian_robinson.christian_robinson(s.T)
    print(r)
    if df.shape[1] == 71:
        df[71] = 0
        df[72] = 0
        df[73] = 0
    df.iloc[ID-1,71:74] = r

    # r = estimate_r.estimate_r(a)
    # print(r)
    # R = 1737.4
    # n = sun_dlp/np.linalg.norm(sun_dlp)
    # a1 = np.array([1,0,0])
    # c1 = a1 - np.dot(a1,n)*n
    # p = c1/np.linalg.norm(c1)
    # s = np.cross(n,p)

    # Rct = np.array([[p[0],s[0],n[0]],
    #                 [p[1],s[1],n[1]],
    #                 [p[2],s[2],n[2]]],dtype=np.float64) 

    
    # Q = np.array([[1/R**2,0,0],
    #                 [0,1/R**2,0],
    #                 [0,0,1/R**2]],dtype=np.float64) 
    
    # M = (1/np.dot(p,np.dot(Q,p)))**0.5
    # m = (1/np.dot(s,np.dot(Q,s)))**0.5

    # T = np.array([[1/M**2,0,0],
    #                 [0,1/m**2,0],
    #                 [0,0,-1]],dtype=np.float64) 
    
    # H1 = np.array([[Rct[0][0], Rct[0][1], r[0]],
    #                 [Rct[1][0], Rct[1][1], r[1]],
    #                 [Rct[2][0], Rct[2][1], r[2]]],dtype=np.float64)
    
    # # K = np.diag([f*1000/7.4,f*1000/7.4,1])
    # focus = 50
    # K = np.array([[focus*1000/7.4,0,0],
    #                 [0,focus*1000/7.4,0],
    #                 [0,0,1]],dtype=np.float64)
    # H = K@H1
    # T_dash = np.linalg.inv(H.T)@T@np.linalg.inv(H)

    # y0 = 247
    # x0 = 329.5
    # a_t = T_dash[0][0]
    # b_t = 2*T_dash[0][1]
    # c_t = T_dash[1][1]
    # d_t = 2*T_dash[0][2]
    # f_t = 2*T_dash[1][2]
    # g_t = T_dash[2][2]
    # g_t = a_t*x0**2 + b_t*x0*y0 + c_t*y0**2 -d_t*x0 -f_t*y0 +g_t
    # f_t = -b_t*x0 - 2*c_t*y0 + f_t
    # d_t= -2*a_t*x0 -b_t*y0 + d_t

    # # b = np.array([a_t,b_t,c_t,d_t,f_t,g_t]) 
    # b = np.array([a_t,b_t,c_t,d_t,f_t,g_t]) 
    # x = np.arange(0,width,1)
    # y = np.arange(0,height,1)
    # x,y = np.meshgrid(x,y)
    # z = b[0]*x**2 + b[1]*x*y + b[2]*y**2 + b[3]*x + b[4]*y + b[5]
    # plt = matplotlib.pyplot.contour(z,[0])
    # x = plt.collections[0].get_paths()[0].vertices[:,0]
    # y = plt.collections[0].get_paths()[0].vertices[:,1]
    # edge_t = np.array([x,y]).T
    # edge_r = np.zeros((np.size(x),2))
    # mat_r = np.array([[sun_dlp[1],-sun_dlp[0]],[sun_dlp[0],sun_dlp[1]]])
    # xmax = -659
    # xmin = 659
    # for i in range(np.size(x)):
    #     # print(mat_r)
    #     edge_r[i] = np.dot(edge_t[i],mat_r.T)
    #     # print(edge_r[i])
    #     if xmax < edge_r[i,1]:
    #         xmax = edge_r[i,1]
    #         # print(xmax)
    #         end = i
    #     if xmin > edge_r[i,1]:
    #         xmin = edge_r[i,1]
    #         # print(xmin)
    #         start = i
    
    # x_3 = np.array([edge_t[start][0], edge_t[start][1]])
    # print("x_3")
    # print(x_3)
    # for i in range(0,len(x)):
    #     if int(y[i]) >= height-1 or int(y[i]) <= 0 or \
    #         int(x[i]) >= width-1 or int(x[i]) <= 0:
    #         continue
    #     original_image = insert_pixel(original_image,x[i],y[i],(0,0,255),line_thickness)
    # original_image = insert_pixel(original_image,x_3[0],x_3[1],(255,255,0),3)

    # original_image.show()

    df.to_csv(csv_file_path, index=False, header=False)     
    # estimate_r.estimate_r(f)
    return r,f