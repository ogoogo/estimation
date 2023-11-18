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

    # image settings
    IMAGE_FILE_NAME = '../output/' +str(date)+'/images/raw/' + str(ID) +'.png'

    # number of illumination scan lines
    m = 100 # number of illumination scan lines across image

    # processing settings
    step_size = 2 # image scan step size
    window_radius = 4 # Sobel window size
    zernike_radius = 4 # Zernike moment size (NOTE: cannot be larger than window_size) 
    sigma_psf = 0.2 # pixel point spread function
    d_min = 100 # usually 5*radius_px_body/100
    N_T = 500 # maximum RANSAC tests
    line_thickness = 1 # pixel line fit thinkness

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
    edge_r_t = np.zeros((pos_t,2))
    edge_r = np.zeros((pos,2))
    xmax = -659
    xmin = 659
    mat_r = np.array([[u2[1],-u2[0]],u2])
    for i in range(pos_t):
        
        # print(mat_r)
        edge_r_t[i] = np.dot(edge_t[i],mat_r)
        # print(edge_r_t[i])
        if xmax < edge_r_t[i,1]:
            xmax = edge_r_t[i,1]
            # print(xmax)
            end_t = i
        if xmin > edge_r_t[i,1]:
            xmin = edge_r_t[i,1]
            # print(xmin)
            start_t = i

    xmax = -659
    xmin = 659
    for i in range(pos):
        
        # print(mat_r)
        edge_r[i] = np.dot(edge[i],mat_r)
        # print(edge_r[i])
        if xmax < edge_r[i,0]:
            xmax = edge_r[i,0]
            # print(xmax)
            end = i
        if xmin > edge_r[i,0]:
            xmin = edge_r[i,0]
            # print(xmin)
            start = i
    
    x_1 = np.array([edge_t[start_t,1],edge_t[start_t,0]])
    x_2 = np.array([edge_t[end_t,1],edge_t[end_t,0]])
    x_4 = np.array([edge[end,1],edge[end,0]])

    def calculate_circle(point1, point2, point3):
        x1, y1 = point1
        x2, y2 = point2
        x3, y3 = point3

        D = 2 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))

        Ux = ((x1**2 + y1**2)*(y2 - y3) + (x2**2 + y2**2)*(y3 - y1) + (x3**2 + y3**2)*(y1 - y2)) / D
        Uy = ((x1**2 + y1**2)*(x3 - x2) + (x2**2 + y2**2)*(x1 - x3) + (x3**2 + y3**2)*(x2 - x1)) / D

        r = math.sqrt((x1 - Ux)**2 + (y1 - Uy)**2)

        return (np.array([Ux, Uy]), r)
    x_c, radius = calculate_circle(x_1,x_2,x_4)

    a_i = 1
    b_i = 0
    c_i = 1
    d_i = -2*x_c[0]
    f_i = -2*x_c[1]
    g_i = x_c[0]**2 + x_c[1]**2 -radius**2
    F_ini = np.array([a_i,b_i,c_i,d_i,f_i,g_i])

    # perform fitzgibbon fit using RANSAC test
    N_test = 0
    N_min = pos
    while N_test < N_T:
        
        # increase the test count
        N_test = N_test + 1

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
    
    y0 = 247
    x0 = 329.5


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


    df.to_csv(csv_file_path, index=False, header=False)     
    # estimate_r.estimate_r(f)
    return r,f