import main_christian_function
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
import r2e
import main_christian_function2
import main_christian_function3
import main_christian_function3_zernike
import vector_estimation
import closest_point
import r_c2r_w


date = "1121"
ID =1
focus = 50
for i in range(800):
    # r_c_e = vector_estimation.vector_estimation(ID,date,focus,"earth")
    # r_c_m = vector_estimation.vector_estimation(ID,date,focus,"moon")
    
    # line1 = r_c2r_w.r_c2r_w(ID,date,r_c_e,"earth")
    # print(line1)
    # line2 = r_c2r_w.r_c2r_w(ID,date,r_c_m,"moon")

    # r_1 = closest_point.find_point_with_min_distance(line1,line2)
    # print(r_1)
    
    # csv_file_path = "../output/" + str(date) + "/information.csv"

    # # CSVファイルを読み込む
    # df = pd.read_csv(csv_file_path, header=None)

    # # 6行目のデータを取得
    # row = df.iloc[ID-1]  # インデックスは0から始まるため
    # if df.shape[1] == 95:
    #     df[95] = 0
    #     df[96] = 0
    #     df[97] = 0
    # df.iloc[ID-1,95:98] = r_c_m
    # if df.shape[1] == 98:
    #     df[98] = 0
    #     df[99] = 0
    #     df[100] = 0
    # df.iloc[ID-1,98:101] = r_1
    try :
        r_c_e = vector_estimation.vector_estimation(ID,date,focus,"earth")
        r_c_m = vector_estimation.vector_estimation(ID,date,focus,"moon")
        
        line1 = r_c2r_w.r_c2r_w(ID,date,r_c_e,"earth")
        print(line1)
        line2 = r_c2r_w.r_c2r_w(ID,date,r_c_m,"moon")

        r_1 = closest_point.find_point_with_min_distance(line1,line2)
        print(r_1)
        
        csv_file_path = "../output/" + str(date) + "/information.csv"

        # CSVファイルを読み込む
        df = pd.read_csv(csv_file_path, header=None)

        # 6行目のデータを取得
        row = df.iloc[ID-1]  # インデックスは0から始まるため
        if df.shape[1] == 95:
            df[95] = 0
            df[96] = 0
            df[97] = 0
        df.iloc[ID-1,95:98] = r_c_m
        if df.shape[1] == 98:
            df[98] = 0
            df[99] = 0
            df[100] = 0
        df.iloc[ID-1,98:101] = r_1
        if df.shape[1] == 101:
            df[101] = 0
            df[102] = 0
            df[103] = 0
        df.iloc[ID-1,101:104] = r_c_e
        df.to_csv(csv_file_path, index=False, header=False)    
    except Exception as e:
        print(e)
    print(str(ID)+"枚目")
        
    
    ID += 1
    