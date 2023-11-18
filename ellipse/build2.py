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
import main_christian_function_shou

date = "1102"
ID =1
focus = 50
for i in range(1):
    # F = np.array([0.1,0.1,1,1,1,1])
    # r,f = main_christian_function3.main_christian(ID,date,"slsqp",0,np.zeros(6))
    r,f = main_christian_function3.main_christian(ID,date,"fitz",0,focus,np.zeros(6))
    # print(i)
    # try :
        # F = np.array([0.1,0.1,1,1,1,1])
        # r,f = main_christian_function3.main_christian(ID,date,"fitz",0,focus,np.zeros(6))
        # r2e.r2e(r)
        
    # except Exception as e:
    #     print(e)
    print(str(ID)+"枚目")
        
    
    ID += 1
    