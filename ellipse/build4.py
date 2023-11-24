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

date = "1123"
ID =1
focus = 50
for i in range(3800):
    # r,f = main_christian_function3.main_christian(ID,date,"fitz",0,np.zeros(6))
    # e = r2e.r2e(r)
    # r = main_christian_function3.main_christian(ID,date,"term",e,f)
    try :
        r,f = main_christian_function3.main_christian(ID,date,"slsqp",0,focus,np.zeros(6),"moon")
        e = r2e.r2e(r)
        r,f = main_christian_function3.main_christian(ID,date,"slsqp2",e,focus, f,"moon")
        
    except Exception as e:
        print(e)
    print(str(ID)+"枚目")
        
    
    ID += 1
    