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

date = "1102"
ID =3
for i in range(1):
    # r = main_christian_function2.main_christian(ID,date,"slsqp",0.8)
    # r2e.r2e(r)
    # print(i)
    try :

        r = main_christian_function2.main_christian(ID,date,"fitz",0)
        e = r2e.r2e(r)
        r = main_christian_function2.main_christian(ID,date,"slsqp",e)
        
    except Exception as e:
        print(e)
    print(str(ID)+"枚目")
        
    
    ID += 1
    