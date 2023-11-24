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
import main_christian_function3_wearth


def vector_estimation(ID,date,focus,celestial):

    r,f = main_christian_function3_wearth.main_christian(ID,date,"slsqp",0,focus,np.zeros(6),celestial)
    e = r2e.r2e(r,celestial)
    r_c,f = main_christian_function3_wearth.main_christian(ID,date,"slsqp",e,focus,np.zeros(6),celestial)
    

    return r_c