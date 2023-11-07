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

date = "1107"
ID =330
for i in range(400):
    # print(i)
    main_christian_function.main_christian(ID,date)
    print(str(ID)+"枚目")
    ID += 1
    