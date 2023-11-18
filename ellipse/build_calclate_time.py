import main_christian_function
import main_christian_function_build
import time
import numpy as np
import r2e

date = "1118_mac"
ID =2
focus = 50
start = time.time()
# fitz
# r,f = main_christian_function_build.main_christian(ID,date,"fitz",0,focus,np.zeros(6))

# a
# r,f = main_christian_function_build.main_christian(ID,date,"slsqp",0,focus,np.zeros(6))

# b
# r,f = main_christian_function_build.main_christian(ID,date,"slsqp",0,focus,np.zeros(6))
# e = r2e.r2e(r)
# r,f = main_christian_function_build.main_christian(ID,date,"slsqp",e,focus,np.zeros(6))


end = time.time()
time_diff = end - start  # 処理完了後の時刻から処理開始前の時刻を減算する
print(time_diff) 
        
# ID += 1
    