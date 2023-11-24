import numpy as np
import pandas as pd

def r_c2r_w(ID,date,r_c,celestial):
    # CSVファイルのパスを指定
    csv_file_path = "../output/" + str(date) + "/information.csv"

    # CSVファイルを読み込む
    df = pd.read_csv(csv_file_path, header=None)

    # 6行目のデータを取得
    row = df.iloc[ID-1] 
    
    if celestial == "moon":
        dcm = np.array([[row[56],row[57],row[58]],[row[59],row[60],row[61]],[row[62],row[63],row[64]]])
        l_cele = np.array(row[5:8])
    else:
        dcm = np.array([[row[74],row[75],row[76]],[row[77],row[78],row[79]],[row[80],row[81],row[82]]])
        l_cele = np.zeros(3)

    # r_c = np.array([1,1,1])
    # dcm = np.array([[np.cos(np.pi/6),np.sin(np.pi/6),0],[-np.sin(np.pi/6),np.cos(np.pi/6),0],[0,0,-1]])
    r_w = np.dot(dcm.T,r_c)
    print(l_cele - r_w)
    return l_cele, -r_w


# if __name__ == "__main__":
#     r_w = r_c2r_w(1,"1121",np.array([0,0,0]),"moon")
#     print(r_w)