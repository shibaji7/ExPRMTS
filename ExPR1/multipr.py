import multiprocessing as mpr
import os
import pandas as pd
import shutil
import sza_calc as szac
import numpy as np
import datetime as dt

np.random.seed(0)
X = []
L = 101
au, al = 4.0, 0.5
if not os.path.exists("design.csv"):
    Bins = np.linspace(au,al,L)
    for i in range(L-1):
        X.append(np.random.uniform(Bins[i],Bins[i+1],1)[0])
        pass
    df = pd.DataFrame()
    df["alpha"] = X
    df.to_csv("design.csv",index=False)
    pass

riom = "ott"
RM = pd.read_csv("rio.csv")
rm = RM[RM.code==riom]
lat = rm["lat"]
lon = rm["lon"]

time = dt.datetime(2015,03,11,16,10)

sza = szac.get_sza_as_deg(time, lat, np.mod((lon+180),360)-180) 
riom = riom + "_" + str(np.round(sza,2))
dirc = "out/" + riom 
if os.path.exists(dirc): shutil.rmtree(dirc)
os.mkdir(dirc)

def run_code(l):
    cmd = "python model.py %.15f %.2f %.2f %s" % (float(l[0]),lat, lon, riom)
    print cmd
    os.system(cmd)
    return

D = pd.read_csv("design.csv")
d = D.as_matrix()

pool = mpr.Pool(10)
pool.map(run_code, d)
