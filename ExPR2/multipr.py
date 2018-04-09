import multiprocessing as mpr
import os
import pandas as pd
import shutil
import sza_calc as szac
import numpy as np
import datetime as dt


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

def run_code(i):
    cmd = "python model.py %d %.2f %.2f %s %s" % (i ,lat, lon, riom, rm["code"].tolist()[0])
    print cmd
    #os.system(cmd)
    return

D = range(100)

pool = mpr.Pool(10)
pool.map(run_code, D)
