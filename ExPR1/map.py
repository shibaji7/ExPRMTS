import pandas as pd
import datetime as dt
import sza_calc as szac


e = dt.datetime(2015,03,11,16,10)
D = pd.read_csv("R/out/design_sza.csv")

for I, row in D.iterrows():
    f = "out_sza/attn_%.15f_%.15f_%.15f_%.15f.csv"%(row["eta"],row["alpha"],row["lat"],row["lon"])
    fn = "out_sza_tag/attn_%.15f_%.15f_%.15f.csv"%(row["eta"],row["alpha"],row["sza"])
    d = pd.read_csv(f)
    d.to_csv(fn,header=True,index=False)
    pass


