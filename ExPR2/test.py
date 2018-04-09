import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
from goes_euvac import EuvacModel as EM
import goes_euvac as GE
import alpha as alp
import pandas as pd

def map_list(em,stime,etime):
    df = em.goes
    prod, T = em.exe()
    for i in range(prod.shape[1]):
        df["bin_"+str(i+1)] = prod[:,i]
        pass
    df = df.sort_values(by="dt")
    df = df[(df.dt>=stime) & (df.dt<etime)]
    df.indx = [i*60 for i in range(len(df))]
    return df

em = EM("Data/2015-03-11.csv")
stime = dt.datetime(2015,03,11,15,50)
etime = dt.datetime(2015,03,11,15,52)
etime = dt.datetime(2015,03,11,17,20)
df = map_list(em,stime,etime)
L = len(df)*60
d = (etime-stime).total_seconds()
Tp = np.array([stime+dt.timedelta(seconds=i) for i in range(int(d))])
fx_aavg = alp.get_extrp_func(df.indx,df["A_AVG"])
fx_bavg = alp.get_extrp_func(df.indx,df["B_AVG"])

fmt = matplotlib.dates.DateFormatter("%H:%M")
fig,axes = plt.subplots(figsize=(10,6),nrows=2)
fig.subplots_adjust(hspace=0.3)
I = 0
ax = axes[I]
ax.xaxis.set_major_formatter(fmt)
#ax.semilogy(df.dt,df["A_AVG"],"r")
ax.semilogy(Tp,fx_aavg(range(L)),"k-.")
#ax.semilogy(df.dt,df["B_AVG"],"b")
ax.semilogy(Tp,fx_bavg(range(L)),"r-.")
ax.set_xlim(stime,etime)


fx_bin1 = alp.get_extrp_func(df.indx,df["bin_1"])
fx_bin2 = alp.get_extrp_func(df.indx,df["bin_2"])
fx_bin3 = alp.get_extrp_func(df.indx,df["bin_3"])
fx_bin4 = alp.get_extrp_func(df.indx,df["bin_4"])
fx_bin5 = alp.get_extrp_func(df.indx,df["bin_5"])
I = 1
ax = axes[I]
ax.xaxis.set_major_formatter(fmt)
#ax.semilogy(df.dt,df["A_AVG"],"r")
ax.semilogy(Tp,fx_bin1(range(L)),"r-.")
#ax.semilogy(df.dt,df["B_AVG"],"b")
ax.semilogy(Tp,fx_bin2(range(L)),"k-.")
ax.semilogy(Tp,fx_bin3(range(L)),"b-.")
ax.semilogy(Tp,fx_bin4(range(L)),"g-.")
ax.semilogy(Tp,fx_bin5(range(L)),"m-.")
ax.set_xlim(stime,etime)

fig.savefig("Production.png")
plt.close()




### Other
riom = "ott"
RM = pd.read_csv("rio.csv")
rm = RM[RM.code==riom]
lat = rm["lat"]
lon = rm["lon"]
lon = np.mod((lon+180),360)-180
EIR = GE.estimate_delta_Ne(lat, lon, stime, etime)

fmt = matplotlib.dates.DateFormatter("%H:%M")
fig,axes = plt.subplots(figsize=(10,6),nrows=3)
fig.subplots_adjust(hspace=0.3)

I = 0
ax = axes[I]
ax.xaxis.set_major_formatter(fmt)
ax.semilogy(EIR.df_high_res.dt,EIR.df_high_res.A_AVG,"k-.")
ax.semilogy(EIR.df_high_res.dt,EIR.df_high_res.B_AVG,"r-.")
ax.set_xlim(stime,etime)


I = 1
ax = axes[I]
ax.xaxis.set_major_formatter(fmt)
ax.semilogy(EIR.out.dt,EIR.out.alt_70,"k-.")
ax.semilogy(EIR.out.dt,EIR.out.alt_100,"b-.")
ax.semilogy(EIR.out.dt,EIR.out.alt_110,"g-.")
ax.semilogy(EIR.out.dt,EIR.out.alt_120,"m-.")
ax.semilogy(EIR.out.dt,EIR.out.alt_145,"r-.")
ax.set_xlim(stime,etime)

I = 2
ax = axes[I]
ax.xaxis.set_major_formatter(fmt)
ax.plot(EIR.out.dt,EIR.out.grd_70,"k-.")
#ax.plot(EIR.out.dt,GE.plotable_logscale(EIR.out.grd_100),"k-.")
#ax.plot(EIR.out.dt,GE.plotable_logscale(EIR.out.grd_110),"k-.")
#ax.plot(EIR.out.dt,GE.plotable_logscale(EIR.out.grd_120),"k-.")
#ax.plot(EIR.out.dt,GE.plotable_logscale(EIR.out.grd_145),"k-.")
#ax.set_yscale("log")
ax.set_xlim(stime,etime)


fig.savefig("Ionization.png")
plt.close()
