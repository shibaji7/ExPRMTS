import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import os

path = "out/ott_51.05/"
code = "ott"
x = pd.read_csv("../Data/csv/20150311_%s_abs.csv"%code)
x.times = pd.to_datetime(x.times)

fmt = matplotlib.dates.DateFormatter("%H:%M")
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(fmt)
ax.plot(x.times,x.absv,"r-.")
for i in range(100):
    if os.path.exists(path+"/attn_%s.csv"%str(i)):
        y = pd.read_csv(path+"/attn_%s.csv"%str(i))
        y.times = pd.to_datetime(y.times)
        ax.plot(y.times,y.euvac,"ko",markersize="0.5")


ax.set_xlabel("Time [UT]")
ax.set_ylabel("Attenuation [dB]")


#ax.legend(["Observed","Modeled"])
ax.set_ylim(0,5)
ax.set_xlim(x.times.tolist()[0],x.times.tolist()[-1])
fig.savefig("out_put_%s.png"%code)
