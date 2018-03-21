import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd

code = "ott"
x = pd.read_csv("../Data/csv/20150311_%s_abs.csv"%code)
x.times = pd.to_datetime(x.times)

y = pd.read_csv("out/attn_0.3980456_1.595794.csv")
y.times = pd.to_datetime(y.times)

fmt = matplotlib.dates.DateFormatter("%H:%M")
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(fmt)
ax.plot(x.times,x.attn,"r-.")
ax.set_xlabel("Time [UT]")
ax.set_ylabel("Attenuation [dB]")

ax.plot(y.times,y.euvac,"ko")

ax.legend(["Observed","Modeled"])
ax.set_ylim(0,3)
ax.set_xlim(x.times.tolist()[0],x.times.tolist()[-1])

plt.show()
