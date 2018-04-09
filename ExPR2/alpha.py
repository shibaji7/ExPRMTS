import datetime as dt
import numpy as np
from pyglow import Point as P
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import os

import swf_util as SU
from scipy.interpolate import interp1d
from scipy import arange, array, exp

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y
    
    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)
    
    def ufunclike(xs):
        return array(map(pointwise, array(xs)))
    
    return ufunclike

def get_extrp_func(x,y):
    f_i = interp1d(x, y)
    f_x = extrap1d(f_i)
    return f_x

import swf_util as SU


Species = ["N2+","O2+","NO+","H+(H2On)"]
Ion_Chemistry = {"N2+":{"C":1.8e-13,"D":-0.39},"O2+":{"C":1.6e-13,"D":-0.55},"NO+":{"C":4.5e-13,"D":-0.83}}


def estimate_alphaH_prior(dn=dt.datetime(2015,03,11,16,10),lat=37.0,lon=-88.0):
    alpha = []
    H = np.linspace(80,150,71)
    for h in H:
        a = 0.0
        p = P(dn,lat,lon,h)
        p.run_iri()
        Ne = p.ne * 1e6
        Tn,Te = p.Tn_msis,p.Te
        for s in Species:
            if s in p.ni.keys():
                ai = Ion_Chemistry[s]["C"] * (Te/300)**(Ion_Chemistry[s]["D"])
                a = a + (ai*p.ni[s]*1e6/Ne)
            else: pass
            pass
        alpha.append(a)
        #print h,a
        #break
        pass
    f_i = interp1d(H, alpha)
    f_x = extrap1d(f_i)
    h = np.linspace(70,79,10)
    alpha_p = f_x(h)
    plt.semilogx(alpha,H,"r-.")
    #print alpha_p
    plt.semilogx(alpha_p,h,"k-.")
    plt.savefig("alpha_daily.png")
    return np.append(alpha_p,alpha),np.append(h,H),f_x


def design_alpha(dn,lat,lon,code):
    alt_bins = SU.get_alt_bin()
    alpha_m,h,f_x = estimate_alphaH_prior(dn,lat,lon)
    alpha = f_x(alt_bins)
    hL = len(alt_bins)
    np.random.seed(0)
    X = []
    L = 101
    au, al = 0.25, 0.01
    if not os.path.exists("design_%s.csv"%code):
        Bins = np.linspace(au,al,L)
        print Bins
        for i in range(L-1):
            x = []
            for k in range(hL):
                #I = np.random.randint(0, high=100, size=1)
                x.append(np.random.uniform(Bins[i],Bins[i+1],1)[0])
                pass
            x = np.array(x)
            X.append(x*alpha)
            pass
        X = np.array(X)
        df = pd.DataFrame()
        for J in range(L-1): df["alpha_run_"+str(J)] = X[J,:]
        df.to_csv("design_%s.csv"%(code),index=False)
        pass
    return

if __name__ == "__main__":
    estimate_alphaH_prior()

    riom = "ott"
    RM = pd.read_csv("rio.csv")
    rm = RM[RM.code==riom]
    lat = rm["lat"]
    lon = rm["lon"]
    time = dt.datetime(2015,03,11,16,10)
    design_alpha(time,lat,np.mod((lon+180),360)-180,riom)
    pass
