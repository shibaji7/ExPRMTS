import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from mpl_toolkits.basemap import Basemap
import sys
import numpy as np
from pyglow import Point as P
import datetime as dt
import pandas as pd
import json
import traceback

import swf_util as SU
from swf_util import SCO
from swf_util import Collision
from swf_util import Absorption_AH as ABS


def __non_magnetic_attn_model_1P(ut,lat,lon,freq,code,etap,alphas,hrefs,bins,J,fname):
    alts = SU.get_alts()
    h_var_absorption_iri = np.zeros((len(alts)))
    h_var_absorption_goes_euvac = np.zeros((len(alts)))
    
    u = ut
    sco = SCO(u, lat, lon)
    for bin in bins:
        for I,h in enumerate(alts):
            p = P(u,lat,lon,h)
            p.run_iri()
            p.run_msis()
            Ne = p.ne * 1e6
            Ni = SU.summ(p.ni.values())
            Nn = (p.nn["AR"] + p.nn["H"] + p.nn["HE"] + p.nn["N"] + p.nn["N2"] + p.nn["O"] + p.nn["O2"] + p.nn["O_anomalous"])
            Tn,Te = p.Tn_msis,p.Te
            nu = Collision.friedrich_torkar_electron_neutral_collision_frequency(Nn, Te, Tn)
            h_var_absorption_iri[I] = h_var_absorption_iri[I] + ABS.no_mag(Ne,freq,nu)
            dNe,sza= SU.get_swf_delta_Ne_bins(ut, lat, lon, h, sco,bins=[int(bin)],etap=etap)
            Ne = SU.incorporate_loss(Ne+dNe,sco,h,hrefs[bin],a=alphas[bin])
            h_var_absorption_goes_euvac[I] = h_var_absorption_goes_euvac[I] + ABS.no_mag(Ne,freq,nu)
            pass
        pass

    beta_L_iri = np.round(SU.summ(h_var_absorption_iri),2)
    beta_L_euvac = np.round(SU.summ(h_var_absorption_goes_euvac),2)
    
    df = pd.DataFrame([[ut,beta_L_iri,beta_L_euvac]],columns=["times","iri","euvac"])
    print "end-",df
    #f = "out/attn_%s_%s.csv"%(code,ut.strftime("%Y-%m-%d"))
    f = fname
    print J
    if J==0: 
        if os.path.exists(f): os.remove(f)
        df.to_csv(f,index=False,header=True,mode="w")
    else: df.to_csv(f,index=False,header=False,mode="a")
    return beta_L_iri, beta_L_euvac

modify = True
def __abs_model(args):
    with open("conf.json","r") as code: conf = json.load(code)
    conf["dn"] = dt.datetime.strptime(conf["dn"],"%Y-%m-%d %H:%M:%S")
    jobs = conf["jobs"]
    code = conf["location_name"]
    freq = float(conf["frequency"])
    lat = float(conf["lat"])
    lon = float(conf["lon"])
    etap = conf["etap"]
    alphas = conf["alphas"]
    hrefs = conf["hrefs"]
    bins = conf["bins"]
    lon = np.mod((lon+180),360)-180
    fname = "out/attn_%s_%s.csv"%(args[0],args[1])
    
    if modify:
        etap = 1.
        alphas["0"] = float(args[0]) * alphas["0"]
        alphas["1"] = float(args[0]) * alphas["1"]
        alphas["2"] = float(args[0]) * alphas["2"]
        lat = float(args[1])
        lon = float(args[2])
        lon = np.mod((lon+180),360)-180
    	fname = "out/"+args[3]+"/attn_%s.csv"%(args[0])
        print fname
        pass
    print etap
    dn = conf["dn"] - dt.timedelta(minutes=10)
    #dn = dt.datetime(2015,03,11,16,17)
    dates = [dn + dt.timedelta(minutes=x) for x in range(90)]
    for I,ut in enumerate(dates):
        __non_magnetic_attn_model_1P(ut,lat,lon,freq,code,etap,alphas,hrefs,bins,I,fname)
        pass
    return


if __name__ == "__main__":
    args = sys.argv[1:]
    __abs_model(args)
    pass
