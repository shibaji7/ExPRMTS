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
import time

import swf_util as SU
from swf_util import SCO
from swf_util import GEO
from swf_util import Collision
from swf_util import Absorption_AH as ABS
import goes_euvac as GE

def get_sco(geo, I):
    sco = SCO()
    sco.msis["Tn"] = geo.msis["Tn"][I,:]
    sco.msis["rho"] = geo.msis["rho"][I,:]
    sco.msis["AR"] = geo.msis["AR"][I,:]
    sco.msis["H"] = geo.msis["H"][I,:]
    sco.msis["HE"] = geo.msis["HE"][I,:]
    sco.msis["N"] = geo.msis["N"][I,:]
    sco.msis["N2"] = geo.msis["N2"][I,:]
    sco.msis["O"] = geo.msis["O"][I,:]
    sco.msis["O2"] = geo.msis["O2"][I,:]
    sco.msis["O_anomalous"] = geo.msis["O_anomalous"][I,:]
    sco.msis["nn"] = geo.msis["nn"][I,:]
    return sco

def __non_magnetic_attn_model_1P(geo,stime,etime,lat,lon,freq,code,alphas,fname):
    h_var_absorption_iri = np.zeros((len(geo.alts)))
    h_var_absorption_goes_euvac = np.zeros((len(geo.alts)))
    
    EIR = GE.estimate_delta_Ne(lat, lon, stime, etime, geo)
    geo = SU.estimate_Ne(EIR,geo,alphas)
    geo.beta = ABS.no_mag(geo.Ne,freq,geo.nu)
    SU.plot_summary(stime,etime,geo,EIR,"Production.png",code)
    print alphas
    #for bin in bins:
    #    for I,h in enumerate(alts):
    #        indx = np.argmin(abs(alt_bins - h))
    #        p = P(u,lat,lon,h)
    #        p.run_iri()
    #        p.run_msis()
    #        Ne = p.ne * 1e6
    #        Ni = SU.summ(p.ni.values())
    #        Nn = (p.nn["AR"] + p.nn["H"] + p.nn["HE"] + p.nn["N"] + p.nn["N2"] + p.nn["O"] + p.nn["O2"] + p.nn["O_anomalous"])
    #        Tn,Te = p.Tn_msis,p.Te
    #        nu = Collision.friedrich_torkar_electron_neutral_collision_frequency(Nn, Te, Tn)
    #        h_var_absorption_iri[I] = h_var_absorption_iri[I] + ABS.no_mag(Ne,freq,nu)
    #        dNe,sza= SU.get_swf_delta_Ne_bins(ut, lat, lon, h, sco,bins=[int(bin)],etap=etap)
    #        Ne = SU.incorporate_loss_ini(dNe,Ne,alphas[indx])
    #        #Ne = SU.incorporate_loss_mod(Ne+dNe,alphas[indx])
    #        h_var_absorption_goes_euvac[I] = h_var_absorption_goes_euvac[I] + ABS.no_mag(Ne,freq,nu)
    #        pass
    #    pass

    #beta_L_iri = np.round(SU.summ(h_var_absorption_iri),2)
    #beta_L_euvac = np.round(SU.summ(h_var_absorption_goes_euvac),2)
    
    #df = pd.DataFrame([[ut,beta_L_iri,beta_L_euvac]],columns=["times","iri","euvac"])
    #f = fname
    #if J==0: 
    #    if os.path.exists(f): os.remove(f)
    #    df.to_csv(f,index=False,header=True,mode="w")
    #else: df.to_csv(f,index=False,header=False,mode="a")
    return

def __abs_model(args):
    with open("conf.json","r") as code: conf = json.load(code)
    conf["dn"] = dt.datetime.strptime(conf["dn"],"%Y-%m-%d %H:%M:%S")
    code = conf["location_name"]
    freq = float(conf["frequency"])
    lat = float(conf["lat"])
    lon = float(conf["lon"])
    etap = 1.0
    bins = conf["bins"]
    lon = np.mod((lon+180),360)-180
    
    i = int(args[0])
    lat = float(args[1])
    lon = float(args[2])
    lon = np.mod((lon+180),360)-180

    d = pd.read_csv("design_%s.csv"%args[4])
    alpha_i = np.array(d["alpha_run_"+str(i)].tolist())
    
    fname = "out/"+args[3]+"/attn_%d.csv"%(i)
    stime = conf["dn"] - dt.timedelta(minutes=10)
    etime = stime + dt.timedelta(minutes=90)
    geo = GEO(stime,etime,lat,lon)
    __non_magnetic_attn_model_1P(geo,stime,etime,lat,lon,freq,args[4],alpha_i,fname)
    return


if __name__ == "__main__":
    args = sys.argv[1:]
    start_time = time.time()
    __abs_model(args)
    print("--- Execution Time - %s seconds ---" % (time.time() - start_time))
    pass
