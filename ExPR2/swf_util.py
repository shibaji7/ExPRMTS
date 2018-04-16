import os
import matplotlib
import matplotlib.pyplot as plt
import goes
import datetime as dt
import pandas as pd
import numpy as np
from pyglow import Point as P

import sza_calc as SZAC
import scale_height as SHC
#from goes_euvac import EstimateIonization as EI
#import efield

goes_ds_file_tmp = "Data/%s.csv"

font = {"family": "serif", "color":  "darkred", "weight": "normal", "size": 12}
fonttext = {"family": "serif", "color":  "darkblue", "weight": "normal", "size": 12}
fonttitle = {"family": "serif", "color":  "black", "weight": "normal", "size": 14}
fontsuptitle = {"family": "serif", "color":  "black", "weight": "bold", "size": 20}

def get_efield(lat,lon,dn):
    return None

def smooth(x,window_len=51,window="hanning"):
    if x.ndim != 1: raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len: raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3: return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': w = numpy.ones(window_len,'d')
    else: w = eval('np.'+window+'(window_len)')
    y = np.convolve(w/w.sum(),s,mode='valid')
    d = window_len - 1
    y = y[d/2:-d/2]
    return y

def get_alts():
    l,u = 70, 150
    d = 1
    n = (u-l)/d
    alts = np.linspace(l,u,n+1)
    return alts

def get_alt_bin():
    alts = np.array(get_alts())
    l,u = alts[0],alts[-1]
    d = 10
    n = (u-l)/d
    alt_bins = np.linspace(l,u,n+1)
    alt_bins = alt_bins + 5
    return alt_bins[:-1]

def summ(x):
    if type(x) is list: x = np.array(x)
    x[np.isnan(x)] = 0.
    return np.sum(x)

#def get_swf_delta_Ne_bins(ut, lat, lon, h, sco, bins=[1,2], etap=None):
#    goes_ds_file = goes_ds_file_tmp % ut.strftime("%Y-%m-%d")
#    utn = ut + dt.timedelta(days=1)
#    if not os.path.exists(goes_ds_file):
#        (d,out) = goes.read_goes(dt.datetime(ut.year,ut.month,ut.day,0,0),dt.datetime(utn.year,utn.month,utn.day,0,0))
#        out["xray"].to_csv(goes_ds_file,index_label="dt")
#    sza = SZAC.get_sza_as_deg(ut, lat, lon, h)
#    dNe = 0.
#    for bin in bins:
#        Hp, h0 = SHC.get_h0_Hp(sco.alts,sco.msis,h,bin_num=bin)
#        ei = EI(sco.alts, goes_ds_file, bin_num=bin)
#        dNe = dNe + ei.exe(ut,h,h0,Hp,sza,etap)[0]
#        pass
#    return dNe*1e6,sza

def incorporate_loss(Ne, sco, h, href, a = 7e-14):
    I = np.where(sco.alts == h)[0]
    i = np.where(sco.alts == href)[0]
    a = a * (sco.msis["O2"][I] / sco.msis["N2"][I]) * (sco.msis["N2"][i] / sco.msis["O2"][i])
    #a = (1/(2*Ne*10))
    for i in range(60):
        Ne = Ne - (a*Ne**2)
    return Ne

def redistribute_alpha(alphas,alts):
    alt_bins = get_alt_bin()
    a = np.zeros(len(alts))
    for K,h in enumerate(alts):
        I = np.argmin(np.abs(alt_bins-h))
        a[K] = alphas[I]
        pass
    return a

def estimate_Ne(EIR,geo,alphas):
    Ne = geo.Ne
    dNe = EIR.height_time_dne_df.as_matrix()[:,1:]
    Ne = geo.iri["Ne"] + dNe
    a = redistribute_alpha(alphas,geo.alts)
    for K in range(len(Ne)):
        for i in range(10):
            Ne[K,:] = Ne[K,:] - (a*Ne[K,:]**2)
            pass
        pass
    geo.Ne = Ne
    for K in range(Ne.shape[1]):
        geo.rate[:,K] = np.gradient(geo.Ne[:,K])
    return geo

def incorporate_loss_mod(Ne, a):
    for i in range(60): Ne = Ne - (a*Ne**2)
    return Ne

def incorporate_loss_ini(dNe, No, a):
    for i in range(60): No = No - (a*No**2)
    Ne = dNe+No
    return Ne

class GEO(object):
    def __init__(self, stime, etime, lat, lon, seconds=10):
        self.alts = get_alts()
        self.stime = stime
        self.etime = etime
        self.lat = lat
        self.lon = lon
        self.sec = seconds
        d = int((etime-stime).total_seconds()/seconds)
        self.dn = [stime+dt.timedelta(seconds=self.sec*i) for i in range(d)]
        self.iri = {
                "Ne":np.zeros((len(self.dn),len(self.alts))),
                "Ni":np.zeros((len(self.dn),len(self.alts))),
                "Te":np.zeros((len(self.dn),len(self.alts)))
                }
        self.msis = {
                "Tn":np.zeros((len(self.dn),len(self.alts))),
                "rho":np.zeros((len(self.dn),len(self.alts))),
                "AR":np.zeros((len(self.dn),len(self.alts))),
                "H":np.zeros((len(self.dn),len(self.alts))),
                "HE":np.zeros((len(self.dn),len(self.alts))),
                "N":np.zeros((len(self.dn),len(self.alts))),
                "N2":np.zeros((len(self.dn),len(self.alts))),
                "O":np.zeros((len(self.dn),len(self.alts))),
                "O2":np.zeros((len(self.dn),len(self.alts))),
                "O_anomalous":np.zeros((len(self.dn),len(self.alts))),
                "nn": np.zeros((len(self.dn),len(self.alts)))
                }
        self.nu = np.zeros((len(self.dn),len(self.alts)))
        self.Ne = np.zeros((len(self.dn),len(self.alts)))
        self.rate = np.zeros((len(self.dn),len(self.alts)))
        self.beta = np.zeros((len(self.dn),len(self.alts)))
        for I,u in enumerate(self.dn):
            for J,h in enumerate(self.alts):
                p = P(u,lat,lon,h)
                p.run_iri()
                p.run_msis()
                self.iri["Ne"][I,J] = p.ne
                self.iri["Ni"][I,J] = summ(p.ni.values())
                self.iri["Te"][I,J] = p.Te
                self.msis["Tn"][I,J] = p.Tn_msis
                self.msis["rho"][I,J] = p.rho
                self.msis["AR"][I,J] = p.nn["AR"]
                self.msis["H"][I,J] = p.nn["H"]
                self.msis["HE"][I,J] = p.nn["HE"]
                self.msis["N2"][I,J] = p.nn["N2"]
                self.msis["O"][I,J] = p.nn["O"]
                self.msis["O2"][I,J] = p.nn["O2"]
                self.msis["O_anomalous"][I,J] = p.nn["O_anomalous"]
                self.msis["nn"][I,J] = (p.nn["AR"] + p.nn["H"] + p.nn["HE"] + p.nn["N"] + p.nn["N2"] + 
                        p.nn["O"] + p.nn["O2"] + p.nn["O_anomalous"])
                self.nu[I,J] = Collision.friedrich_torkar_electron_neutral_collision_frequency(self.msis["nn"][I,J], p.Te, p.Tn_msis)
                pass
            pass
        return

class SCO(object):
    def __init__(self):
        self.alts = get_alts()
        self.msis = {
                "Tn":np.zeros(self.alts.shape),
                "rho":np.zeros(self.alts.shape),
                "AR":np.zeros(self.alts.shape),
                "H":np.zeros(self.alts.shape),
                "HE":np.zeros(self.alts.shape),
                "N":np.zeros(self.alts.shape),
                "N2":np.zeros(self.alts.shape),
                "O":np.zeros(self.alts.shape),
                "O2":np.zeros(self.alts.shape),
                "O_anomalous":np.zeros(self.alts.shape),
                "nn": np.zeros(self.alts.shape)
                }
        return

class Collision(object):
    
    @staticmethod
    def schunk_nagy_electron_ion_collision_frequency(Ne,Ni,Te,Ti,gamma=0.5572):
        e = 4.8e-10
        eps0 = 1
        k = 1.38e-16
        
        ki2 = (e**2 / (k * (4 * np.pi * eps0))) * (Ni / Ti)
        ke2 = (e**2 / (k * (4 * np.pi * eps0))) * (Ne / Te)
        ke = np.sqrt(ke2)
        kG = np.log((4*k)*(4*np.pi*eps0)*(Te/ke)/(gamma**2*e**2))-(((ke2+ki2)/ki2)*np.log(np.sqrt((ke2+ki2)/ke2)))
        nu = 3.63e-6 * Ni * Te**(3/2) * kG
        return nu

    
    @staticmethod
    def atmospheric_collision_frequency(ni, nn, T):
        na_profile = lambda T,nn: (1.8*1e-8*nn*np.sqrt(T/300))
        ni_profile = lambda T,ni: (6.1*1e-3*ni*(300/T)*np.sqrt(300/T))
        nu = (ni_profile(T,ni) + na_profile(T,nn))
        return nu
    
    @staticmethod
    def schunk_nagy_electron_nutral_collision_frequency(nn, Te):
        nu_N2 = 2.33e-11 * nn["N2"] * (1 - (1.12e-4 * Te)) * Te
        nu_O2 = 1.82e-10 * nn["O2"] * (1 + (3.6e-2 * np.sqrt(Te))) * np.sqrt(Te)
        nu_O  = 8.9e-11 * nn["O"] * (1 + (5.7e-4 * Te)) * np.sqrt(Te)
        nu_He = 4.6e-10 * nn["HE"] * np.sqrt(Te)
        nu_H  = 4.5e-9 * nn["H"] * (1-(1.35e-4 * Te)) * np.sqrt(Te)
        nu = nu_N2 + nu_O2 + nu_O + nu_He + nu_H
        return nu
    
    @staticmethod
    def friedrich_torkar_electron_neutral_collision_frequency(nn, Te, Tn, avg = False):
        k = 1.38e-23
        p = nn*1e6*k*Tn
        nu = ((2.637e6/np.sqrt(Te)) + 4.945e5) * p
        if avg: nu = (3*nu/2)
        return nu

    @staticmethod
    def atmospheric_ion_neutral_collision_frequency(Nn):
        nu_i = 3.8e-11*Nn
        return nu_i

class Absorption_AH(object):
    
    @staticmethod
    def no_mag(N,f,nu):
        e = 1.6e-19
        m = 9.1e-31
        eps0 = 1e-9/(36*np.pi)
        c = 3e8
        k = (2 * np.pi * f) / c
        omega = 2 * np.pi * f
        X = (N * e**2) / (eps0 * m * omega**2)
        Z = nu / omega
        size = Z.shape
        beta = np.zeros(size)
        for i in range(size[0]):
            for j in range(size[1]):
                n2 = 1 - (X[i,j]/np.complex(1,Z[i,j]))
                n = np.sqrt(n2)
                beta[i,j] = 8.68 * k * 1e3 * n.imag
                pass
            pass
        return beta

    @staticmethod
    def with_no_mag_ql(N,f,nu):
        e = 1.6e-19
        m = 9.1e-31
        eps0 = 1e-9/(36*np.pi)
        c = 3e8
        k = (2 * np.pi * f) / c
        omega = 2 * np.pi * f
        X = (N * e**2) / (eps0 * m * omega**2)
        Z = nu / omega
        YL, YT = (e * B) / (m * omega), 0
        n2L,n2R = 1 - (X/(np.complex(1,Z)+YL)), 1 - (X/(np.complex(1,Z)-YL))
        nL,nR = np.sqrt(n2L),np.sqrt(n2R)
        R,L = 8.68 * k * 1e3 * nR.imag, 8.68 * k * 1e3 * nL.imag
        return R,L

    @staticmethod
    def with_no_mag_qt(N,f,nu):
        e = 1.6e-19
        m = 9.1e-31
        eps0 = 1e-9/(36*np.pi)
        c = 3e8
        k = (2 * np.pi * f) / c
        omega = 2 * np.pi * f
        X = (N * e**2) / (eps0 * m * omega**2)
        Z = nu / omega
        YL, YT = 0, (e * B) / (m * omega)
        n2O,n2X = 1 - (X/np.complex(1,Z)), 1 - ((2*X*np.complex(1-X,Z))/((2*np.complex(1-X,Z)*np.complex(1,Z))-(2*YT**2)))
        nO,nX = np.sqrt(n2O),np.sqrt(n2X)
        O,X = 8.68 * k * 1e3 * nO.imag, 8.68 * k * 1e3 * nX.imag
        return O,X



def plot_summary(stime,etime,geo,EIR,fname,code,alts=[80,83,86],colors=["r-.","b-.","k-."]):
    fmt = matplotlib.dates.DateFormatter("%H:%M")
    fig,axes = plt.subplots(figsize=(15,9),nrows=4)
    fig.subplots_adjust(hspace=0.5)
    
    I = 0
    ax = axes[I]
    ax.xaxis.set_major_formatter(fmt)
    ax.semilogy(EIR.df_high_res.dt,EIR.df_high_res.A_AVG,"k-.")
    ax.semilogy(EIR.df_high_res.dt,EIR.df_high_res.B_AVG,"r-.")
    ax.set_xlim(stime,etime)
    ax.set_xlabel(r"$Time [UT]$",fontdict=font)
    ax.set_ylabel(r"$\Phi_0 [Wm^{-2}]$",fontdict=font)

    I = 1
    ax = axes[I]
    ax.xaxis.set_major_formatter(fmt)
    for k,alt in enumerate(alts):
        ax.semilogy(EIR.out.dt,EIR.out["alt_"+str(alt)],colors[k])
    ax.set_xlim(stime,etime)
    ax.set_xlabel(r"$Time [UT]$",fontdict=font)
    ax.set_ylabel(r"$N_e [m^{-3}]$",fontdict=font)
    
    I = 2
    ax = axes[I]
    ax.xaxis.set_major_formatter(fmt)
    for k,alt in enumerate(alts):
        ax.plot(EIR.out.dt,smooth(geo.rate[:,k]),colors[k])
    ax.set_xlim(stime,etime)
    ax.set_xlabel(r"$Time [UT]$",fontdict=font)
    ax.set_ylabel(r"$\frac{dN_e}{dt} [m^{-3}sec^{-1}]$",fontdict=font)

    I=3
    ax = axes[I]
    ax.xaxis.set_major_formatter(fmt)
    ax.plot(EIR.out.dt,np.sum(geo.beta,axis=1),"ko",markersize=0.5)
    ax.set_xlim(stime,etime)
    ax.set_xlabel(r"$Time [UT]$",fontdict=font)
    ax.set_ylabel(r"$\beta_L [dB]$",fontdict=font)
    
    x = pd.read_csv("../Data/csv/20150311_%s_abs.csv"%code)
    x.times = pd.to_datetime(x.times)
    ax.plot(x.times,x.absv,"r-.")
    fig.autofmt_xdate()
    fig.savefig(fname)
    pass
