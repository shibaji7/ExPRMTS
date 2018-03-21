import os
import goes
import datetime as dt
import pandas as pd
import numpy as np
from pyglow import Point as P

import sza_calc as SZAC
import scale_height as SHC
from goes_euvac import EstimateIonization as EI
#import efield

goes_ds_file_tmp = "Data/%s.csv"

font = {"family": "serif", "color":  "darkred", "weight": "normal", "size": 12}
fonttext = {"family": "serif", "color":  "darkblue", "weight": "normal", "size": 12}
fonttitle = {"family": "serif", "color":  "black", "weight": "normal", "size": 14}
fontsuptitle = {"family": "serif", "color":  "black", "weight": "bold", "size": 20}

def get_efield(lat,lon,dn):
    return None

def get_alts():
    l,u = 50, 150
    n = (u - l)*1
    alts = np.round(np.linspace(l+1,u,n),1)
    return alts

def summ(x):
    if type(x) is list: x = np.array(x)
    x[np.isnan(x)] = 0.
    return np.sum(x)

def get_swf_delta_Ne_bins(ut, lat, lon, h, sco, bins=[1,2], etap=None):
    goes_ds_file = goes_ds_file_tmp % ut.strftime("%Y-%m-%d")
    utn = ut + dt.timedelta(days=1)
    if not os.path.exists(goes_ds_file):
        (d,out) = goes.read_goes(dt.datetime(ut.year,ut.month,ut.day,0,0),dt.datetime(utn.year,utn.month,utn.day,0,0))
        out["xray"].to_csv(goes_ds_file,index_label="dt")
    sza = SZAC.get_sza_as_deg(ut, lat, lon, h)
    dNe = 0.
    for bin in bins:
        Hp, h0 = SHC.get_h0_Hp(sco.alts,sco.msis,h,bin_num=bin)
        ei = EI(sco.alts, goes_ds_file, bin_num=bin)
        dNe = dNe + ei.exe(ut,h,h0,Hp,sza,etap)[0]
        pass
    return dNe*1e6,sza

def incorporate_loss(Ne, sco, h, href, a = 7e-14):
    I = np.where(sco.alts == h)[0]
    i = np.where(sco.alts == href)[0]
    a = a * (sco.msis["O2"][I] / sco.msis["N2"][I]) * (sco.msis["N2"][i] / sco.msis["O2"][i])
    #a = (1/(2*Ne*10))
    for i in range(60):
        Ne = Ne - (a*Ne**2)
    return Ne

class SCO(object):
    def __init__(self, dn, lat, lon):
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
        for I,alt in enumerate(self.alts):
            pt = P(dn, lat, lon, alt)
            pt.run_msis()
            self.msis["Tn"][I] = pt.Tn_msis
            self.msis["rho"][I] = pt.rho
            self.msis["AR"][I] = pt.nn["AR"]
            self.msis["H"][I] = pt.nn["H"]
            self.msis["HE"][I] = pt.nn["HE"]
            self.msis["N"][I] = pt.nn["N"]
            self.msis["N2"][I] = pt.nn["N2"]
            self.msis["O"][I] = pt.nn["O"]
            self.msis["O2"][I] = pt.nn["O2"]
            self.msis["O_anomalous"][I] = pt.nn["O_anomalous"]
            self.msis["nn"][I] = (pt.nn["AR"] + pt.nn["H"] + pt.nn["HE"] + pt.nn["N"] + pt.nn["N2"] + pt.nn["O"] + pt.nn["O2"] + pt.nn["O_anomalous"])
            pass
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
        n2 = 1 - (X/np.complex(1,Z))
        n = np.sqrt(n2)
        beta = 8.68 * k * 1e3 * n.imag
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
