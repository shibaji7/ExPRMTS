import pandas as pd
import datetime as dt
import numpy as np
import json

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import pyglow
import scale_height as SHC

import alpha as alp
import swf_util as SU
from swf_util import SCO
import sza_calc as SZAC


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

class EuvacModel(object):

    def __init__(self, goes_data_file, goes_tags=["A_AVG","B_AVG"]):
        self.goes_data_file = goes_data_file
        self.goes = pd.read_csv(goes_data_file)
        self.goes.dt = pd.to_datetime(self.goes.dt)
        self.goes_tags = goes_tags
        with open("euvac.json","r") as code: self.euvac_config = json.load(code)
        self.lam_l = (1e-9 * np.array(self.euvac_config["lambda_l"]))
        self.lam_u = (1e-9 * np.array(self.euvac_config["lambda_u"]))
        self.lam = ((self.lam_l + self.lam_u) / 2)
        return

    def exe(self):
        x = []
        for J, goes_tag in enumerate(self.goes_tags):
            flux = self.goes[goes_tag]
#            flux = flux[(flux.index >= self.start_time) & (flux.index <= self.end_time)]
            flux = flux.groupby(flux.index).first()
            m_value = np.mean(flux[-100:])
            flux[np.isnan(flux)] = m_value
            h,c = 6.626e-34,3e8
            phi_0 = flux / (h*c/self.lam[J])
            x.append(phi_0.tolist())
            pass
        X = np.transpose(np.array(x))
        E = self.__create_norm_euvac_model(X)
        XE = np.concatenate((X, E), axis=1)
        T = np.array([t.to_pydatetime() for t in self.goes.dt])
        #production = XE[:,self.bin_num]
        production = XE
        return production, T

    def __create_norm_euvac_model(self,xx):
        start_time = dt.datetime(1980,1,1)
        end_time = dt.datetime(2017,8,31)
        days = (end_time - start_time).days
        xe = np.zeros(shape=(days,37))
        for i, day in enumerate(range(days)):
            dn = start_time + dt.timedelta(days=day)
            kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae = pyglow.get_kpap.get_kpap(dn)
            solspec = self.__euvac_model(f107,f107a)
            xe[i,:] = solspec["current"][:]
            pass
        cols = range(37)
        XE = pd.DataFrame(xe,columns=cols)
        XE = XE.dropna(axis=0,how="any")
        X = XE.as_matrix(cols[:2])
        E = XE.as_matrix(cols[2:])
        ee = []
        for m in range(35):
            regr = linear_model.LinearRegression()
            regr.fit(X,E[:,m])
            ee.append(regr.predict(xx).tolist())
            pass
        ee = np.transpose(np.array(ee))
        return ee

    def __euvac_model(self,f107,f107a):
        _n = 37
        solspec= {
                "n_wave" : _n,
                "current" : np.empty(_n),
                "reference" : np.array([3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,\
                        5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,\
                        2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,\
                        8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,\
                        4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,\
                        3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,\
                        1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,\
                        8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,\
                        5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,\
                        5.010e+01]),
                "afac" : np.array([5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,\
                        1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,\
                        2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,\
                        2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,\
                        5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,\
                        4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,\
                        3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,\
                        7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,\
                        1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,\
                        6.240e-01])
                }
        P = (f107 + f107a)/2.
        if P < 80.: solspec["current"] = 0.8*solspec["reference"]*1e4
        else: solspec["current"] = solspec["reference"]*(1. + solspec["afac"]*(P-80.))*1e4                                                                          
        solspec["current"] = 0.5*solspec["current"]
        return solspec



class EstimateIonizationRate(object):

    def __init__(self, alts, goes_data_file, stime, etime, geo, bins=[1,2]):
        self.alts = alts
        self.goes_data_file = goes_data_file
        self.model = EuvacModel(goes_data_file)
        self.stime = stime
        self.etime = etime
        self.__map_list()
        self.geo = geo
        self.bins = bins
        with open("euvac.json","r") as code: self.euvac_config = json.load(code)
        self.height_time_dne_df = pd.DataFrame()
        self.height_time_dne_df["dt"] = self.df_high_res["dt"]
        for h in self.alts: self.height_time_dne_df["alt_"+str(int(h))] = np.zeros(len(self.df_high_res["dt"]))
        return

    def __map_list(self):
        df = self.model.goes
        prod, T = self.model.exe()
        for i in range(prod.shape[1]):
           df["bin_"+str(i+1)] = prod[:,i]
           pass
        df = df.sort_values(by="dt")
        df = df[(df.dt>=self.stime) & (df.dt<self.etime)]
        df.indx = [i*60 for i in range(len(df))]
        self.df = df
        
        self.ds = 10.
        d = (self.etime-self.stime).total_seconds()/self.ds
        Tp = np.array([self.stime+dt.timedelta(seconds=self.ds*i) for i in range(int(d))])
        self.df_high_res = pd.DataFrame()
        self.df_high_res["dt"] = Tp
        L = len(Tp)
        keys = ["A_AVG","B_AVG"]
        for i in range(37): keys.append("bin_"+str(i+1))
        for k in keys:
            f_x = alp.get_extrp_func(df.indx,df[k])
            self.df_high_res[k] = f_x(np.arange(L)*self.ds)
            pass
        return

    def __get_eta(self, bin_num):
        eta = 0.
        for sps in self.euvac_config["species"]: eta = eta + np.array(self.euvac_config[sps]["eta"],dtype=np.float)[bin_num]
        eta = eta / len(self.euvac_config["species"])
        return eta

    def __estimate_photo_ionization(self, e_time, Hp, z, sza, bin):
        xe = self.df_high_res[self.df_high_res.dt == e_time]["bin_"+str(bin)].tolist()[0]
        self.dQ = 0.
        self.dQ = self.__estimate_iono_h0(xe, Hp, sza, z, bin-1)
        return

    def __chapman_function(self, z, sza):
        if sza >= 90.: n_Ne = 0.
        else: n_Ne = np.exp(1-z-(np.exp(-z)/np.cos(np.deg2rad(sza))))
        return n_Ne

    def __estimate_iono_h0(self, xe, Hp, sza, z, bin):
        n_Ne = self.__chapman_function(z, sza)
        eta = self.__get_eta(bin)
        q = (n_Ne[0] * xe * eta  / (Hp[0] * np.exp(1)))
        return q

    def __exe(self, e_time, h, h0, Hp, sza, bin):
        z = (h - h0) / Hp
        self.__estimate_photo_ionization(e_time, Hp, z, sza, bin)
        return self.dQ * 1e-6

    def __execute_T(self,ut,lat,lon):
        deltaNe_h = []
        print ut
        for h in self.alts:
            dNe = 0.
            sza = SZAC.get_sza_as_deg(ut, lat, lon, h)
            for bin in self.bins:
                Hp, h0 = SHC.get_h0_Hp(self.sco.alts,self.sco.msis,h,bin_num=bin)
                dNe = dNe + self.__exe(ut,h,h0,Hp,sza,bin)
                pass
            deltaNe_h.append(dNe)
            pass
        deltaNe_h = np.array(deltaNe_h)
        return deltaNe_h

    def __calc_grad(self):
        print len(self.height_time_dne_df)
        for i in range(len(self.alts)):
            y = self.height_time_dne_df["alt_"+str(int(self.alts[i]))]
            self.height_time_dne_df["grd_"+str(int(self.alts[i]))] = np.gradient(y)/self.ds
        return

    def execu(self,lat,lon):
        deltaNe_T = []
        for I,ut in enumerate(self.df_high_res["dt"]):
            self.sco = get_sco(self.geo,I)
            deltaNe_h = self.__execute_T(ut,lat,lon)
            deltaNe_T.append(deltaNe_h)
            #break
            pass
        deltaNe_T = np.transpose(np.array(deltaNe_T))*1.e6
        #print deltaNe_T.shape,len(deltaNe_T), self.height_time_dne_df.as_matrix().shape
        for i in range(len(deltaNe_T)):
            self.height_time_dne_df["alt_"+str(int(self.alts[i]))] = deltaNe_T[i,:]
        self.out = self.height_time_dne_df
        #self.__calc_grad()
        return self.height_time_dne_df, self.df_high_res


def estimate_delta_Ne(lat, lon, stime, etime, geo, alts=SU.get_alts()):
    goes_data_file = "Data/%s.csv"%stime.strftime("%Y-%m-%d")
    eir = EstimateIonizationRate(alts, goes_data_file, stime, etime, geo)
    eir.execu(lat,lon)
    return eir
