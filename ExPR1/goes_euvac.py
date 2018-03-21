import pandas as pd
import datetime as dt
import numpy as np
import json

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import pyglow

class EuvacModel(object):

    def __init__(self, goes_data_file, bin_num = 1, goes_tags=["A_AVG","B_AVG"]):
        self.goes_data_file = goes_data_file
        self.goes = pd.read_csv(goes_data_file)
        self.goes.dt = pd.to_datetime(self.goes.dt)
        self.bin_num = bin_num
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
        production = XE[:,self.bin_num]
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


class EstimateIonization(object):

    def __init__(self, alts, goes_data_file, bin_num = 1, model = None):
        self.alts = alts
        self.goes_data_file = goes_data_file
        self.bin_num = bin_num
        if not model: model = EuvacModel(goes_data_file, bin_num)
        self.prod, self.T = model.exe()
        with open("euvac.json","r") as code: self.euvac_config = json.load(code)
        return

    def __estimate_photo_ionization(self, e_time, Hp, z, sza, etap=None):
        xe = self.prod[np.where(self.T == e_time)][0]
        if etap is None or etap < 1.:
            eta = 0.
            for sps in self.euvac_config["species"]:
                eta = eta + np.array(self.euvac_config[sps]["eta"],dtype=np.float)[self.bin_num]
                pass
            eta = eta / len(self.euvac_config["species"])
            if etap is not None: eta = eta * etap
        else:
            eta = etap
        self.dQ = self.__estimate_iono_h0(xe, Hp, eta, sza, z)
        return

    def __chapman_function(self, z, sza):
        if sza >= 90.: n_Ne = 0.
        else: n_Ne = np.exp(1-z-(np.exp(-z)/np.cos(np.deg2rad(sza))))
        return n_Ne

    def __estimate_ionization(self, e_time, Hp, z, sza):
        self.dQ = 0.0
        xe = self.prod[np.where(self.T == e_time)][0]
        for sps in self.euvac_config["species"]:
            eta = np.array(self.euvac_config[sps]["eta"],dtype=np.float)[self.bin_num]
            self.dQ = self.dQ + self.__estimate_iono_h0(xe, Hp, eta, sza, z)
            pass
        return

    def __estimate_iono_h0(self, xe, Hp, eta, sza, z):
        n_Ne = self.__chapman_function(z, sza)
        q = (n_Ne * xe * eta  / (Hp * np.exp(1)))
        return q

    def exe(self, e_time, h, h0, Hp, sza, eta=None):
        self.I, = np.where( self.alts==h )
        z = (h - h0) / Hp
        self.__estimate_photo_ionization(e_time, Hp, z, sza, eta)
        #self.__estimate_ionization(e_time, Hp, z, sza)
        return self.dQ * 1e-6
