import numpy as np
import json

class PSH_Estimator(object):

    def __init__(self,alts,msis,h,bin_num=1,Re=6317):
        self.I = np.where(alts==h)[0]
        self.bin_num = bin_num
        self.alts = alts
        self.msis = msis
        self.Re = Re
        self.Hp = None
        self.sigma = 0.
        with open("euvac.json", "r") as code: self.euvac_config = json.load(code)
        for sps in self.euvac_config["species"]: self.sigma = (msis[sps][self.I]/msis["nn"][self.I]) * np.array(self.euvac_config[sps]["sigma"],dtype=np.float)[self.bin_num]
        #self.sigma = 1e-18 * self.sigma / len(self.euvac_config["species"]);
        self.sigma = 1e-18 * self.sigma;
        self.h0 = 0.0
        return

    def __estimate_Hp(self):
        Tn = self.msis["Tn"]
        self.Hp = 0.85 * (1 + (self.alts/self.Re) )**2 * Tn / self.__mole_mass(self.msis)
        return

    def __calculate_ref_height(self):
        if self.sigma > 0:
            P = self.msis["nn"] * self.sigma * self.Hp * 100000
            J = np.argmin(np.abs(P-1))
            self.h0 = self.alts[J]
            pass
        return

    def __mole_mass(self,density):
        m_O = 15.994
        m_O2 = 2*15.994
        m_N = 14.0067
        m_N2 = 2*14.0067
        m_H = 1.00794
        m_He = 4.002602
        m_Ar = 39.948
        rho = density["nn"]
        M = (density["O"]*m_O + density["O2"]*m_O2 + density["N"]*m_N + density["N2"]*m_N2 + density["H"]*m_H + 
                density["HE"]*m_He + density["AR"]*m_Ar + density["O_anomalous"]*m_O) / rho
        return M

    def exe(self):
        self.__estimate_Hp()
        self.__calculate_ref_height()
        return self.Hp[self.I],self.h0

    @staticmethod
    def convert_to_ref_alts(alts,Hp,h0):
        z = (alts-h0)/Hp
        return z

def get_h0_Hp(alts,msis,h,bin_num=1):
    PSHE = PSH_Estimator(alts,msis,h,bin_num=bin_num)
    Hp, h0 = PSHE.exe()
    return Hp, h0
