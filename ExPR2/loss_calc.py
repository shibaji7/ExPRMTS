import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


font = {"family": "serif", "color":  "darkred", "weight": "normal", "size": 12}
fonttext = {"family": "serif", "color":  "darkblue", "weight": "normal", "size": 12}
fonttitle = {"family": "serif", "color":  "black", "weight": "normal", "size": 14}
fontsuptitle = {"family": "serif", "color":  "black", "weight": "bold", "size": 20}


class Recombination(object):

    def __init__(self, ds):
        self.ds = ds
        self.dr = self.ds.dr
        self.rr = self.ds.rr
        self.dt = 1e-4
        self.sim_time = 2/self.dt
        self.dQ = np.transpose(self.ds.dQ)
        self.Ne = np.zeros(self.dQ.shape)
        self.__initialization()
        return

    def __initialization(self):
        self.Ne[0,:] = self.ds.Ne
        return

    def __apply_dissociative_loss(self,i):
        sim_time = self.sim_time
        dt = self.dt
        while sim_time >= 0:
            self.Ne[i,:] = self.Ne[i,:] - (dt*np.power(self.Ne[i,:],2)*self.dr)
            sim_time -= 1
            pass
        return

    def __apply_radiative_loss(self,i):
        sim_time = self.sim_time
        dt = self.dt
        while sim_time >= 0:
            self.Ne[i,:] = self.Ne[i,:] - (dt*np.power(self.Ne[i,:],2)*self.rr)
            sim_time -= 1
            pass
        return

    def __plot_ionization(self,I):
        if np.mod(I,100) == 0:
            fig,ax = plt.subplots(figsize=(15,10),nrows=1,ncols=1, dpi=80)
            ax.semilogx(self.Ne[I,:],self.ds.alts,label=r"$n_{e}$")
            ax.set_ylabel("Heights [in km]", fontdict=font)
            ax.set_xlabel(r"$n_e$ [in $cm^{-3}$]", fontdict=font)
            ax.set_title("Electron density vs. hight profile", fontdict=fonttitle)
            fig.savefig(self.ds.this_plasma.replace("{}",str(I)))
            plt.close()
        return

    def estimate_density_evolution_with_loss(self):
        for I,dq in enumerate(self.dQ):
            self.Ne[I,:] = self.Ne[I,:] + dq
            self.__apply_dissociative_loss(I)
            self.__apply_radiative_loss(I)
            print I, np.max(dq)/1e10
            #self.__plot_ionization(I)
            pass
        self.ds.Ne_evol = self.Ne
        return self.ds


class Absorption(object):

    def __init__(self, ds):
        self.ds = ds
        return

    def __estimate_nr2(self, Ne):
        e = 1.6e-19
        m = 9.1e-31
        eps0 = 1e-9/(36*np.pi)
        nr2 = 1 - ((Ne * 1e6 * e**2 / (eps0 * m)) / self.ds.f0**2) 
        return nr2

    def __estimate_absorption(self):
        self.ds.absorptions["abs"] = np.zeros(self.ds.Ne_evol.shape)
        for I in range(len(self.ds.time)):
            Ne = self.ds.Ne_evol[I,:]
            nu = self.ds.collision_freq["nu"]
            w = 2 * np.pi * self.ds.f0
            self.ds.nr2[:,I] = self.__estimate_nr2(Ne)
            self.ds.absorptions["abs"][I,:] = 4.6e-2 * Ne * nu / (w**2 + nu**2)
        return

    def __plot_absorption_evolution(self):
        fig, ax = plt.subplots(1, figsize=(10,6))
        ax.set_xlabel(r"Time [in $UT$]",fontdict=font)
        ax.set_ylabel(r"Height [in $km$]",fontdict=font)
        ax.set_title(r"Time evolution of absorption coefficient $\beta$",fontdict=fonttitle)
        vmin = self.ds.absorptions["abs"].min()
        vmax = self.ds.absorptions["abs"].max()
        print len(self.ds.time),len(self.ds.alts), self.ds.absorptions["abs"].shape
        p = ax.pcolormesh(self.ds.time, self.ds.alts, self.ds.absorptions["abs"].transpose(), cmap="RdBu",vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(p)
        cbar.ax.set_ylabel(r"Absorption coefficient $\beta$ [in $dB/km$]")
        fig.savefig(self.ds.this_absorption_coeff_evolution)
        #plt.show()
        plt.close()
        return

    def exe(self):
        self.__estimate_absorption()
        self.__plot_absorption_evolution()
        return
