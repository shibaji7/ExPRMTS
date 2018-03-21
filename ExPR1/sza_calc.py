"""
Created on Jun 9, 2016

@author: Shibaji Chakraborty
"""

from datetime import datetime
from dateutil import parser
import pytz
from math import atan, atan2, degrees, floor, pi, sqrt, radians
from skyfield.api import earth, JulianDate, sun
import ephem
import numpy as np

class SZAProperties():
    """
        This class is used to set the properties needed to calculate SZA; it consist all setter functions.  
    """
    def __init__(self,time,lat,lon,alt):
        self.time           = time
        self.lat            = lat
        self.colat          = None
        self.colon          = None
        self.lon            = lon
        self.alt            = alt
        return
    
    def set_time(self):
        if type(self.time) is str:
            self.time = parser.parse(self.time)
            pass
        elif type(self.time) is datetime:
            pass
        else:
            pass
        self.time = self.time.replace(tzinfo=pytz.UTC)
        return
    
    def modify_to_deg(self):
        self.lat            = degrees(self.lat)
        self.colat          = 90-self.lat
        self.lon            = degrees(self.lon)
        self.colon          = 180-self.lon
        return
    
    def modify_to_rad(self):
        self.lat            = radians(self.lat)
        self.colat          = (pi/2)-self.lat
        self.lon            = radians(self.lon)
        self.colon          = pi-self.lon
        return
    

class SolarZenithAngleCalculator(object):
    """
    This class is used to find solar zenith angle calculation with respect to any geographical location on Earth, at any instance of 
    UT time.
    """

    def __init__(self,time,lat,lon,alt,loccation_type="degrees"):
        """
        Constructor takes location parameter & time  
        """
        self.prop           = SZAProperties(time,lat,lon,alt)
        self.prop.set_time()
        if loccation_type   == "degrees":
            self.prop.modify_to_rad()
        self.lat            = self.prop.lat
        self.colat          = self.prop.colat
        self.lon            = self.prop.lon
        self.colon          = self.prop.colon
        self.alt            = self.prop.alt
        self.time           = self.prop.time
        self.sun_lat        = None
        self.sun_lon        = None
        self.R              = 1.0
        return

    def find_subsolar_location_based_on_skyfield_locator(self):
        x, y, z             = earth(JulianDate(utc=self.time)).observe(sun).apparent().position.km
        julian_date         = JulianDate(utc=self.time).tt
        dublin_julian_date  = julian_date - 2415020
        sidereal_solar      = 1.0027379093
        sid_day             = floor(dublin_julian_date)
        t                   = (sid_day - 0.5) / 36525
        sid_reference       = (6.6460656 + (2400.051262 * t) + (0.00002581 * (t**2))) / 24
        sid_reference       -= floor(sid_reference)
        lon                 = 2 * pi * ((dublin_julian_date - sid_day) *  sidereal_solar + sid_reference) - atan2(y, x)
        lon                 = lon % (2 * pi)
        lon                 -= pi
        self.sun_lon        = -lon
        self.sun_lat        = atan(z / sqrt(x**2 + y**2))
        self.sun_colat      = (pi/2)-self.sun_lat
        self.sun_colon      = pi-self.sun_lon
        return
    
    def find_subsolar_location_based_on_ephem_locator(self):
        sun                 = ephem.Sun(epoch="date");
        str_time            = str(self.time).split(".")[0]
        sun.compute(str_time)
        self.sun_lat        = sun.dec
        self.sun_lon        = sun.ra
        self.sun_colat      = (pi/2)-self.sun_lat
        self.sun_colon      = pi-self.sun_lon
        return
    
    def get_sun_location_rc(self):
        return {"sun_lat":self.sun_lat, "sun_lon":self.sun_lon}
    
    def get_sun_location_degree(self):
        return {"sun_lat":degrees(self.sun_lat), "sun_lon":degrees(self.sun_lon)}
        
    def calc_sza(self, locator):
        if locator.lower()    == "skyfield":
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_skyfield_locator()
        else:
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_ephem_locator()
        n_hat               = np.array([(np.cos(self.lat)*np.cos(self.lon)),(np.cos(self.lat)*np.sin(self.lon)),np.sin(self.lat)])
        s_hat               = np.array([(np.cos(self.sun_lat)*np.cos(self.sun_lon)),(np.cos(self.sun_lat)*np.sin(self.sun_lon)),np.sin(self.sun_lat)])
        ndots               = np.dot(n_hat,s_hat)
        #n_mod               = np.sqrt((n_hat*n_hat).sum())
        #s_mod               = np.sqrt((s_hat*s_hat).sum())
        cos_angle           = ndots#/(n_mod*s_mod)
        self.sza            = np.arccos(cos_angle) # in radian
        return
    
    def calc_sza_from_cordlen(self, locator):
        if locator.lower()    == "skyfield":
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_skyfield_locator()
        else:
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_ephem_locator()
        dX                  = (np.cos(self.lat)*np.cos(self.lon))-(np.cos(self.sun_lat)*np.cos(self.sun_lon))
        dY                  = (np.cos(self.lat)*np.sin(self.lon))-(np.cos(self.sun_lat)*np.sin(self.sun_lon))
        dZ                  = (np.sin(self.lat)-np.sin(self.sun_lat))
        C                   = np.sqrt(dX**2 + dY**2 + dZ**2)
        self.sza            = 2 * np.arcsin(C/2) # in radian
        return
    
    def calc_dsza(self, locator):
        if locator.lower()    == "skyfield":
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_skyfield_locator()   
        else:
            pass
            if self.sun_lat is None or self.sun_lon is None:
                self.find_subsolar_location_based_on_ephem_locator()
        self.diff_lat       = abs(self.lat - self.sun_lat)
        self.diff_lon       = abs(self.lon - self.sun_lon)
        return
    
    def get_dsza_rc(self):
        return {"diff_lat":self.diff_lat, "diff_lon":self.diff_lon}
    
    def get_dsza_degree(self):
        return {"diff_lat": degrees(self.diff_lat), "diff_lon": degrees(self.diff_lon)}
    
    def get_sza_rc(self):
        return self.sza
    
    def get_sza_degree(self):
        return degrees(self.sza)
    
    def get_sun_colocation_rc(self):
        return {"sun_colat":self.sun_colat, "sun_colon":self.sun_colon}
    
    def get_sun_colocation_degree(self):
        return {"sun_colat": degrees(self.sun_colat), "sun_colon": degrees(self.sun_colon)}
    
    def get_obs_colocation_rc(self):
        return {"obs_colat":self.colat, "obs_colon":self.colon}
    
    def get_obs_colocation_degree(self):
        return {"obs_colat": degrees(self.colat), "obs_colon": degrees(self.colon)}
    
    def get_details_rc(self):
        return {"sun_colat":self.sun_colat, "sun_colon":self.sun_colon, "obs_colat":self.colat, "obs_colon":self.colon, "diff_lat":self.diff_lat, "diff_lon":self.diff_lon,
                "sun_lat":self.sun_lat, "sun_lon":self.sun_lon, "obs_lat": self.lat, "obs_lon": self.lon, "sza":self.sza}
    
    def get_details_degree(self):
        return {"sun_colat": degrees(self.sun_colat), "sun_colon": degrees(self.sun_colon), "obs_colat": degrees(self.colat), "obs_colon": degrees(self.colon),
                "diff_lat": degrees(self.diff_lat), "diff_lon": degrees(self.diff_lon), "sun_lat":degrees(self.sun_lat), "sun_lon":degrees(self.sun_lon), "sza":degrees(self.sza),
                "obs_lat": degrees(self.lat), "obs_lon": degrees(self.lon)}

def get_sza_from_cordlen_as_deg(time, lat, lon, alt, locator="skyfield"):
    calc = SolarZenithAngleCalculator(time, lat, lon, alt)
    calc.calc_sza_from_cordlen(locator)
    calc.calc_dsza(locator)
    return calc.get_details_degree()

def get_sza_as_deg(time, lat, lon, alt=100, locator="skyfield"):
    calc = SolarZenithAngleCalculator(time, lat, lon, alt)
    calc.calc_sza(locator)
    calc.calc_dsza(locator)
    return calc.get_details_degree()["sza"]

def get_subsolar_point_deg(time, locator="skyfield"):
    calc = SolarZenithAngleCalculator(time, 0, 0, 0)
    if locator == "skyfield":
        calc.find_subsolar_location_based_on_skyfield_locator()
    else:
        calc.find_subsolar_location_based_on_ephem_locator()
    return calc.get_sun_location_degree()
