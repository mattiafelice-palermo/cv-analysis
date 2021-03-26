#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:51:31 2021

@author: mattia
"""
import utils

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy import integrate
from sys import exit


class SVAnalysis():
    def __init__(self, data):
        self.voltage, self.current = data
        self.ox_voltage = None
        self.ox_current = None
        self.red_voltage = None
        self.red_current = None

        self._split()

        self.anode_data = self.peak_current(self.ox_voltage, self.ox_current,
                                            'anode')
        self.cathode_data = self.peak_current(self.red_voltage, self.red_current,
                                            'cathode')
        
        self.cv_data = {}

        self._charge()

        self.plot(self.anode_data)
        self.plot(self.cathode_data)
        
    def ip_ratio(self):
        ip_a = self.anode_data['peak current']
        ip_c = abs(self.cathode_data['peak current'])
        
        self.cv_data['peak current ratio'] = ip_c/ip_a
        
        return ip_c/ip_a
    
    def _charge(self):
        self.anode_data['charge'] = a = integrate.trapezoid(self.ox_current,
                                                        x=self.ox_voltage)
                                                        
        self.cathode_data['charge'] = c = integrate.trapezoid(self.red_current,
                                                        x=-self.red_voltage)
        
        self.cv_data['charge ratio'] = c/a
        self.cv_data['charge difference'] = a-c
        
        
        
    def plot(self, data):
        mode = data['mode']
        peak = data['current maxval']
        fit = data['capacitive fit']
        peak_pos= data['peak voltage']
        
        
        plt.rcParams.update({'font.size': 12})
        plt.rc('legend',fontsize=9) 
        plt.xlabel('V vs ref [V]')
        plt.ylabel('i [mA]')            
        
        n_points = len(self.ox_voltage)
        rounded = n_points - (n_points % 5)
        
        if mode == 'anode':
            volt = self.ox_voltage[:rounded]
            curr = self.ox_current[:rounded]
        else:
            volt = self.red_voltage[:rounded]
            curr = self.red_current[:rounded]
            
        plt.plot(volt, curr,label = mode)

        # ranges
        if mode == 'anode':
            volt_min = np.amin(volt)
        else:
            volt_min = np.amax(volt)
        volt_max = peak_pos

        def line(x, a, b): return a*x + b
        x = np.linspace(volt_min, volt_max, 2)
        y = line(x, fit[0], fit[1])

        y_vert = volt_max*fit[0]+fit[1]

        plt.plot(x,y, color='green', label='capacitive '+mode)
        plt.plot([volt_max, volt_max], [y_vert, peak], 
                 color='orange', label = r'$i_{p,c}$ '+mode)

        plt.legend()
        plt.show()

    def _split(self):
        split = np.argmax(self.voltage) +1
        self.ox_voltage = self.voltage[:split]
        self.ox_current = self.current[:split]
        self.red_voltage = self.voltage[split:]
        self.red_current = self.current[split:]        
        
    def _maximas(self, arr, ave = 1, gauss = None):
        n_points = len(arr)
        rounded = n_points - (n_points % ave)
        arr = arr[:rounded]
        
        arr = np.mean(arr.reshape(-1, ave), axis=1)
        
        if gauss != None:
            arr = gaussian_filter1d(arr, gauss)
        
        peaks, _ = find_peaks(arr, prominence=0.0001)
        
        return peaks*ave
        
    def peak_current(self, in_volt, in_curr, mode, 
                       average = 10, d2_limit = 10, group_thr = 0.15):
        
        # reduce array lenght to divisible number by average
        n_points = len(in_volt)
        trim = n_points - (n_points % average)
        
        volt = in_volt[0:trim]
        curr = in_curr[0:trim]
        
        # block average on source
        volt = np.mean(volt.reshape(-1,average), axis=1)
        curr = np.mean(curr.reshape(-1,average), axis=1)
        
        curr = gaussian_filter1d(curr, 1) # smoothing
        
        # find CV peaks and consider just the first one
        # index is relative to the original, non averaged array
        if mode == 'anode':
            peaks = self._maximas(in_curr, 5, 1)
            peak_idx = peaks[0]
            peak_xpos = in_volt[peak_idx]
        else:
            peaks = self._maximas(-in_curr, 5, 1) # 3,3 were found empirically
            peak_idx = peaks[-1]
            peak_xpos = in_volt[peak_idx] # position of the first peak on x axis
            
        
        # Compute derivatives
        curr_d1 = np.gradient(curr, volt) # first derivative
        curr_d1 = gaussian_filter1d(curr_d1, 3)
        curr_d2 = np.gradient(curr_d1, volt) # second derivative

        # filter values outside of the vertical range & following the 1st peak
        curr_d2_max = np.amax(curr_d2)
        curr_d2_min = np.amin(curr_d2)
        max_thr = curr_d2_max/d2_limit
        min_thr = curr_d2_min/d2_limit

        y_mask = np.logical_and(curr_d2 < max_thr, curr_d2 > min_thr)
        if mode == 'anode':
            x_mask = volt < peak_xpos # before first peak!
        else:
            x_mask = volt > peak_xpos # before first peak!
        
        mask = np.logical_and(y_mask, x_mask) # indicates elements within thrs
        
        # Masked array with only points within thresolds
        volt_base = volt[mask]
        curr_base = curr[mask]
        curr_d2_base = curr_d2[mask]

        # thresholds for grouping (on second derivative!)
        volt_thr = 3*(np.amax(volt)-np.amin(volt))/n_points*average
        curr_thr = (np.amax(curr_d2)-np.amin(curr_d2))/10.0


        # Find clusters of data where each point is closer than group_thr
        clusters = utils.grouping(volt_base, curr_d2_base, volt_thr, curr_thr,
                                  order = 'size-decremental')
                        
        cap_idx = clusters[0]
        
        # fit!
        def line(x, a, b): return a*x + b
        params, _ = curve_fit(line, volt_base[cap_idx], curr_base[cap_idx])
        
        # Intersection between capacitive fit and peak position
        ip_ground = peak_xpos*params[0]+params[1]
        
        curr_max = in_curr[peaks[0]]
        
        ip = curr_max - ip_ground
        
        data = {'mode': mode,
                'peak current': ip,
                'peak ground': ip_ground,
                'current maxval': curr_max,
                'peak voltage': peak_xpos,
                'peak x index': peak_idx,
                'capacitive fit': params,
                'fit data indexes': cap_idx,
                'current 1st der': curr_d1,
                'current 2nd der': curr_d2,
                }
        
        
        return data
    
    def integrate(self):
        pass
    

