# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 13:35:19 2020

@author: ilail
"""

import pandas as pd
import re
from os import path

#%% Reads MPT or DTA to common data structure

class CyclicVoltammetry:
    def __init__(self, path):
        self.filepath = path # path of input file
        self.settings = {} # info from input header
        self.data = None # contains Pandas dataframe of CV
        
        self._load_cv() # read input data
        
        
    def _load_cv(self):
        # choose input parser based on file extension
        extension = path.splitext(self.filepath)[1]
        
        if extension.lower() == '.dta':
            self.settings['format'] = 'Gamry'
            self._read_DTA()    
        if extension.lower() == '.mpt':
            self._read_MPT()
            self.settings['format'] = 'Biologic'
    
    def _read_DTA(self):
        # only consider data after label CURVE\d and ignore CURVEOCV
        is_it_curve = re.compile('CURVE\d')
        with open(self.filepath, "r",  encoding="utf8", errors='ignore') as f:        
            curves=[]
            curveN=-1
            
            iterator = iter(f) # iterator, avoid checking settings after header
            
            # headers for self.settings dict
            for line in iterator:
                if 'SCANRATE' in line:
                    scanrate = line.replace(',','.').split('\t')[2]
                    self.settings['scan rate'] = float(scanrate)
                    break
            
            # read CV data
            for line in iterator:
                if is_it_curve.match(line):
                    curveN+=1
                    next(f)
                    next(f)
                    continue
                if curveN>=0:
                    line_float = [float(ele) for ele in line.strip().replace(',','.').split('\t')[0:7]]
                    curves.append([curveN] + line_float)
        
            self.settings['n cycles'] = curveN            
            header=['Cycle', 'Pt', 'T', 'Vf', 'Im', 'Vu', 'Sig0', 'Ach']
            # just keep the relevant ones and use cycle num as vertical key
            self.data = pd.DataFrame(curves, columns=header)[['Cycle', 'T', 'Vf', 'Im' ]]
            self.data.set_index('Cycle', inplace=True)
            
    def _read_MPT(self):
        pass
            
    def __getitem__(self, key):
        if key in self.data.columns:
            return self.data[key]
        elif isinstance(key, int):
            return self.data.loc[key]
        #else:
        #    return self.data.filter(like = key, axis = 0)
    
    def __repr__(self):
        return self.data.__repr__()
    
    def _current_dens(self,area):
        self.data['Idens'] = self.data['Im']*1000/area
    


#%% 


