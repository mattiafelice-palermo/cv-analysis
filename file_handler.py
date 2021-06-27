#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 11:26:48 2021

@author: mattia
"""

from read_input import CyclicVoltammetry
from cvanalysis import SCVAnalysis

import os
import glob


class FileHandler:
    def __init__(self):
        self.cvs = {}  # filepath: cyclic voltammetry objects
        self.analyses = {}  # filepath: SCVAnalysis objects

    def open_folder(self, folder):
        globber = glob.glob(folder + "/*")

        self.open_files(globber) 
    
    def open_files(self, filepaths):
        for filepath in filepaths:
            if not os.path.isfile(filepath):
                continue
            # Fixes mixed slashes/blackslashes, so that filepath can be unique key
            # Thanks Windows.
            filepath = os.path.normpath(filepath)
            
            # Try reading cyclicvoltammetry file
            try:
                cv = CyclicVoltammetry(filepath)
            # TODO: the following should be fixed with a custom FileTypeError
            except Exception:
                print(f"Could not read file {filepath}")
                continue
            self.cvs[filepath] = cv

            n_cycles = cv.settings["n_cycles"]

            self.analyses[filepath] = []

            for cycle in range(n_cycles):
                try:
                    analysis = SCVAnalysis(cv, cycle)
                # If analysis fails for any reason, discard cycle (e.g. single point cycle)
                except Exception:
                    cv.settings["n_cycles"] -= 1
                    print(f"Invalid cycle {cycle} in file {filepath}")
                    continue

                analysis.compute_ip()
                self.analyses[filepath].append(analysis)

    def __iter__(self):
        for key in self.cvs:
            yield key, self.cvs[key]
