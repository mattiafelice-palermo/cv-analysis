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

    def add_folder(self, folder):
        globber = glob.glob(folder + "/*")

        for filename in globber:
            if not os.path.isfile(filename):
                continue
            filepath = os.path.join(folder, filename)
            try:
                cv = CyclicVoltammetry(filepath)
            # the following should be fixed with a custom FileTypeError
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

    # Likely not used: to be deleted soon
    #def add_cvs(self, filepath):
    #    self.cvs[filepath] = CyclicVoltammetry(filepath)

    def __iter__(self):
        for key in self.cvs:
            yield key, self.cvs[key]
