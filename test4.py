#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:25:35 2021

@author: mattia
"""
import numpy as np
from read_input import CyclicVoltammetry
from cvanalysis import SCVAnalysis, CCVAnalysis
from graphical_tools import MaxMinFinder, LineFit
import os
import glob
from sys import exit

# folder = r"/home/mattia/Downloads/DTAs/"

# files = []
# glob = glob.glob(folder + "*")
# glob.sort(key=os.path.getmtime)

# for i, filename in enumerate(glob):
#     print(i, filename)
#     path = os.path.join(folder, filename)
#     files.append(path)


#path = files[10]
path = r"/home/mattia/Downloads/DTAs/PtvsPt_Mn(IV)-Mn(II)_CV_10mVs_06-16_C03.mpt"
data = CyclicVoltammetry(path)


exit()

data = CyclicVoltammetry(path)
analysis = SCVAnalysis(data, 2)
# analysis2 = SCVAnalysis(data, 0)

analysis.compute_ip(automatic=False, work_mode="both")
# analysis.compute_ip(work_mode="cathode", automatic=False)

exit()
