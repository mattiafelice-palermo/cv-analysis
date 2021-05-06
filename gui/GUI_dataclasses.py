#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 15:49:14 2021

@author: mattia
"""
from dataclasses import dataclass


@dataclass
class MplChekboxStates:
    anode_curr: bool
    anode_ip: bool
    anode_capfit: bool
    cathode_curr: bool
    cathode_ip: bool
    cathode_capfit: bool
