#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 08:40:37 2026

@author: shill
"""

#def dummy():
from scipy.signal import convolve2d
import sys
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Services/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Maps/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/HST/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Winds/")
import get_spice_ephem
import convert_system3_to_I_II_spice
import PlanetMapper_Spice_Furnish
import HST_Analysis_Script

