#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 11:06:15 2018

@author: Wei Wei

Main script which inclues major workflow.

"""

import time
import numpy as np
from parameters import *
from network_check import *
from log import log
#from sampling import *
from milestones import *
#from analysis import analysis_kernel

# run free trajectories without sampling
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--skipSampling', action='store_true', help='skip sampling',
                    required=False)
args = parser.parse_args()  
status = 1 if args.skipSampling else 0

# initialize environment
MFPT_temp = 1
MFPT_converged = False
work_path = '/work/xy3866/ScMiles_tinker/test_Alanine_Dipeptide_12'
parameter = parameters(work_dir=work_path)
parameter.initialize()
#jobs = run(parameter)
#samples = sampling(parameter, jobs)

# initialize with reading anchor info and identifying milestones 
MS = milestones(parameter)
MS.initialize(status=status)
print(MS.parameter.MS_list)

# initialize iteration number
MS.parameter.iteration = 0

MS.sampling()

