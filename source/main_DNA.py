#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
work_path = '/work/xy3866/work2/ScMiles_tinker/DNA'
parameter = parameters(work_dir=work_path)
parameter.initialize()
#jobs = run(parameter)
#samples = sampling(parameter, jobs)

# initialize with reading anchor info and identifying milestones 
MS = milestones(parameter)
MS.initialize(status=1)
print(MS.parameter.MS_list)

# initialize iteration number
MS.parameter.iteration = 0

MS.sampling()

