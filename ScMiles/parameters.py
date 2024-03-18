#/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:21:21 2018

@author: Wei Wei

A class that stores all of the information from input 
"""

import os,sys
import pandas as pd

class parameters:
    def __init__(self, work_dir=None, path_jobsubmit=None, path_tinker=None, forcefield=None,
                 temp=298.15, pressure=1., ensemble=None, timeFactor=1.,
                 key_free=None, seek_traj=0, seekTime=1., seek_saveFrequency=1., RestartFrequency=None,
                 anchors_txt=None, path_anchors=None, AnchorNum=None,
                 reactant=None,product=None, key_restrain=None, restrain=None, restrain_grp=None, trajWidths=13,
                 bincoordinates=None, binvelocities=None, coor=None, saveFrequency=1., interval=1,
                 restrain_md_time=1.0, restrain_eq_time=0.5, res_saveFrequency=1.,
                 startTraj=None, trajPerLaunch=None, max_md_time=10., stepPerLaunch=1000,
                 coordinates=None, outputname=None, tolerance=1e-5, err_sampling=None) -> None:

        # Directory of your work
        self.work_dir = work_dir

        # Path of submission directory
        self.path_jobsubmit = path_jobsubmit  if path_jobsubmit else "/home/xy3866/bin_Miles/JobPool/"

        # Path of tinker executives
        self.path_tinker = path_tinker

        # Name of forcefield prm
        self.forcefield = forcefield

        # Path of forcefield prm
        self.path_forcefield = None

        # Path of program
        tmp_ = os.getcwd()
        self.path_program = tmp_

        # Path of log file
        self.path_log = None

        # Path of record file
        # Record the necessary intermediate math data
        self.path_math = None

        # The current iteration number of the complete milestoning processes
        self.iteration = 0

        # Parallelize the submission
        self.num_sub = 1

        # Type of milestoning
        # 0 for classic milestoning; 1 for exact milestoning: iteration
        self.method = 0
        
        # Max num. of iterations
        self.maxIteration = 100
        
        # Time step size used by all MD jobs
        # unit: fs, default: 1.0 fs
        self.timeFactor = float(timeFactor)

        # Temperature of the system
        # unit: K, default: 298.15 K
        self.temp = float(temp)

        # Pressure of the system
        # unit: atm, default: 1.0 atm
        self.pressure = float(pressure)
        
        # Ensemble: NVT or NPT
        self.ensemble = ensemble

        # Network for free trajectories
        self.network = {}

        # Search method for possible milestones
        # 0 for traverse; 1 for seek
        self.milestone_search = 0
        
        # Milestone list
        self.MS_list = []

        # Milestone property dict
        self.MS_property = {}

        # Periodic boundary for 1D case
        self.pbc = []
        
        # Current iteration time
        self.current_iteration_time = {}
        
        ## Seek setting
        # Name of tinker key file for unbiased traj (both in seeking process & production MD)
        self.key_free = key_free

        # Number of trajs starting from each anchor in seeking
        # Default: 0
        self.seek_traj = int(seek_traj)
        
        # Total MD running time for each seeking traj
        # unit: ns, default: 1.0
        self.seekTime = float(seekTime)

        # Saving frequency (to arc) for each seeking trajs
        # unit: ps
        self.seek_saveFrequency = float(seek_saveFrequency)

        # Restarting frequency
        self.RestartFrequency = RestartFrequency

        # Threshold to tell if any config hits the specified milestone
        self.MS_threshold = 0.

        # ignore new milestones found by free trajs or not
        self.ignorNewMS = False

        # Name of the description txt file of anchors
        self.anchors_txt = anchors_txt

        # Path of directory of anchor tinker xyz
        self.path_anchors = path_anchors

        # Total num of anchors
        self.AnchorNum = AnchorNum

        # Find new anchors
        self.new_anchor = False

        # New anchor separation distance
        self.anchor_dist = 100.0

        # New MS reached
        self.MS_new = []

        # MS num of reactant and product
        self.boundary = [-1, -1]

        # MS list for reactant
        self.reactant_milestone = []

        # MS list for product
        self.product_milestone = []

        # Voronoi cell for reactant
        self.reactant = reactant

        # Voronoi cell for product
        self.product = product

        ## Restrained MD setting
        # Name of tinker key file for restrained traj
        self.key_restrain = key_restrain

        # Type of restraint
        # type: dihedral, 2D dihedral, RMSD, gyration
        self.restrain = restrain

        # Restrain group
        self.restrain_grp = restrain_grp

        # For RMSD milestoing only
        # RMSD max
        self.RMSD_max = 0.
        # the number of divisions between starting and final config along RMSD traj
        self.RMSD_num = 0
        # RMSD starting config path
        self.path_RMSD_start = ''
        # RMSD final config path
        self.path_RMSD_final = ''
        # atom number of the solute in RMSD system
        self.RMSD_solute_num = 0

        # Force constant of the harmonic restraint
        self.forceConst = 1.

        # Width for each traj output, default = 13
        self.trajWidths = int(trajWidths)

        # Coordinates file name
        self.bincoordinates = bincoordinates

        # Velocity file name
        self.binvelocities = binvelocities
 
        # Coordinates
        self.coor = coor

        # Interval between each frame (choosing config in restrained MD for production MD)
        self.interval = int(interval)

        # Saving frequency in restrained MD
        # unit: ps, default: 1.0 ps
        self.res_saveFrequency = float(res_saveFrequency)

        # Saving frequency in production MD
        # unit: ps, default: 1.0 ps
        self.saveFrequency = float(saveFrequency)

        # Max time for restrained MD
        # unit: ns, default: 1.0 ns
        self.restrain_md_time = float(restrain_md_time)

        # Equilibrium time for restrained MD
        # unit: ns, default: 0.5 ns
        self.restrain_eq_time = float(restrain_eq_time)

        ## Production MD setting
        # Initial sampling frame
        self.startTraj = startTraj

        # Num of the free trajs to launch each time
        self.trajPerLaunch = trajPerLaunch

        # Num of the free trajs to launch each time (real-time)
        self.real_trajPerLaunch = trajPerLaunch

        # Max num of the free trajs to launch each time
        self.max_trajPerLaunch = 100

        # 
        self.sampling_time = 100 # sampling time for a single MD in final sampling (ps)
        
        # Max time of the final production MD
        # unit: ns, default: 10.0 ns
        self.max_md_time = float(max_md_time)

        # Structure file name
        self.coordinates = coordinates

        # Output name
        self.outputname = outputname if outputname else "output"

        # Path of output
        self.path_output = None

        # Tolerance of MFPT convergence
        self.tolerance = float(tolerance)

        # Num of error sampling
        self.err_sampling = 5000

        # MS whose free traj is completed
        self.Finished = []

        # MS whose sampling is completed
        self.finished_restrain = []
        self.skip_restrain = False

        # k matrix singularity
        self.sing = True

        # Smooth list limit
        self.smooth_list_limit = 5

        # kij matrix
        self.kij = []

        # MS index
        self.index = []

        # Flux
        self.flux = []

        # Probability
        self.prob = []

        # Free energy
        self.free_energy = []

        # Energy error
        self.energy_err = []

        # MFPT
        self.MFPT_1 = 0.0
        self.MFPT_2 = 0.0
        self.MFPT_err = 0.0
        self.MFPT_err2 = 0.0
        self.MFPT_1_rev = 0.0
        self.MFPT_2_rev = 0.0
        self.MFPT_err_rev = 0.0
        self.MFPT_err2_rev = 0.0

        # Committor vector
        self.c = []

        self.__read_input()
        
    def __read_input(self):
        from log import log
        scriptPath = os.path.dirname(os.path.abspath(__file__)) 
        inputfolder = os.path.abspath(os.path.join(scriptPath, os.pardir)) + '/my_project'
        if not self.work_dir:
            self.work_dir = inputfolder
        inputFile = os.path.join(self.work_dir, 'input.txt')
        self.path_output = os.path.join(self.work_dir, 'output')
        self.path_log = os.path.join(self.work_dir, 'work.log')
        self.path_math = os.path.join(self.work_dir, 'math.log')
        with open(file=inputFile) as f:
            for line in f:
                line = line.rstrip()
                info = line.split(" ")
                
                if(len(line) > 0 and line[0] == '#'):
                    continue
                
                if "path_jobsubmission" in info:
                    self.path_jobsubmit = str(info[1])
                if "path_tinker" in info:
                    self.path_tinker = str(info[1])
                if "forcefield" in info:
                    self.forcefield = str(info[1])
                if "method" in info:
                    self.method = int(info[1])
                if "num_sub" in info:
                    self.num_sub = int(info[1])

                if "initial_iteration" in info:
                    self.iteration = int(info[1]) - 1
                if "max_iteration" in info:
                    self.maxIteration = int(info[1])
                
                if "milestoneSearch" in info:
                    self.milestone_search = int(info[1])
                if "pbc" in info:
                    rm = line.replace(","," ").replace("  "," ").split(" ")
                    rm.pop(0)
                    self.pbc = list(map(int, rm))

                if "coordinates" in info:
                    self.coordinates = info[1]

                if "outputname" in info:
                    self.outputname = info[1]
                    self.path_output = os.path.join(self.work_dir, self.outputname)

                if "key_free" in info:
                    self.key_free = str(info[1])
                if "key_restrain" in info:
                    self.key_restrain = str(info[1])

                if "ensemble" in info:
                    if str(info[1]).lower() == 'nvt':
                        self.ensemble = 'NVT'
                    elif str(info[1]).lower() == 'npt':
                        self.ensemble = 'NPT'
                if "temperature" in info:
                    self.temp = float(info[1])
                if "pressure" in info:
                    self.pressure = float(info[1])

                if "time_step" in info:
                    self.timeFactor = float(info[1])    

                if "seek_traj" in info:
                    self.seek_traj = int(info[1])
                if "seek_time" in info:
                    self.seekTime = float(info[1])                  
                if "seek_save_frequency" in info:
                    self.seek_saveFrequency = float(info[1])                  
                if "ignore_new_ms" in info:
                    self.ignorNewMS = True if info[1] == 'yes' or info[1] == 'on' else False   

                if "anchorsNum" in info:
                    self.AnchorNum = int(info[1])    
                if "find_new_anchor" in info:
                    if str(info[1]).lower() == 'true' or 'yes' or 'on':
                        self.new_anchor = True
                if "new_anchor_dist" in info:
                    self.anchor_dist = float(info[1])         
                if "anchors_txt" in info:
                    self.anchors_txt = info[1]
                if "anchors_dir" in info:
                    self.path_anchors = info[1]

                if "force_const" in info:
                    self.forceConst = float(info[1])
                if "restrain_save_frequency" in info:
                    self.res_saveFrequency = float(info[1])
                if "save_frequency" in info:
                    self.saveFrequency = float(info[1])
                if "skip_restrain" in info:
                    self.skip_restrain = True
                if "restrain_type" in info:
                    self.restrain = str(info[1])
                    if(self.restrain == "dihedral"):
                        self.restrain_grp = [int(info[2]),int(info[3]),int(info[4]),int(info[5])]
                        self.MS_threshold = 0.1
                    elif(self.restrain == "2Ddihedral"):
                        self.restrain_grp = [int(info[2]),int(info[3]),int(info[4]),int(info[5]),
                                             int(info[6]),int(info[7]),int(info[8]),int(info[9])]
                        self.MS_threshold = 0.1
                    elif(self.restrain == "RMSD"):
                        self.restrain_grp = [0]
                        self.RMSD_max = float(info[2])
                        self.RMSD_num = int(info[3])
                        self.RMSD_solute_num = int(info[4])
                        if(self.RMSD_num == 0 or self.RMSD_max == 0. or self.RMSD_solute_num == 0):
                            log("Error in command of RMSD restraint!\n", self.path_log)
                            log("check RMSD_max, RMSD_num or RMSD_solute_num\n", self.path_log)
                            sys.exit(0)
                    else:
                        log("Error in command of restraint!\n", self.path_log)
                        sys.exit(0)
                if "restrain_md_time" in info:
                    self.restrain_md_time = float(info[1])
                if "restrain_eq_time" in info:
                    self.restrain_eq_time = float(info[1])
                if "milestone_threshold" in info:
                    self.MS_threshold = float(info[1])
                if "max_md_time" in info:
                    self.max_md_time = float(info[1])
                if "interval" in info:
                    self.interval = int(info[1])

                # reactant_milestone
                if "reactant" in info:
                    rm = line.replace(","," ").replace("  "," ").split(" ")
                    rm.pop(0)
                    if len(rm) == 2:
                        self.reactant = list(map(int, rm))
                    elif len(rm) == 1:
                        self.reactant = [int(rm[0])]
                # product_milestone
                if "product" in info:
                    pm = line.replace(","," ").replace("  "," ").split(" ")
                    pm.pop(0)
                    if len(pm) == 2:
                        self.product = list(map(int, pm))
                    elif len(pm) == 1:
                        self.product = [int(pm[0])]                       

                if "start_traj" in info:
                    self.startTraj = int(info[1])                  
                if "traj_per_launch" in info:
                    self.trajPerLaunch = int(info[1])  
                    self.real_trajPerLaunch = int(info[1])  
                if "max_traj_per_launch" in info:
                    self.max_trajPerLaunch = int(info[1])  
                
                if "tolerance" in info:
                    self.tolerance = float(info[1])   
                if "error_sampling" in info:    
                    self.err_sampling = int(info[1])   
                if "sampling_time" in info:
                    self.sampling_time = float(info[1])
                if "finished_restrain" in info:
                    i1 = int(info[1].split(',')[0])
                    i2 = int(info[1].split(',')[1])
                    self.finished_restrain.append((i1, i2))
                      

    def print_properties(self):
        from log import log_properties, log_properties_2
        log_properties(self.flux, self.prob, self.free_energy, self.c, self.energy_err, self.MS_list, self.path_log)
        log_properties_2(self.MFPT_1, self.MFPT_2, self.MFPT_err, self.MFPT_err2, self.path_log)
        log_properties_2(self.MFPT_1_rev, self.MFPT_2_rev, self.MFPT_err_rev, self.MFPT_err2_rev, self.path_log)

    def initialize(self):
        from log import log
        if os.path.exists(self.path_log):
            os.remove(self.path_log) 
        if not self.anchors_txt:
            self.anchors_txt = os.path.join(self.work_dir, 'anchors.txt')
        else:
            self.anchors_txt = os.path.join(self.work_dir, self.anchors_txt)
        if not self.path_anchors:
            self.path_anchors = os.path.join(self.work_dir, 'anchors')
        else:
            self.path_anchors = os.path.join(self.work_dir, self.path_anchors)
        if(self.restrain != "RMSD" and not os.path.isfile(self.anchors_txt)):
            log("Error: No anchor description is found in " + self.anchors_txt + '\n', self.path_log)
            sys.exit(0)
        if not os.path.isdir(self.path_anchors):
            log("Error: No anchor file is found in " + self.path_anchors + '\n', self.path_log)
            sys.exit(0)

        if(self.restrain != "RMSD"):
            self.anchors = pd.read_fwf(self.anchors_txt, header=None).values
        else:
            self.path_RMSD_start = os.path.join(self.work_dir, "start.xyz")
            self.path_RMSD_final = os.path.join(self.work_dir, "final.xyz")
            if(not os.path.isfile(self.path_RMSD_start) or not os.path.isfile(self.path_RMSD_final)):
                log("Error: No starting or final config is found in " + self.work_dir + " (start.xyz and final.xyz)", self.path_log)
                sys.exit(0)

        if not self.forcefield:
            self.forcefield = 'amoeba.prm'
        self.path_forcefield = os.path.join(self.work_dir, self.forcefield)
        if not os.path.isfile(self.path_forcefield):
            log("Error: No tinker prm file is found in " + self.path_forcefield + '\n', self.path_log)
            sys.exit(0)

        if not self.key_free:
            self.key_free = os.path.join(self.work_dir, 'free.key')
        else:
            self.key_free = os.path.join(self.work_dir, self.key_free)
        if not os.path.isfile(self.key_free):
            log("Error: No tinker key file for free MD is found in " + self.key_free + '\n', self.path_log)
            sys.exit(0)

        if not self.key_restrain:
            self.key_restrain = os.path.join(self.work_dir, 'restrain.key')
        else:
            self.key_restrain = os.path.join(self.work_dir, self.key_restrain)
        if not os.path.isfile(self.key_restrain):
            log("Error: No tinker key file for restrained MD is found in " + self.key_restrain + '\n', self.path_log)
            sys.exit(0)

        crdfolder = os.path.join(self.work_dir,'crd')
        if not os.path.exists(crdfolder):
            os.makedirs(crdfolder)
        log("Initialized with {} anchors.\n".format(self.AnchorNum), self.path_log)

if __name__ == '__main__':
    new = parameters()
    new.initialize()
    print(new.timeFactor)
