#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 17:25:27 2018

@author: Wei Wei

This subroutine generates NAMD configuration based on templete and submits jobs.
"""

import os, time
from shutil import copy
import subprocess
from shutil import copyfile
#from find_milestone import *
#from milestones import *
from log import log


class run:
    
    def __init__(self, parameter, work_path, job_name, coord_name, steps, save_frequency, additional_path=None, more_cmd='') -> None:
        self.parameter = parameter
        self.work_path = work_path
        self.job_name = job_name
        self.coord_name = coord_name
        self.steps = steps
        self.save_frequency = save_frequency
        self.additional_path = additional_path # For reference txyz
        self.more_cmd = more_cmd
        
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        return 
            
    def __repr__(self) -> str:
        return ('Submit jobs.')

    def submit(self):
        '''job submission'''
        hoststr = os.getenv('HOSTNAME').split('.')[0]
        timestr = str(time.time()).replace('.', '')
        scriptfile = os.path.join(self.parameter.path_jobsubmit, f"{hoststr}-{timestr}.sh")
        shfile = os.path.join(self.work_path, f"{self.job_name}.sh")
        temperature = self.parameter.temp
        pressure = self.parameter.pressure
        with open(shfile, 'w') as f0:
            c = self.coord_name
            log(f"Attention: submit job {self.job_name}", self.parameter.path_log)
            f0.write("source ~/.forcebalanceOrganic\n")
            if(self.parameter.ensemble == 'NVT'):
                f0.write("$DYNAMIC %s -k %s.key %i %f %f 2 %f > %s-results.log \n" %(c+'.xyz',self.job_name,self.steps,self.parameter.timeFactor,self.save_frequency,temperature,c))
            elif(self.parameter.ensemble == 'NPT'):
                f0.write("$DYNAMIC %s -k %s.key %i %f %f 4 %f %f > %s-results.log \n" %(c+'.xyz',self.job_name,self.steps,self.parameter.timeFactor,self.save_frequency,temperature,pressure,c))
            if(self.parameter.restrain == "RMSD" and self.job_name == c):
                path_ref = self.parameter.path_RMSD_start
                n = self.parameter.RMSD_solute_num
                f0.write("superpose %s -k %s.key %s 1 1 %d N N M N 0.0 > RMSD-%s-start.log \n" %(path_ref,c,c+'.arc',n,c))
                path_ref = self.parameter.path_RMSD_final
                cmdstr = f"grep 'IMPOSE  --  After Rotation' RMSD-{c}-start.log > RMSD-{c}-start.dat\n"
                f0.write(cmdstr)
                f0.write("superpose %s -k %s.key %s 1 1 %d N N M N 0.0 > RMSD-%s-final.log \n" %(path_ref,c,c+'.arc',n,c))
                cmdstr = f"grep 'IMPOSE  --  After Rotation' RMSD-{c}-final.log > RMSD-{c}-final.dat"
                f0.write(cmdstr)
            f0.write(self.more_cmd)

        shstr = f"python /home/xy3866/bin_Miles/TinkerGPU2022/submitTinker.py -x {self.job_name}.sh -t GPU -p {self.work_path}"
        with open(scriptfile, 'w') as f0:
            f0.write(shstr)

    def check(self, num_snapshot):
        os.chdir(self.work_path)
        c = self.coord_name
        if not os.path.isfile(f"{c}-results.log"):
            return False
        else:
            cmdstr = f"grep 'Picosecond' {c}-results.log > {c}-time.dat"
            subprocess.run(cmdstr, shell=True)
            nLines = sum(1 for line in open(f"{c}-time.dat"))
            if(nLines == num_snapshot):
                return True
            else:
                return False

    def RMSD_check(self, num_snapshot):
        os.chdir(self.work_path)
        c = self.coord_name
        if not os.path.isfile(f"RMSD-{c}-start.log"):
            return False
        elif not os.path.isfile(f"RMSD-{c}-final.log"):
            return False
        else:
            if not os.path.isfile(f"RMSD-{c}-start.dat"):
                return False
            nLines_start = sum(1 for line in open(f"RMSD-{c}-start.dat"))
            if(nLines_start != num_snapshot):
                return False

            if not os.path.isfile(f"RMSD-{c}-final.dat"):
                return False
            nLines_final = sum(1 for line in open(f"RMSD-{c}-final.dat"))
            if(nLines_final != num_snapshot):
                return False

            if(nLines_start == num_snapshot and nLines_final == num_snapshot):
                return True
            else:
                return False

