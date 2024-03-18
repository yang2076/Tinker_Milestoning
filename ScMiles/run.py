#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
                path_ref_0 = self.parameter.path_RMSD_start
                path_ref_1 = self.parameter.path_RMSD_final
                mol = c.split('-')[0] + ".xyz"
                n = self.parameter.RMSD_solute_num
                path_analyze = os.path.join(self.parameter.path_program, "analyze.py")
                f0.write("python %s %s %s RMSD %s %s %d > RMSD-%s.dat\n" %(path_analyze,c+'.arc',mol,path_ref_0,path_ref_1,n,c))
            elif(self.parameter.restrain == "dihedral" and self.job_name == c):
                mol = c.split('-')[0] + ".xyz"
                path_analyze = os.path.join(self.parameter.path_program, "analyze.py")
                cmd_ = "python %s %s %s dihedral " %(path_analyze,c+'.arc',mol)
                grp_ = self.parameter.restrain_grp
                cmd_ += "%d %d %d %d > dih-%s.dat\n" %(grp_[0], grp_[1], grp_[2], grp_[3], c)
                f0.write(cmd_)
            elif(self.parameter.restrain == "2Ddihedral" and self.job_name == c):
                mol = c.split('-')[0] + ".xyz"
                path_analyze = os.path.join(self.parameter.path_program, "analyze.py")
                cmd_ = "python %s %s %s 2Ddihedral " %(path_analyze,c+'.arc',mol)
                grp_ = self.parameter.restrain_grp
                cmd_ += "%d %d %d %d %d %d %d %d > dih2D-%s.dat\n" %(grp_[0], grp_[1], grp_[2], grp_[3], grp_[4], grp_[5], grp_[6], grp_[7], c)
                f0.write(cmd_)
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
            cmdstr = f"tail -n1 {c}-time.dat"
            f_ = subprocess.Popen(cmdstr,shell=True, encoding="utf8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            f_ = f_.stdout.read()
            if(nLines == num_snapshot):
                return True
            elif("Binary" in f_):
                cmdstr = f"tail -n1 {c}-results.log"
                f_ = subprocess.Popen(cmdstr,shell=True, encoding="utf8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_ = f_.stdout.read().split()
                if(len(f_) >= 1 and f_[0] == "Atoms"):
                    log(f"Binary Problem in {self.job_name}", self.parameter.path_log)
                    err_path = job.work_path + f"/{c}.err"
                    with open(err_path, 'w') as f0:
                        f0.write("Binary problem")
                return False
            else:
                return False

    def structure_check(self, num_snapshot, restrn):
        os.chdir(self.work_path)
        c = self.coord_name
        structure_file = ""
        if(restrn == "RMSD"):
           structure_file = f"RMSD-{c}.dat" 
        elif(restrn == "dihedral"):
           structure_file = f"dih-{c}.dat" 
        elif(restrn == "2Ddihedral"):
           structure_file = f"dih2D-{c}.dat" 
        elif(restrn == "distance"):
           structure_file = f"dist-{c}.dat" 
        elif(restrn == "angle"):
           structure_file = f"angle-{c}.dat" 

        if not os.path.isfile(structure_file):
            return False
        else:
            RMSD_file = open(structure_file)
            lines = RMSD_file.readlines()
            RMSD_file.close()
            nLines_ = len(lines)
            if(nLines_ == num_snapshot):
                return True
            else:
                return False
