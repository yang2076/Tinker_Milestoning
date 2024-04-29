#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Running free trajectories.
"""

import os, time, re, sys
from run import *
from parameters import *
from datetime import datetime
from log import log
from shutil import copyfile
from compute import *
import inspect 
from itertools import combinations
import numpy as np
import math

class traj:
    def __init__(self, traj_name, start_MS, random_num, MS_dir):
        self.traj_name = traj_name
        self.start_MS = start_MS
        self.random_num = random_num
        self.MS_dir = MS_dir
        self.MS_visit = []
        self.MS_visit.append(start_MS)
        f = open(self.MS_dir + '/' + self.traj_name + '.xyz')
        self.length = len(f.readlines())
        self.time = 0 # (ps)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return
    
    def __repr__(self) -> str:
        return ('Free trajectories.')    

    def time_collect(self, time):
        self.time += time

    def create_job(self, parameter):
        '''launch free trajectories'''
        steps = int(parameter.sampling_time*1000 / parameter.timeFactor)
        job = run(parameter, self.MS_dir, self.traj_name, self.traj_name, steps, parameter.saveFrequency)
        return job

    def restart(self, job, is_continue=True):
        if(is_continue):
            self.time_collect(job.steps*job.parameter.timeFactor/1000.)
        if(self.time < job.parameter.max_md_time*1000):
            job.submit()
        else:
            log(f"Error: Exceed the max md time for job {self.traj_name} in {self.MS_dir}", job.parameter.path_log)
            #sys.exit(0)

    def finish_iter(self, MS_name, index, ratio, save_frequency):
        self.time_collect(save_frequency*(index+2.-ratio))
        self.MS_visit.append(MS_name)

    def prep_for_next_iter(self):
        self.time = 0
        self.start_MS = self.MS_visit[1]
        self.MS_visit = []
        self.MS_visit.append(self.start_MS)

class trajPool: 
    def __init__(self, parameter, traj_dict):
        self.parameter = parameter
        self.traj_dict = traj_dict
        self.smooth_K_list = []
        self.smooth_list_limit = self.parameter.smooth_list_limit
        
    def __enter__(self):
        return self
               
    def __exit__(self, exc_type, exc_value, traceback):
        return 
    
    def __repr__(self) -> str:
        return ('Current reached by {} milestones'
                .format(self.count_ms()))
        
    def launch(self):
        log("Iteration # {}".format(self.parameter.iteration), self.parameter.path_log)
        job_list = []

        for name_, traj_ in self.traj_dict.items():
            log(f"Start job {name_}", self.parameter.path_log)
            job = traj_.create_job(self.parameter)
            job_list.append(job)

        return job_list

    def count_ms(self):
        MS_list = self.parameter.MS_list
        dim = len(MS_list)
        k_count = np.zeros((dim, dim))
        for traj_ in self.traj_dict.values():
            i = MS_list.index(traj_.start_MS)
            j = MS_list.index(traj_.MS_visit[1])
            if k_count[i, j] == 0:
                k_count[i, j] = 1
            else:
                k_count[i, j] += 1
        return k_count

    def count_time(self):
        MS_list = self.parameter.MS_list
        dim = len(MS_list)
        time_list = []
        for i in range(dim):
            time_list.append([])
        time_vec = []
        time_std = []
        for traj_ in self.traj_dict.values():
            i = MS_list.index(traj_.start_MS)
            if(type(traj_.time) == np.ndarray):
                time_list[i].append(traj_.time[0])
                traj_.time = traj_.time[0]
            else:
                time_list[i].append(traj_.time)
        for i in range(dim):
            time_list_ = np.array(time_list[i])
            time_vec.append(np.average(time_list_))
            time_std.append(np.std(time_list_))
        return np.array(time_vec), np.array(time_std)

    def smooth_list(self, data):
        if(len(self.smooth_K_list) == self.smooth_list_limit):
            self.smooth_K_list.pop(0)
        new_K = data.copy()
        dim = len(new_K)
        self.smooth_K_list.append(new_K)
        tmp = np.zeros((dim,dim))
        for j in self.smooth_K_list:
            tmp += j
        return tmp/len(self.smooth_K_list)

    def compute_(self):
        k_count = self.count_ms()
        k_n = k_average(k_count)
        self.parameter.kij = self.smooth_list(k_n)
        np.set_printoptions(threshold=np.inf)
        print("K:")
        print(self.parameter.kij)

        self.parameter.flux = flux(self.parameter.kij)
        print("flux")
        print(self.parameter.flux)

        t, t_std = self.count_time()
        t[np.isnan(t)] = 1e-20
        t_std[np.isnan(t_std)] = 1e-20
        print("t")
        print(t)
        print(t_std)

        self.parameter.prob = prob(self.parameter.flux, t)
        print("prob")
        print(self.parameter.prob)

        self.parameter.free_energy = free_energy(self.parameter.prob)
        self.parameter.MFPT_1 = MFPT(self.parameter.kij, t, self.parameter.reactant_milestone, self.parameter.product_milestone)
        self.parameter.MFPT_1_rev = MFPT(self.parameter.kij, t, self.parameter.product_milestone, self.parameter.reactant_milestone)
        self.parameter.MFPT_2 = MFPT2(self.parameter.kij, t, self.parameter.reactant_milestone, self.parameter.product_milestone, self.parameter)
        self.parameter.MFPT_2_rev = MFPT2(self.parameter.kij, t, self.parameter.product_milestone, self.parameter.reactant_milestone, self.parameter)
        dim = len(t)

        energy_samples = []
        MFPT_samples = []
        MFPT2_samples = []
        MFPT_samples_rev = []
        MFPT2_samples_rev = []
        for i in range(self.parameter.err_sampling):
            k_err = k_error(k_count)
            t_err = t_error(t, t_std)
            q_temp = flux(k_err)
            p_temp = prob(q_temp, t)
            energy_samples.append(free_energy(p_temp))
            MFPT_er = MFPT(k_err, t_err, self.parameter.reactant_milestone, self.parameter.product_milestone)
            MFPT_samples.append(MFPT_er)
            MFPT_er = MFPT(k_err, t_err, self.parameter.product_milestone, self.parameter.reactant_milestone)
            MFPT_samples_rev.append(MFPT_er)
            MFPT_er2 = MFPT2(k_err, t_err, self.parameter.reactant_milestone, self.parameter.product_milestone, self.parameter)
            MFPT2_samples.append(MFPT_er2)
            MFPT_er2 = MFPT2(k_err, t_err, self.parameter.product_milestone, self.parameter.reactant_milestone, self.parameter)
            MFPT2_samples_rev.append(MFPT_er2)

        self.parameter.energy_err = []
        for i in range(dim):
            self.parameter.energy_err.append(np.std(np.array(energy_samples)[:,i], ddof=1))
        self.parameter.MFPT_err = float(np.std(MFPT_samples, ddof=1))
        self.parameter.MFPT_err_rev = float(np.std(MFPT_samples_rev, ddof=1))
        self.parameter.MFPT_err2 = float(np.std(MFPT2_samples, ddof=1))
        self.parameter.MFPT_err2_rev = float(np.std(MFPT2_samples_rev, ddof=1))

        self.parameter.c = committor(self.parameter, self.parameter.kij)

    def prep_for_next_iter(self):
        for name_, traj_ in self.traj_dict.items():
            traj_.prep_for_next_iter()

