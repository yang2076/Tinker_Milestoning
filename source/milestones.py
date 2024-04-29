#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:35:54 2019

@author: Wei Wei

This subroutine stores the milestone information.
It initializes the milestone list, also contains the function that provides initial and final milestone.

Modified on Mon Feb 13 2023

@author: Xudong Yang

Designed for Tinker-compatible interface

"""

import os, re, time
import pandas as pd
import numpy as np
import sympy
from run import *
from keyfile import *
from log import log
from parameters import *
from network_check import *
from traj import *


class milestones: 

    def __init__(self, parameter) -> None:
        self.parameter = parameter
        self.possible_MS = []
        self.possible_MS_property = {}
        self.seek_MS = []
        self.anchor_nLines = []
        self.traj_dict = {}
        self.trajPool_ = None

    def __enter__(self):
        return self
               
    def __exit__(self, exc_type, exc_value, traceback):
        return 
            
    def __repr__(self) -> str:
        return ('miletstones details.'
                .format(self.anchor_orig, self.anchor_dest, self.lifetime))

    def __index_modifying(self, index, max_digit):
        index_ = str(index)
        length = len(index_)
        z = '0'
        zero_num = max_digit - length if max_digit >= length else 0
        z = z*zero_num
        return z + index_


    def __random_num_assign(self, random_num_exist):
        i = 123456789
        while True:
            i = sympy.randprime(10000000, 999999999)
            if(i not in random_num_exist):
                break
        return i

    def __get_next_frame_num(self, struPath):
        import subprocess
        next_frame = 1
        random_num_exist = []
        while True:
            tmp = self.__index_modifying(next_frame, 4)
            pdbPath = struPath + '/' + 'traj' + tmp + '.key'
            if os.path.exists(pdbPath):
                next_frame += 1
                f = subprocess.Popen(f"grep 'randomseed' {pdbPath}",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
                f_ = f.stdout.read()
                if(next_frame != 2 and f_ == b''):
                    log("Error: No random seed is given in existing key files", self.parameter.path_log)
                    sys.exit(0)
                out = int(f_.split()[1])
                random_num_exist.append(out)
            else:
                break
        return next_frame, random_num_exist
    
    def __check_dihedral_pbc(self, point):
        point_ = point
        if(point > 180.0):
            point_ -= 360.0
        elif(point < -180.0):
            point_ += 360.0
        return point_

    def __cross_point(self, p1, p2, p3, p4):
        x = y = 0
        x1,y1,x2,y2 = p1[0], p1[1], p2[0], p2[1] 
        x3,y3,x4,y4 = p3[0], p3[1], p4[0], p4[1]
        p31 = np.array(p1) - np.array(p3)
        p34 = np.array(p4) - np.array(p3)
        p32 = np.array(p2) - np.array(p3)
        p23 = - p32
        p21 = np.array(p1) - np.array(p2)
        p24 = np.array(p4) - np.array(p2)
    
        if(max(x1,x2)>=min(x3,x4)
        and max(x3,x4)>=min(x1,x2)
        and max(y1,y2)>=min(y3,y4)
        and max(y3,y4)>=min(y1,y2)
        and np.dot(np.cross(p31,p34),np.cross(p32,p34))<=0
        and np.dot(np.cross(p23,p21),np.cross(p24,p21))<=0):
            if (x2 - x1 == 0.):
                k1 = None
                b1 = 0
            else:
                k1 = (y2 - y1) * 1.0 / (x2 - x1)
                b1 = y1 * 1.0 - x1 * k1 * 1.0
        
            if (x4 - x3) == 0:
                k2 = None
                b2 = 0
            else:
                k2 = (y4 - y3) * 1.0 / (x4 - x3)
                b2 = y3 * 1.0 - x3 * k2 * 1.0
        
            if k1 is None:
                if not k2 is None:
                    x = x1
                    y = k2 * x1 + b2
            elif k2 is None:
                x = x3
                y = k1 * x3 + b1
            elif not k2 == k1:
                x = (b2 - b1) * 1.0 / (k1 - k2)
                y = k1 * x * 1.0 + b1 * 1.0
        return np.array([x, y])

    def __MS_property_analysis(self, anchor1, anchor2, pbc):
        # For 1D dihedral & RMSD MS pre-analysis
        restrn = self.parameter.restrain
        if(restrn == "dihedral"):
            # Any dihedral value should be inside -180~180
            an_1 = self.__check_dihedral_pbc(self.parameter.anchors[anchor1-1])
            an_2 = self.__check_dihedral_pbc(self.parameter.anchors[anchor2-1])
            if not pbc:
                return (an_1 + an_2)/2.0
            else:
                tmp = (an_1 + an_2 + 360.0)/2.0
                return self.__check_dihedral_pbc(tmp)
        elif(restrn == "RMSD"):
            return [anchor1, anchor2]

    def __RMSD_reactant_product(self, MS_list_):
        if(len(self.parameter.reactant) == 6):
            if((self.parameter.reactant[0], self.parameter.reactant[1]) in MS_list_):
                self.parameter.reactant = [self.parameter.reactant[0], self.parameter.reactant[1]]
            elif((self.parameter.reactant[2], self.parameter.reactant[3]) in MS_list_):
                self.parameter.reactant = [self.parameter.reactant[2], self.parameter.reactant[3]]
            elif((self.parameter.reactant[4], self.parameter.reactant[5]) in MS_list_):
                self.parameter.reactant = [self.parameter.reactant[4], self.parameter.reactant[5]]
            else:
                log("No reactant is found!\n", self.parameter.path_log)
                log(f"Check if there is MS{self.parameter.reactant[0]}_{self.parameter.reactant[1]}", self.parameter.path_log)
                log(f" MS{self.parameter.reactant[2]}_{self.parameter.reactant[3]}", self.parameter.path_log)
                log(f" MS{self.parameter.reactant[4]}_{self.parameter.reactant[5]}\n", self.parameter.path_log)

        if(len(self.parameter.product) == 6):
            if((self.parameter.product[0], self.parameter.product[1]) in MS_list_):
                self.parameter.product = [self.parameter.product[0], self.parameter.product[1]]
            elif((self.parameter.product[2], self.parameter.product[3]) in MS_list_):
                self.parameter.product = [self.parameter.product[2], self.parameter.product[3]]
            elif((self.parameter.product[4], self.parameter.product[5]) in MS_list_):
                self.parameter.product = [self.parameter.product[4], self.parameter.product[5]]
            else:
                log("No product is found!\n", self.parameter.path_log)
                log(f"Check if there is MS{self.parameter.product[0]}_{self.parameter.product[1]}", self.parameter.path_log)
                log(f" MS{self.parameter.product[2]}_{self.parameter.product[3]}", self.parameter.path_log)
                log(f" MS{self.parameter.product[4]}_{self.parameter.product[5]}\n", self.parameter.path_log)

    def __Voronoi_tessellation(self):
        # For 2D dihedral voronoi tessellation only
        from scipy.spatial import Voronoi
        import numpy as np
        right_ext = np.hstack((np.ones((self.parameter.AnchorNum,1))*360.0,np.zeros((self.parameter.AnchorNum,1))))
        upper_ext = np.hstack((np.zeros((self.parameter.AnchorNum,1)),np.ones((self.parameter.AnchorNum,1))*360.0))
        points = self.parameter.anchors
        points_copy1 = right_ext + points
        points_copy2 = upper_ext + points
        points_copy3 = right_ext + upper_ext + points
        points_final = np.vstack((points, points_copy1, points_copy2, points_copy3))
        vor = Voronoi(points_final)
        vertices = list(vor.vertices)
        ridge_points = list(vor.ridge_points)
        ridge_vertices = list(vor.ridge_vertices)
        for k, i in enumerate(ridge_points):
            pt_pbc_1 = [0., 0.]
            pt_pbc_2 = [0., 0.]
            anchor_1 = min(i[0]%self.parameter.AnchorNum, i[1]%self.parameter.AnchorNum)+1
            anchor_2 = max(i[0]%self.parameter.AnchorNum, i[1]%self.parameter.AnchorNum)+1
            index_1 = ridge_vertices[k][0]
            index_2 = ridge_vertices[k][1]
            if(index_1 == -1 or index_2 == -1):
                continue
            elif((anchor_1, anchor_2) not in self.possible_MS):
                ver_1 = vertices[index_1]
                ver_2 = vertices[index_2]
                self.possible_MS.append((anchor_1,anchor_2))
                pt_pbc_1 = self.__cross_point(ver_1, ver_2, [180.,-180.],[180.,180.])
                pt_pbc_2 = self.__cross_point(ver_1, ver_2, [-180.,180.],[180.,180.])
                if(ver_1[0] < ver_2[0]):
                    ver_1[0] = self.__check_dihedral_pbc(ver_1[0])
                    ver_1[1] = self.__check_dihedral_pbc(ver_1[1])
                    ver_2[0] = self.__check_dihedral_pbc(ver_2[0])
                    ver_2[1] = self.__check_dihedral_pbc(ver_2[1])
                    self.possible_MS_property[(anchor_1,anchor_2)] = [ver_1, ver_2, pt_pbc_1, pt_pbc_2]
                elif(ver_1[0] > ver_2[0]):
                    ver_1[0] = self.__check_dihedral_pbc(ver_1[0])
                    ver_1[1] = self.__check_dihedral_pbc(ver_1[1])
                    ver_2[0] = self.__check_dihedral_pbc(ver_2[0])
                    ver_2[1] = self.__check_dihedral_pbc(ver_2[1])
                    self.possible_MS_property[(anchor_1,anchor_2)] = [ver_2, ver_1, pt_pbc_1, pt_pbc_2]
                elif(ver_1[0] == ver_2[0] and ver_1[1] <= ver_2[1]):
                    ver_1[0] = self.__check_dihedral_pbc(ver_1[0])
                    ver_1[1] = self.__check_dihedral_pbc(ver_1[1])
                    ver_2[0] = self.__check_dihedral_pbc(ver_2[0])
                    ver_2[1] = self.__check_dihedral_pbc(ver_2[1])
                    self.possible_MS_property[(anchor_1,anchor_2)] = [ver_1, ver_2, pt_pbc_1, pt_pbc_2]
                elif(ver_1[0] == ver_2[0] and ver_1[1] > ver_2[1]):
                    ver_1[0] = self.__check_dihedral_pbc(ver_1[0])
                    ver_1[1] = self.__check_dihedral_pbc(ver_1[1])
                    ver_2[0] = self.__check_dihedral_pbc(ver_2[0])
                    ver_2[1] = self.__check_dihedral_pbc(ver_2[1])
                    self.possible_MS_property[(anchor_1,anchor_2)] = [ver_2, ver_1, pt_pbc_1, pt_pbc_2]

    def __read_tinker_arc(self, arc_name, frame_id):
        arc = arc_name.split('.')[0] + '.arc'
        with open(arc) as f:
            frame_0 = f.readline()
            num_of_atoms = int(frame_0)
            for i, l in enumerate(f):
                frame_0 += l
                if i == num_of_atoms:
                    break
            if frame_id == 0:
                return frame_0
            frame_byte_size = len(frame_0.encode())
            f.seek(frame_byte_size * frame_id, 0)
            frame_i = f.read(frame_byte_size)
            frame_i = frame_i.strip().splitlines()
            tmp_ = ['\n']
            tmp_ = tmp_ * (num_of_atoms+2)
            frame_i = list(map(lambda x, y: x + y, frame_i, tmp_))
        return frame_i

    def __hit_MS_dih(self, dih):
        for i in self.possible_MS_property.keys():
            if(np.abs(self.possible_MS_property[i] - dih) <= self.parameter.MS_threshold):
                return i
            elif(self.parameter.MS_threshold + self.possible_MS_property[i] > 180.):
                if(np.abs(self.possible_MS_property[i] - dih - 360.0) <= self.parameter.MS_threshold):
                    return i
            elif(self.possible_MS_property[i] - self.parameter.MS_threshold < -180.):
                if(np.abs(self.possible_MS_property[i] - dih + 360.0) <= self.parameter.MS_threshold):
                    return i
        return None

    def __avg_config(self, config_1, config_2, ratio):
        config = []
        for i, j in enumerate(config_1):
            if(i == 0 or i == 1):
                config.append(j)
            else:
                terms_1 = j.split()
                if(len(terms_1) == 1):
                    print(terms_1)
                    continue
                terms_2 = config_2[i].split()
                i_ = int(terms_1[0])
                a_ = terms_1[1]
                x_ = float(terms_1[2]) * ratio + float(terms_2[2]) * (1-ratio)
                y_ = float(terms_1[3]) * ratio + float(terms_2[3]) * (1-ratio)
                z_ = float(terms_1[4]) * ratio + float(terms_2[4]) * (1-ratio)
                t_ = int(terms_1[5])
                if(len(terms_1) == 7):
                    c_1 = int(terms_1[6])
                    config.append(" %6d %3s %12.6f %12.6f %12.6f %6d %6d\n" %(i_,a_,x_,y_,z_,t_,c_1))
                elif(len(terms_1) == 8):
                    c_1 = int(terms_1[6])
                    c_2 = int(terms_1[7])
                    config.append(" %6d %3s %12.6f %12.6f %12.6f %6d %6d %6d\n" %(i_,a_,x_,y_,z_,t_,c_1,c_2))
                elif(len(terms_1) == 9):
                    c_1 = int(terms_1[6])
                    c_2 = int(terms_1[7])
                    c_3 = int(terms_1[8])
                    config.append(" %6d %3s %12.6f %12.6f %12.6f %6d %6d %6d %6d\n" %(i_,a_,x_,y_,z_,t_,c_1,c_2,c_3))
                elif(len(terms_1) == 10):
                    c_1 = int(terms_1[6])
                    c_2 = int(terms_1[7])
                    c_3 = int(terms_1[8])
                    c_4 = int(terms_1[9])
                    config.append(" %6d %3s %12.6f %12.6f %12.6f %6d %6d %6d %6d %6d\n" %(i_,a_,x_,y_,z_,t_,c_1,c_2,c_3, c_4))
        return config

    def __write2txyz(self, txyz_path, config_lines):
        with open(txyz_path, 'w') as f:
            for line in config_lines:
                f.write(line)

    def __hit_MS_dih_2(self, dih_list):
        for i in range(len(dih_list)-1):
            first_index = min(dih_list[i], dih_list[i+1])
            second_index = max(dih_list[i], dih_list[i+1])
            if(second_index - first_index < 180.0):
                for j in self.possible_MS_property.keys():
                    if(self.possible_MS_property[j] >= first_index and self.possible_MS_property[j] <= second_index):
                        dist_1 = np.abs(self.possible_MS_property[j] - dih_list[i])
                        dist_2 = np.abs(self.possible_MS_property[j] - dih_list[i+1])
                        ratio = dist_1 / (dist_1 + dist_2)
                        return i, ratio, j
            else:
                for j in self.possible_MS_property.keys():
                    if(self.possible_MS_property[j] <= first_index or self.possible_MS_property[j] >= second_index):
                        dist_1 = np.abs(self.__check_dihedral_pbc(self.possible_MS_property[j] - dih_list[i]))
                        dist_2 = np.abs(self.__check_dihedral_pbc(self.possible_MS_property[j] - dih_list[i+1]))
                        ratio = dist_1 / (dist_1 + dist_2)
                        return i, ratio, j
        return 0, 0, None

    def __hit_MS_2Ddih(self, dih_list):
        diff = {}
        dih = np.array(dih)
        for j, an in enumerate(self.parameter.anchors):
            diff[j+1] = np.linalg.norm(dih-an)
        tmp = sorted(diff.items(), key=lambda i: i[1])
        if(np.abs(tmp[0][1]-tmp[1][1]) <= self.parameter.MS_threshold):
            if(tmp[0][0] < tmp[1][0]):
                return (tmp[0][0], tmp[1][0])
            else:
                return (tmp[1][0], tmp[0][0])
        return None

    def __hit_MS_RMSD(self, RMSD_start, RMSD_final):
        for i in range(1, self.parameter.RMSD_num*2+1):
            for j in range(1, self.parameter.RMSD_num*2+1):
                tmp_1 = i+(j-1)*(self.parameter.RMSD_num*2+1)
                tmp_2 = i+j*(self.parameter.RMSD_num*2+1)
                hit_1 = hit_2 = hit_3 = hit_4 = False
                if((tmp_1,tmp_1+1) in self.possible_MS):
                    if(RMSD_start >= self.possible_MS_property[(tmp_1,tmp_1+1)][0][0] and
                       RMSD_start <= self.possible_MS_property[(tmp_1,tmp_1+1)][0][1] and
                       RMSD_final >= self.possible_MS_property[(tmp_1,tmp_1+1)][1][0] and
                       RMSD_final <= self.possible_MS_property[(tmp_1,tmp_1+1)][1][1]):
                        hit_1 = True
                if((tmp_1,tmp_2) in self.possible_MS):
                    if(RMSD_start >= self.possible_MS_property[(tmp_1,tmp_2)][0][0] and
                       RMSD_start <= self.possible_MS_property[(tmp_1,tmp_2)][0][1] and
                       RMSD_final >= self.possible_MS_property[(tmp_1,tmp_2)][1][0] and
                       RMSD_final <= self.possible_MS_property[(tmp_1,tmp_2)][1][1]):
                        hit_2 = True
                if((tmp_1+1,tmp_2+1) in self.possible_MS):
                    if(RMSD_start >= self.possible_MS_property[(tmp_1+1,tmp_2+1)][0][0] and
                       RMSD_start <= self.possible_MS_property[(tmp_1+1,tmp_2+1)][0][1] and
                       RMSD_final >= self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][0] and
                       RMSD_final <= self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][1]):
                        hit_3 = True
                if((tmp_2,tmp_2+1) in self.possible_MS):
                    if(RMSD_start >= self.possible_MS_property[(tmp_2,tmp_2+1)][0][0] and
                       RMSD_start <= self.possible_MS_property[(tmp_2,tmp_2+1)][0][1] and
                       RMSD_final >= self.possible_MS_property[(tmp_2,tmp_2+1)][1][0] and
                       RMSD_final <= self.possible_MS_property[(tmp_2,tmp_2+1)][1][1]):
                        hit_4 = True

                if(not hit_1 and not hit_2 and not hit_3 and not hit_4):
                    continue
                elif(hit_1 and not hit_2 and not hit_3 and not hit_4):
                    return (tmp_1,tmp_1+1)
                elif(not hit_1 and hit_2 and not hit_3 and not hit_4):
                    return (tmp_1,tmp_2)
                elif(not hit_1 and not hit_2 and hit_3 and not hit_4):
                    return (tmp_1+1,tmp_2+1)
                elif(not hit_1 and not hit_2 and not hit_3 and hit_4):
                    return (tmp_2,tmp_2+1)
                elif(hit_1 and hit_2 and hit_3 and hit_4):
                    return (tmp_1,tmp_1+1)
                elif(hit_1 and hit_2):
                    x_diff = (self.possible_MS_property[(tmp_1,tmp_1+1)][0][0]+self.possible_MS_property[(tmp_1,tmp_1+1)][0][1])/2.0
                    x_diff = np.abs(RMSD_start - x_diff)
                    y_diff = (self.possible_MS_property[(tmp_1,tmp_2)][1][0]+self.possible_MS_property[(tmp_1,tmp_2)][1][1])/2.0
                    y_diff = np.abs(RMSD_final - y_diff)
                    if(x_diff <= y_diff):
                        return (tmp_1, tmp_1+1)
                    else:
                        return (tmp_1,tmp_2)
                elif(hit_1 and hit_3):
                    x_diff = (self.possible_MS_property[(tmp_1,tmp_1+1)][0][0]+self.possible_MS_property[(tmp_1,tmp_1+1)][0][1])/2.0
                    x_diff = np.abs(RMSD_start - x_diff)
                    y_diff = (self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][0]+self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][1])/2.0
                    y_diff = np.abs(RMSD_final - y_diff)
                    if(x_diff <= y_diff):
                        return (tmp_1, tmp_1+1)
                    else:
                        return (tmp_1+1,tmp_2+1)
                elif(hit_4 and hit_2):
                    x_diff = (self.possible_MS_property[(tmp_2,tmp_2+1)][0][0]+self.possible_MS_property[(tmp_2,tmp_2+1)][0][1])/2.0
                    x_diff = np.abs(RMSD_start - x_diff)
                    y_diff = (self.possible_MS_property[(tmp_1,tmp_2)][1][0]+self.possible_MS_property[(tmp_1,tmp_2)][1][1])/2.0
                    y_diff = np.abs(RMSD_final - y_diff)
                    if(x_diff < y_diff):
                        return (tmp_2,tmp_2+1)
                    else:
                        return (tmp_1,tmp_2)
                elif(hit_4 and hit_3):
                    x_diff = (self.possible_MS_property[(tmp_2,tmp_2+1)][0][0]+self.possible_MS_property[(tmp_2,tmp_2+1)][0][1])/2.0
                    x_diff = np.abs(RMSD_start - x_diff)
                    y_diff = (self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][0]+self.possible_MS_property[(tmp_1+1,tmp_2+1)][1][1])/2.0
                    y_diff = np.abs(RMSD_final - y_diff)
                    if(x_diff < y_diff):
                        return (tmp_2,tmp_2+1)
                    else:
                        return (tmp_1+1,tmp_2+1)
        return None          

    def __whether_hit_new_MS(self, work_path, coord_name, num_snapshot, index, nLines, MS_list_, seek_):
        filePath = self.parameter.work_dir
        crddir = filePath + '/crd'
        os.chdir(work_path)
        restrn = self.parameter.restrain
        fine_hit = True if seek_ else False
        hit_MS = None
        hit_MS_seek = []
        val_1 = []
        val_2 = []
        structure_job = run(self.parameter, work_path, coord_name, coord_name, 0, 0.)
        while True:
            rFinish = structure_job.structure_check(num_snapshot, restrn)
            if rFinish:
                break
            else:
                time.sleep(0.5)

        if(restrn == "dihedral"):
            val_1 = np.loadtxt(f"dih-{coord_name}.dat",usecols=(0,))
        elif(restrn == "2Ddihedral"):
            val_1 = np.loadtxt(f"dih-{coord_name}.dat",usecols=(0,))
            val_2 = np.loadtxt(f"dih-{coord_name}.dat",usecols=(1,))
        elif(restrn == "RMSD"):
            val_1 = np.loadtxt(f"RMSD-{coord_name}.dat",usecols=(0,))
            val_2 = np.loadtxt(f"RMSD-{coord_name}.dat",usecols=(1,))

        for i in range(num_snapshot):
            if(restrn == "dihedral"):
                hit_MS = self.__hit_MS_dih(val_1[i])
            elif(restrn == "2Ddihedral"):
                hit_MS = self.__hit_MS_2Ddih([val_1[i], val_2[i]])
            elif(restrn == "RMSD"):
                hit_MS = self.__hit_MS_RMSD(val_1[i], val_2[i])

            if(hit_MS and hit_MS not in MS_list_):
                MS_list_.append(hit_MS)
                if not seek_:
                    return hit_MS, i, 1
                new_MS_dir = crddir + '/MS' + str(hit_MS[0]) + '_' + str(hit_MS[1])
                os.makedirs(new_MS_dir)
                lines = self.__read_tinker_arc(coord_name, i)
                MS_config_path = new_MS_dir + '/MS' + str(hit_MS[0]) + '_' + str(hit_MS[1]) + '.xyz'
                self.__write2txyz(MS_config_path, lines)
                fine_hit = True
                hit_MS_seek.append(hit_MS)

            # avg for 1D only right now
            if (restrn == "dihedral" and not fine_hit):
                n, ratio, hit_MS = self.__hit_MS_dih_2(val_1)
                if(hit_MS and hit_MS not in MS_list_):
                    MS_list_.append(hit_MS)
                    if not seek_:
                        return hit_MS, n, ratio
                    config_1 = self.__read_tinker_arc(coord_name, n)
                    config_2 = self.__read_tinker_arc(coord_name, n+1)
                    config_ = self.__avg_config(config_1, config_2, ratio)
                    new_MS_dir = crddir + '/MS' + str(hit_MS[0]) + '_' + str(hit_MS[1])
                    os.makedirs(new_MS_dir)
                    MS_config_path = new_MS_dir + '/MS' + str(hit_MS[0]) + '_' + str(hit_MS[1]) + '.xyz'
                    self.__write2txyz(MS_config_path, config_)
                    hit_MS_seek.append(hit_MS)

        if not seek_:
            return None, 0, 0
        return hit_MS_seek

    def __seek_milestones(self):
        from shutil import copy
        filePath = self.parameter.work_dir
        seekdir = filePath + '/crd/seek'
        path_keyfile = filePath + '/free.key'
        log("Ready to seek milestones from the anchors given", self.parameter.path_log)

        if not os.path.exists(seekdir):
            os.makedirs(seekdir)
        
        for an in range(1, self.parameter.AnchorNum+1):
            random_int_list = []
            an_ = str(an)
            path_anchor_txyz = os.path.join(self.parameter.path_anchors, '%s.xyz' %(an_))
            self.anchor_nLines.append(sum(1 for line in open(path_anchor_txyz)))
            an_ = self.__index_modifying(an_, 3)
            log(f"Create the directory of Anchor {an_}", self.parameter.path_log)
            if not os.path.exists(seekdir + '/anchor' + an_):
                os.makedirs(seekdir + '/anchor' + an_)
            seekNum, random_int_list = self.__get_next_frame_num(seekdir + '/anchor' + an_)
            for i in range(self.parameter.seek_traj):
                i_ = str(i + seekNum)
                i_ = self.__index_modifying(i_, 4)
                key_path = seekdir + '/anchor' + an_ + '/anchor' + an_  + '_traj_%s.key' %(i_)
                copy(path_keyfile, key_path)
                copy(path_anchor_txyz, seekdir + '/anchor' + an_ + '/anchor' + an_ + '_traj_%s.xyz' %(i_))
                r = self.__random_num_assign(random_int_list)
                random_int_list.append(r)
                keyfile(self.parameter, free=True, initial=True, key_path=key_path, random_num=r).generate()
            log(f"Finish creating {self.parameter.seek_traj} key files for Anchor {an_}", self.parameter.path_log)

        for i in range(1, self.parameter.seek_traj+1):
            i_ = str(i)
            i_ = self.__index_modifying(i_, 4)
            job_list = []
            num_snapshot = int(self.parameter.seekTime * 1000 / self.parameter.seek_saveFrequency)
            for an in range(1, self.parameter.AnchorNum+1):
                an_ = str(an)
                an_ = self.__index_modifying(an_, 3)
                work_path = seekdir + '/anchor' + an_
                log(f"Begin to run free MD on traj {i_} of Anchor {an_}", self.parameter.path_log)
                steps=int(self.parameter.seekTime*1000*1000 / self.parameter.timeFactor)
                job=run(self.parameter, work_path, f"anchor{an_}_traj_{i_}", 
                        f"anchor{an_}_traj_{i_}", steps, self.parameter.seek_saveFrequency)
                job_list.append(job)
            for job in job_list:
                job.submit()
            Finish_num = 0
            len_job_list = len(job_list)
            Finish_job_list = []
            while True:
                for job in job_list:
                    if (job.job_name in Finish_job_list):
                        continue
                    iFinish = job.check(num_snapshot)
                    if (iFinish):
                        Finish_num += 1
                        Finish_job_list.append(job.job_name)
                        index = int(job.job_name[6:9]) - 1
                        hit_MS_list = self.__whether_hit_new_MS(job.work_path, job.coord_name, num_snapshot, index, self.anchor_nLines, self.seek_MS, True)
                        log(f"Job {job.job_name} has been finished.", self.parameter.path_log)
                        log(f"There are {Finish_num}/{len_job_list} jobs completed.", self.parameter.path_log)
                        if(len(hit_MS_list) > 0):
                            for ii in hit_MS_list:
                                log(f"MS{ii[0]}_{ii[1]} is found.", self.parameter.path_log)
                        else:
                            log(f"No new MS is found.", self.parameter.path_log)
                        tmp_1 = job.work_path + '/' + job.job_name + '.arc'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + '/' + job.job_name + '.dyn'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + '/' + job.job_name + '-results.log'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        if(self.parameter.restrain == "RMSD"):
                            tmp_1 = job.work_path + '/' + 'RMSD-' + job.job_name + '.dat'
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                    else:
                        time.sleep(1.0)
                if(Finish_num == len_job_list):
                    log(f"Round {i_} of seeking completed", self.parameter.path_log)
                    break
            if(len(self.seek_MS) == len(self.possible_MS)):
                break

    def __read_milestone_folder(self):
        filePath = self.parameter.work_dir 
        crddir = filePath + '/crd'
        restrn = self.parameter.restrain
        if(restrn == "dihedral"):
            range_1 = self.parameter.AnchorNum
            range_2 = self.parameter.AnchorNum + 1
            for i in range(1, range_1):
                for j in range(i, range_2):
                    name = 'MS' + str(i) + '_' + str(j)
                    if os.path.exists(crddir + '/' + name):
                        self.parameter.MS_list.append((i,j))
                        self.parameter.MS_property[(i,j)] = self.possible_MS_property[(i,j)]
        elif(restrn == "RMSD"):
            range_1 = (self.parameter.RMSD_num*2+1)**2
            for i in range(1, range_1):
                name = 'MS' + str(i) + '_' + str(i+1)
                if os.path.exists(crddir + '/' + name):
                    j = i + 1
                    self.parameter.MS_list.append((i,j))
                    self.parameter.MS_property[(i,j)] = self.possible_MS_property[(i,j)]
                name = 'MS' + str(i) + '_' + str(i+self.parameter.RMSD_num*2+1)
                if os.path.exists(crddir + '/' + name):
                    j = i + self.parameter.RMSD_num * 2 + 1
                    self.parameter.MS_list.append((i,j))
                    self.parameter.MS_property[(i,j)] = self.possible_MS_property[(i,j)]

    def __restrain_MS_sampling(self):
        from shutil import copy
        import subprocess
        job_list = []
        filePath = self.parameter.work_dir 
        crddir = filePath + '/crd'
        MS_list = self.parameter.MS_list.copy()
        finished = self.parameter.Finished.copy()
        steps=int(self.parameter.restrain_md_time*1000*1000 / self.parameter.timeFactor)
        num_snapshot = int(self.parameter.restrain_md_time*1000/self.parameter.res_saveFrequency)
        eq_num = int(self.parameter.restrain_eq_time*1000/self.parameter.res_saveFrequency)
        possible_config_Num = int((num_snapshot - eq_num) / self.parameter.interval)
        for i in MS_list:
            if(i in finished or i in self.parameter.finished_restrain):
                log(f"Restrained MD on MS{i[0]}_{i[1]} already completed", self.parameter.path_log)
                continue
            MS_dir = crddir + '/MS' + str(i[0]) + '_' + str(i[1])
            key_path = MS_dir + '/r_MD-MS%d_%d.key' %(i[0],i[1])
            path_keyfile = filePath + '/restrain.key'
            path_config_py = os.path.join(self.parameter.path_program, 'read_config.py')
            copy(path_keyfile, key_path)
            copy(path_config_py, MS_dir)
            r = self.__random_num_assign([])
            keyfile(self.parameter, anchor=i, free=False, initial=True, key_path=key_path, random_num=r).generate()
            if(possible_config_Num < self.parameter.trajPerLaunch):
                log(f"Error: No enough configurations generated by restrained MD for unbiased sampling!", self.parameter.path_log)
                sys.exit(0)
            job=run(self.parameter, MS_dir, f"r_MD-MS{i[0]}_{i[1]}", 
                    f"MS{i[0]}_{i[1]}", steps, self.parameter.res_saveFrequency)
            log(f"Begin to run Restrained MD on MS{i[0]}_{i[1]}", self.parameter.path_log)
            job_list.append(job)

        for job in job_list:
            job.submit()

        len_job_list = len(job_list)
        Finish_Num = 0
        Config_Num = 0

        while True:
            for job in job_list:
                iFinish = job.check(num_snapshot)
                i1 = int(job.coord_name.split('_')[0][2:])
                i2 = int(job.coord_name.split('_')[1])
                anchor_name = (i1, i2)
                if(iFinish and anchor_name not in self.parameter.finished_restrain):
                    self.parameter.finished_restrain.append(anchor_name)
                    MS_dir = crddir + '/' + job.coord_name
                    Finish_Num += 1
                    index = 1
                    timestr = str(time.time()).replace('.', '')
                    scriptfile = os.path.join(self.parameter.path_jobsubmit, f"read_{timestr}.sh")
                    shfile = os.path.join(job.work_path, f"read_traj.sh")
                    end_config = eq_num + self.parameter.trajPerLaunch * self.parameter.interval
                    with open(shfile, 'w') as f0:
                        f0.write(f"python read_config.py {job.coord_name} {crddir} {eq_num} {end_config} {self.parameter.interval}")
                    shstr = f"python /home/xy3866/bin_Miles/TinkerGPU2022/submitTinker.py -x read_traj.sh -t CPU -p {MS_dir}"
                    with open(scriptfile, 'w') as f0:
                        f0.write(shstr)
                    for j in range(eq_num, num_snapshot, self.parameter.interval):
                        random_int_list = []
                        index_ = self.__index_modifying(index, 4)
                        key_path = MS_dir + f"/{job.coord_name}-traj_{index_}.key"
                        path_keyfile = filePath + '/restrain.key'
                        copy(path_keyfile, key_path)
                        r = self.__random_num_assign(random_int_list)
                        random_int_list.append(r)
                        keyfile(self.parameter, anchor=anchor_name, free=True, initial=True, key_path=key_path, random_num=r).generate()
                        index += 1
                        if(index > self.parameter.trajPerLaunch):
                            break
                    Config_Num = index-1
                    log(f"Complete Restrained MD {job.coord_name}", self.parameter.path_log)
                else:
                    time.sleep(1.0)
            if(Finish_Num == len_job_list):
                log("Restrained MD completed!", self.parameter.path_log)
                break

        while True:
            Finish_reading = True
            for i in MS_list:
                MS_dir = crddir + '/MS' + str(i[0]) + '_' + str(i[1])
                os.chdir(MS_dir)
                for j in range(Config_Num):
                    index_ = self.__index_modifying(j+1, 4)
                    path_traj = os.path.join(MS_dir, f"MS{str(i[0])}_{str(i[1])}-traj_{index_}.xyz")
                    if(os.path.exists(path_traj)):
                        f = subprocess.Popen(f"wc -l MS{str(i[0])}_{str(i[1])}-traj_{index_}.xyz",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
                        f_ = int(f.stdout.read().split()[0])
                        f = subprocess.Popen(f"wc -l MS{str(i[0])}_{str(i[1])}.xyz",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
                        nLines = int(f.stdout.read().split()[0])
                        if(f_ != nLines):
                            Finish_reading = False
                            break
                    else:
                        Finish_reading = False
                        break
            if(Finish_reading):
                log(f"Configurations for all MS have been chosen. Ready for production MD", self.parameter.path_log)
                break
            else:
                time.sleep(1.0)
        
        for i in MS_list:
            MS_dir = crddir + '/MS' + str(i[0]) + '_' + str(i[1])
            os.chdir(MS_dir)
            for j in range(Config_Num):
                index_ = self.__index_modifying(j+1, 4)
                traj_ = traj(f"MS{i[0]}_{i[1]}-traj_{index_}", i, r, MS_dir)
                self.traj_dict[f"MS{i[0]}_{i[1]}-traj_{index_}"] = traj_

            if(os.path.isfile(MS_dir + '/MS' + str(i[0]) + '_' + str(i[1]) + '.dyn')):
                os.remove(MS_dir + '/MS' + str(i[0]) + '_' + str(i[1]) + '.dyn')
            if(os.path.isfile(MS_dir + '/MS' + str(i[0]) + '_' + str(i[1]) + '.arc')):
                os.remove(MS_dir + '/MS' + str(i[0]) + '_' + str(i[1]) + '.arc')
            if(os.path.isfile(MS_dir + '/r_MD-MS' + str(i[0]) + '_' + str(i[1]) + '.key')):
                os.remove(MS_dir + '/r_MD-MS' + str(i[0]) + '_' + str(i[1]) + '.key')

    def __final_sampling(self):
        from shutil import copy
        job_list = self.trajPool_.launch()
        for job in job_list:
            job.submit()

        steps = int(self.parameter.sampling_time*1000 / self.parameter.timeFactor)
        num_snapshot = int(self.parameter.sampling_time / self.parameter.saveFrequency)
        len_job_list = len(job_list)
        remain_job_num = len_job_list
        Finish_num = 0
        Finish_job_list = []
        while True:
            iFinish = False
            for job in job_list:
                if (job.job_name in Finish_job_list):
                    continue
                iFinish = job.check(num_snapshot)

                i = job.job_name
                err_path = job.work_path + f"/{self.traj_dict[i].traj_name}.err"
                if os.path.isfile(err_path):
                    log(f"Problem in {job.work_path}/{job.coord_name}", self.parameter.path_log)
                    os.remove(err_path)
                    tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-results.log'
                    if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                    tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-time.dat'
                    if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                    tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.arc'
                    if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                    tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.dyn'
                    if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                    key_path_ = job.work_path + f'/{self.traj_dict[i].traj_name}.key'
                    copy(self.parameter.work_dir + '/free.key', job.work_path + f'/{self.traj_dict[i].traj_name}.key')
                    r = self.__random_num_assign([])
                    keyfile(self.parameter, anchor=i, free=True, initial=True, key_path=key_path_, random_num=r).generate()
                    log(f"Restart at {job.work_path}/{job.coord_name}", self.parameter.path_log)
                    self.traj_dict[i].restart(job, False)
                    continue

                if (iFinish):
                    MS_list_ = [self.traj_dict[i].start_MS]
                    hit_MS, index, ratio = self.__whether_hit_new_MS(job.work_path, job.coord_name, num_snapshot, 0, [self.traj_dict[i].length], MS_list_, False)
                    os.chdir(job.work_path)
                    txyz_path = self.traj_dict[i].MS_dir + f"/{self.traj_dict[i].traj_name}.xyz"
                    length = self.traj_dict[i].length
                    if(hit_MS):
                        if(ratio != 1):
                            #config_1 = self.__read_tinker_arc(self.traj_dict[i].traj_name, index)
                            #config_2 = self.__read_tinker_arc(self.traj_dict[i].traj_name, index+1)
                            #config_ = self.__avg_config(config_1, config_2, ratio)
                            #self.__write2txyz(txyz_path, config_)
                            config_ = self.__read_tinker_arc(self.traj_dict[i].traj_name, index)
                            self.__write2txyz(txyz_path, config_)
                        else:
                            config_ = self.__read_tinker_arc(self.traj_dict[i].traj_name, index)
                            self.__write2txyz(txyz_path, config_)

                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-results.log'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-time.dat'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)

                        while True:
                            if(hit_MS in self.parameter.MS_list or hit_MS == None):
                                break
                            if(hit_MS):
                                MS_list_.append(hit_MS)
                            new_MS_dir = self.parameter.work_dir + f'/crd/MS{hit_MS[0]}_{hit_MS[1]}'
                            if(not os.path.isdir(new_MS_dir)):
                                log(f"MS{hit_MS[0]}_{hit_MS[1]} is not in the list!!!", self.parameter.path_log)
                                os.makedirs(new_MS_dir)
                                new_MS_coord = new_MS_dir + f'/MS{hit_MS[0]}_{hit_MS[1]}.xyz'
                                copy(txyz_path, new_MS_coord)
                            hit_MS, index, ratio = self.__whether_hit_new_MS(job.work_path, job.coord_name, num_snapshot, 0, [self.traj_dict[i].length], MS_list_, False)

                        if(self.parameter.restrain == "RMSD"):
                            tmp_1 = job.work_path + f"/RMSD-{job.coord_name}.dat"
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        elif(self.parameter.restrain == "dihedral"):
                            tmp_1 = job.work_path + f"/dih-{job.coord_name}.dat"
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)

                        if(not hit_MS):
                            log(f"Restart at {job.work_path}/{job.coord_name}", self.parameter.path_log)
                            config_ = self.__read_tinker_arc(self.traj_dict[i].traj_name, num_snapshot-1)
                            self.__write2txyz(txyz_path, config_)
                            tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.arc'
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                            self.traj_dict[i].restart(job)
                            continue
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.arc'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.dyn'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        Finish_job_list.append(job.job_name)
                        Finish_num += 1
                        remain_job_num -= 1
                        MS_ = self.traj_dict[i].start_MS
                        self.traj_dict[i].finish_iter(hit_MS, index, ratio, self.parameter.saveFrequency)
                        time_ = self.traj_dict[i].time
                        msg_ = f"{self.traj_dict[i].traj_name} "
                        msg_ += f"MS{MS_[0]}_{MS_[1]} to {hit_MS[0]}_{hit_MS[1]}."
                        msg_ += " Time: %7.3f ps" %(time_)
                        msg_ += " Remaining: %5d" %(remain_job_num)
                        log(msg_, self.parameter.path_log)
                    else:
                        log(f"Restart at {job.work_path}/{job.coord_name}", self.parameter.path_log)
                        config_ = self.__read_tinker_arc(self.traj_dict[i].traj_name, num_snapshot-1)
                        self.__write2txyz(txyz_path, config_)
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}.arc'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-time.dat'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        tmp_1 = job.work_path + f'/{self.traj_dict[i].traj_name}-results.log'
                        if(os.path.isfile(tmp_1)): os.remove(tmp_1)

                        if(self.parameter.restrain == "RMSD"):
                            tmp_1 = job.work_path + f"/RMSD-{job.coord_name}.dat"
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)
                        elif(self.parameter.restrain == "dihedral"):
                            tmp_1 = job.work_path + f"/dih-{job.coord_name}.dat"
                            if(os.path.isfile(tmp_1)): os.remove(tmp_1)

                        self.traj_dict[i].restart(job)
            if(Finish_num == len_job_list):
                log(f"Iteration {self.parameter.iteration} completed", self.parameter.path_log)
                break

    def __create_more_traj(self):
        filePath = self.parameter.work_dir 
        crddir = filePath + '/crd'
        more_new = 1
        for i in self.parameter.MS_list:
            more_new = 1
            MS_dir = crddir + '/MS' + str(i[0]) + '_' + str(i[1])
            os.chdir(MS_dir)
            if(self.parameter.real_trajPerLaunch > self.parameter.max_trajPerLaunch):
                break
            for j in range(1, self.parameter.real_trajPerLaunch+1):
                random_int_list = []
                index_ = self.__index_modifying(j, 4)
                txyz_path = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{index_}.xyz"
                key_path = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{index_}.key"
                i_new = self.traj_dict[f"MS{str(i[0])}_{str(i[1])}-traj_{index_}"].start_MS

                lines = open(key_path).readlines()
                for line in lines:
                    terms = line.split()
                    if(len(terms) > 0 and terms[0].upper() == "RANDOMSEED"):
                        random_int_list.append(int(terms[1]))
                if(self.parameter.real_trajPerLaunch + more_new > self.parameter.max_trajPerLaunch):
                    break
                more_new += 1

                path_keyfile = filePath + '/free.key'
                copy(path_keyfile, key_path)
                r = self.__random_num_assign(random_int_list)
                random_int_list.append(r)
                keyfile(self.parameter, anchor=i, free=True, initial=True, key_path=key_path, random_num=r).generate()

                index_ = self.__index_modifying(j+self.parameter.real_trajPerLaunch, 4)
                txyz_path_2 = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{index_}.xyz"
                key_path = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{index_}.key"
                copy(path_keyfile, key_path)
                copy(txyz_path, txyz_path_2)
                r = self.__random_num_assign(random_int_list)
                random_int_list.append(r)
                keyfile(self.parameter, anchor=i, free=True, initial=True, key_path=key_path, random_num=r).generate()
                traj_ = traj(f"MS{str(i[0])}_{str(i[1])}-traj_{index_}", i_new, r, MS_dir)
                self.traj_dict[f"MS{str(i[0])}_{str(i[1])}-traj_{index_}"] = traj_
        self.parameter.real_trajPerLaunch += more_new

    def initialize(self, status=0):
        import copy
        restrn = self.parameter.restrain
        # Pre-analysis
        if(restrn == "2Ddihedral"):
            msg = "Begin Voronoi tessellation on 2D dihedral space"
            log(msg, self.parameter.path_log)
            self.__Voronoi_tesselletion()
        elif(restrn == "dihedral"):
            msg = "Make sure that the dihedral value of each anchor in the list follows the increasing or decreasing sequence\n"
            log(msg, self.parameter.path_log)
            # Examples of increasing sequence: -160, -120, -80, -40, 0, 40, 80, 120, 160
            # Examples of decreasing sequence: 160, 120, 80, 40, 0, -40, -80, -120, -160
            # Do not put the pbc (jumping point) in the intermediate point: -80, -40, 0, 40, 80, 120, 160, -160, -120 is not acceptable
            msg = "Default dihedral range: -180.0~180.0\n"
            log(msg, self.parameter.path_log)
            for i in range(1, self.parameter.AnchorNum):
                self.possible_MS.append((i,i+1))
                self.possible_MS_property[(i,i+1)] = self.__MS_property_analysis(i, i+1, False)
            if self.parameter.pbc:
                if(self.parameter.pbc == []):
                    log("Error: no pbc is given for dihedral restraint!\n", self.parameter.path_log)
                    sys.exit(0)
                pbc_1 = min(self.parameter.pbc[0], self.parameter.pbc[1])
                pbc_2 = max(self.parameter.pbc[0], self.parameter.pbc[1])
                self.possible_MS.append((pbc_1,pbc_2))
                self.possible_MS_property[(pbc_1,pbc_2)] =  self.__MS_property_analysis(pbc_1, pbc_2, True)
        elif(restrn == "RMSD"):
            msg = f"The RMSD between starting config and final config is {self.parameter.RMSD_max}.\n"
            diff = self.parameter.RMSD_max / self.parameter.RMSD_num
            start_interval = diff / 2.
            tmp_num = self.parameter.RMSD_num*2+1
            if(not self.parameter.reactant):
                self.parameter.reactant = [tmp_num*(self.parameter.RMSD_num-1)+1, tmp_num*self.parameter.RMSD_num+1]
                self.parameter.reactant += [tmp_num*self.parameter.RMSD_num+1, tmp_num*self.parameter.RMSD_num+2]
                self.parameter.reactant += [tmp_num*self.parameter.RMSD_num+1, tmp_num*(self.parameter.RMSD_num+1)+1]
            if(not self.parameter.product):
                self.parameter.product = [self.parameter.RMSD_num, self.parameter.RMSD_num+1]
                self.parameter.product += [self.parameter.RMSD_num+1, self.parameter.RMSD_num+tmp_num+1]
                self.parameter.product += [self.parameter.RMSD_num+1, self.parameter.RMSD_num+2]
            msg += f"The number of MS on one axis will be increased to {2*self.parameter.RMSD_num}\n"
            msg += f"The total number of possible MS: {4*self.parameter.RMSD_num*self.parameter.RMSD_num}\n"
            msg += f"Range: {start_interval} to {-start_interval+diff*self.parameter.RMSD_num*2} A for RMSD from each config.\n"
            log(msg, self.parameter.path_log)
            for i in range(1, self.parameter.RMSD_num*2+1):
                for j in range(1, self.parameter.RMSD_num*2+1):
                    tmp_1 = i+(j-1)*(self.parameter.RMSD_num*2+1)
                    tmp_2 = i+j*(self.parameter.RMSD_num*2+1)
                    self.possible_MS.append((tmp_1, tmp_1+1))
                    self.possible_MS.append((tmp_1, tmp_2))
                    tmp_3 = [start_interval+(i-1)*diff-self.parameter.MS_threshold, start_interval+(i-1)*diff+self.parameter.MS_threshold]
                    tmp_4 = [max(0., start_interval+diff*(j-2)), start_interval+diff*(j-1)]
                    self.possible_MS_property[(tmp_1, tmp_1+1)] = [tmp_3, tmp_4]
                    tmp_5 = [max(0., start_interval+diff*(i-2)), start_interval+diff*(i-1)]
                    tmp_6 = [start_interval+(j-1)*diff-self.parameter.MS_threshold, start_interval+(j-1)*diff+self.parameter.MS_threshold]
                    self.possible_MS_property[(tmp_1, tmp_2)] = [tmp_5, tmp_6]

        if self.parameter.milestone_search == 0 and status == 0:
            self.parameter.MS_list = copy.deepcopy(self.possible_MS)
            self.parameter.MS_property = copy.deepcopy(self.possible_MS_property)
        elif self.parameter.milestone_search == 0 and status == 1:
            self.__read_milestone_folder()
            if(restrn == "RMSD"):
                self.__RMSD_reactant_product(self.parameter.MS_list)
        else:
            while True:    
        ## Start here
            #   free runs from each anchors, markdown once it reaches another cell (i.e. closer to another anchor ).
                self.__seek_milestones()
            #   Determine reactant and product for RMSD Milestoning
                if(restrn == "RMSD"):
                    self.__RMSD_reactant_product(self.seek_MS)
            #   check if reactant and product are connected.
                if network_check(self.parameter, MS_list=self.seek_MS) == True:
                    break
            log("Seeking is completed\n", self.parameter.path_log)
            # read folders to get the milestones list 
            self.__read_milestone_folder()

                
        log("MS list completed\n", self.parameter.path_log)
        reactant_MS_name = (self.parameter.reactant[0],  self.parameter.reactant[1])
        product_MS_name = (self.parameter.product[0], self.parameter.product[1])
        log(f" Reactant: MS{self.parameter.reactant[0]}_{self.parameter.reactant[1]}", self.parameter.path_log)
        log(f" Product: MS{self.parameter.product[0]}_{self.parameter.product[1]}\n", self.parameter.path_log)
        self.parameter.reactant_milestone = [self.parameter.MS_list.index(reactant_MS_name)]
        self.parameter.product_milestone = [self.parameter.MS_list.index(product_MS_name)]

    def sampling(self):
        # Sampling begins here: restrained MD on each milestone and then unbiased MD to sample the whole configuration space
        if(self.parameter.skip_restrain == False and set(self.parameter.MS_list) != set(self.parameter.finished_restrain)):
            log("Begin restrained MD on each milestone\n", self.parameter.path_log)
            self.__restrain_MS_sampling()
        else:
            log("Restrained MD is skipped!\n", self.parameter.path_log)
            filePath = self.parameter.work_dir 
            crddir = filePath + '/crd'
            for i in self.parameter.MS_list:
                MS_dir = crddir + '/MS' + str(i[0]) + '_' + str(i[1])
                config_Num = 0
                for j in range(self.parameter.trajPerLaunch):
                    j_ = self.__index_modifying(j+1, 4)
                    txyz_path = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{j_}.xyz"
                    key_path = MS_dir + f"/MS{str(i[0])}_{str(i[1])}-traj_{j_}.key"
                    if(os.path.exists(txyz_path) and os.path.exists(key_path)):
                        r = 0
                        with open(key_path) as f:
                            lines = f.readlines()
                            for line in lines:
                                terms = line.split()
                                if(len(terms) > 0 and terms[0].upper() == "RANDOMSEED"):
                                    r = int(terms[1])
                        if(r == 0):
                            log(f"Error: No randomseed for traj{j_}", self.parameter.path_log)
                            sys.exit(0)
                        traj_ = traj(f"MS{i[0]}_{i[1]}-traj_{j_}", i, r, MS_dir)
                        self.traj_dict[f"MS{i[0]}_{i[1]}-traj_{j_}"] = traj_
                    else:
                        log(f"Error: No xyz file or key file for MS{i[0]}_{i[1]}-traj_{j_}", self.parameter.path_log)
                        sys.exit(0)
        log("Begin unbiased MD for the final sampling\n", self.parameter.path_log)
        self.trajPool_ = trajPool(self.parameter, self.traj_dict)
        MFPT_temp = 1e-19
        MFPT_temp_list = []
        while True:
            if(self.parameter.iteration >= 1):
                self.__create_more_traj()
            self.__final_sampling()
            self.trajPool_.compute_()
            if(self.parameter.iteration >= self.parameter.maxIteration):
                log("Reach max iteration", self.parameter.path_log)
                self.parameter.print_properties()
                break
            elif(np.isnan(self.parameter.MFPT_1) or self.parameter.MFPT_1 < 0.):
                log("Prepare for more free trajectories", self.parameter.path_log)
                self.parameter.print_properties()
                self.parameter.MFPT_1 = 0.
                self.parameter.MFPT_2 = 0.
                self.parameter.MFPT_1_rev = 0.
                self.parameter.MFPT_2_rev = 0.
            elif(np.abs(self.parameter.MFPT_1 - MFPT_temp) / MFPT_temp > self.parameter.tolerance):
                print("Comparison:", self.parameter.MFPT_1, MFPT_temp)
                log("Prepare for the next iteration", self.parameter.path_log)
                if(len(MFPT_temp_list) >= 5):
                    MFPT_temp_list.pop(0)

                if(self.parameter.iteration == 0):
                    MFPT_temp = self.parameter.MFPT_1
                else:
                    MFPT_temp_list.append(self.parameter.MFPT_1)
                    MFPT_temp = sum(MFPT_temp_list)/len(MFPT_temp_list)

                self.parameter.print_properties()
                self.parameter.MFPT_1 = 0.
                self.parameter.MFPT_2 = 0.
                self.parameter.MFPT_1_rev = 0.
                self.parameter.MFPT_2_rev = 0.
            else:
                self.parameter.print_properties()
                log("MFPT converged", self.parameter.path_log)
                break
            self.parameter.iteration += 1
            self.trajPool_.prep_for_next_iter()

