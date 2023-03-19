#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:37:37 2018

@author: Wei Wei


This subroutine takes the k and t, and calculate flux, probability,
free energy, committor, and MFPT.
"""

import os
import pandas as pd
import numpy as np
from network_check import pathway
#from voronoi_plot import voronoi_plot

#__all__ = ['k_average','compute']


def find_ms_index(ms, ms_index):
    ms = sorted(ms)
    return(int(list(ms_index.keys())[list(ms_index.values()).index(ms)]))


def get_ms_of_cell(cell, ms_index):
    '''if reactant/product is a cell, return all milestones associated with this cell'''
    ms_of_cell = []
    for item in ms_index.values():
        if int(cell) in item:
            ms_of_cell.append(int(list(ms_index.keys())[list(ms_index.values()).index(item)]))
    return ms_of_cell


def k_average(k_count):
    '''convert count matrix to probability matrix'''
    dim = len(k_count)
    k_ave = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if i != j:
                if np.sum(k_count[i, :]) != 0.0:
                    k_ave[i, j] = k_count[i, j] / np.sum(k_count[i, :]) 
    return k_ave


def k_error(k_count):
    '''return a new random k based on beta function'''
    dim = len(k_count)
    k = np.zeros((dim, dim))
    for i in range(dim):
        total = np.sum(k_count[i]) 
        for j in range(dim):
            a = k_count[i,j]
            if i == j or a == 0:
                continue
            if total == a:
                k[i, j] = 1.0
            b = total - a
            if b > 0: 
                k[i, j] = np.random.beta(a, b)
        if sum(k[i, :]) != 0:
            k[i, :] = k[i, :] / sum(k[i, :])
    return k

 
def t_error(t, t_std):
    '''return a new set of life time based on the std'''
    return np.abs(np.random.normal(t, t_std, len(t))).tolist()


def committor(parameter, k):
    '''committor function'''
    kk = k.copy()
    for i in parameter.reactant_milestone:
        kk[i] = [0 for j in k[i]]
    for i in parameter.product_milestone:
        kk[i] = [0 for j in k[i]]
        kk[i][i] = 1.0
    c = np.linalg.matrix_power(kk,1000000000)
    A = np.ones((len(c),1))
    return np.matmul(c,A)


def flux(k):
    '''flux calculation'''
    kk = k.copy()
    kk_trans = np.transpose(kk)
    e_v, e_f = np.linalg.eig(kk_trans)
    idx = np.abs(e_v - 1).argmin()  
    q = [i.real for i in e_f[:, idx]]
    q = np.array(q)
    if np.all(q < 0):
        q = -1 * q
    return q


def prob(q, t):
    '''probability calculation'''
    p = np.transpose(q) * np.squeeze(t)
    p_norm = p / np.sum(p)
    p_norm[p_norm <= 0.] = 1e-50
    return p_norm


def free_energy(p):
    '''free energy from probability'''
    return -1.0 * np.log(p)
    

def get_boundary(parameter, ms_index):
    '''
    get the index number for reactant and product state

    If a K matrix is like this, reactant is 1_2, this function will return 0 as the first row/column indicates 1_2
    1_2 2_3 3_4
      0   1   0
    0.5   0 0.5
      0   1   0
    '''
    if len(parameter.reactant) == 2:
        bc1 = sorted(parameter.reactant)
        if bc1 in ms_index.values():
            parameter.reactant_milestone.append(int(list(ms_index.keys())[list(ms_index.values()).index(bc1)]))
        else:
            parameter.reactant_milestone.append(-1)
    else:
        for item in ms_index.values():
            if int(parameter.reactant[0]) in item and int(list(ms_index.keys())[list(ms_index.values()).index(item)]) \
                    not in parameter.reactant_milestone:
                parameter.reactant_milestone.append(int(list(ms_index.keys())[list(ms_index.values()).index(item)]))
            
    if len(parameter.product) == 2:
        bc2 = sorted(parameter.product)
        if bc2 in ms_index.values():   
            parameter.product_milestone.append(int(list(ms_index.keys())[list(ms_index.values()).index(bc2)]))
        else:
            parameter.product_milestone.append(-1)
    else:
        for item in ms_index.values():
            if int(parameter.product[0]) in item and int(list(ms_index.keys())[list(ms_index.values()).index(item)]) \
                    not in parameter.product_milestone:
                parameter.product_milestone.append(int(list(ms_index.keys())[list(ms_index.values()).index(item)]))            


def MFPT(kk, t, start, end):
    '''MFPT based on flux'''
    k = kk.copy()
    for i in end:
        k[i] = [0 for j in k[i]]
        for j in start:
            k[i][j] = 1.0 / len(start)
    q = flux(k)
    qf = 0
    for i in end:
        qf += q[i]
    tau = np.dot(q, t) / qf
    return float(tau)
    
def MFPT2(k, t, start, end, parameter):
    '''MFPT based on inverse of K'''
    dim = len(k)
    I = np.identity(dim)
    k2 = k.copy()
    for i in end:
        if i == -1:
            return -1
        k2[i] = [0 for i in k2[i]]
    
    p0 = np.zeros(dim)
    for i in start:
        if i == -1:
            return -1
        p0[i] = 1 / (len(start))
        
    if np.linalg.det(np.mat(I) - np.mat(k2)) == 0.0:
        parameter.sing = True
        return -1
    else:
        parameter.sing = False
        tau = p0 * np.linalg.inv(I - np.mat(k2)) * np.transpose(np.mat(t))
    return float(tau)

