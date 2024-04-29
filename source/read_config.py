#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, time, sys

def read_tinker_arc(coord_name):
    f = open(coord_name + '.arc')
    lines = f.readlines()
    f.close()
    return lines

def index_modifying(index, max_digit):
    index_ = str(index)
    length = len(index_)
    z = '0'
    zero_num = max_digit - length if max_digit >= length else 0
    z = z*zero_num
    return z + index_

def write2txyz(txyz_path, config_lines):
    with open(txyz_path, 'w') as f:
        for line in config_lines:
            f.write(line)

MS_name = sys.argv[1]
crddir = sys.argv[2]
eq_num = int(sys.argv[3])
num_snapshot = int(sys.argv[4])
interval = int(sys.argv[5])
MS_dir = crddir + '/' + MS_name
lines = read_tinker_arc(MS_name)
index = 1
nLines = int(lines[0].split()[0]) + 2
for j in range(eq_num, num_snapshot, interval):
    index_ = index_modifying(index, 4)
    txyz_path = MS_dir + f"/{MS_name}-traj_{index_}.xyz"
    write2txyz(txyz_path, lines[j*nLines:(j+1)*nLines])
    index += 1
