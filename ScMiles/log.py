# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 15:49:08 2018

@author: Wei Wei

This subroutine writes running informations to log file.

"""

__all__ = ['log']


def log(msg, log_path):
    import os
    from datetime import datetime
    with open(log_path, 'a+') as f1:
        loginfo = str(datetime.now()).split('.')[0] + "    " + msg
        print(loginfo, file=f1)

def log_properties(q, p, energy, c, energy_err, ms_list, log_path):
    from datetime import datetime
    import numpy as np
    with open(log_path, 'a+') as f1:
        loginfo = str(datetime.now()).split('.')[0] + "    results:\n"
        loginfo += "{:>4} {:>4} {:>10} {:>8} {:>13} {:>10}".format('a1', 'a2', 'q', 'p', 'freeE(kT)', 'freeE_err')
        print(loginfo, file=f1)
        for i in range(len(q)):
            print('{:4d} {:4d} {:10.5f} {:8.5f} {:13.5f} {:10.5f}'.format(ms_list[i][0], ms_list[i][1], q[i], p[i],
                  energy[i], energy_err[i]), file=f1)
        print("committor func:",file=f1)
        c_ = np.squeeze(c)
        for i in range(len(ms_list)):
            print('{:4d} {:4d} {:15.8f}'.format(ms_list[i][0], ms_list[i][1], c_[i]), file=f1)

def log_properties_2(tau1, tau2, MFPT_err, MFPT_err2, log_path):
    with open(log_path, 'a+') as f2:
        print("MFPT is {:15.8e} ps, with an error of {:15.8e}, from eigenvalue method.".format(tau1, MFPT_err),file=f2)  
        print("MFPT is {:15.8e} ps, with an error of {:15.8e}, from inverse method.".format(tau2, MFPT_err2),file=f2)

def log_properties_rev(tau1, tau2, MFPT_err, MFPT_err2, log_path):
    with open(log_path, 'a+') as f2:
        print("Reverse the reactant and product state:")
        print("MFPT is {:15.8e} ps, with an error of {:15.8e}, from eigenvalue method.".format(tau1, MFPT_err),file=f2)  
        print("MFPT is {:15.8e} ps, with an error of {:15.8e}, from inverse method.".format(tau2, MFPT_err2),file=f2)

