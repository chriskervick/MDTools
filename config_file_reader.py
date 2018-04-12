# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

def read_config(filename):
    file1 = open(filename,"r")
    lines = file1.readlines()
    file1.close()

    atom_groups = []
    dens_exp = []

    for i in lines:
        if i.split()[0] == 'u':
            scaling = [float(j) for j in i.split()[1:]]
        if i.split()[0] == 'n':
            atom_groups.append([int(i.split()[1]), int(i.split()[2])])
        if i.split()[0] == 'p':
            zmin = float(i.split()[1])
            zmax = float(i.split()[2])
            zstep = float(i.split()[3])
        if i.split()[0] == 'd':
            tmp = i.split()[1:]
            for j in tmp:
                dens_exp.append(float(j))

    return scaling, atom_groups, np.array([zmin + i*zstep for i in range(int(np.floor((zmax-zmin)/zstep)+1))]), np.array(dens_exp)