# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:49:24 2016

@author: btreece
"""

import numpy as np
from matplotlib import pyplot as plt
import MDAnalysis
import config_file_reader
import Utility_Functions

def Plot_MD_Specify(Target=False, COM_Target=False, Bilayer=False, Protein=False, Cholesterol=False, config_file = None, sim_info = {}, time_range = None, lipid_list = None, gmx = True):
    
    if not sim_info.has_key('struct_file'):
        sim_info.update({'struct_file':None})
    if not sim_info.has_key('traj_file'):
        sim_info.update({'traj_file':None})
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    if Target:
        if config_file == None:
            raise ReferenceError('No configuration file specified.')

        scaling, atom_groups, z_target, target_dens = config_file_reader.read_config(config_file)
        # GROMACs uses nm not A
        if gmx:
            z_target = 10.*z_target
            target_dens = target_dens/10.
        # The protein analysis in this script doesn't normalize for number of atoms, so this needs scaled
        target_dens = target_dens*(atom_groups[0][1]-atom_groups[0][0]+1)
            
        ax.plot(z_target, target_dens, 'r',linewidth = 3.0, label = 'target')
        
    if Bilayer:
        if not sim_info.has_key('universe'):
            sim_info.update({'universe':None})
            
        if all([sim_info[i]==None for i in ['struct_file', 'traj_file', 'universe']]):
            raise ReferenceError('Failure to specify atomselection "universe" or the trajectory and structure.')
            return

        z_lipid,head_dens,tail_dens,methyl_dens = Utility_Functions.Density_MD_Bilayer(sim_info, lipid_list,time_range)
        ax.plot(z_lipid,sum([head_dens[lpd] for lpd in lipid_list]),'c',linewidth = 3.0, label = 'headgroups')
        ax.plot(z_lipid,sum([tail_dens[lpd] for lpd in lipid_list])+sum([methyl_dens[lpd] for lpd in lipid_list]),'b',linewidth = 3.0, label = 'tails')

    if Protein:
        if not sim_info.has_key('protein'):
            sim_info.update({'protein':None})
            
        if all([sim_info[i]==None for i in ['struct_file', 'traj_file', 'protein']]):
            raise ReferenceError('Failure to specify atomselection "protein" or the trajectory and structure.')
            return
 
        z_prot,prot_dens = Utility_Functions.Density_MD_Protein(time_range, sim_info)
        ax.plot(z_prot,prot_dens,'k',linewidth=3.0,label= 'protein')
    
    if Cholesterol:
        if not sim_info.has_key('cholesterol'):
            sim_info.update({'cholestererol':None})
            
        if all([sim_info[i]==None for i in ['struct_file', 'traj_file', 'cholesterol']]):
            raise ReferenceError('Failure to specify atomselection "cholesterol" or the trajectory and structure.')
            return

        z_chol, chol_dens = Utility_Functions.Density_Chol(sim_info, time_range)
        ax.plot(z_chol, chol_dens, 'magenta', linewidth=3.0, label= 'cholesterol')
    
    if COM_Target:
        if not Protein:
            raise ReferenceError('No protein specified to find COM.')
            return

        if config_file==None:
            raise ReferenceError('No configuration file specified.')
            return

        if not Target:
            scaling, atom_groups, z_target, target_dens = config_file_reader.read_config(config_file)
            # GROMACs uses nm not A
            if gmx:
                z_target = 10.*z_target
                target_dens = target_dens/10.
            # The protein analysis in this script doesn't normalize for number of atoms, so this needs scaled
            target_dens = target_dens*(atom_groups[0][1]-atom_groups[0][0]+1)

        p_COM = sum([z_prot[i]*prot_dens[i]*(z_prot[1]-z_prot[0]) for i in range(len(z_prot)) if z_prot[i] > z_target[0]])
        boundary = (z_prot[-1] + 0.05*z_prot[0])/1.05 # This comes from how the z_prot was defined
        p_COM += sum([(z_prot[i]+boundary)*prot_dens[i]*(z_prot[1]-z_prot[0]) for i in range(len(z_prot)) if z_prot[i] < z_target[0]])

        t_COM = sum([z_target[i]*target_dens[i]*(z_target[1]-z_target[0]) for i in range(len(z_target))])
        
        ax.plot(z_target+(p_COM-t_COM)/(atom_groups[0][1]-atom_groups[0][0]+1), target_dens, 'g',linewidth = 3.0, label = 'target_COM')
        
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.show()
    fig.tight_layout(rect=[0,0,0.8,1])

def Plot_Potential(config_file, struct_file, traj_file, time, plot_densities=False, remove_scaling = False, target_offset = True):
    
    results = Utility_Functions.Calculate_Potential(config_file, struct_file, traj_file, time)
    scaling = results['scale']
    atom_groups = results['atom_groups']
    z_pot = results['z_pot']
    potential = ['potential']
    
    if plot_densities:
        z_prot = results['z_prot']
        prot_dens = results['prot_dens']
        z_target = results['z_target']
        target_dens = results['target_dens']
    if remove_scaling:
        potential = potential/scaling[0]
    if target_offset:
        offset = results['offset']
        z_target = z_target + offset
        
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(z_pot, potential, 'purple',linewidth = 3.0, label = 'Potential')

    if plot_densities:
        ax.plot(z_prot,prot_dens,'k',linewidth=3.0,label= 'protein')
        ax.plot(z_target, target_dens, 'r',linewidth = 3.0, label = 'target')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig.show()
    fig.tight_layout(rect=[0,0,0.8,1])    