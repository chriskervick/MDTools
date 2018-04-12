#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:27:00 2018

@author: btreece
"""

import numpy as np
import MDAnalysis
import config_file_reader

#VDW_TABLE = {'C': 1.7, 'F': 1.47, 'H': 1.2, 'N': 1.55, 'O': 1.52, 'P': 1.8, 'S': 1.8}

""" Need to reevaluate methyls POPC and POPG"""
POPC_RANGE={'headgroups':range(0,24),'tails':range(24,130),'methyls':range(130,134)}
POPG_RANGE={'headgroups':range(0,17),'tails':range(17,123),'methyls':range(123,127)}
DOPC_RANGE={'headgroups':range(0,24),'tails':range(24,87)+range(91,134),'methyls':range(87,91)+range(134,138)}
DOPS_RANGE={'headgroups':range(0,17),'tails':range(17,80)+range(84,127),'methyls':range(80,84)+range(127,131)}
lipid_grp_dict = {'POPC':POPC_RANGE,'POPG':POPG_RANGE,'DOPC':DOPC_RANGE,'DOPS':DOPS_RANGE}

def Lip_Grp_Select(univ,resname,name_range):
    NAMES=univ.select_atoms("resname "+resname).names
    str_tmp = "name "+NAMES[name_range[0]]
    for i in name_range[1:]:
        str_tmp += " or name "+NAMES[i]
    GROUP = univ.select_atoms("resname "+resname+" and ("+str_tmp+")")
    return GROUP

def Density_Chol(time_range, chol_info):
    if not chol_info.has_key('struct_file'):
        chol_info.update({'struct_file':None})
    if not chol_info.has_key('traj_file'):
        chol_info.update({'traj_file':None})
    if not chol_info.has_key('cholesterol'):
        chol_info.update({'cholesterol':None})
    
    if all([chol_info[i]==None for i in ['struct_file', 'traj_file', 'cholesterol']]):
        raise ReferenceError('Failure to specify atomselection "cholesterol" or the trajectory and structure.')
        return
    
    if chol_info['cholesterol'] == None:
        print "Using Files Supplied In Info."
        u = MDAnalysis.Universe(chol_info['struct_file'],chol_info['traj_file'])
        chol = u.select_atoms("resname CHL1")
        if len(chol.resids) == 0:
            raise ReferenceError('Failed to find cholesterol under the resname CHL1')
            return
    else:
        chol = chol_info['cholesterol']
    
    nframes = len(time_range)
    zmin = np.min(chol.universe.atoms.positions[:,2])
    zmax = np.max(chol.universe.atoms.positions[:,2])
    
    z = np.arange(zmin, zmax, (zmax - zmin)/500.)
    dens = 0*z
    
    sigma = 1.0
    norm = 1. / ( (2.0*np.pi)**0.5 * sigma ) # There are 74 atoms in cholesterol
    for ts in time_range:
        chol.universe.trajectory[ts]
        
        for i in chol.positions[:,2]:
            dens += norm * np.exp(-(z-i)**2.0 / (2.0*sigma**2.0))
    dens = dens/nframes
    
    return z, dens
            
        
    
def Density_MD_Bilayer(sim_info,lipids,time_range):
    """ lipids is a list of strings containing the resnames of the lipid components
    'POPC', 'POPG', and 'DOPC' are the current resnames available. """
    
    if not sim_info.has_key('struct_file'):
        sim_info.update({'struct_file':None})
    if not sim_info.has_key('traj_file'):
        sim_info.update({'traj_file':None})
    if not sim_info.has_key('universe'):
        sim_info.update({'universe':None})
    
    if all([sim_info[i]==None for i in ['struct_file', 'traj_file', 'universe']]):
        raise ReferenceError('Failure to specify atomselection "universe" or the trajectory and structure.')
        return
    
    if sim_info['universe'] == None:
        print "Using Files Supplied In Info."
        u = MDAnalysis.Universe(sim_info['struct_file'],sim_info['traj_file'])
    else:
        u = sim_info['universe']
    
    nframes=len(time_range)
        
#    for name in VDW_TABLE:
#        u.select_atoms("type "+name).set_radius(VDW_TABLE[name])
    zmin = np.min(u.trajectory.ts.positions[:,2])
    zmax = np.max(u.trajectory.ts.positions[:,2])
    
    z = np.arange(zmin,zmax,(zmax-zmin)/500.)
    head_dens={}
    tail_dens={}
    methyl_dens={}

    for lpd in lipids:
        head_dens[lpd]=0*z
        tail_dens[lpd]=0*z
        methyl_dens[lpd]=0*z
        
        heads = Lip_Grp_Select(u,lpd,lipid_grp_dict[lpd]['headgroups'])
        tails = Lip_Grp_Select(u,lpd,lipid_grp_dict[lpd]['tails'])
        methyls = Lip_Grp_Select(u,lpd,lipid_grp_dict[lpd]['methyls'])
        
        sigma = 1.0
        for ts in time_range:
            u.trajectory[ts]

#            for i,rad in zip(heads.positions,heads.radii):
#                head_dens[lpd] += (1/(2*np.pi*(rad*sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(rad*sigma)**2.0))
#            for i,rad in zip(tails.positions,tails.radii):
#                tail_dens[lpd] += (1/(2*np.pi*(rad*sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(rad*sigma)**2.0))
#            for i,rad in zip(methyls.positions,methyls.radii):
#                methyl_dens[lpd] += (1/(2*np.pi*(rad*sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(rad*sigma)**2.0))
            for i in heads.positions:
                head_dens[lpd] += (1/(2*np.pi*(sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(sigma)**2.0))
            for i in tails.positions:
                tail_dens[lpd] += (1/(2*np.pi*(sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(sigma)**2.0))
            for i in methyls.positions:
                methyl_dens[lpd] += (1/(2*np.pi*(sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(sigma)**2.0))
        
        head_dens[lpd] = head_dens[lpd]/nframes
        tail_dens[lpd] = tail_dens[lpd]/nframes
        methyl_dens[lpd] = methyl_dens[lpd]/nframes

    return z,head_dens,tail_dens,methyl_dens
    
def Density_MD_Protein(time_range, protein_info = {}):
    
    # Check for input files
    if not protein_info.has_key('struct_file'):
        protein_info.update({'struct_file':None})
    if not protein_info.has_key('traj_file'):
        protein_info.update({'traj_file':None})
    # Check for protein selection (preprocessed files)
    if not protein_info.has_key('protein'):
        protein_info.update({'protein':None})
    
    # If nothing is provided, raise error
    if all([protein_info[i]==None for i in ['struct_file', 'traj_file', 'protein']]):
        raise ReferenceError('Failure to specify protein for trajectory or structure.')
        return
    # If no selection is provided, import data. Else, use the selection
    if protein_info['protein'] == None:
        print "Using Files Supplied In Info."
        u = MDAnalysis.Universe(protein_info['struct_file'],protein_info['traj_file'])
        protein = u.select_atoms("protein")
    else:
        protein = protein_info['protein']
    
    # Number of frames
    nframes=len(time_range)
    
    # For modifying sigma by VdW radius
#    for name in VDW_TABLE:
#        u.select_atoms("type "+name).set_radius(VDW_TABLE[name])
    # Array based on extent of protein
    zmin = np.min(protein.ts.positions[:,2])
    zmax = np.max(protein.ts.positions[:,2])
    # Expanding for Gaussian tails
    zmin = zmin - 0.5*(zmax-zmin)
    zmax = zmax + 0.5*(zmax-zmin)
    
    # z and density array initialize
    z = np.arange(zmin,zmax,(zmax-zmin)/500.)
    dens = 0*z
    
    # Gaussian width
    sigma = 1.0
    
    for ts in time_range:
        # Update frame
        protein.universe.trajectory[ts]
        z_bbox = protein.universe.dimensions[2]
        
        # VdW radius
#        for i,rad in zip(protein.positions,protein.radii):
#            dens += (1/(2*np.pi*(rad*sigma)**2)**0.5)*np.exp(-(z-i[2])**2.0/(2*(rad*sigma)**2.0))
        for i in protein.positions[:,2]:
            if i < 0.2*z_bbox:
                i += z_bbox
            dens += (1/(2*np.pi*(sigma)**2)**0.5)*np.exp(-(z-i)**2.0/(2*(sigma)**2.0))

    # Normalize frames
    dens = dens/nframes

    return z,dens

def Calculate_Potential_and_Force(config_file, time, protein_info = {}, gmx = True, New_Potential = True):
    
    if not protein_info.has_key('struct_file'):
        protein_info.update({'struct_file':None})
    if not protein_info.has_key('traj_file'):
        protein_info.update({'traj_file':None})
    if not protein_info.has_key('protein'):
        protein_info.update({'protein':None})
    
    #### Import Target Density
    scaling, atom_groups, z_target, target_dens = config_file_reader.read_config(config_file)
    # MDAnalysis uses A not nm (gmx)
    # z = [A]
    # rho = [A^{-1}]
    # scaling shape = [kJ-A/mol]
    # scaling offset = [kJ/mol-A]*[1/A]
    if gmx:
        z_target = 10.*z_target
        target_dens = target_dens/10.
        scaling = [10.*scaling[0], 0.01*scaling[1]]
    
    #### Calculate Simulation Density as it is done in GROMACS
    
    zmin = z_target[0]
    zstep = z_target[1] - z_target[0]
    if all([protein_info[i]==None for i in protein_info.keys()]):
        raise ReferenceError('Failure to specify protein for trajectory or structure.')
        return
    if protein_info['protein'] == None:
        print "Using Files Supplied In Info."
        u = MDAnalysis.Universe(protein_info['struct_file'],protein_info['traj_file'])
        protein = u.select_atoms("protein")
    else:
        protein = protein_info['protein']
        
    protein.universe.trajectory[time]
    z_bbox = protein.universe.dimensions[2]
    
    prot_dens = 0*np.arange(0.,10000.)
    sigma = 1.0
    nrm   = (atom_groups[0][1] - atom_groups[0][0] + 1)*sigma*(2*np.pi)**0.5
   
    for i in protein.positions:
        z = i[2]
        if z < zmin:
            z += z_bbox
        z_low = np.floor( (z - 3.0 - zmin) / zstep)
        z_hi  = np.ceil( (z + 3.0 - zmin) / zstep)
        for j in range(int(z_low), int(z_hi)):
                prot_dens[j] += np.exp(-0.5*(z - (zmin + zstep*j))**2.0/((sigma)**2.0)) / nrm
    
    potential = 0*prot_dens

    #### Old Potential
    if not New_Potential:
        for i in range(len(potential)):
            potential[i] += prot_dens[i]
            if i < len(target_dens):
                potential[i] += -target_dens[i]
        
        potential = potential*scaling[0]
        offset = 0
        
    #### New Potential
    if New_Potential:
        sim_mean = sum([(zmin + zstep*i)*prot_dens[i] for i in range(len(prot_dens))])*zstep
        exp_mean = sum([(zmin + zstep*i)*target_dens[i] for i in range(len(target_dens))])*zstep
        offset = sim_mean - exp_mean
        
        for i in range(len(potential)):
            potential[i] += prot_dens[i]
            indx_e = i - int(np.floor(offset/zstep))
            if indx_e >=0 and indx_e < len(target_dens):
                potential[i] += - target_dens[indx_e]
#        for i in range(len(target_dens)):
#            indx_offset = int(np.floor(((zmin + i*zstep) - offset - zmin)/zstep))
#            potential[indx_offset] += -target_dens[i]

        potential = potential*scaling[0]
        
    #### Forces
    z_f_by_atom = []
    CoM = protein.center_of_mass()
    tau = np.array([0.,0.,0.])

    for r in protein.positions:
        z = r[2]
        indx = int(np.floor( (z - zmin)/ zstep ))
        f1s = scaling[0]*( ( prot_dens[indx] - prot_dens[indx+1] )/ zstep )
        indx = int(np.floor( (z - offset - zmin)/zstep))
        if indx >= 0 and indx < len(target_dens)-1:
            f1e = -scaling[0]*( ( target_dens[indx] - target_dens[indx+1] )/ zstep )
        f2 = - scaling[1]*offset/(atom_groups[0][1] - atom_groups[0][0] + 1)
        z_f_by_atom.append([z, f1s, f1e, f2])
        tau += np.cross(r-CoM, np.array([0.,0.,f1s+f1e+f2]))

    return {'scale':scaling, 'atom_groups':atom_groups, 'zmin':zmin, 'zstep':zstep, 'prot_dens':prot_dens,
            'z_target':z_target, 'target_dens':target_dens, 'potential':potential,
            'offset':offset, 'z_f_by_atom':z_f_by_atom, 'torque':tau}
