# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:51:25 2016

@author: btreece
"""

import numpy as np
from matplotlib import pyplot as plt

# Dictionaries for location of data in the data file
bilayer_dict={'z':0, 'lipid2':1, 'substrate':2, 'bME':3, 'protein':4, 'lipid1':5, 'headgroup2':6, 'headgroup2_2':7, 'headgroup2_3':8, 'normarea':9, 'tetherg':10, 'headgroup1_3':11, 'headgroup1_2':12, 'defect_hg':13, 'methyl1':14, 'tether':15, 'defect_hc':16, 'methyl2':17, 'headgroup1':18, 'summol':19, 'water':20, 'waterperc':21}
protein_dict={'median_area':0, 'p2sigma_nsld':1, 'm2sigma_nsld':2, 'p2sigma_area':3, 'median_nsl':4, 'm2sigma_nsl':5, 'm2sigma_area':6, 'psigma_area':7, 'msigma_nsl':8, 'msigma_nsld':9, 'p2sigma_nsl':10, 'zaxis':11, 'psigma_nsld':12, 'psigma_nsl':13, 'median_nsld':14, 'msigma_area':15}
group_dict = {'Outer_headgroup':['headgroup2','headgroup2_2','headgroup2_3'],'Inner_headgroup':['headgroup1','headgroup1_2','headgroup1_3'],'Tether':['bME','tetherg','tether'],'Outer_lipid':['lipid2','methyl2'],'Inner_lipid':['lipid1','methyl1']}


def Extract_Data_Bilayer(filepath):
    global bilayer_dict
    
    # Read the file
    file1=open(filepath,"r")
    lines=file1.readlines()
    file1.close()
    
    # Verify the dictionaries match
    for i in range(len(lines[0].split())):
        if i != bilayer_dict[lines[0].split()[i]]:
            print "Incorrect dictionary preset, defining new dictionary."
            temp_dict = {}
            for j in range(len(lines[0].split())):
                temp_dict[lines[0].split()[j]]=j
                bilayer_dict = temp_dict

    # Process the data
    dat=[]
    for i in range(1,len(lines)):
        dat.append([float(j) for j in lines[i].split()])
        
    return np.array(dat)
 
def Extract_Data_Protein(filepath):
    global protein_dict
    
    # Read the file
    file1=open(filepath,"r")
    lines=file1.readlines()
    file1.close()
    
    # Verify the dictionaries match
    for i in range(len(lines[0].split())):
        if i != protein_dict[lines[0].split()[i]]:
            print "Incorrect dictionary preset, defining new dictionary."
            temp_dict = {}
            for j in range(len(lines[0].split())):
                temp_dict[lines[0].split()[j]]=j
                protein_dict = temp_dict

    # Process the data
    dat=[]
    for i in range(1,len(lines)):
        dat.append([float(j) for j in lines[i].split()])
        
    return np.array(dat)

   
def Plotting(groups,data,colors=None):
    """ Plotting(groups): groups is a list of lists denoting which data sets to plot. Elements of a single list are summed prior to plotting.
    colors is a list of strings to assign colors to the elements to plot, i.e. ['b','k','r']. Other plotting strings can be used ['.']
    
    groups = [['headgroup1','headgroup1_2','headgroup1_3'],['substrate']]
    
    groups = [['lipid1','methyl1']]
    """
    Z = np.array([i[bilayer_dict['z']] for i in data])

    plot_items=[]
    lines=[]

    for group in groups:
        plot_items.append(np.array([sum([i[bilayer_dict[element]] for element in group]) for i in data]))

    fig = plt.figure()
    ax=fig.add_subplot(1,1,1)
    if colors==None:
        for item in plot_items:
            lines.append(ax.plot(Z,item))
            plt.show()
    else:        
        for i in range(len(plot_items)):
            lines.append(ax.plot(Z,plot_items[i],colors[i]))
            plt.show()
    
    return lines,ax
    
def Generate_Formatted_Plot(bilayer_path,protein_path=None):
    
    data=Extract_Data_Bilayer(bilayer_path)
    lines,ax=Plotting([['substrate'],group_dict['Tether'],group_dict['Outer_headgroup']+group_dict['Inner_headgroup'],group_dict['Outer_lipid'],group_dict['Inner_lipid']],data,['y','g','c','b','b'])
    for i in range(len(lines)):    
        plt.setp(lines[i],linewidth=3.0)

    if protein_path != None:
        pdata = Extract_Data_Protein(protein_path)
        lines.append(ax.plot(np.array([i[protein_dict['zaxis']] for i in pdata]),np.array([i[protein_dict['median_area']] for i in pdata]),'r',linewidth = 3.0,label='protein with 68%\nconfidence interval'))
        lines.append(ax.fill_between(np.array([i[protein_dict['zaxis']] for i in pdata]),np.array([i[protein_dict['msigma_area']] for i in pdata]),np.array([i[protein_dict['psigma_area']] for i in pdata]),color='r',alpha = 0.2,linewidth=0.))
    ax.set_axisbelow(True)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')
    
    plt.tick_params(axis='both',labelsize=20)
    plt.xlabel(r'distance $\AA$',fontsize=20)
    plt.ylabel(r'cross-sectional area $\AA^2$',fontsize=20)
    
    plt.setp(lines[0],label='substrate')
    plt.setp(lines[1],label='tether')
    plt.setp(lines[2],label='headgroups')
    plt.setp(lines[3],label='hydrocarbons')
    ax.legend(loc='upper right',fontsize=18)
