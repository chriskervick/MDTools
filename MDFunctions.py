#import MDAnalysis
#import numpy as np
#import sys

###### SOME PARSING




def GetDensity(u,sel="protein",N=10000,fname="output",start=0,end=1,z_start=0,z_end=500,timestep=0):
    import numpy as np
    if timestep > len(u.trajectory):
        print("Given timestep is outside range of trajectory")
    
    def gaussian(x,mu,sig):
        return 1/np.sqrt(2*3.14159265359)*np.exp(-np.power(x-mu,2.) / (2.0*np.power(sig,2.)))


    fname = fname + ".dat"
    pro = u.select_atoms(sel)
    pospro = pro.positions
    u.trajectory[timestep]
    LEN = np.size(pospro[:,2])
    box = u.universe.dimensions
    zlim = box[2]

    if sel=="protein":
        for i in range(0,LEN):
            if pospro[i,2] < zlim/2:
                pospro[i,2] = pospro[i,2] + zlim

    array = np.zeros((N,2))
    array[:,0] = np.linspace(z_start,z_end,N)
            
    for i in range(0,LEN):
        array[:,1] = array[:,1] + gaussian(np.linspace(z_start,z_end,N),pospro[i,2],1.0)
        
        
    #array[:,1] = array[:,1] / len(pro)
    np.savetxt(fname,array)
    
    print("Density saved to " + fname)
    return(array)


def GetZforce(u,sel="protein",start=0,end=1):
    import numpy as np
    pro=u.select_atoms(sel)
    forces = pro.forces
    zforce = np.sum(forces[:,2])
    return(zforce)



def GetTorque(u,sel="protein",start=0,end=1):
    import numpy as np
    pro = u.select_atoms(sel)
    com = pro.center_of_mass()
    forces = pro.forces
    pos_wrt_com = pro.positions - pro.center_of_mass()
    torques = np.cross(pos_wrt_com,forces)
    torque = np.array([np.sum(torques[:,0]),np.sum(torques[:,1]),np.sum(torques[:,2])])
    return(torque)


def DensityFromExp(fname):
    import numpy as np
    data = open(fname)
    lines = data.readlines()
    data.close()
    ###The below are hardcoded for Brad's naming conventions of the experimental data file
    l0 = lines[0].split()
    l1 = lines[1].split()
    l2 = lines[2].split()
    l3 = lines[3].split() 
    scale = float(l0[1])
    atom_start = float(l1[1])
    atom_end = float(l1[2])
    z_start = float(l2[1])
    z_end = float(l2[2])
    z_step = float(l2[3])
    del(l3[0])
    density = np.array([float(i) for i in l3])
    density = density/10
    num = (z_end - z_start)/z_step + 2 
    zarray = np.arange(z_start,z_start+5000*z_step,z_step)
    zarray = 10*zarray
    zeros = np.zeros(np.size(zarray) - np.size(density))
    density = np.append(density,zeros)
    
    return(zarray,density)


def Forcing(sim,exp,scale):
    import numpy as np
    import math

    ###Make sure exp and sim are the same size
    if (np.size(exp[:,1]) != np.size(sim[:,1])):
        raise NameError('The experimental and simulation datasets are not the same size')
    from scipy.integrate import simps
    exp_mean = simps(exp[:,0]*exp[:,1],exp[:,0])
    print(exp_mean)
    
    sim_mean = simps(sim[:,0]*sim[:,1],sim[:,0])
    print(sim_mean)
    te = exp[:,1]
    ts = sim[:,1]


    spacing = sim[2,0] - sim[1,0]

    offset = sim_mean - exp_mean

    movebyidx = int(math.floor(offset/spacing))

    temp_te = te

    temp_zeros = np.zeros(movebyidx)

    te_padded = np.append(temp_zeros,temp_te)
    te_cut = te_padded[0:np.size(ts)]

    if (np.size(te_cut) != np.size(ts)):
        raise NameError('The shifted experimental and simulation datasets are not the same size')

    
    offset_force = -1.0*scale*(sim_mean - exp_mean)

    difference = ts-te_cut
    gradient = np.zeros((np.size(sim[:,1]),2))
    for i in range(0,np.size(sim[:,1])-1):
        gradient[i,1] = difference[i+1] - difference[i]
    gradient[:,1] = gradient[:,1] * -1.0 * scale / (sim[2,0] - sim[1,0])
    gradient[:,0] = sim[:,0]
    gradient[:,1] = gradient[:,1] + offset_force

    return(gradient)
