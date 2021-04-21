###############################################################################
# msys.py is created for rebuild model simulation test of Mechsys             #
# In dynamic simulation, especially for shaking table tests simulation, the   # 
# result should be redrew via relative position between SPH particles and DEM #
# particles.                                                                  #
# Created by Zheng Zhou 'zhougeohhu@foxmail.com' April 2021                   #
# Please follow the explanation and edit it for your simulation result        #
# Under editing                                                               #
###############################################################################

import numpy as np
import h5py
import os
import shutil
import matplotlib.pyplot as plt

# Python libraries above should be checked bedore your beginning
# Python 3.X is suggested
# If you use Spyder, you could have an interface similar to MatLAB

############################## files operation ################################
os.mkdir('/root/msys_sandbox/BM/')
os.mkdir('/root/msys_sandbox/BM/Deformation/')
os.mkdir('/root/msys_sandbox/BM/Deformation/pics/')
os.mkdir('/root/msys_sandbox/BM/Deformation/video/')
# Attention please, the new folder should be created in the fold including 
# result files and this py.file      

filePath = '/root/msys_sandbox/'
# Edit this variable definiation with aforementioned definition
for i in os.listdir(filePath):
    if i.startswith('Bigmodel2_'):
        # edit '...' string with the result name you set
        if i.endswith('.h5'):
            synpo   = i.index('_')
            ii      = i[synpo+1]
            # select the result files named by number
            if ii.isnumeric() == True :
                ## print (i,ii) ## for testing, skip this please
                shutil.copy2(i,'/root/msys_sandbox/BM/')
# copy h5.file to new folder for next step
                
fin = h5py.File('/root/msys_sandbox/Bigmodel2_Initial.h5',mode = 'r')
# print(fin.keys())
# read the initial model, and record some information
iniSPHposition = fin.get('Position')
dsp            = len(iniSPHposition)
rowdp          = int(dsp/3)
iniSPHposition = np.array(iniSPHposition).reshape(rowdp,3)
spxin          = np.min(iniSPHposition[:,0])
spxinmax       = np.max(iniSPHposition[:,0])
spyin          = np.min(iniSPHposition[:,1])
spyinmax       = np.max(iniSPHposition[:,1])
iniDEMposition = fin.get('DEMPosition')  
dp             = len(iniDEMposition)
rowdp          = int(dp/3)
iniDEMposition = np.array(iniDEMposition).reshape(rowdp,3)
dexin          = np.min(iniDEMposition[:,0])
dexinmax       = np.max(iniDEMposition[:,0])
deyin          = np.min(iniDEMposition[:,1])
deyinmax       = np.max(iniDEMposition[:,1])
dezin          = np.min(iniDEMposition[:,2])

fnn = (deyinmax-deyin)/(dexinmax-dexin)

print ('finish file operation')


########################## redraw deformatin results ##########################
filePath    = '/root/msys_sandbox/BM/'
for i in os.listdir(filePath):
    f           = h5py.File(i,mode='r')
    SPHposition = f.get('Position')
    sp          = len(SPHposition)
    rowsp       = int(sp/3)
    SPHposition = np.array(SPHposition).reshape(rowsp,3)
    # read, reshape and save SPH particles position data
    
    DEMposition = f.get('DEMPosition')  
    dp          = len(DEMposition)
    rowdp       = int(dp/3)
    DEMposition = np.array(DEMposition).reshape(rowdp,3)
    dexm        = DEMposition[1,0]
    deym        = DEMposition[0,1]
    dezm        = 0
    # read, reshape and save DEM particles position data
    
    deltm       = np.ones((rowsp,3))
    deltm[:,0]  = deltm[:,0]*dexm
    deltm[:,1]  = deltm[:,1]*deym
    deltm[:,2]  = deltm[:,2]*dezm
    # create a new array to save DEM axia limitation
    
    AA = np.zeros((rowsp,3))
    AA[:,0] = SPHposition[:,0] - deltm[:,0]
    AA[:,1] = SPHposition[:,1] - deltm[:,1]
    AA[:,2] = SPHposition[:,2] - deltm[:,2]
    # record the relative position to DEM axia limitation of SPH particles
    
    filename = i.rstrip('h5')
    filename = filename.rstrip('.')
    name = filename+'.txt'
    np.savetxt(filePath+'Deformation/'+name,AA)
    # name and save txt.files in the folder

    picPath = '/root/msys_sandbox/BM/Deformation/pics/'
    plt.subplots(constrained_layout=True);
    x = AA[:,0];
    y = AA[:,1];
    plt.figure(figsize=(7,7*fnn),dpi=300); 
    # plot setting, unit inch
    plt.axis([0,dexinmax-dexin,0,deyinmax-dexin]);
    # axia setting
    plt.xlabel('Distance');
    plt.ylabel('Elevation');
    plt.scatter(x,y,s=1.0);
    filename = filename.lstrip('Bigmodel2_')
    # Attention please, pictures are suggested to name by number for other
    # function could read pictures files in the actually orders.
    plt.savefig(picPath+filename+'.png');
    # jpg format couldn't be applied here
    
print ('finish drawing')

############################ redraw strain results ############################





