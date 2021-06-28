#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:41:17 2021

@author: zheng
"""
import matplotlib.ticker as ticker
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import shutil
import matplotlib as mpl
mpl.use('Agg')
# Shut down console graphic display

resultpath = '/root/msys_sandbox/'
newpath = '/home/zheng/Msys_result/' + 'Buiseismicfi/'
resultstart = 'buiseismicfi_'
os.mkdir(newpath)
os.mkdir(newpath + 'Position/')
os.mkdir(newpath + 'Deformation/')
os.mkdir(newpath + 'Deformation/pics/')
os.mkdir(newpath + 'Deformation/video/')
os.mkdir(newpath + 'Strain/')
os.mkdir(newpath + 'Strain/' + 'pics/')
os.mkdir(newpath + 'Strain/' + 'video/')

filePath = resultpath
for i in os.listdir(filePath):
    if i.startswith(resultstart):
        if i.endswith('.h5'):
            synpo = i.index('_')
            ii = i[synpo + 1]
            if ii.isnumeric() is True:
                shutil.copy2(i, newpath)
# Copy corresponding h5 files to new folder

fin = h5py.File(resultpath+resultstart+'Initial.h5', mode='r')
iniSPHposition = fin.get('Position')
dsp = len(iniSPHposition)
rowsp = int(dsp/3)
iniSPHposition = np.array(iniSPHposition).reshape(rowsp, 3)
np.savetxt(newpath+'iniSPHposition.txt', iniSPHposition)
# Plot setting start
X = iniSPHposition[:, 0]
X = X.flatten()
Y = iniSPHposition[:, 1]
Y = Y.flatten()
triang = tri.Triangulation(X, Y)

iniDEMposition = fin.get('DEMPosition')
dp = len(iniDEMposition)
rowdp = int(dp/3)
iniDEMposition = np.array(iniDEMposition).reshape(rowdp, 3)
np.savetxt(newpath+'iniDEMposition.txt', iniDEMposition)
dexin = np.min(iniDEMposition[:, 0])
dexinmax = np.max(iniDEMposition[:, 0])
deyin = np.min(iniDEMposition[:, 1])
deyinmax = np.max(iniDEMposition[:, 1])
dezin = np.min(iniDEMposition[:, 2])
fnn = (deyinmax-deyin)/(dexinmax-dexin)
# PLot setting end

# Max and min strain
filePath = newpath
stmax = 1*1e-10
stmin = 1*1e-10
for i in os.listdir(filePath):
    if i.endswith('.h5'):
        f = h5py.File(i, mode='r')
        SPHstrain = f.get('Strain')
        # For 2D simulation, 6 elements for each particle.
        # For 3D simulation, 9 elements for each particle.
        dst = len(SPHstrain)
        rowst = int(dst/6)
        SPHstrain = np.array(SPHstrain).reshape(rowst, 6)
        SPHstring = 0.5*(SPHstrain[:, 0]+SPHstrain[:, 4])
        fstmax = np.max(SPHstring)
        fstmin = np.min(SPHstring)
        if fstmax > stmax:
            stmax = fstmax
        if fstmin < stmin:
            stmin = fstmin

# Final step for plotting
filePath = newpath
for i in os.listdir(filePath):
    if i.endswith('.h5'):
        f = h5py.File(i, mode='r')
        SPHstrain = f.get('Strain')
        # For 2D simulation, 6 elements for each particle.
        # For 3D simulation, 9 elements for each particle.
        dst = len(SPHstrain)
        rowst = int(dst/6)
        SPHstrain = np.array(SPHstrain).reshape(rowst, 6)
        # read, reshape and save SPH particles strain data

        SPHposition = f.get('Position')
        sp = len(SPHposition)
        rowsp = int(sp/3)
        SPHposition = np.array(SPHposition).reshape(rowsp, 3)
        # read, reshape and save SPH particles position data

        DEMposition = f.get('DEMPosition')
        dp = len(DEMposition)
        rowdp = int(dp/3)
        DEMposition = np.array(DEMposition).reshape(rowdp, 3)
        dexm = DEMposition[1, 0]
        deym = DEMposition[0, 1]
        dezm = 0
        # read, reshape and save DEM particles position data

        deltm = np.ones((rowsp, 3))
        deltm[:, 0] = deltm[:, 0]*dexm
        deltm[:, 1] = deltm[:, 1]*deym
        deltm[:, 2] = deltm[:, 2]*dezm
        # create a new array to save DEM axia limitation

        AA = np.zeros((rowsp, 3))
        AA[:, 0] = SPHposition[:, 0] - deltm[:, 0]
        AA[:, 1] = SPHposition[:, 1] - deltm[:, 1]
        AA[:, 2] = SPHposition[:, 2] - deltm[:, 2]
        # record the relative position to DEM axia limitation of SPH particles

# =============================================================================
#         Please do not delete
#         name        = filename+'.txt'
#         namesp      = filename.lstrip(resultstart)+'.txt'
#         np.savetxt(filePath+'Strain/'+name,SPHstrain)
#         np.savetxt(filePath+'Position/'+namesp,SPHposition)
# =============================================================================

        filename = i.rstrip('h5')
        filename = filename.rstrip('.')
        # Attention please, pictures are suggested to name by number for other
        # function could read pictures files in the actually orders.
        # ===== Strain Contour Plotting
        name = filename+'.npy'
        namesp = filename.lstrip(resultstart)+'.npy'
        np.save(filePath+'Strain/'+name, SPHstrain)
        np.save(filePath+'Position/'+namesp, SPHposition)
        np.save(filePath+'Deformation/'+namesp, AA)
        # Strain Plot start
        Z = 0.5*SPHstrain[:, 0]+0.5*SPHstrain[:, 4]
        # 0.5*(Epsilon11 + 22) for 2D
        Z = Z.flatten()
        # Z :accumulated deviatoric strain
        fig1, ax1 = plt.subplots(
            constrained_layout=True, figsize=(3.5, 3.5), dpi=300)
        # ax1 = plt.figure(figsize=(7,7*fnn),dpi=300)
        ax1.set_aspect('equal')
        ax1.set_title('SPH particle number = ' + str(rowsp) +
                      ' Step = '+filename.lstrip(resultstart), fontsize=8)
        tt = ax1.tricontourf(triang, Z, cmap="jet",
                             levels=np.linspace(stmin, stmax, 100))
        plt.xticks(np.linspace(0, 0.9, 10), fontsize=8)
        plt.yticks(np.linspace(0, 0.5, 6), fontsize=8)
        plt.xlabel('Distance', fontsize=8)
        plt.ylabel('Elevation', fontsize=8)
        cb = fig1.colorbar(tt, orientation='horizontal', fraction=0.1)
        ticker_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = ticker_locator
        cb.set_ticks([stmin, 0.25*(stmin+stmax), 0.5 *
                      (stmin+stmax), 0.75*(stmin+stmax), stmax])
        cb.ax.tick_params(labelsize=8)
        cb.update_ticks()
        filename = filename.lstrip(resultstart)
        plt.savefig(newpath+'Strain/'+'pics/'+filename+'.png')
        plt.clf()
        plt.cla()
        plt.close()
        # ===== Strain Plot end

        picPath = newpath + 'Deformation/pics/'
        plt.plot(constrained_layout=True)
        x = AA[:, 0]
        y = AA[:, 1]
        plt.figure(figsize=(7, 5), dpi=300)
        # plot setting, unit inch
        plt.axis([0, dexinmax-dexin, 0, deyinmax-dexin])
        # axia setting
        plt.title('SPH particle number = ' + str(rowsp) +
                  ' Step = '+filename, fontsize=8)
        plt.xticks(np.linspace(0, 2, 11), fontsize=8)
        plt.yticks(np.linspace(0, 0.5, 6), fontsize=8)          
        plt.xlabel('Distance')
        plt.ylabel('Elevation')
        plt.scatter(x, y, s=1.0)
        axd = plt.gca()
        axd.set_aspect(1)
        plt.savefig(picPath+filename+'.png')
        plt.clf()
        plt.cla()
        plt.close()
        # jpg format couldn't be applied here
