import numpy as np
# import h5py
import os
import pandas as pd
# import shutil
# import matplotlib.pyplot as plt
# import matplotlib.tri as tri
# import matplotlib.ticker as ticker

#         yset2=np.append(yset2,[[xi,yi]],axis=0)



workpath = "/root/msys_sandbox/Buislopeseven/Position/"
# os.mkdir('/root/msys_sandbox/Buislope/PointTrack/')
fildn = len(os.listdir(workpath))
finre = np.zeros([fildn, 14])
for i in os.listdir(workpath):
    f = np.load(workpath+i)
    fsize = np.shape(f)
    flpnb = fsize[0]
    newname = i.rstrip('.npy')
    numod = int(newname)
    sfixxset = [0.025, 0.225, 0.425]
    # Monitor horizontal displacement, fix y, set name x
    sfixyset = [0.200, 0.350, 0.650, 0.850]
    # Monitor vertical displacement, fix x, set name y
    yset1 = np.zeros([1, 2])
    yset2 = np.zeros([1, 2])
    yset3 = np.zeros([1, 2])
    xset1 = np.zeros([1, 2])
    xset2 = np.zeros([1, 2])
    xset3 = np.zeros([1, 2])
    xset4 = np.zeros([1, 2])
    sresult = np.zeros([7, 4])
    for ii in range(0, 3):
        sresult[ii, 0] = 0
        # Define horizontal direction by 0, fix y
        sresult[ii, 1] = sfixxset[ii]
    for ii in range(3, 7):
        sresult[ii, 0] = 1
        # Define vertical direction by 1, fix x
        sresult[ii, 1] = sfixyset[ii-len(sfixxset)]
    for j in range(0, flpnb):
        xi = f[j, 0]
        yi = f[j, 1]
        if yi >= sfixxset[0] - 0.025 and yi <= sfixxset[0] + 0.025:
            if xi >= yset1[0, 0]:
                yset1[0, 0] = xi
                yset1[0, 1] = yi
        if yi >= sfixxset[1] - 0.025 and yi <= sfixxset[1] + 0.025:
            if xi >= yset2[0, 0]:
                yset2[0, 0] = xi
                yset2[0, 1] = yi
        if yi >= sfixxset[2] - 0.025 and yi <= sfixxset[2] + 0.025:
            if xi >= yset3[0, 0]:
                yset3[0, 0] = xi
                yset3[0, 1] = yi
        if xi >= sfixyset[0] - 0.025 and xi <= sfixyset[0] + 0.025:
            if yi >= xset1[0, 1]:
                xset1[0, 0] = xi
                xset1[0, 1] = yi
        if xi >= sfixyset[1] - 0.025 and xi <= sfixyset[1] + 0.025:
            if yi >= xset2[0, 1]:
                xset2[0, 0] = xi
                xset2[0, 1] = yi
        if xi >= sfixyset[2] - 0.025 and xi <= sfixyset[2] + 0.025:
            if yi >= xset3[0, 1]:
                xset3[0, 0] = xi
                xset3[0, 1] = yi
        if xi >= sfixyset[3] - 0.025 and xi <= sfixyset[3] + 0.025:
            if yi >= xset4[0, 1]:
                xset4[0, 0] = xi
                xset4[0, 1] = yi

    sresult[0, 2] = yset1[0, 0]
    sresult[0, 3] = yset1[0, 1]
    sresult[1, 2] = yset2[0, 0]
    sresult[1, 3] = yset2[0, 1]
    sresult[2, 2] = yset3[0, 0]
    sresult[2, 3] = yset3[0, 1]
    sresult[3, 2] = xset1[0, 0]
    sresult[3, 3] = xset1[0, 1]
    sresult[4, 2] = xset2[0, 0]
    sresult[4, 3] = xset2[0, 1]
    sresult[5, 2] = xset3[0, 0]
    sresult[5, 3] = xset3[0, 1]
    sresult[6, 2] = xset4[0, 0]
    sresult[6, 3] = xset4[0, 1]
    
    finre[numod-1, :] = [yset1[0, 0], yset1[0, 1], yset2[0, 0], yset2[0, 1],
                         yset3[0, 0], yset3[0, 1], xset1[0, 0], xset1[0, 1],
                         xset2[0, 0], xset2[0, 1], xset3[0, 0], xset3[0, 1],
                         xset4[0, 0], xset4[0, 1]]


filename = workpath.lstrip('/root/msys_sandbox/')
filename = filename.rstrip('/Position/')
filename = '/home/zheng/Desktop/MonitorData/' + filename + '.xls'
writer = pd.ExcelWriter(filename)
dataframe = pd.DataFrame(finre)
dataframe.to_excel(writer)
writer.save()
