import mumpy as np
import h5py
# h5 file is a type of multiple dimensions file. In this file, data would be saved in 
# different subsets. All the data are format as 1 dimensional array. If that data are multiple
# dimensions, use numpy.array.reshape to reshape the array.

f = h5py.File('filename', mode = 'r')
# read h5 file
print(f.key())
# output the labels of h5 file

ns = f.get('Sigma')
# write data from label 'Sigma' in ns

