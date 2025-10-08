import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import shutil
from test_callback import *
from functions_double import open_h5, normalization

#%%File import

#Unpack the data
measurement_number='26805'
file='/0000'+measurement_number+'-WavelengthScanBaseExperiment.h5'
path = os.path.join(r'/home/benedict/Documents/Studium/FOPRA/Photonics/data/2025-10-08/15'+file)
# path = os.path.join(r'/Volumes/eqn/quantumnetworks/Data/2024/RT1/2024-07-18/16/'+file) #for mac
save_path = "Plots/"
WL, RR, sample_nr = open_h5(path)

directory = save_path + sample_nr[2:-1]
if os.path.exists(directory):
    shutil.rmtree(directory)
    
#os.makedirs(directory)


#%%figure    
fig1=plt.figure(num=1, figsize=(10,4))
fig1.clf()
ax1=fig1.add_subplot(111)

ax1.set_title('Structure')
graph, = ax1.plot(WL, RR, '.')
ax1.set_ylabel('Relative reflection [-]')   
ax1.set_xlabel('Wavelength [nm]')

# selectivefitter= SelectiveFitter(graph, 0.7)

timestamp = assign_filename(directory, '_ParamFile_54.txt')
selectivefitter= SelectiveFitter(graph, 0.7, timestamp)



fig1.tight_layout()
#ax1.grid()
#fig1.savefig(File_to_save,dpi=1200)