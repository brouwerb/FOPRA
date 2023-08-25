import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy import sparse
from scipy.sparse import linalg
import numpy as np
from numpy.linalg import norm
import peakutils
from peakutils.plot import plot as pplot


# import data from xlsx file

data = pd.read_excel('NV/Berta/ 8peaks.xlsx', sheet_name='Sheet1').to_numpy().T

#cast 2 and 3 colums of pandas dataframe to numpy array
x = np.array(data[1][1:])
y = -np.array(data[2][1:])

indexes = peakutils.indexes(y, thres=0.5)


ybase = peakutils.baseline(y)

#plot data
fig, ax = plt.subplots()
ax.plot(x, y, color='red', label='data')
ax.plot(x, ybase, '-', color='blue', label='baseline')
ax.plot(x[indexes], y[indexes], 'x', color='black', label='peaks')

ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Intensity (a.u.)')
ax.set_title('NV spectrum')
ax.legend()

plt.show()

