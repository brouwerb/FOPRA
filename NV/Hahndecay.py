#import stufff
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy import sparse
from scipy.sparse import linalg
import numpy as np
from numpy.linalg import norm
import peakutils
from peakutils.plot import plot as pplot
from matplotlib.backend_bases import MouseButton
from scipy.optimize import curve_fit

data = [[], [], []]

data[0] = pd.read_csv('NV/Berta/FoPra-HahnEcho-Revival-600t1-100-1100-data.csv', sep="       ").to_numpy().T
data[1] = pd.read_csv('NV/Berta/HahnEcho-t1000-500-1500-data.csv', sep="       ").to_numpy().T
data[2] = pd.read_csv('NV/Berta/FoPra-HahnEcho-Revival-1500t1-1000-2000-data.csv', sep="       ").to_numpy().T
print(data)
def fitfunc(x, a, b, x0, h):
    return a*np.exp(-b*(x-x0)**2) + h

for i in range(len(data)):
#cast 2 and 3 colums of pandas dataframe to numpy array
    x = np.array(data[i][0][1:])
    y = np.array(data[i][1][1:])

    print(x, y)



    hp = np.min(y)
    bp = 1/600
    ap= np.max(y)
    x0p = 600


    popt, pconv = curve_fit(fitfunc, x, y, p0 = [ap, bp, x0p, hp])
    print([ap, bp, x0p, hp])
    print(popt)
    err = np.sqrt(np.diag(pconv))

    print("Maximum bei", popt[2],"pm", err[2], "ns")


fig, ax = plt.subplots()
ax.plot(x, y, label="data")
ax.plot(x, fitfunc(x, *popt), label="fit")

plt.show()