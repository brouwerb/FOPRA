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

data = pd.read_excel('NV/Berta/ FoPra-Rabi-No1.xlsx', sheet_name='Sheet1').to_numpy().T

#cast 2 and 3 colums of pandas dataframe to numpy array
x = np.array(data[1][1:])
y = -np.array(data[2][1:])


#delete first two data points
x = x[10:]
y = -y[10:]

print(x, y)

def fitfunc(x, a, b, x0, h, e):
    return a*np.sin(b*(x)+x0)*np.exp(-e*x) + h

hp = np.mean(y)
bp = 2*np.pi/350
ap= (np.max(y)-np.min(y))/2


popt, pconv = curve_fit(fitfunc, x, y, p0 = [ap, bp, 0, hp, 0], bounds = ([ap*0.2, bp*0.2, 0, hp*0.2, 0], [ap*2, bp*2, 2*np.pi, hp*2, 1]))
print(popt)
freq = popt[1]*1000/np.pi
halb = np.log(2)/popt[4]
print("Frequenz in Mhz:", freq)
print("Halbwertszeit in ns:", halb)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, fitfunc(x, *popt))
ax.plot(x, [popt[3] for i in x])

plt.show()