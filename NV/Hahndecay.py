#import stufff
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy import sparse
from scipy.sparse import linalg
import numpy as np
from numpy.linalg import norm
from matplotlib.backend_bases import MouseButton
from scipy.optimize import curve_fit

data = [[], [], []]

data[0] = pd.read_csv('NV/Berta/FoPra-HahnEcho-Revival-600t1-100-1100-data.csv', sep="       ").to_numpy().T
data[1] = pd.read_csv('NV/Berta/HahnEcho-t1000-500-1500-data.csv', sep="       ").to_numpy().T
data[2] = pd.read_csv('NV/Berta/FoPra-HahnEcho-Revival-1500t1-1000-2000-data.csv', sep="       ").to_numpy().T

def fitfunc(x, a, b, x0, h):
    return a*np.exp(-b*(x-x0)**2) + h

times = np.array([600, 1000, 1500])

#create x and y data
x = [[], [], []]
y = [[], [], []]
popt = [[], [], []]
err = [[], [], []]

for i in range(len(data)):
#cast 2 and 3 colums of pandas dataframe to numpy array
    x[i] = np.array(data[i][0][1:])
    y[i] = np.array(data[i][1][1:])


    hp = np.min(y[i])
    bp = 1/times[i]
    ap= np.max(y[i])
    x0p = times[i]


    popt[i], pconv = curve_fit(fitfunc, x[i], y[i], p0 = [ap, bp, x0p, hp])
    print([ap, bp, x0p, hp])
    print(popt)
    err[i] = np.sqrt(np.diag(pconv))

    print("Maximum bei", popt[i][2],"pm", err[i][2], "ns")


#plot all
fig, axs = plt.subplot_mosaic([["a", "b", "c"], ["hight", "hight", "hight"]])
#share y axis
axs["b"].sharey(axs["a"])
axs["c"].sharey(axs["a"])

for i in range(len(axs.items())-1):
    axs[list(axs.keys())[i]].plot(x[i], y[i], label='Messwerte')
    axs[list(axs.keys())[i]].plot(x[i], fitfunc(x[i], *popt[i]), label='Fit')
    axs[list(axs.keys())[i]].set_title(f'Hahn Echo bei t1={times[i]}ns')
    axs[list(axs.keys())[i]].set_xlabel('Zeit in ns')
    axs[list(axs.keys())[i]].set_ylabel('Intensität')

hight = np.array([popt[0][0], popt[1][0], popt[2][0]])
errhight = np.array([err[0][0], err[1][0], err[2][0]])

#fit exponential
def fitexp(x, a, b):
    return a*np.exp(-b*x)

popt, pconv = curve_fit(fitexp, times, hight, p0=[np.max(hight), 1/np.max(times)])
errexp = np.sqrt(np.diag(pconv))
T2 = 1/popt[1]
axs["hight"].errorbar(times, hight, yerr=errhight, fmt='o')
axs["hight"].plot(np.arange(np.min(times)-50, np.max(times)+50, 2), fitexp(np.arange(np.min(times)-50, np.max(times)+50, 2), *popt), label='Fit with T2 = '+str(round(T2, 2))+'ns')
axs["hight"].set_xlabel('t1 in ns')
axs["hight"].set_ylabel('Stärke des Echos')
axs["hight"].set_title('Stärke des Echos in Abhängigkeit von t1')
axs["hight"].legend()

plt.tight_layout()

plt.show()


