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


def baseline_arPLS(y, ratio=1e-6, lam=100, niter=10, full_output=False):
    L = len(y)

    diag = np.ones(L - 2)
    D = sparse.spdiags([diag, -2*diag, diag], [0, -1, -2], L, L - 2)

    H = lam * D.dot(D.T)  # The transposes are flipped w.r.t the Algorithm on pg. 252

    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)

    crit = 1
    count = 0

    while crit > ratio:
        z = linalg.spsolve(W + H, W * y)
        d = y - z
        dn = d[d < 0]

        m = np.mean(dn)
        s = np.std(dn)

        w_new = 1 / (1 + np.exp(2 * (d - (2*s - m))/s))

        crit = norm(w_new - w) / norm(w)

        w = w_new
        W.setdiag(w)  # Do not create a new matrix, just update diagonal values

        count += 1

        if count > niter:
            print('Maximum number of iterations exceeded')
            break

    if full_output:
        info = {'num_iter': count, 'stop_criterion': crit}
        return z, d, info
    else:
        return z



# import data from xlsx file

data = pd.read_excel('NV/Berta/ Test1.xlsx', sheet_name='Sheet1').to_numpy().T

#cast 2 and 3 colums of pandas dataframe to numpy array
x = np.array(data[1][1:])
y = -np.array(data[2][1:])

#baseline correction
ybase = baseline_arPLS(y)


#plot data
fig, axs = plt.subplots(2)
axs[0].plot(x, -y, color='red', label='data')
axs[0].plot(x, -ybase, '-', color='blue', label='baseline')
#axs[0].plot(x[indexes], y[indexes], 'x', color='black', label='peaks')

axs[0].set_xlabel('Frequency (GHz)')
axs[0].set_ylabel('Intensity (a.u.)')
axs[0].set_title('NV spectrum')
axs[0].legend()

#plot baseline corrected data
axs[1].plot(x, -y+ybase, color='red', label='data')
axs[1].set_xlabel('Frequency (GHz)')
axs[1].set_ylabel('Intensity (a.u.)')
axs[1].set_title('NV spectrum - baseline')
axs[1].legend()

# select points by mouse click and mark them with a cross





marked = []
def onclick(event):
    if event.button == MouseButton.LEFT:
        print('x = %f, y = %f' % (event.xdata, event.ydata))
        marked.append([event.xdata, event.ydata])
        axs[1].plot(event.xdata, event.ydata, 'x', color='black', label='peaks')
        

        fig.canvas.draw()



cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()

numberofpeaks = len(marked)

#marked = [[2.6136060087564794, -6531.950040046049], [2.691105882944995, -17860.102775250758], [2.830605656484324, -21636.153686985665], [2.8895055608675966, -10150.665497125334], [2.942205475315787, -10622.671861092196], [2.993355392280208, -20692.14095905194], [3.1034052136279007, -15500.070955416446], [3.143705148205929, -5745.272766767946]]

def gaussianfit(x, *args):
    y=0
    for i in range(numberofpeaks):
        y = y + args[3*i]*np.exp(-(((x-args[3*i+1])/args[3*i+2])**2))
    return y

def gaussian(x, a, x0, sig):
    return a*np.exp(-((x-x0)*sig)**2)


pO = []
for i in marked:
    pO.append(i[1])
    pO.append(i[0])
    pO.append(0.01)

popt, pconv = curve_fit(gaussianfit, x, -y+ybase, p0 = pO)

print(popt)

#sort data into hight and position
height = []
position = []
for i in range(int(len(popt)/3)):
    height.append(popt[3*i])
    position.append(popt[3*i+1])

print("Hight:", height)
print("Position:", position)

#calculate middle position
middle = np.mean(position)
print("Splitting:", position-middle)

    


fig, axs = plt.subplots(2)
axs[0].plot(x, -y, color='red', label='data')
axs[0].plot(x, -ybase, '-', color='blue', label='baseline')
axs[0].plot(x, gaussianfit(x, *popt)-ybase, '-', color='green', label='fit')
axs[0].legend()

axs[1].plot(x, -y+ybase, color='red', label='data, corrected')
axs[1].plot(x, gaussianfit(x, *popt), '-', color='green', label=f'fit, x0 = {round(position[0], 5)} GHz')
axs[1].set_xlabel('Frequency (GHz)')
axs[1].set_ylabel('Intensity (a.u.)')
axs[1].legend()

plt.show()

