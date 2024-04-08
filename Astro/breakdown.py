import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backend_bases import MouseButton
from uncertainties import unumpy as unp
from inspect import getsourcefile
import os.path as path, sys
current_dir = path.dirname(path.abspath(getsourcefile(lambda:0)))
sys.path.insert(0, current_dir[:current_dir.rfind(path.sep)])
from AP import *

def gaussianfit(x, *args):
    y=0
    for i in range(numberofpeaks):
        y = y + args[3*i]*np.exp(-(((x-args[3*i+1])/args[3*i+2])**2))
    return y

def gaussian(x, a, x0, sig):
    return a*np.exp(-((x-x0)*sig)**2)




#import data for all voltages
voltages = [53, 54, 55, 56, 57, 58, 59, 60]

data = []
for i in voltages:
    data.append(pd.read_csv('Astro/data/Bertha-{}V_histo.txt'.format(i), header=None, sep="\t" , decimal=",").to_numpy().T)

fits = []
uncertainties = []
inputs = []

# for i, voltage in enumerate(voltages):
#     x = data[i][0]
#     y = data[i][1]
#     x = x[y > 5]
#     y = y[y > 5]
#     # Plot histogram
#     fig, ax = plt.subplots()
#     ax.plot(x, y, color='red', label='data')
#     ax.set_xlabel('Height (a.u.)')
#     ax.set_ylabel('Counts)')
#     ax.set_title('Histogram for Voltage {}V'.format(voltage))
#     ax.legend()

#     #  set x limits, exclude values with y = 0
 



    
    
#     # Wait for click to select starting values for Gaussian fit
#     marked = []
#     def onclick(event):
#         if event.button == MouseButton.LEFT:
#             print('x = %f, y = %f' % (event.xdata, event.ydata))
#             marked.append([event.xdata, event.ydata])
#             ax.plot(event.xdata, event.ydata, 'x', color='black', label='peaks')
#             fig.canvas.draw()
#     cid = fig.canvas.mpl_connect('button_press_event', onclick)
#     plt.show()
    
#     # Perform Gaussian fit
#     numberofpeaks = len(marked)
#     p0 = []
#     for point in marked:
#         p0.append(point[1])  # Amplitude
#         p0.append(point[0])  # Mean
#         p0.append(5)      # Standard deviation
#     popt, pconv = curve_fit(gaussianfit, x, y, p0=p0)
#     fits.append(popt)
#     uncertainties.append(np.sqrt(np.diag(pconv)))
#     inputs.append(marked)
#     # Plot fit
#     fig, ax = plt.subplots()
#     ax.plot(x, y, color='red', label='data')
#     ax.plot(x, gaussianfit(x, *popt), color='blue', linestyle='--', label='fit')
#     ax.set_xlabel('Height (a.u.)')
#     ax.set_ylabel('Counts')
#     ax.set_title('Fit for Voltage {}V'.format(voltage))
#     ax.legend()
#     plt.show()

# DATA:
Inputs = [[[196.86001169354836, 4093.7702380952396], [377.83152379032276, 1600.3059523809534], [12.269069354838848, 6190.547023809527]], [[290.0109266129032, 4307.507142857144], [571.5812576612902, 2149.0723086734693], [24.763513306451614, 3531.699604591837], [849.0708592741935, 689.6323852040817]], [[387.8222383064517, 5154.737500000001], [782.1096159274196, 2748.300000000001], [-0.9117959677419094, 4012.112500000001], [1165.290306854839, 965.1125000000004]], [[478.70688689516123, 6069.476190476192], [934.417891935484, 3786.380952380953], [-51.57500987903222, 3866.023809523811], [1456.4141340725805, 1556.380952380953]], [[580.275222580645, 1882.3642857142866], [1135.0770596774191, 1152.0428571428577], [16.069964516128948, 1142.6797619047625], [1708.685738709677, 515.3523809523813]], [[682.6560893145161, 2891.800000000001], [1367.6864695564514, 1951.6125000000006], [2052.716849798387, 934.6750000000006], [-15.060038709677428, 1740.5500000000006]], [[783.0213709677419, 3531.853571428572], [1527.9320010080646, 2527.8416666666667], [2372.164048387097, 1139.9428571428575], [-28.10353729838721, 1789.5976190476194]], [[910.5545580645162, 3180.2000000000003], [1752.6398290322577, 2171.4500000000003], [1.9888709677420593, 1835.2000000000003], [2483.9244064516124, 490.2000000000003]]]
fits =  [np.array([3813.78900442,  187.53223788,   67.12136834, 1399.29801287,
        363.03333467,  135.50577776, 6200.9179252 ,   14.33238178,
         63.24411993]), np.array([4118.31040398,  291.55191175,   70.55369348, 2039.80768392,
        566.34534206,   90.89413624, 2503.06364647,   23.50346038,
         60.10093737,  599.48911036,  843.92574814,  108.95102143]), np.array([5058.38258345,  391.48228539,   74.16150106, 2599.97267577,
        761.10976505,  103.89903548, 3075.75692168,   24.99432444,
         62.27008512,  804.82165244, 1125.38079234,  153.44730215]), np.array([5648.09142699,  488.31239114,   67.36580316, 2980.71917726,
        958.65015786,   73.46629035, 3625.6677001 ,   24.85204985,
         56.59092096,  973.63352052, 1037.4910027 ,  852.27482179]), np.array([1684.58272952,  584.14263038,   69.59487749,  877.23046839,
       1147.15795782,   78.00948042, 1177.82478768,   22.34412463,
         53.67546759,  299.89061866, 1252.24445768, 1078.96416648]), np.array([2514.64399763,  683.31716671,   74.82908164, 1385.15642539,
       1341.89535276,   86.20088201,  585.70148545, 1560.44423814,
       1328.03196357, 2081.13745491,   23.75478727,   53.71147959]), np.array([3125.39361572,  785.47056234,   78.22808552, 1676.65009342,
       1537.06949296,   92.71356352,  799.57692347, 1846.35903181,
       1628.44129176, 2919.42510912,   23.53745732,   53.25869069]), np.array([2.42963369e+03, 8.87647622e+02, 7.83705693e+01, 1.33973555e+03,
       1.73652726e+03, 1.00703103e+02, 1.09628095e+04, 3.87831292e+00,
       4.84540476e+00, 7.67569925e+02, 1.95990578e+03, 2.42568354e+03])]
uncertainties = [np.array([76.14846509,  0.9786311 ,  1.57722218, 33.14902337,  5.60930448,
        7.26243448, 46.94900225,  0.42407885,  0.63955886]), np.array([66.16748683,  0.92754134,  1.32348023, 58.85671418,  2.23837062,
        3.35079753, 71.55123395,  1.4027238 ,  1.98431568, 54.00417363,
        8.27295648, 12.46430725]), np.array([ 96.96242269,   1.16071649,   1.64187743,  82.21637856,
         2.79625784,   4.06909991, 105.81196649,   1.74901806,
         2.47348636,  68.50562981,  10.85378723,  16.29700338]), np.array([136.17625357,   1.28443417,   1.98729026, 130.79906738,
         2.54305289,   3.94522233, 146.63664941,   1.82816347,
         2.75143728,  49.03275474,  33.07865791,  47.21236409]), np.array([44.94983505,  1.48055522,  2.24239227, 42.5952545 ,  3.01119298,
        4.58496109, 50.75393365,  1.85621515,  2.75857135, 14.08651651,
       39.60216105, 58.00641726]), np.array([78.18807796,  1.86163316,  2.79051011, 73.00253169,  3.63011628,
        5.45903894, 21.86278348, 39.74600156, 57.2125033 , 91.48092571,
        1.90389144,  2.79252715]), np.array([104.08891385,   2.09229422,   3.10456101,  95.7977406 ,
         4.24854569,   6.326685  ,  26.20820924,  43.40120207,
        62.60746171, 125.29110697,   1.84711879,   2.69301481]), np.array([6.87722748e+01, 1.79342657e+00, 2.61220004e+00, 6.09005686e+01,
       3.68708804e+00, 5.42328865e+00, 2.19097153e+02, 2.88585859e-01,
       1.83252306e-01, 1.44415986e+01, 3.79051764e+01, 5.90721678e+01])]

distances = []
# Print fits and uncertainties
for i, voltage in enumerate(voltages):
    first_peak_mean = unp.uarray(fits[i][1], uncertainties[i][1])
    second_peak_mean = unp.uarray(fits[i][4], uncertainties[i][4])
    distance = second_peak_mean - first_peak_mean
    print("Distance between first and second peak:", round_errtexU(distance))
    distances.append([voltage, distance])
print('Inputs:', inputs)
print('Fits:', fits)
print('Uncertainties:', uncertainties)

colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'black']

# Plot voltages 53, 56, and 59 in one plot with fits
fig, ax = plt.subplots()
for i, voltage in enumerate(voltages):
    if voltage in [53, 56, 59]:
        x = data[i][0]
        y = data[i][1]
        ax.plot(x, y, color=colors[i], label=' {}V'.format(voltage))
        popt = fits[i]
        numberofpeaks = int(len(popt)/3)
        ax.plot(x[:-300], gaussianfit(x, *popt)[:-300], color=colors[i], linestyle='--', label='Fit for {}V'.format(voltage))
ax.set_xlabel('Height (a.u.)')
ax.set_ylabel('Counts')
ax.set_title('Histograms for Voltages 53V, 56V, and 59V')
ax.set_xlim(0, 2000)
ax.legend()
plt.show()
#save
fig.savefig('Astro/figures/histograms.pdf')

# Extract voltages and distances
voltages = [distance[0] for distance in distances]
distances = [distance[1] for distance in distances]

# Fit a linear function to the data
popt, pconv = np.polyfit(voltages, unp.nominal_values(distances), deg=1, cov=True, w=1/unp.std_devs(distances))

# Calculate x for y=0
x1 = ufloat(popt[0], np.sqrt(pconv[0, 0]))
x2 = ufloat(popt[1], np.sqrt(pconv[1, 1]))
x_intercept = -x2/x1


# Print the result
print("Breakdown Voltage:", round_errtexU(x_intercept))

# Plot distances with error bars
fig, ax = plt.subplots()
ax.errorbar(voltages, unp.nominal_values(distances), yerr=unp.std_devs(distances), fmt=".", color='red', label='Data', linestyle='None', elinewidth=10)
x = np.linspace(51, 61, 1000)
ax.plot(x, np.polyval(popt, x), color='blue', linestyle='--', label='Fit with Breakdown Voltage at {}V'.format(round_errtexU(x_intercept)))
ax.axhline(0, color='grey', linestyle='--')
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Distance')
ax.set_title('Distance between first and second peak')
ax.legend()
plt.show()
#save
fig.savefig('Astro/figures/distance.pdf')
