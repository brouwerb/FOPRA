import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import curve_fit

# Load the images
#images = [Image.open(f'/home/benedict/Studium/FOPRA/Microfluidics/Analysierte_Pics/Marked_Image_{i}.jpg') for i in range(31)]

plt.rcParams['text.usetex'] = True
slices = []
for i in range(31):
    print(f'Loading and slicing image {i+1} of 30...')
    img = Image.open(f'./Microfluidics/Pics/Daten_t{str(i+1).zfill(2)}_c02.jpg')
    slices.append(np.array(img)[105:1310, img.width//2-20:img.width//2+20])
    img.close()

#slices[0] = slices[0][:, :20]
#slices[30] = slices[30][:, :20]


positions = [0.17275747508305647, 0.3289036544850498, 0.4285714285714286, 0.5049833887043189, 0.5730897009966778, 0.6345514950166113, 0.6926910299003323, 0.7441860465116279, 0.7940199335548173, 0.8372093023255814, 0.8803986710963455, 0.9235880398671097, 0.9667774086378738, 0.9983388704318937, 1.0315614617940199, 1.064784053156146, 1.0980066445182723, 1.1295681063122924, 1.159468438538206, 1.1926910299003322, 1.2192691029900333, 1.2508305647840532, 1.282392026578073, 1.3106312292358804, 1.3388704318936877, 1.3687707641196014, 1.4019933554817277, 1.4335548172757475, 1.4651162790697674, 1.5016611295681064, 1.5348837209302326]

x = np.array(positions)
t = np.array(range(0, 61, 2))*60

# fit to diffusion curve

def func(x, d, v):
    return x**2/(2*d)-v

popt, pcov = curve_fit(func, x, t, p0=[1e-10, 0.01])

#inverse

def infunc(t, d, v):
    return np.sqrt(2*d*(t+v))



print(f'Diffusion coefficient: {popt[0]*1e-6}', f'Error: {np.sqrt(pcov[0][0])*1e-6}')



# Display the slices as a graph
fig, ax = plt.subplots()
ax.imshow(np.hstack(slices), extent=[-1, 61, 2, 0], aspect='auto')
ax.set_xticks(range(0, 61, 2))
ax.scatter(range(0, 61, 2), positions, color='red')

# add label with fitted diffusion coefficient rounded with error
ax.text(0.05, 0.95, f'Diffusion coefficient: ${round(popt[0]*1e4, 4)} \pm {round(np.sqrt(pcov[0][0])*1e4, 4)} \cdot 10^{{-10}} \, \\textrm{{m}}^{{2}} \, \\textrm{{s}}^{{-1}}$', transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=dict(boxstyle="square", facecolor='white', alpha=0.5))
#ax.set_xticklabels(range(0, 31))

#plot the fitted curve
ax.plot(func(np.linspace(0, infunc(61*60, *popt), 100), *popt)/60, np.linspace(0, infunc(60*60, *popt), 100), color='red')

ax.set_xlim(-1, 61)

ax.set_xlabel('Time (minutes)')
ax.set_ylabel('Distance (mm)')
plt.show()

fig.savefig('./Microfluidics/Finalplot.pdf', bbox_inches='tight')
