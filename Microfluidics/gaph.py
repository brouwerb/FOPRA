import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Load the images
#images = [Image.open(f'/home/benedict/Studium/FOPRA/Microfluidics/Analysierte_Pics/Marked_Image_{i}.jpg') for i in range(31)]

slices = []
for i in range(30):
    print(f'Loading and slicing image {i+1} of 30...')
    img = Image.open(f'/home/benedict/Studium/FOPRA/Microfluidics/Analysierte_Pics/Marked_Image_{i}.jpg')
    slices.append(np.array(img)[105:1310, img.width//2-20:img.width//2+20])
    img.close()


# Display the slices as a graph
fig, ax = plt.subplots()
ax.imshow(np.hstack(slices), extent=[0, 60, 2, 0], aspect='auto')
ax.set_xticks(range(0, 61, 2))
#ax.set_xticklabels(range(0, 31))
ax.set_xlabel('Time (minutes)')
ax.set_ylabel('Distance (mm)')
plt.show()
