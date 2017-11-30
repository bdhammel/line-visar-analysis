import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import sys
sys.path.append("/Users/bdhammel/Documents/programming/research_programming/ellipse-fitting")
from ellipses import LSqEllipse
import numpy as np


plt.ion()
plt.close('all')
img = plt.imread("/Users/bdhammel/Desktop/SimulatedVISARImage.png")
plt.imshow(img)

plt.figure()
plt.plot(img[:,100])

plt.figure()
plt.plot(img[350,:])
plt.plot(img[357,:])

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.plot(img[350,:], img[357,:])
lsqe = LSqEllipse()
lsqe.fit((img[350,:], img[357,:]))
center, width, height, phi = lsqe.parameters()
ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                       edgecolor='r', fc='None', lw=2, label='Fit', zorder = 2)
ax.add_patch(ellipse)
ax.axis('equal')
