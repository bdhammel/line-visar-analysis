import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import sys
sys.path.append("/Users/bdhammel/Documents/programming/research_programming/ellipse-fitting")
from ellipses import LSqEllipse
import numpy as np


plt.ion()
plt.close('all')
img = plt.imread("/Users/bdhammel/Desktop/SimulatedVISARImage.png")
#img = plt.imread("/Users/bdhammel/Documents/research_programing/visar_scripts/media/smith.png")
plt.imshow(img)

imsin = img[:-30]
imcos = img[7:-23]
imnegsin = img[15:-15]
imnegcos = img[23:-7]

plt.figure()
plt.plot(imsin[300,:])
plt.plot(imcos[300,:])
plt.plot(imnegsin[300,:])
plt.plot(imnegcos[300,:])

subsin = (imsin - imnegsin)/2
subcos = (imcos - imnegcos)/2

plt.figure()
plt.plot(subsin[300,:])
plt.plot(subcos[300,:])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(subsin[300,:], subcos[300,:])
ax.axis('equal')
lsqe = LSqEllipse()
lsqe.fit((subsin[300,:], subcos[300,:]))
center, width, height, phi = lsqe.parameters()
ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                       edgecolor='r', fc='None', lw=2, label='Fit', zorder=20)
ax.add_patch(ellipse)

wrapped = np.zeros_like(subsin)
for row in range(subsin.shape[0]):
    wrapped[row] = np.arctan(subsin[row,:]/subcos[row,:])

plt.figure()
plt.plot(wrapped[300,:])

dif = np.diff(wrapped, axis=1)
dif[np.where(dif > 1.5)] -= np.pi
dif[np.where(dif < -1.5)] += np.pi

unwrapped_phase = np.cumsum(dif, axis=1)
plt.figure()
plt.plot(unwrapped_phase[300,:])

plt.figure()
plt.imshow(unwrapped_phase)
