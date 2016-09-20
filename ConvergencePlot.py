import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

from numpy import *


data = loadtxt("Check_Conver.txt")
data = transpose(data)
m = data[0]
print m.size

plt.scatter(np.log10(data[0]), data[1], s=1, facecolor='0.8', lw = 1.0)
#plt.plot((0.0,np.amax(dataX)),(0.0,0.0),'-k')
plt.ylabel('Lyapunov')
plt.xlabel('steps')
#plt.axis((3.0,10.0,6.35,6.45))
plt.show()
