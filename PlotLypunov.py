import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

from numpy import *

'''
dataY1 = loadtxt("LyapunovSlice0.txt")
dataX = np.linspace(21,24, dataY1[:,0].size)

plt.plot(dataX, dataY1[:,0], 'b', linewidth=1.0)
plt.plot(dataX, dataY1[:,1], 'g', linewidth=1.0)
'''

Param = loadtxt("Parameters.txt")

xlow  = Param[6]
xhigh = Param[7]

dataY = loadtxt("LyapunovSlice.txt")

dataX = np.linspace(xlow, xhigh, dataY.size)

plt.plot(dataX, dataY, 'b', linewidth=1.0)

plt.plot((xlow,np.amax(dataX)),(0.0,0.0),'-k', linewidth = 0.5)
plt.ylabel('Lyapunov')
plt.xlabel('K')
#plt.axis((21.0,24.0,-5,12))
#plt.text(50, -1.5, 'Delta = 0.0', fontsize=15)
#plt.axis((0.0,100.0,-0.1,0.1))
plt.show()
