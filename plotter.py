import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

from numpy import *

data = loadtxt("RungeKutta.txt")
data = transpose(data)

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[0], data[1], data[2], label='K = 20.0')
ax.set_xlabel('Ex')
ax.set_ylabel('Ey')
ax.set_zlabel('N')
ax.legend()

plt.show()

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[0][1500:5000], data[1][1500:5000], data[2][1500:5000], label='K = 20.0')
ax.legend()
ax.set_xlabel('Ex')
ax.set_ylabel('Ey')
ax.set_zlabel('N')
plt.show()

data1 = loadtxt("MagVsTime.txt")
plt.scatter(data[0], data1[1])
#plt.axis((0.0,1.0,-3.1,3.1))
plt.ylabel('Mag')
plt.xlabel('Time')
plt.show()


data1 = loadtxt("MagVsTime.txt")
data1 = transpose(data1)
plt.scatter(data1[0], data1[1])
plt.axis((0.0,1.0,-3.1,3.1))
plt.ylabel('Mag')
plt.xlabel('Time')
plt.show()


data2 = loadtxt("LocalMax.txt")
data2 = transpose(data2)
plt.scatter(data2[0], data2[1], s=2, facecolor='0.0', lw = 0)
plt.ylabel('Max')
plt.xlabel('Time')
plt.show()
