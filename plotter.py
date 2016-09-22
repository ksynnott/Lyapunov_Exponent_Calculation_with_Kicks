import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math as ma
import matplotlib.pyplot as plt

from numpy import *

Lam = 10
Dp  = 2

Am = ma.sqrt(Lam - 1 - Dp*Dp)
Nn = 1 + Dp*Dp

data = loadtxt("AddedValues.txt")

print data.size
print data[0].size

Nd = data[0].size
Nd_2 = 1000

x  = np.linspace(0, 2*np.pi, Nd_2)
zz = np.linspace(Nn, Nn, Nd_2)
xx = Am*np.sin(x)
yy = Am*np.cos(x)

#data = transpose(data)

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[0], data[1], data[2], '-g')
ax.plot(xx, yy, zz, '-r')
#ax.plot(data[3], data[4], data[5])
ax.set_xlabel('Ex')
ax.set_ylabel('Ey')
ax.set_zlabel('N')

plt.show()

plt.plot(data[0], data[1], '-g')
plt.plot(data[0][0], data[1][0], '+g', mew=5, ms=10)
plt.plot(xx, yy, '-r')
plt.show()

"""
data = loadtxt("bois.txt")
#data = transpose(data)

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[0], data[1], data[2])
ax.set_xlabel('Ex')
ax.set_ylabel('Ey')
ax.set_zlabel('N')

plt.show()

plt.plot(data[0], data[1])
plt.show()

"""

"""
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
"""
