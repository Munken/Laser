from numpy import *
import os
import matplotlib.pyplot as plt

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
# ax1.semilogx(data[:,1],data[:,2])

ax2 = plt.axes([.65, .6, .2, .2], axisbg='y')
# ax2.semilogx(data[3:8,1],data[3:8,2])
plt.setp(ax2, xticks=[], yticks=[])

plt.show()