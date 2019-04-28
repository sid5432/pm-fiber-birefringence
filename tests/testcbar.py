#!/usr/bin/python
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy
import sys

fig = plt.figure()
plot = fig.add_subplot(111)
scatter = plot.scatter([1, 2], [3, 4], cmap=cm.spring, color="red")
scatter.set_array(numpy.array([5, 6]))

# first time....
# cb = fig.colorbar(scatter, use_gridspec=True)
cb = fig.colorbar(scatter, use_gridspec=False)
fig.delaxes(cb.ax)
fig.subplots_adjust(right=0.90)

for i in range(10):
    cb = fig.colorbar(scatter, use_gridspec=False)
    # cb = fig.colorbar(scatter)
    fig.delaxes(cb.ax)
    fig.subplots_adjust(right=0.90)

# last time keep
cb = fig.colorbar(scatter, use_gridspec=False)
# cb = fig.colorbar(scatter)
# fig.subplots_adjust(right=0.80)

plt.draw()
plt.show()
