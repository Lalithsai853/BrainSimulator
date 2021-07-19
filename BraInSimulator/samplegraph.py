import os
import random
import sys
import numpy as np
import scipy as scpy
import scipy.integrate as spin
from matplotlib import pyplot as plt
import random
import matplotlib.patches as mpatches
import pandas
import math

regularfile = open('regular.txt', 'r')
seizurefile = open('seizure.txt', 'r')
num = 1
totalnum = 1000
t = np.linspace(0, totalnum, totalnum)
count = 0
data1 = np.ones(totalnum)
data2 = np.ones(totalnum)
for line in regularfile:
    if line.startswith("-"):
        num = -1
        line = line.lstrip("-")
        line = line.rstrip('\n')
        num = num * int(line)
    else:
        num = num * int(line.strip('\n'))
    data1[count] = num
    count += 1
    num = 1

count = 0
for line in seizurefile:
    if line.startswith("-"):
        num = -1
        line = line.lstrip("-")
        line = line.rstrip('\n')
        num = num * int(line)
    else:
        num = num * int(line.strip('\n'))
    data2[count] = num
    count += 1
    num = 1

plt.plot(t, data1, 'b.-')
plt.axis([0,1000,-150,150])
plt.grid(True)
plt.title("Graph of regular data")
plt.show()

plt.plot(t, data2, 'r.-')
plt.axis([0,1000,-150,150])
plt.grid(True)
plt.title("Graph of seizure data")
plt.show()
