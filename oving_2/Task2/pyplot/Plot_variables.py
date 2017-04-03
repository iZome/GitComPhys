from matplotlib import pyplot as plt
import numpy as np
import sys

sys.path.append("/home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/data")

data = np.loadtxt("data/P_inf_of_p.csv", delimiter = ",")
data2 = np.loadtxt("data/s_of_p.csv", delimiter = ",")
data3 = np.loadtxt("data/P_infSquared_of_p.csv", delimiter = ",")
data4 = np.loadtxt("data/giantComponents.csv", delimiter = ",")
data5 = np.loadtxt("data/weightedAvgClusterSize.csv", delimiter = ",")

x = np.linspace(0, 1, len(data))
x3 = np.linspace(0, 1, len(data2))
x2 = np.linspace(0, 1, len(data4))

plt.figure(1)
ax = plt.subplot(221)
ax.set_title("P inf before convo")
plt.plot(x2,data4)
ax = plt.subplot(222)
ax.set_title("P inf after convo")
plt.plot(x,data)
ax = plt.subplot(223)
ax.set_title("s before convo")
plt.plot(x2, data5)
ax = plt.subplot(224)
ax.set_title("s after convo")
plt.plot(x3,data2)
plt.show()
