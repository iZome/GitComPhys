from matplotlib import pyplot as plt
import numpy as np
import sys

sys.path.append("/home/karsten/Documents/c-/6th-Semester/Computational_Physics/oving_2/Task2/data")

data = np.loadtxt("data/test1.csv", delimiter = ",")



N = int(data[0])
print(N)

#x = []
#y = []

def getXY():
    lattice = np.zeros((N,N))
    print(lattice)
    for i in range(1,len(data)):
        x = (np.floor(data[i]/N)) #maafikses
        y = (data[i]%N)
        lattice[x,y] = 1
    return lattice

lattice = getXY()

#grid = temp.reshape


plt.matshow(lattice, cmap = 'plasma')
plt.show()
