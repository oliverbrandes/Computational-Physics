import numpy as np
import matplotlib.pyplot as plt

class GridFunction:
    
    def __init__(self, x, name):
        # print("Creating grid function", name)
        self.x = x
        self.name = name
        self.delta = x[1] - x[0]
        self.n = np.size(x)
        self.data = np.zeros((self.n))
        self.sec_deriv = np.zeros(self.n)
        
    def __getitem__(self, index):
        return self.data[index]
    
    def __setitem__(self, index, value):
        self.data[index] = value
        
    def equals_const(self, value):
        for i in range(self.n):
            self.data[i] = value
            
    def equals(self, gf):
        for i in range(self.n):
            self.data[i] = gf.data[i]
            
    def plus(self, factor, gf):
        for i in range(self.n):
            self.data[i] += factor * gf.data[i]
    
    def dx(self, i):
        return (self.data[self.ind(i+1)] - self.data[self.ind(i-1)]) / (2.0 * self.delta)

    def ind(self,i):
        # implements periodic boundary conditions
        if (i == -1):
            return self.n - 1
        elif (i == self.n):
            return 0
        else:
            return i
        
    def plot(self, t):
        time = "t = " + "{:.2f}".format(t)
        plt.title(time)
        plt.plot(self.x, self.data, marker = "o")
        plt.xlabel("x", size = 14)
        plt.ylabel(self.name, size = 14)
        plt.grid(True)
        plt.show()

    def name():
        return self.name    
