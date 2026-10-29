import gridfunction

class State:
    
    def __init__(self, x, name):
        print("Setting up state", name)
        self.name = name
        self.x = x
        
        self.rho = gridfunction.GridFunction(x, "rho")
        
    def equals(self, state):
        self.rho.equals(state.rho)
        
    def plus(self, factor, state):
        self.rho.plus(factor, state.rho)
        
    def plot(self, t):
        self.rho.plot(t)
        # self.psi.plot()
