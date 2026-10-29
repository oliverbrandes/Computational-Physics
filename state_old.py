import gridfunction

class State:
    
    def __init__(self, x, name):
        print("Setting up state", name)
        self.name = name
        self.x = x
        
        self.phi = gridfunction.GridFunction(x, "phi")
        self.psi = gridfunction.GridFunction(x, "psi")
        
        
    def equals(self, state):
        self.phi.equals(state.phi)
        self.psi.equals(state.psi)
        
    def plus(self, factor, state):
        self.phi.plus(factor, state.phi)
        self.psi.plus(factor, state.psi)
        
    def plot(self):
        self.phi.plot()
        # self.psi.plot()
