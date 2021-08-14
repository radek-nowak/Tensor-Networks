import numpy as np
from ncon import ncon

class Hamiltonian:

    def __init__(self, L, J, g):
        self.L = L
        self.J = J
        self.g = g
        self.sigma_x = np.array([[0.,1.],[1.,0,]])
        self.sigma_z = np.array([[1.,0.],[0.,-1.]])
        self.id = np.eye(2)
        self.ising
        
    def ising(self,d):
        self.d = d
        sx, sz, id = self.sigma_x, self.sigma_z, self.id
        nbonds = self.L
        h = []
        for i in range(nbonds):
            bonds = -self.J*np.kron(sz,sz) - self.g*np.kron(sx,id) - self.g*np.kron(id,sx)
            h.append(np.reshape(bonds,[d,d,d,d]))
            #h.append(bonds)
        return h
    