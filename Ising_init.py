import numpy as np
from ncon import ncon

class Hamiltonian:

    def __init__(self, L, J, g):
        self.L = L
        self.J = J
        self.g = g
        #self.sigma_x = np.array([[0.,1.],[1.,0,]])
        self.sigma_z = np.array([[1.,0.],[0.,-1.]])
        #self.one = np.eye(2)
        self.ising
        
    def ising(self,d):
        self.d = d
        #sx, sz, one = self.sigma_x, self.sigma_z, self.one
        sz = self.sigma_z
        nbonds = self.L
        h = []
        for i in range(nbonds):
            bonds = -self.J*np.kron(sz,sz)# - self.g*np.kron(sx,one) - self.g*np.kron(one,sx)
            h.append(bonds.reshape(d,d,d,d))
        return h
    
