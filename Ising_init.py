import numpy as np
from ncon import ncon

class Hamiltonian:

    def __init__(self, L, J):
        self.L = L
        self.J = J
        self.sigma_z = np.array([[1.,0.],[0.,-1.]])
        self.one = np.eye(2)
        #self.tensor_L    
        #self.ising

    # tensor L should have 3 indices
    def tensor_L(self, beta):
        sz, one = self.sigma_z, self.one
        return np.array([one*np.sqrt(np.cosh(beta)), sz*np.sqrt(np.sinh(beta))])
    
    # needs testing
    # Ls commute
    # this is for corner tensor
    def tensor_T(self, beta):
        tensor_L = self.tensor_L(beta)
        index_order_1 = [[-1,-2,1],[1,-3,-4]]
        index_order_2 = [[-1,-2,1,-3],[1,-4,-5,-6]]
        L1, L2 = tensor_L, tensor_L
        L3, L4 = tensor_L, tensor_L
        L12 = ncon([L1,L2],index_order_1)
        L34 = ncon([L3,L4],index_order_1)
        T = ncon([L12,L34], index_order_2)
        return T

    def tensor_O(self, beta):
        T = self.tensor_T(beta)
        return ncon(T,[1,-1,-2,-3,-4,1])
        

    def ising(self, beta):
        n_sites = self.L
        H = []
        for i in range(n_sites):
            H.append(self.tensor_O(beta))
        return np.array(H)
