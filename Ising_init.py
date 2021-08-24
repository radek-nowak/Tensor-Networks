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
    def tensor_T(self, beta, x):
        if x == 'edge':
            tensor_L = self.tensor_L(beta)
            index_order_1 = [[-1,-2,1],[1,-3,-4]]
            index_order_2 = [[-1,-2,1],[1,-3]]
            index_order_3= [[-1,-2,-3,1],[1,-4,-5]]
            L1, L2, L3 = tensor_L, tensor_L, tensor_L
            L4 = tensor_L[0]
            L12 = ncon([L1,L2],index_order_1)
            L34 = ncon([L3,L4],index_order_2)
            T = ncon([L12,L34], index_order_3)
            return T
        elif x == 'corner':
            tensor_L = self.tensor_L(beta)
            index_order_1 = [[-1,-2,1],[1,-3,-4]]
            index_order_2 = [[-1,1],[1,-2]]
            index_order_3= [[-1,-2,-3,1],[1,-4]]
            L1, L2 = tensor_L, tensor_L
            L3, L4 = tensor_L[0], tensor_L[0]
            L12 = ncon([L1,L2],index_order_1)
            L34 = ncon([L3,L4],index_order_2)
            T = ncon([L12,L34], index_order_3)
            return T
        else:
            tensor_L = self.tensor_L(beta)
            index_order_1 = [[-1,-2,1],[1,-3,-4]]
            index_order_2 = [[-1,-2,1,-3],[1,-4,-5,-6]]
            L1, L2 = tensor_L, tensor_L
            L3, L4 = tensor_L, tensor_L
            L12 = ncon([L1,L2],index_order_1)
            L34 = ncon([L3,L4],index_order_1)
            T = ncon([L12,L34], index_order_2)
            return T
            
        
    def tensor_O(self, beta, x):
        T = self.tensor_T(beta, x)
        index_order = [-x for x in range(1,len(T.shape)-1)]
        index_order.insert(0,1)
        index_order.insert(len(T.shape),1)
        return ncon(T,index_order)
        

    def ising(self, beta):
        n_sites = self.L
        H = []
        for i in range(n_sites):
            H.append(self.tensor_O(beta))
        return np.array(H)
