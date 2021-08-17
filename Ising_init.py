import numpy as np
import tensorflow as tf

class Hamiltonian:

    def __init__(self, L, J, g):
        self.L = L
        self.J = J
        self.sigma_z = np.array([[1.,0.],[0.,-1.]])
        self.one = np.eye(2)
        self.tensor_L        
        self.ising

    # tensor L should have 3 indices
    def tensor_L(self, beta):
        sz, one = self.sigma_z, self.one
        return [one*np.sqrt(np.cosh(beta)), sz*np.sqrt(np.sinh(beta))]
    
    # needs testing
    # Ls commute
    # this is for corner tensor
    def tensor_T(beta):
        L1, L2 = tensor_L(beta), tensor_L(beta)
        L3, L4 = tensor_L(beta)[0], tensor_L(beta)[0]
        L12 = tf.einsum('ijk, lkm -> ijlm', L1, L2) # i and l correspond do alpha
        L123 = tf.einsum('ijlm, mb -> ijlb', L12, L3) 
        T = tf.einsum('ijlb, bc -> ijlc', L123, L4)
        return T

    def tensor_O(beta):
        O = tf.einsum('ijlc')
        

    def ising(self,d):
        self.d = d
        sz, one = self.sigma_z, self.one
        sz = self.sigma_z
        nbonds = self.L
        
        for i in range(nbonds):
            bonds = -self.J*np.kron(sz,sz)
            h.append(np.reshape(bonds,[d,d,d,d]))
        return h
        #test change
    

