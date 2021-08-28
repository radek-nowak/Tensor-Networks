import numpy as np
from ncon import ncon
import matplotlib.pyplot as plt

from Ising_init import Hamiltonian as H


tc = lambda x: H(1,1).tensor_O(x,'corner')
te = lambda x: H(1,1).tensor_O(x,'edge')
to = lambda x: H(1,1).tensor_O(x,'o')

tcz = lambda x: H(1,1).tensor_Z(x,'corner')
tez = lambda x: H(1,1).tensor_Z(x,'edge')
toz = lambda x: H(1,1).tensor_Z(x,'o')

def sum_edge(beta,n,q):
    chi = 2
    if q == 'O':
        e = te(beta)
        o = to(beta)
        indices = [[-1,1,-2],[-3,-5,-4,1]]
        if n == 1:
            m1 = ncon([e,o],indices).reshape(chi**2,chi,chi**2)
            frob_norm_m1 = np.sqrt(ncon(m1.reshape(chi**2,chi**3)@m1.reshape(chi**2,chi**3).transpose(1,0),[1,1]))
            return m1/frob_norm_m1
        else:
            mn = ncon([sum_edge(beta,n-1,'O'),o],indices).reshape(chi**(n+1),chi,chi**(n+1))
            frob_norm_mn = np.sqrt(ncon(
                            mn.reshape(chi**(n+1),chi**(n+2))@mn.reshape(chi**(n+1),chi**(n+2)).transpose(1,0),[1,1]))
            return mn/frob_norm_mn
    else:
        e = tez(beta)
        o = toz(beta)
        indices = [[-1,1,-2],[-3,-5,-4,1]]
        if n == 1:
            m1 = ncon([e,o],indices).reshape(chi**2,chi,chi**2)
            return m1
        else:
            mn = ncon([sum_edge(beta,n-1,'Z'),o],indices).reshape(chi**(n+1),chi,chi**(n+1))
            frob_norm_mn = np.sqrt(ncon(
                            mn.reshape(chi**(n+1),chi**(n+2))@mn.reshape(chi**(n+1),chi**(n+2)).transpose(1,0),[1,1]))
            return mn/frob_norm_mn
    
def sum_corner(beta,n,q):
    chi = 2
    if q == 'O':
        c = tc(beta)
        e1,e2 = te(beta),te(beta)
        o = to(beta)
        indices = [[1,2],[-1,3,1],[2,4,-4],[3,-2,-3,4]]
        if n == 1:
            m1 = ncon([c,e1,e2,o],indices).reshape(chi**2,chi**2)
            frob_norm_m1 = np.sqrt(ncon([m1 @ m1.transpose(1,0)],[1,1]))
            return m1/frob_norm_m1      
        else:
            mn = ncon(
                [sum_corner(beta,n-1,'O'),sum_edge(beta,n-1,'O'),sum_edge(beta,n-1,'O'),o],
                indices).reshape(chi**(n+1),chi**(n+1))
            frob_norm_mn = np.sqrt(ncon([mn @ mn.transpose(1,0)],[1,1]))
            return mn/frob_norm_mn
    else:
        c = tcz(beta)
        e1,e2 = tez(beta),tez(beta)
        z = toz(beta)
        indices = [[1,2],[-1,3,1],[2,4,-4],[3,-2,-3,4]]
        if n == 1:
            m1 = ncon([c,e1,e2,z],indices).reshape(chi**2,chi**2)
            frob_norm_m1 = np.sqrt(ncon([m1 @ m1.transpose(1,0)],[1,1]))
            return m1/frob_norm_m1      
        else:
            mn = ncon(
                    [sum_corner(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),z],
                        indices).reshape(chi**(n+1),chi**(n+1))
            frob_norm_mn = np.sqrt(ncon([mn @ mn.transpose(1,0)],[1,1]))
            return mn/frob_norm_mn