import numpy as np
from ncon import ncon
import matplotlib.pyplot as plt

from Ising_init import Hamiltonian as H


def sum_edge(beta,n,q):
    chi = 2
    if q == 'O':
        e = te(beta)
        o = to(beta)
        indices = [[-1,1,-2],[-3,-5,-4,1]]
        if n == 1:
            m1 = ncon([e,o],indices).reshape(chi**2,chi,chi**2)
            return m1
        else:
            mn = ncon([sum_edge(beta,n-1,'O'),o],indices).reshape(chi**(n+1),chi,chi**(n+1))
            return mn
    else:
        e = tez(beta)
        o = toz(beta)
        indices = [[-1,1,-2],[-3,-5,-4,1]]
        if n == 1:
            m1 = ncon([e,o],indices).reshape(chi**2,chi,chi**2)
            return m1
        else:
            mn = ncon([sum_edge(beta,n-1,'Z'),o],indices).reshape(chi**(n+1),chi,chi**(n+1))
            return mn
    
def sum_corner(beta,n,q):
    chi = 2
    if q == 'O':
        c = tc(beta)
        e1,e2 = te(beta),te(beta)
        o = to(beta)
        indices = [[1,2],[-1,3,1],[2,4,-4],[3,-2,-3,4]]
        if n == 1:
            m1 = ncon([c,e1,e2,o],indices).reshape(chi**2,chi**2)
            return m1      
        else:
            mn = ncon(
                [sum_corner(beta,n-1,'O'),sum_edge(beta,n-1,'O'),sum_edge(beta,n-1,'O'),o],
                indices).reshape(chi**(n+1),chi**(n+1))
            return mn
    else:
        c = tcz(beta)
        e1,e2 = tez(beta),tez(beta)
        z = toz(beta)
        indices = [[1,2],[-1,3,1],[2,4,-4],[3,-2,-3,4]]
        if n == 1:
            m1 = ncon([c,e1,e2,z],indices).reshape(chi**2,chi**2)
            return m1      
        else:
            mn = ncon(
                [sum_corner(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),z],
                indices).reshape(chi**(n+1),chi**(n+1))
            return mn