import numpy as np
from ncon import ncon
from numpy.linalg import svd
from numpy import linalg as LA
#from functools import lru_cache

from Ising_init import Hamiltonian as H


tc = lambda x: H().tensor_O(x,'corner')
te = lambda x: H().tensor_O(x,'edge')
to = lambda x: H().tensor_O(x,'o')

tcz = lambda x: H().tensor_Z(x,'corner')
tez = lambda x: H().tensor_Z(x,'edge')
toz = lambda x: H().tensor_Z(x,'o')

#@njit(nopython=True)
#@lru_cache
def sum_edge(beta,n,q):
    chi = 2
    eo = te(beta)
    o = to(beta)
    ez = tez(beta)
    z = toz(beta)
    indices = [[-1,1,-2],[-3,-5,-4,1]]
    if q == 'O':
        m1 = ncon([eo,o],indices).reshape(chi**2,chi,chi**2)
        #frob_norm_m1 = np.sqrt(ncon([m1,m1],[[1,2,3],[3,2,1]]))
        norm = LA.norm(m1)#abs(np.max(m1))
        if n == 1:
            return m1/norm#/frob_norm_m1
        while n>1:
            mn = ncon([sum_edge(beta,n-1,'O'),o],indices).reshape(chi**(n+1),chi,chi**(n+1))
            #frob_norm_mn = np.sqrt(ncon([mn,mn],[[1,2,3],[3,2,1]]))
            norm = LA.norm(mn)#abs(np.max(mn))
            return mn/norm#/frob_norm_mn
    else:
        m1 = ncon([ez,z],indices).reshape(chi**2,chi,chi**2)
        #frob_norm_m1 = np.sqrt(ncon([m1,m1],[[1,2,3],[3,2,1]]))
        norm = LA.norm(m1)#abs(np.max(m1))
        if n == 1:
            return m1/norm#/frob_norm_m1
        while n>1:
            mn = ncon([sum_edge(beta,n-1,'Z'),z],indices).reshape(chi**(n+1),chi,chi**(n+1))
            #frob_norm_mn = np.sqrt(ncon([mn,mn],[[1,2,3],[3,2,1]]))
            norm = LA.norm(mn)#abs(np.max(mn))
            return mn/norm#/frob_norm_mn

#@njit(nopython=True) 
#@lru_cache
def sum_corner(beta,n,q):
    chi = 2
    co = tc(beta)
    eo = te(beta)
    o = to(beta)
    cz = tcz(beta)
    ez = tez(beta)
    z = toz(beta)
    indices = [[1,2],[-1,3,1],[2,4,-4],[3,-2,-3,4]]
    if q == 'O':
        m1 = ncon([co,eo,eo,o],indices).reshape(chi**2,chi**2)
        #frob_norm_m1 = np.sqrt(ncon([m1,m1],[[1,2],[2,1]]))
        norm = LA.norm(m1)#abs(np.max(m1))
        if n == 1:
            return m1/norm     
        while n>1:
            mn = ncon(
                [sum_corner(beta,n-1,'O'),sum_edge(beta,n-1,'O'),sum_edge(beta,n-1,'O'),o],
                indices).reshape(chi**(n+1),chi**(n+1))
            #frob_norm_mn = np.sqrt(ncon([mn,mn],[[1,2],[2,1]]))
            norm = LA.norm(mn)#abs(np.max(mn))
            return mn/norm
    else:
        m1 = ncon([cz,ez,ez,z],indices).reshape(chi**2,chi**2)
        #frob_norm_m1 = np.sqrt(ncon([m1,m1],[[1,2],[2,1]]))
        norm = LA.norm(m1)#abs(np.max(m1)) 
        if n == 1:
            return m1/norm#/frob_norm_m1      
        while n>1:
            mn = ncon(
                [sum_corner(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),sum_edge(beta,n-1,'Z'),z],
                indices).reshape(chi**(n+1),chi**(n+1))
            #frob_norm_m1frob_norm_mn = np.sqrt(ncon([mn,mn],[[1,2],[2,1]]))
            norm = LA.norm(mn)#abs(np.max(mn))
            return mn/norm#/frob_norm_mn
#@lru_cache
def ctmrg(beta,n, q):
    chi = 2
    if q == 'O':
        c = sum_corner(beta,n,'O')
        e = sum_edge(beta,n,'O')
        o = to(beta)
        
        tensor_list1 = [e,c,e,c,e]
        index_list1 = [[-1,-2,1],[1,2],[2,-5,3],[3,4],[-4,-3,4]]
        t1 = ncon(tensor_list1,index_list1).reshape(4*chi**(2*n),chi**3)
        
        #tensor_list2 = [e,c,e,c,e]
        #index_list2 = [[1,-2,-1],[2,1],[2,-3,3],[3,4],[4,-4,-5]]
        #t2 = ncon(tensor_list2, index_list2).reshape(chi**3,4*chi**(2*n))
        
        tensor_list3 = [o,o]
        index_list3 = [[-1,-2,-3,1],[-5,1,-6,-7]]
        t3 = ncon(tensor_list3,index_list3).reshape(chi**3,chi**3)
        
        w1 = ncon([t1,t3,t1.transpose(1,0)],[[1,2],[2,3],[3,1]])
        #trw1 = ncon(w1.reshape(4*chi**(2*n),4*chi**(2*n)),[1,1])
        return w1#/trw1
        #tensor_list = [t1,t3,t2]
        #index_list = [[1,2,3,4,5],[2,3,4,6,7,8],[1,6,7,8,5]]
        #return ncon(tensor_list,index_list)
    else:
        c = sum_corner(beta,n,'Z')
        e = sum_edge(beta,n,'Z')
        o = toz(beta)
        
        tensor_list1 = [e,c,e,c,e]
        index_list1 = [[-1,-2,1],[1,2],[2,-5,3],[3,4],[-4,-3,4]]
        t1 = ncon(tensor_list1,index_list1).reshape(4*chi**(2*n),chi**3)
        
        #tensor_list2 = [e,c,e,c,e]
        #index_list2 = [[1,-2,-1],[2,1],[2,-3,3],[3,4],[4,-4,-5]]
        #t2 = ncon(tensor_list2, index_list2).reshape(chi**3,4*chi**(2*n))
        
        tensor_list3 = [o,o]
        index_list3 = [[-1,-2,-3,1],[-5,1,-6,-7]]
        t3 = ncon(tensor_list3,index_list3).reshape(chi**3,chi**3)
        
        return ncon([t1,t3,t1.transpose(1,0)],[[1,2],[2,3],[3,1]])

def new_function(beta,n):
    P1, C_prime, P2 = svd(sum_corner(beta,n,'O'),full_matrices=False)
    P1_matrix = P1.reshape(2,2**n,2**(n+1))
    C_diag = np.diag(C_prime)
    P2_matrix = P2.reshape(2**(n+1),2**n,2)
    tensor_list = [P1_matrix, C_diag, P2_matrix]
    index_list = [[-1,-2,1],[1,2],[2,-3,-4]]
    return ncon(tensor_list,index_list).reshape(2**(n+1),2**(n+1))
    #return P1_matrix, C_diag, P2_matrix


def new_function2(beta,n):
    P1, C_prime, P2 = svd(sum_corner(beta,n,'O'),full_matrices=False)
    chi = int(len(C_prime))
    new_chi = int(len(C_prime)/2)
    P1_matrix = P1[0:new_chi,0:new_chi]
    C_diag = np.diag(C_prime[0:new_chi])
    P2_matrix = P2[0:new_chi,0:new_chi]
    tensor_list = [P1_matrix, C_diag, P2_matrix]
    index_list = [[-1,1],[1,2],[2,-4]]
    return ncon(tensor_list,index_list)
    #return P1_matrix, C_diag, P2_matrix
    #return new_chi
