import numba as nb
import numpy as np
import re

@nb.jit(nopython=True, fastmath=True)
def make_nei_ele(m, n, l, nei_ele):
    """
    Computes nei_ele array

    Only for hexahedral mesh

    parameter
    ---------
    m : x direction elements
    n : y direction elements
    l : z direction elements
    """
    
    idx = 0
    ml = m*l
    nml = n*m*l

    for j in range(n):
        for i in range(m):
            for k in range(l):
                # Next or previous neighbor cell count
                ny = idx-ml
                py = idx+ml
                nx = idx-l
                px = idx+l
                nz = idx-1
                pz = idx+1

                # Handling grid boundary
                # -y direction
                if ny < 0:
                    nei_ele[0, idx] = idx
                else:
                    nei_ele[0, idx] = ny

                # -x direction
                if nx < j*ml:
                    nei_ele[1, idx] = idx
                else:
                    nei_ele[1, idx] = nx

                # +z direction                
                if pz >= j*ml + (i+1)*l:
                    nei_ele[2, idx] = idx
                else:
                    nei_ele[2, idx] = pz
                
                # +x direction                
                if px >= (j+1)*ml:
                    nei_ele[3, idx] = idx
                else:
                    nei_ele[3, idx] = px

                # -z direction                
                if nz < j*ml + i*l:
                    nei_ele[4, idx] = idx
                else:
                    nei_ele[4, idx] = nz
                
                # +y direction
                if py >= nml:
                    nei_ele[5, idx] = idx
                else:
                    nei_ele[5, idx] = py
                
                idx += 1


@nb.jit(nopython=True, fastmath=True)
def make_coloring(m, n, l, icolor, lcolor):
    """
    2-Color Algorithm for hexahedral mesh
    """
    neles = m*n*l
    denom = neles//2
    cb = denom
    ml = m*l
    
    ca = 0
    ele = 0

    while ele < neles:
        for j in range(n):
            ist = (j*ml)//2
            iend = (j+1)*ml//2
            
            for ind in range(ist, iend):
                icolor[ca+ind] = ele
                lcolor[ele] = ca//denom
                ele += 1
                icolor[cb+ind] = ele
                lcolor[ele] = cb//denom
                ele += 1
            
            tmp = ca
            ca = cb
            cb = tmp


def read_input(fname):
    isnum = re.compile('\d+')

    try:
        src = open(fname).read()
    except:
        print("[Error] Cannot open the input file")
        return None
    
    nums = re.findall(isnum, src)
    if len(nums) != 7:
        print("[Error] Missing input parameter(s)")
        return None
    
    return np.array([eval(num) for num in nums])
