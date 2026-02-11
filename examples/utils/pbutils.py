import numpy as np
import re

def make_nei_ele(m, n, l, dtype=np.int32):
    """
    Build adjacency for a regular m*n*l grid (x,y,z) with linear index:
        idx = x + m*y + (m*n)*z,   x∈[0,m), y∈[0,n), z∈[0,l)
    Output:
        adj_arr shape (6, m*n*l), order = [-x, +x, -y, +y, -z, +z]
        Boundary faces point to self.
    """
    N = m*n*l
    idx3d = np.arange(N, dtype=dtype).reshape(l, n, m)  # axes: [z, y, x]
    adj = np.empty((6, l, n, m), dtype=dtype)

    # Start by pointing every face to self; then overwrite interior neighbors.
    adj[0] = idx3d  # -x
    adj[1] = idx3d  # +x
    adj[2] = idx3d  # -y
    adj[3] = idx3d  # +y
    adj[4] = idx3d  # -z
    adj[5] = idx3d  # +z

    # x-neighbors (axis=-1)
    adj[0][..., 1:] = idx3d[..., :-1]   # -x for x>=1
    adj[1][..., :-1] = idx3d[..., 1:]   # +x for x<=m-2

    # y-neighbors (axis=1)
    adj[2][:, 1:, :] = idx3d[:, :-1, :] # -y for y>=1
    adj[3][:, :-1, :] = idx3d[:, 1:, :] # +y for y<=n-2

    # z-neighbors (axis=0)
    adj[4][1:, :, :] = idx3d[:-1, :, :] # -z for z>=1
    adj[5][:-1, :, :] = idx3d[1:, :, :] # +z for z<=l-2

    return adj.reshape(6, N)


def make_coloring(m, n, l):
    """
    Linearization (x-fastest): idx = x + m*y + (m*n)*z
    Returns:
      A, B  -> 1D arrays of intrinsic indices where (x+y+z) is even/odd.
    Order matches the intrinsic numbering.
    """
    N = m * n * l
    i = np.arange(N, dtype=np.int32)

    # Invert x-fastest mapping
    x =  i              % m
    y = (i // m)        % n
    z =  i // (m * n)

    even = ((x + y + z) & 1) == 0

    return np.concatenate((i[even], i[~even]))


def read_input(fname):
    isnum = re.compile(r'\d+')

    try:
        src = open(fname).read()
    except:
        print("[Error] Cannot open the input file")
        return None
    
    nums = re.findall(isnum, src)
    if len(nums) != 6:
        print("[Error] Missing input parameter(s)")
        return None
    
    return np.array([eval(num) for num in nums])
