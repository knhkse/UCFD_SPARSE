from ctypes import CDLL, c_double, c_void_p, c_int
from mpi4py import MPI
import numpy as np
import sys

from utils.pbutils import make_nei_ele, read_input

"""
doxygen?
"""
def run(params, comm, lib):
    # Constants
    nvars, nface = 7, 6
    gamma = 1.4
    nx, ny, nz, npx, npy, npz = params
    m, n, l = nx//npx, ny//npy, nz//npz
    nele = m*n*l
    rank = comm.Get_rank()

    # Ctypes functions
    update = lib.lusgs_serial_ns_update
    pre_lusgs = lib.serial_pre_lusgs
    lower_sweep = lib.ns_serial_lower_sweep
    upper_sweep = lib.ns_serial_upper_sweep

    update.argtypes = [c_int, c_void_p, c_void_p]
    pre_lusgs.argtypes = [c_int, c_int, c_double,
                          c_void_p, c_void_p, c_void_p, c_void_p]
    lower_sweep.argtypes = [c_int, c_int,
                            c_void_p, c_void_p, c_void_p, c_void_p,
                            c_void_p, c_void_p, c_void_p, c_void_p]
    upper_sweep.argtypes = [c_int, c_int,
                            c_void_p, c_void_p, c_void_p, c_void_p,
                            c_void_p, c_void_p, c_void_p, c_void_p]
    
    if rank == 0:
        print("[Run] Allocating arrays...", end='')
    
    # Flow variables & arrays
    rho, u, v, w, p = gamma, 1.5, 0, 0, 1.0
    et = p / (gamma - 1) + 0.5*rho*u**2
    upts = np.empty((nvars, nele), dtype=np.float64)
    rhs = np.empty_like(upts)
    dub = np.empty_like(upts)
    diag = np.ones((nele,), dtype=np.float64)
    fspr = np.ones((6, nele), dtype=np.float64)
    dt = np.ones((nele,), dtype=np.float64)*0.1
    fnorm_vol = np.ones((nface, nele), dtype=np.float64)
    vfi = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1], [1, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=np.float64)
    vec_fnorm = np.repeat(vfi[:, :, np.newaxis], nele, axis=2)
    comm.Barrier()
    if rank == 0:
        print(" Done")

    # Initial condition
    if rank == 0:
        print("[Run] Initializing arrays...", end='')
    upts[0] = rhs[0] = dub[0] = rho
    upts[1] = rhs[1] = dub[1] = rho*u
    upts[2] = rhs[2] = dub[2] = rho*v
    upts[3] = rhs[3] = dub[3] = rho*w
    upts[4] = rhs[4] = dub[4] = et
    nei_ele = make_nei_ele(m, n, l)
    comm.Barrier()
    if rank == 0:
        print(" Done")

    # Iteration
    if rank == 0:
        print("[Run] LU-SGS computation...", end='')
    pre_lusgs(nele, nface, 1.0, fnorm_vol.ctypes.data, dt.ctypes.data,
                diag.ctypes.data, fspr.ctypes.data)
    lower_sweep(nele, nface, nei_ele.ctypes.data,
                fnorm_vol.ctypes.data, vec_fnorm.ctypes.data,
                upts.ctypes.data, rhs.ctypes.data, dub.ctypes.data,
                diag.ctypes.data, fspr.ctypes.data)
    upper_sweep(nele, nface, nei_ele.ctypes.data,
                fnorm_vol.ctypes.data, vec_fnorm.ctypes.data,
                upts.ctypes.data, rhs.ctypes.data, dub.ctypes.data,
                diag.ctypes.data, fspr.ctypes.data)
    update(nele, upts.ctypes.data, rhs.ctypes.data)
    comm.Barrier()
    if rank == 0:
        print(" Done")

def write_output(nprocs, params):
    m, n, l, dm, dn, dl = params
    neles = m*n*l
    
    print("========= LU-SGS example info =========\n")
    print("Number of cores: {}\n".format(nprocs))
    print("Number of elements: {} = {} x {} x {}\n".format(neles, m, n, l))
    print("Partitioning with {}: {} x {} x {}\n".format(dm*dn*dl, dm, dn, dl))
    print("Elements per core: {} = {} x {} x {}\n\n".format(neles//nprocs, m//dm, n//dn, l//dl))


if __name__ == "__main__":
    # MPI initialize
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()

    if rank == 0:
        if (len(sys.argv)) < 3:
            print("[Error] Input file, dynamic library argument is needed")
            err = 1
        else:
            inputname = sys.argv[1]
            params = read_input(inputname)
            if params is None:
                err = 1
            else:
                dm = params[3]
                dn = params[4]
                dl = params[5]

                if nproc != dm*dn*dl:
                    print("[Error] Partitioning number mismatch")
                    print("{} processors, but total {} partitioning".format(nproc, dm*dn*dl))
                    err = 1
                else:
                    print("[Main] 3D hexahedral example starts...")
                    err = 0
    else:
        err = 0
        params = np.empty(6, dtype=np.int64)

    err = comm.bcast(err, root=0)
    if err == 1:
        sys.exit()
    
    params = comm.bcast(params, root=0)

    # Dynamic Linking
    try:
        lib = CDLL(sys.argv[2])
    except FileNotFoundError:
        print("[Error] Dynamic library file not found")
        sys.exit()
    
    run(params, comm, lib)
