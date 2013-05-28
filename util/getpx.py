from ctypes import *
import os
os.environ['PATH'] = os.path.dirname(__file__) + ';' + os.environ['PATH']
api = CDLL('getpx.dll')

def get_px(area, max_j, all_c, lattice):
    lattice_arr = [p[i] for p in lattice for i in xrange(2)]
    lattice_arr = (c_float * len(lattice_arr))(*lattice_arr)
    c_arr = [c  for j in xrange(max_j+1) \
                for kx in xrange(2**j) for ky in xrange(2**j) \
                for c in all_c[(j, kx,ky)]]
    c_arr = (c_float * len(c_arr))(*c_arr)
    area = c_float(area)
    px_num, max_j = c_int(len(lattice)), c_int(max_j)
    px_arr = (c_float * len(lattice))(*([0] * len(lattice)))
    api.get_px(area, max_j, px_num, c_arr, lattice_arr, px_arr)
    return [px for px in px_arr]