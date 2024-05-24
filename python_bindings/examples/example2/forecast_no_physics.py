from .pytorch_fv3lm import FV3LMModule
from pyfv3jedilm import runtime_mod

import isodate
from mpi4py import MPI
import netCDF4
import numpy as np
import torch


# MPI Setup
ROOT = 0
COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()
assert SIZE == 6

# Initialize the runtime (FMS, MPP, Field Manager, etc)
fms_nml_file = '../data/input_gfs_c12.nml'
field_table_file = '../data/field_table_gmao'
runtime_mod.runtime_init(COMM.py2f(), fms_nml_file, field_table_file)

# Model parameters
dt = isodate.parse_duration('PT15M').total_seconds()
npx = 13
npy = 13
npz = 127
with netCDF4.Dataset('../data/akbk127.nc4') as ds:
    ak = ds['ak'][:]
    bk = ds['bk'][:]
ptop = ak[0]
ntracers = 4
do_phy_trb = False
do_phy_mst = False
inputpert_filename = '../data/inputpert_4dvar.nml'

mod = FV3LMModule(dt, npx, npy, npz, ptop, ak, bk, ntracers, do_phy_mst, do_phy_trb, inputpert_filename)

isc = mod.lm.conf.isc
iec = mod.lm.conf.iec
jsc = mod.lm.conf.jsc
jec = mod.lm.conf.jec

nx = iec - isc + 1
ny = jec - jsc + 1

u = torch.zeros((nx, ny, npz), requires_grad=True)
v = torch.zeros((nx, ny, npz), requires_grad=True)
t = torch.zeros((nx, ny, npz), requires_grad=True)
delp = torch.zeros((nx, ny, npz), requires_grad=True)
tracers = torch.zeros((nx, ny, npz, ntracers), requires_grad=True)

if not mod.lm.conf.hydrostatic:
    w = torch.zeros((nx, ny, npz), requires_grad=True)
    delz = torch.zeros((nx, ny, npz), requires_grad=True)
else:
    w = None
    delz = None

# if do_phy_mst:
#     cfcn = torch.zeros((nx, ny, npz), requires_grad=False)
#     qls = torch.zeros((nx, ny, npz), requires_grad=False)
#     qcn = torch.zeros((nx, ny, npz), requires_grad=False)
# else:
#     cfcn = None
#     qls = None
#     qcn = None
# phis = torch.zeros((nx, ny), requires_grad=False)
# ps = torch.zeros((nx, ny), requires_grad=False)
# frocean = torch.zeros((nx, ny), requires_grad=False)
# frland = torch.zeros((nx, ny), requires_grad=False)
# varflt = torch.zeros((nx, ny), requires_grad=False)
# ustar = torch.zeros((nx, ny), requires_grad=False)
# bstar = torch.zeros((nx, ny), requires_grad=False)
# zpbl = torch.zeros((nx, ny), requires_grad=False)
# cm = torch.zeros((nx, ny), requires_grad=False)
# ct = torch.zeros((nx, ny), requires_grad=False)
# cq = torch.zeros((nx, ny), requires_grad=False)
# kcbl = torch.zeros((nx, ny), requires_grad=False)
# ts = torch.zeros((nx, ny), requires_grad=False)
# khl = torch.zeros((nx, ny), requires_grad=False)
# khu = torch.zeros((nx, ny), requires_grad=False)

z = mod.forward(u, v, t, delp, tracers, w=w, delz=delz)

# Shutdown the runtime
runtime_mod.runtime_exit()