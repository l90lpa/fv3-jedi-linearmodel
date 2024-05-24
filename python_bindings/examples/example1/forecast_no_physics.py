from pyfv3jedilm import fv3jedi_lm_mod, runtime_mod

import isodate
from mpi4py import MPI
import netCDF4
import numpy as np

# MPI Setup
ROOT = 0
COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()
assert SIZE == 6

# FMS Parameters
fms_nml_file = '../data/input_gfs_c12.nml'
field_table_file = '../data/field_table_gmao'

# Initialize the runtime (FMS, MPP, Field Manager, etc)
runtime_mod.runtime_init(COMM.py2f(), fms_nml_file, field_table_file)

# Create the model
lm = fv3jedi_lm_mod.fv3jedi_lm_type()

# Disable physics (it is enabled by default)
lm.conf.do_phy = 0
lm.conf.do_phy_trb = 0
lm.conf.do_phy_mst = 0

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
include_pert_tracers = False

# Initialize the model
lm.create(dt, npx, npy, npz, ptop, ak, bk)
if ntracers > 0:
    lm.allocate_tracers(lm.conf.isc, lm.conf.iec, lm.conf.jsc, lm.conf.jec, npz, ntracers, include_pert_tracers)
lm.init_nl()

# Set dynamical core variables
fp = ['../data/20201215.000000.fv_core.res.tile1.nc',
      '../data/20201215.000000.fv_core.res.tile2.nc',
      '../data/20201215.000000.fv_core.res.tile3.nc',
      '../data/20201215.000000.fv_core.res.tile4.nc',
      '../data/20201215.000000.fv_core.res.tile5.nc',
      '../data/20201215.000000.fv_core.res.tile6.nc']
with netCDF4.Dataset(fp[RANK]) as ds:
    lm.traj.u[...] = np.transpose(np.squeeze(ds['ua']))
    lm.traj.v[...] = np.transpose(np.squeeze(ds['va']))
    lm.traj.t[...] = np.transpose(np.squeeze(ds['T']))
    lm.traj.delp[...] = np.transpose(np.squeeze(ds['delp']))
    lm.traj.phis[...] = np.transpose(np.squeeze(ds['phis']))

# Set tracers
fp = ['../data/20201215.000000.fv_tracer.res.tile1.nc',
      '../data/20201215.000000.fv_tracer.res.tile2.nc',
      '../data/20201215.000000.fv_tracer.res.tile3.nc',
      '../data/20201215.000000.fv_tracer.res.tile4.nc',
      '../data/20201215.000000.fv_tracer.res.tile5.nc',
      '../data/20201215.000000.fv_tracer.res.tile6.nc']

with netCDF4.Dataset(fp[RANK]) as ds:
    encoding = 'ascii'

    # Initialize tracer names with 'space' character 
    lm.traj.tracer_names[...] = 32

    tracer_names = ['sphum', 'liq_wat', 'ice_wat', 'o3mr']
    for index, name in enumerate(tracer_names):

        # Set tracer data
        lm.traj.tracers[:,:,:,index] = np.transpose(np.squeeze(ds[name]))

        # Set tracer name
        encoded_name = np.frombuffer("sphum".encode(encoding), dtype=np.uint8)
        lm.traj.tracer_names[:np.size(encoded_name),index] = encoded_name

lm.step_nl()

lm.final_nl()

lm.delete()

# Shutdown the runtime
runtime_mod.runtime_exit()