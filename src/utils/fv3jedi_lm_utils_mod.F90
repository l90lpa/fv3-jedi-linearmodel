module fv3jedi_lm_utils_mod

use fv3jedi_lm_kinds_mod

implicit none
private

public :: fv3jedi_lm_conf, fv3jedi_lm_pert, fv3jedi_lm_traj

!> Fortran derived type to hold the linearized model configuration
type :: fv3jedi_lm_conf
  real(kind_real) :: dt                  !<Model time step
  real(kind_real) :: ptop                !<Pressure of top level
  integer         :: isc,iec,jsc,jec     !<Grid, compute region
  integer         :: isd,ied,jsd,jed     !<Grid, with halo
  integer         :: npz                 !<Number of vertical levels
  integer         :: do_dyn
  integer         :: do_phy_trb
  integer         :: do_phy_mst
  real(kind_real), allocatable, dimension(:) :: ak, bk
end type fv3jedi_lm_conf

!> Fortran derived type to hold the linearized model increment
type :: fv3jedi_lm_pert
  real(kind_real), allocatable, dimension(:,:,:) :: u, v, t, delp
  real(kind_real), allocatable, dimension(:,:,:) :: ua, va
  real(kind_real), allocatable, dimension(:,:,:) :: qv, ql, qi, o3, cfcn
  real(kind_real), allocatable, dimension(:,:,:) :: w, delz
end type fv3jedi_lm_pert

!> Fortran derived type to hold the linearized model trajectory
type :: fv3jedi_lm_traj
  real(kind_real), allocatable, dimension(:,:,:) :: u, v, ua, va, t, delp
  real(kind_real), allocatable, dimension(:,:,:) :: w, delz
  real(kind_real), allocatable, dimension(:,:,:) :: qv, ql, qi, o3
  real(kind_real), allocatable, dimension(:,:,:) :: qls, qcn, cfcn
  real(kind_real), allocatable, dimension(:,:)   :: phis, ps
  real(kind_real), allocatable, dimension(:,:)   :: frocean, frland
  real(kind_real), allocatable, dimension(:,:)   :: varflt, ustar, bstar
  real(kind_real), allocatable, dimension(:,:)   :: zpbl, cm, ct, cq
  real(kind_real), allocatable, dimension(:,:)   :: kcbl, ts, khl, khu
end type fv3jedi_lm_traj

end module fv3jedi_lm_utils_mod
