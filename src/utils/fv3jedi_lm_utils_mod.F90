module fv3jedi_lm_utils_mod

use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod, only: kappa

implicit none
private

public :: fv3jedi_lm_conf, fv3jedi_lm_pert, fv3jedi_lm_traj
public :: allocate_pert, deallocate_pert
public :: allocate_traj, deallocate_traj, copy_traj

!> Fortran derived type to hold the linearized model configuration
type :: fv3jedi_lm_conf
  real(kind_real) :: dt                                !<Model time step
  logical         :: saveltraj = .false.               !<Option to save local trajectories (physics)
  integer         :: n = 1                             !<Current time step in window
  integer         :: nt = 1                            !<Number of timesteps in window
  real(kind_real) :: ptop                              !<Pressure of top level
  integer         :: isc,iec,jsc,jec                   !<Cube grid, compute region
  integer         :: isd,ied,jsd,jed                   !<Cube grid, with halo
  integer         :: npx, npy, npz                     !<Number of grid points, dynamics
  integer         :: im, jm, lm                        !<Number of grid points, physics, 1:im etc
  integer         :: do_dyn = 1                        !<Dynamics switch
  integer         :: do_phy_trb = 1                    !<Physics switch for BL turb
  integer         :: do_phy_mst = 1                    !<Physics switch for convection and cloud
  real(kind_real), allocatable, dimension(:) :: ak, bk !<Vertical grid
  logical         :: hydrostatic                       !<Hydrostatic dy core
  logical         :: rpe                               !<True if root process
end type fv3jedi_lm_conf

!> Fortran derived type to hold the linearized model increment
type :: fv3jedi_lm_pert
  real(kind_real), allocatable, dimension(:,:,:) :: u, v, t, delp  !Dynamics
  real(kind_real), allocatable, dimension(:,:,:) :: qv, ql, qi, o3 !Tracers
  real(kind_real), allocatable, dimension(:,:,:) :: w, delz        !nh vars
  real(kind_real), allocatable, dimension(:,:,:) :: ua, va, cfcn   !Internal not part of increment
end type fv3jedi_lm_pert

!> Fortran derived type to hold the linearized model trajectory
type :: fv3jedi_lm_traj
  real(kind_real), allocatable, dimension(:,:,:) :: u, v, t, delp  !Dynamics
  real(kind_real), allocatable, dimension(:,:,:) :: qv, ql, qi, o3 !Tracers
  real(kind_real), allocatable, dimension(:,:,:) :: w, delz        !nh vars
  real(kind_real), allocatable, dimension(:,:,:) :: ua, va, cfcn   !Internal not part of increment
  real(kind_real), allocatable, dimension(:,:,:) :: qls, qcn
  real(kind_real), allocatable, dimension(:,:)   :: phis, ps
  real(kind_real), allocatable, dimension(:,:)   :: frocean, frland
  real(kind_real), allocatable, dimension(:,:)   :: varflt, ustar, bstar
  real(kind_real), allocatable, dimension(:,:)   :: zpbl, cm, ct, cq
  real(kind_real), allocatable, dimension(:,:)   :: kcbl, ts, khl, khu
end type fv3jedi_lm_traj

!> Compute ice fraction from temperature
public IceFraction
interface IceFraction
 module procedure IceFraction_r4
 module procedure IceFraction_r8
end interface

!> Compute pressures from delp
public compute_pressures
interface compute_pressures
 module procedure compute_pressures_r4
 module procedure compute_pressures_r8
end interface


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine allocate_pert(pert,isc,iec,jsc,jec,npz,hydrostatic)

 implicit none
 type(fv3jedi_lm_pert), intent(inout) :: pert
 logical,               intent(in   ) :: hydrostatic
 integer,               intent(in   ) :: isc, iec, jsc, jec, npz

 allocate(pert%u      (isc:iec, jsc:jec, npz))
 allocate(pert%v      (isc:iec, jsc:jec, npz))
 allocate(pert%ua     (isc:iec, jsc:jec, npz))
 allocate(pert%va     (isc:iec, jsc:jec, npz))
 allocate(pert%t      (isc:iec, jsc:jec, npz))
 allocate(pert%delp   (isc:iec, jsc:jec, npz))
 allocate(pert%qv     (isc:iec, jsc:jec, npz))
 allocate(pert%ql     (isc:iec, jsc:jec, npz))
 allocate(pert%qi     (isc:iec, jsc:jec, npz))
 allocate(pert%o3     (isc:iec, jsc:jec, npz))

 allocate(pert%cfcn     (isc:iec, jsc:jec, npz))
 
 if (.not. hydrostatic) then
   allocate(pert%w      (isc:iec, jsc:jec, npz))
   allocate(pert%delz   (isc:iec, jsc:jec, npz))
 endif

end subroutine allocate_pert

! ------------------------------------------------------------------------------

subroutine deallocate_pert(pert)

 implicit none
 type(fv3jedi_lm_pert), intent(inout) :: pert

 if(allocated(pert%u)   ) deallocate(pert%u)
 if(allocated(pert%v)   ) deallocate(pert%v)
 if(allocated(pert%ua)  ) deallocate(pert%ua)
 if(allocated(pert%va)  ) deallocate(pert%va)
 if(allocated(pert%t)   ) deallocate(pert%t)
 if(allocated(pert%delp)) deallocate(pert%delp)
 if(allocated(pert%qv)  ) deallocate(pert%qv)
 if(allocated(pert%ql)  ) deallocate(pert%ql)
 if(allocated(pert%qi)  ) deallocate(pert%qi)
 if(allocated(pert%o3)  ) deallocate(pert%o3)
 if(allocated(pert%w)   ) deallocate(pert%w)
 if(allocated(pert%delz)) deallocate(pert%delz)

end subroutine deallocate_pert

! ------------------------------------------------------------------------------

subroutine allocate_traj(traj,isc,iec,jsc,jec,npz,hydrostatic,dpm)

 implicit none
 type(fv3jedi_lm_traj), intent(inout) :: traj
 logical,               intent(in   ) :: hydrostatic
 integer,               intent(in   ) :: dpm
 integer,               intent(in   ) :: isc, iec, jsc, jec, npz

 allocate(traj%u      (isc:iec, jsc:jec, npz))
 allocate(traj%v      (isc:iec, jsc:jec, npz))
 allocate(traj%ua     (isc:iec, jsc:jec, npz))
 allocate(traj%va     (isc:iec, jsc:jec, npz))
 allocate(traj%t      (isc:iec, jsc:jec, npz))
 allocate(traj%delp   (isc:iec, jsc:jec, npz))
 allocate(traj%qv     (isc:iec, jsc:jec, npz))
 allocate(traj%ql     (isc:iec, jsc:jec, npz))
 allocate(traj%qi     (isc:iec, jsc:jec, npz))
 allocate(traj%o3     (isc:iec, jsc:jec, npz))
 
 if (.not. hydrostatic) then
   allocate(traj%w      (isc:iec, jsc:jec, npz))
   allocate(traj%delz   (isc:iec, jsc:jec, npz))
 endif
 
 if (dpm .ne. 0) then
   allocate(traj%qls    (isc:iec, jsc:jec, npz))
   allocate(traj%qcn    (isc:iec, jsc:jec, npz))
   allocate(traj%cfcn   (isc:iec, jsc:jec, npz))
 endif
 
 allocate(traj%phis   (isc:iec, jsc:jec))
 allocate(traj%ps     (isc:iec, jsc:jec))
 allocate(traj%frocean(isc:iec, jsc:jec))
 allocate(traj%frland (isc:iec, jsc:jec))
 allocate(traj%varflt (isc:iec, jsc:jec))
 allocate(traj%ustar  (isc:iec, jsc:jec))
 allocate(traj%bstar  (isc:iec, jsc:jec))
 allocate(traj%zpbl   (isc:iec, jsc:jec))
 allocate(traj%cm     (isc:iec, jsc:jec))
 allocate(traj%ct     (isc:iec, jsc:jec))
 allocate(traj%cq     (isc:iec, jsc:jec))
 allocate(traj%kcbl   (isc:iec, jsc:jec))
 allocate(traj%ts     (isc:iec, jsc:jec))
 allocate(traj%khl    (isc:iec, jsc:jec))
 allocate(traj%khu    (isc:iec, jsc:jec))

end subroutine allocate_traj

! ------------------------------------------------------------------------------

subroutine deallocate_traj(traj)

 implicit none
 type(fv3jedi_lm_traj), intent(inout) :: traj
 
 if (allocated(traj%u      )) deallocate(traj%u      )
 if (allocated(traj%v      )) deallocate(traj%v      )
 if (allocated(traj%ua     )) deallocate(traj%ua     )
 if (allocated(traj%va     )) deallocate(traj%va     )
 if (allocated(traj%t      )) deallocate(traj%t      )
 if (allocated(traj%delp   )) deallocate(traj%delp   )
 if (allocated(traj%qv     )) deallocate(traj%qv     )
 if (allocated(traj%ql     )) deallocate(traj%ql     )
 if (allocated(traj%qi     )) deallocate(traj%qi     )
 if (allocated(traj%o3     )) deallocate(traj%o3     )
 if (allocated(traj%w      )) deallocate(traj%w      )
 if (allocated(traj%delz   )) deallocate(traj%delz   )
 if (allocated(traj%qls    )) deallocate(traj%qls    )
 if (allocated(traj%qcn    )) deallocate(traj%qcn    )
 if (allocated(traj%cfcn   )) deallocate(traj%cfcn   )
 if (allocated(traj%phis   )) deallocate(traj%phis   )
 if (allocated(traj%ps     )) deallocate(traj%ps     )
 if (allocated(traj%frocean)) deallocate(traj%frocean)
 if (allocated(traj%frland )) deallocate(traj%frland )
 if (allocated(traj%varflt )) deallocate(traj%varflt )
 if (allocated(traj%ustar  )) deallocate(traj%ustar  )
 if (allocated(traj%bstar  )) deallocate(traj%bstar  )
 if (allocated(traj%zpbl   )) deallocate(traj%zpbl   )
 if (allocated(traj%cm     )) deallocate(traj%cm     )
 if (allocated(traj%ct     )) deallocate(traj%ct     )
 if (allocated(traj%cq     )) deallocate(traj%cq     )
 if (allocated(traj%kcbl   )) deallocate(traj%kcbl   )
 if (allocated(traj%ts     )) deallocate(traj%ts     )
 if (allocated(traj%khl    )) deallocate(traj%khl    )
 if (allocated(traj%khu    )) deallocate(traj%khu    )

end subroutine deallocate_traj

! ------------------------------------------------------------------------------

subroutine copy_traj( traj_in, traj_out, hydrostatic, dpm )

 implicit none
 type(fv3jedi_lm_traj), intent(in)    :: traj_in
 type(fv3jedi_lm_traj), intent(inout) :: traj_out
 logical,               intent(in)    :: hydrostatic
 integer,               intent(in)    :: dpm
 
 traj_out%u    = traj_in%u
 traj_out%v    = traj_in%v
 traj_out%ua   = traj_in%ua
 traj_out%va   = traj_in%va
 traj_out%t    = traj_in%t
 traj_out%delp = traj_in%delp
 traj_out%qv   = traj_in%qv
 traj_out%ql   = traj_in%ql
 traj_out%qi   = traj_in%qi
 traj_out%o3   = traj_in%o3
 
 if (.not. hydrostatic) then
 traj_out%w    = traj_in%w
 traj_out%delz = traj_in%delz
 endif
 
 if (dpm /= 0) then
 traj_out%qls  = traj_in%qls
 traj_out%qcn  = traj_in%qcn
 traj_out%cfcn = traj_in%cfcn
 endif
 
 traj_out%phis    = traj_in%phis
 traj_out%ps      = traj_in%ps
 traj_out%frocean = traj_in%frocean
 traj_out%frland  = traj_in%frland
 traj_out%varflt  = traj_in%varflt
 traj_out%ustar   = traj_in%ustar
 traj_out%bstar   = traj_in%bstar
 traj_out%zpbl    = traj_in%zpbl
 traj_out%cm      = traj_in%cm
 traj_out%ct      = traj_in%ct
 traj_out%cq      = traj_in%cq
 traj_out%kcbl    = traj_in%kcbl
 traj_out%ts      = traj_in%ts
 traj_out%khl     = traj_in%khl
 traj_out%khu     = traj_in%khu

end subroutine copy_traj

! ------------------------------------------------------------------------------

subroutine icefraction_r4(temp, icefrct)

 implicit none

 !arguments
 real(4), intent(in) :: temp
 real(4), intent(out) :: icefrct

 !locals
 real(4), parameter :: t_ice_all = 233.16_4, t_ice_max = 273.16_4
 integer, parameter :: icefrpwr = 4

 icefrct  = 0.00_4
 if ( temp <= t_ice_all ) then
    icefrct = 1.000_4
 else if ( (temp > t_ice_all) .and. (temp <= t_ice_max) ) then
    icefrct = 1.00_4 -  ( temp - t_ice_all ) / ( t_ice_max - t_ice_all )
 end if

 icefrct = min(icefrct,1.00_4)
 icefrct = max(icefrct,0.00_4)

 icefrct = icefrct**icefrpwr

end subroutine icefraction_r4

! ------------------------------------------------------------------------------

subroutine icefraction_r8(temp, icefrct)

 implicit none

 !arguments
 real(8), intent(in) :: temp
 real(8), intent(out) :: icefrct

 !locals
 real(8), parameter :: t_ice_all = 233.16_8, t_ice_max = 273.16_8
 integer, parameter :: icefrpwr = 4

 icefrct  = 0.0_8
 if ( temp <= t_ice_all ) then
    icefrct = 1.000_8
 else if ( (temp > t_ice_all) .and. (temp <= t_ice_max) ) then
    icefrct = 1.00_8 -  ( temp - t_ice_all ) / ( t_ice_max - t_ice_all )
 end if

 icefrct = min(icefrct,1.00_8)
 icefrct = max(icefrct,0.00_8)

 icefrct = icefrct**icefrpwr

end subroutine icefraction_r8

! ------------------------------------------------------------------------------

subroutine compute_pressures_r4(im,jm,lm,ptop,delp,pe,p,pk)

 integer, intent(in)  :: im,jm,lm
 real(4), intent(in)  :: ptop
 real(4), intent(in)  :: delp(1:im,1:jm,1:lm) !Pressure thickness
 real(4), intent(out) :: pe  (1:im,1:jm,0:lm) !Pressure at edges
 real(4), intent(out) :: p   (1:im,1:jm,1:lm) !Pressure at mid-point
 real(4), intent(out) :: pk  (1:im,1:jm,1:lm) !Pressure to the kappa at mid-point
 
 real(4) :: lpe(1:im,1:jm,0:lm)
 real(4) :: pek(1:im,1:jm,0:lm)

 integer :: l

 !Pressure at the level edge
 pe(:,:,0) = ptop
 do l=1,lm
    pe(:,:,l) = pe(:,:,l-1) + delp(:,:,l)
 end do

 !Pressure mid-point
 p(:,:,1:lm) = 0.5_4 * (pe(:,:,1:lm) + pe(:,:,0:lm-1))

 !Log of edge pressure
 lpe = log(pe)

 !pe to the kappa at edge
 pek = pe**kappa

 !p to the kappa
 pk(:,:,1:lm) = (pek(:,:,1:lm)-pek(:,:,0:lm-1))/(kappa*(lpe(:,:,1:lm)-lpe(:,:,0:lm-1)))

endsubroutine compute_pressures_r4

! ------------------------------------------------------------------------------

subroutine compute_pressures_r8(im,jm,lm,ptop,delp,pe,p,pk)

 integer, intent(in)  :: im,jm,lm
 real(8), intent(in)  :: ptop
 real(8), intent(in)  :: delp(1:im,1:jm,1:lm) !Pressure thickness
 real(8), intent(out) :: pe  (1:im,1:jm,0:lm) !Pressure at edges
 real(8), intent(out) :: p   (1:im,1:jm,1:lm) !Pressure at mid-point
 real(8), intent(out) :: pk  (1:im,1:jm,1:lm) !Pressure to the kappa at mid-point
 
 real(8) :: lpe(1:im,1:jm,0:lm)
 real(8) :: pek(1:im,1:jm,0:lm)

 integer :: l

 !Pressure at the level edge
 pe(:,:,0) = ptop
 do l=1,lm
    pe(:,:,l) = pe(:,:,l-1) + delp(:,:,l)
 end do

 !Pressure mid-point
 p(:,:,1:lm) = 0.5_8 * (pe(:,:,1:lm) + pe(:,:,0:lm-1))

 !Log of edge pressure
 lpe = log(pe)

 !pe to the kappa at edge
 pek = pe**kappa

 !p to the kappa
 pk(:,:,1:lm) = (pek(:,:,1:lm)-pek(:,:,0:lm-1))/(kappa*(lpe(:,:,1:lm)-lpe(:,:,0:lm-1)))

endsubroutine compute_pressures_r8

end module fv3jedi_lm_utils_mod
