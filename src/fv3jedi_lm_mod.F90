module fv3jedi_lm_mod

use fv3jedi_lm_kinds_mod
use fv3jedi_lm_utils_mod
use fv3jedi_lm_dynamics_mod, only: fv3jedi_lm_dynamics_type

!> Top level for fv3jedi linearized model

implicit none
private
public :: fv3jedi_lm_type

type fv3jedi_lm_type
 type(fv3jedi_lm_conf) :: conf
 type(fv3jedi_lm_pert) :: pert
 type(fv3jedi_lm_traj) :: traj
 type(fv3jedi_lm_dynamics_type) :: fv3jedi_lm_dynamics
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type

contains

! ------------------------------------------------------------------------------

subroutine create(self,dt,npx,npy,npz,ptop,ak,bk)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 real(kind=kind_real), intent(in) :: dt, ptop
 integer, intent(in) :: npx,npy,npz
 real(kind=kind_real), intent(in) :: ak(npz+1), bk(npz+1)

 self%conf%dt = dt
 self%conf%ptop = ptop

 self%conf%do_phy_mst = 1

 allocate(self%conf%ak(npz+1))
 allocate(self%conf%bk(npz+1))
 self%conf%ak = ak
 self%conf%bk = bk

 call self%fv3jedi_lm_dynamics%create(self%conf)

 !Make grid available to all components
 self%conf%isc = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%isc
 self%conf%iec = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%iec
 self%conf%jsc = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jsc
 self%conf%jec = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jec
 self%conf%isd = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%isd
 self%conf%ied = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%ied
 self%conf%jsd = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jsd
 self%conf%jed = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jed
 self%conf%npz = self%fv3jedi_lm_dynamics%FV_Atm(1)%npz

 !Convenience
 self%conf%hydrostatic = self%fv3jedi_lm_dynamics%FV_Atm(1)%flagstruct%hydrostatic

 !Allocate main traj and pert structures
 call allocate_traj(self%traj,self%conf%isc,self%conf%iec,self%conf%jsc,self%conf%jec,&
                    self%conf%npz,self%conf%hydrostatic,self%conf%do_phy_mst)
 call allocate_pert(self%pert,self%conf%isc,self%conf%iec,self%conf%jsc,self%conf%jec,self%conf%npz,self%conf%hydrostatic)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%init_nl(self%pert,self%traj)

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%init_tl(self%pert,self%traj)

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%init_ad(self%pert,self%traj)

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%step_nl(self%conf,self%traj)

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%step_tl(self%conf,self%traj,self%pert)

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%step_ad(self%conf,self%traj,self%pert)

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self)

 implicit none

 class(fv3jedi_lm_type), intent(inout) :: self

 call self%fv3jedi_lm_dynamics%delete()

 deallocate(self%conf%ak)
 deallocate(self%conf%bk)

 call deallocate_traj(self%traj)
 call deallocate_pert(self%pert)

endsubroutine delete

! ------------------------------------------------------------------------------

end module fv3jedi_lm_mod
