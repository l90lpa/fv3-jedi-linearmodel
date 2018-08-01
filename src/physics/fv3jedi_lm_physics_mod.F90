module fv3jedi_lm_physics_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

use fv3jedi_lm_moist_mod, only: fv3jedi_lm_moist_type
use fv3jedi_lm_turbulence_mod, only: fv3jedi_lm_turbulence_type

!> Physics driver for fv3-jedi linearized model
!> Just calls its children in turn if turned on

implicit none
private
public :: fv3jedi_lm_physics_type

type fv3jedi_lm_physics_type
 type(fv3jedi_lm_moist_type) :: fv3jedi_lm_moist
 type(fv3jedi_lm_turbulence_type) :: fv3jedi_lm_turbulence
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type fv3jedi_lm_physics_type

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_physics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%create(conf)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%create(conf)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_physics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%init_nl(pert,traj)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%init_nl(pert,traj)

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_physics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%init_tl(pert,traj)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%init_tl(pert,traj)

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_physics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%init_ad(pert,traj)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%init_ad(pert,traj)

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 class(fv3jedi_lm_physics_type), intent(inout), target :: self
 type(fv3jedi_lm_traj), intent(inout) :: traj
 type(fv3jedi_lm_conf), intent(in) :: conf

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%step_nl(conf,traj)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%step_nl(conf,traj)

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_physics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%step_tl(conf,traj,pert)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%step_tl(conf,traj,pert)

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_physics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%step_ad(conf,traj,pert)
 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%step_ad(conf,traj,pert)

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self,conf)

 implicit none
 class(fv3jedi_lm_physics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 if (conf%do_phy_mst.ne.0) call self%fv3jedi_lm_moist%delete(conf)
 if (conf%do_phy_trb.ne.0) call self%fv3jedi_lm_turbulence%delete(conf)

endsubroutine delete

! ------------------------------------------------------------------------------

end module fv3jedi_lm_physics_mod
