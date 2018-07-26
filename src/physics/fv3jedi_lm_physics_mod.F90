module fv3jedi_lm_physics_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

use fv3jedi_lm_moist_mod

!> Physics driver for fv3jedi linearized model
!> Just calls its children in turn

implicit none
private
public :: fv3jedi_lm_dynamics_type

type fv3jedi_lm_dynamics_type
  type(fv3jedi_lm_moist_type) :: moist
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type fv3jedi_lm_dynamics_type

contains

! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_physics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 call self%fv3jedi_lm_moist_type%create()

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 call self%fv3jedi_lm_moist_type%init_nl()

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 call self%fv3jedi_lm_moist_type%init_tl()

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

 call self%fv3jedi_lm_moist_type%init_ad()

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout), target :: self
 type(fv3jedi_lm_traj), intent(inout) :: traj
 type(fv3jedi_lm_conf), intent(in) :: conf

 call self%fv3jedi_lm_moist_type%step_nl()

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 call self%fv3jedi_lm_moist_type%step_tl()

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 call self%fv3jedi_lm_moist_type%step_ad()

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self)

 implicit none
 class(fv3jedi_lm_dynamics_type), intent(inout) :: self

 call self%fv3jedi_lm_moist_type%delete()

endsubroutine delete

! ------------------------------------------------------------------------------

end module fv3jedi_lm_physics_mod
