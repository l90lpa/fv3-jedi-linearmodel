!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
! $Id: fv_control_nlm.F90,v 1.3 2017/11/13 21:58:43 drholdaw Exp $
!
!----------------
! FV contro panel
!----------------

module fv_control_nlm_mod

 use fv_arrays_mod,     only: fv_atmos_type 
 use fv_arrays_nlm_mod, only: fv_atmos_pert_type, allocate_fv_atmos_pert_type, deallocate_fv_atmos_pert_type 
 use fms_mod,           only: open_namelist_file, check_nml_error, close_file
 use mpp_mod,           only: stdlog

 implicit none
 private

 logical, pointer :: split_hord 
 integer, pointer :: hord_mt_pert 
 integer, pointer :: hord_vt_pert 
 integer, pointer :: hord_tm_pert 
 integer, pointer :: hord_dp_pert 
 integer, pointer :: hord_tr_pert 
 logical, pointer :: split_kord 
 integer, pointer :: kord_mt_pert 
 integer, pointer :: kord_wz_pert 
 integer, pointer :: kord_tm_pert 
 integer, pointer :: kord_tr_pert 
 logical, pointer :: split_damp 
 logical, pointer :: do_vort_damp_pert 
 integer, pointer :: nord_pert 
 real,    pointer :: dddmp_pert 
 real,    pointer :: d2_bg_pert 
 real,    pointer :: d4_bg_pert 
 real,    pointer :: vtdm4_pert 
 logical, pointer :: split_damp_tr 
 integer, pointer :: nord_tr_pert 
 real,    pointer :: trdm2_pert

 public fv_init_pert, fv_end_pert

 contains

!-------------------------------------------------------------------------------

 subroutine fv_init_pert(Atm, AtmP) 

 type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
 type(fv_atmos_pert_type), allocatable, intent(inout), target :: AtmP(:)

 integer :: n, ntilesMe

  call run_setup_pert(AtmP,Atm)

  ntilesMe = size(AtmP(:)) 

  do n=1,ntilesMe

     call allocate_fv_atmos_pert_type( AtmP(n), Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, &
                                                Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                                                Atm(n)%flagstruct%npz,  Atm(n)%flagstruct%ncnst )

  enddo

 end subroutine fv_init_pert

!-------------------------------------------------------------------------------

 subroutine fv_end_pert(AtmP)

  type(fv_atmos_pert_type), intent(inout) :: AtmP(:)

  integer :: n, ntilesMe

  ntilesMe = size(AtmP(:)) 

  do n=1,ntilesMe

   call deallocate_fv_atmos_pert_type(AtmP(n))

  enddo

 end subroutine fv_end_pert

!-------------------------------------------------------------------------------

 subroutine setup_pointers_pert(AtmP)

  type(fv_atmos_pert_type), intent(INOUT), target :: AtmP

   !Linearized model pointers
   split_hord        => AtmP%flagstruct%split_hord
   hord_mt_pert      => AtmP%flagstruct%hord_mt_pert
   hord_vt_pert      => AtmP%flagstruct%hord_vt_pert
   hord_tm_pert      => AtmP%flagstruct%hord_tm_pert
   hord_dp_pert      => AtmP%flagstruct%hord_dp_pert
   hord_tr_pert      => AtmP%flagstruct%hord_tr_pert
   split_kord        => AtmP%flagstruct%split_kord
   kord_mt_pert      => AtmP%flagstruct%kord_mt_pert
   kord_wz_pert      => AtmP%flagstruct%kord_wz_pert
   kord_tm_pert      => AtmP%flagstruct%kord_tm_pert
   kord_tr_pert      => AtmP%flagstruct%kord_tr_pert
   split_damp        => AtmP%flagstruct%split_damp
   nord_pert         => AtmP%flagstruct%nord_pert
   dddmp_pert        => AtmP%flagstruct%dddmp_pert
   d2_bg_pert        => AtmP%flagstruct%d2_bg_pert
   d4_bg_pert        => AtmP%flagstruct%d4_bg_pert
   do_vort_damp_pert => AtmP%flagstruct%do_vort_damp_pert
   vtdm4_pert        => AtmP%flagstruct%vtdm4_pert
   split_damp_tr     => AtmP%flagstruct%split_damp_tr
   nord_tr_pert      => AtmP%flagstruct%nord_tr_pert
   trdm2_pert        => AtmP%flagstruct%trdm2_pert

  end subroutine setup_pointers_pert

!-------------------------------------------------------------------------------

 subroutine run_setup_pert(AtmP,Atm)

  type(fv_atmos_pert_type), intent(inout), target :: AtmP(:)
  type(fv_atmos_type), intent(inout), target :: Atm(:)

  integer :: f_unit, n, ierr, ios, unit
  character(len=80) :: nested_grid_filename
 
  namelist /fv_core_pert_nml/split_hord, hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, hord_tr_pert, &
                             split_kord, kord_mt_pert, kord_wz_pert, kord_tm_pert, kord_tr_pert, split_damp, do_vort_damp_pert, &
                             nord_pert, dddmp_pert, d2_bg_pert, d4_bg_pert, vtdm4_pert, &
                             split_damp_tr, nord_tr_pert, trdm2_pert

  do n=1,size(AtmP)

     call setup_pointers_pert(AtmP(n))

     if (size(AtmP) == 1) then
        f_unit = open_namelist_file('inputpert.nml')
     else if (n == 1) then
        f_unit = open_namelist_file('inputpert.nml')
     else 
        write(nested_grid_filename,'(A10, I2.2, A4)') 'input_nest_pert', n, '.nml'
        f_unit = open_namelist_file(nested_grid_filename)
     endif

     !Read linearized FVCORE namelist
     rewind (f_unit)
     read (f_unit,fv_core_pert_nml,iostat=ios)
     ierr = check_nml_error(ios,'fv_core_pert_nml')

     call close_file(f_unit)

     unit = stdlog()
     write(unit, nml=fv_core_pert_nml)

     !Unless specfied the trajectory uses the coeffs suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_damp) then
        Atm(n)%flagstruct%nord         = AtmP(n)%flagstruct%nord_pert
        Atm(n)%flagstruct%dddmp        = AtmP(n)%flagstruct%dddmp_pert
        Atm(n)%flagstruct%d2_bg        = AtmP(n)%flagstruct%d2_bg_pert
        Atm(n)%flagstruct%d4_bg        = AtmP(n)%flagstruct%d4_bg_pert
        Atm(n)%flagstruct%do_vort_damp = AtmP(n)%flagstruct%do_vort_damp_pert
        Atm(n)%flagstruct%vtdm4        = AtmP(n)%flagstruct%vtdm4_pert
     endif

     !Unless specfied the trajectory uses the coeffs suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_damp_tr) then
        Atm(n)%flagstruct%nord_tr = AtmP(n)%flagstruct%nord_tr_pert
        Atm(n)%flagstruct%trdm2   = AtmP(n)%flagstruct%trdm2_pert
     endif

     !Unless specfied the trajectory uses hord suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_hord) then
        Atm(n)%flagstruct%hord_mt = AtmP(n)%flagstruct%hord_mt_pert
        Atm(n)%flagstruct%hord_vt = AtmP(n)%flagstruct%hord_vt_pert
        Atm(n)%flagstruct%hord_tm = AtmP(n)%flagstruct%hord_tm_pert
        Atm(n)%flagstruct%hord_dp = AtmP(n)%flagstruct%hord_dp_pert
        Atm(n)%flagstruct%hord_tr = AtmP(n)%flagstruct%hord_tr_pert
     endif

     !Unless specfied the trajectory uses hord suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_kord) then
        Atm(n)%flagstruct%kord_mt = AtmP(n)%flagstruct%kord_mt_pert
        Atm(n)%flagstruct%kord_wz = AtmP(n)%flagstruct%kord_wz_pert
        Atm(n)%flagstruct%kord_tm = AtmP(n)%flagstruct%kord_tm_pert
        Atm(n)%flagstruct%kord_tr = AtmP(n)%flagstruct%kord_tr_pert
     endif

   enddo

  end subroutine run_setup_pert
       
end module fv_control_nlm_mod
