! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_lm_kinds_mod

  use, intrinsic :: iso_fortran_env, only: REAL64

  implicit none
  private
  public kind_real
  
  integer, parameter :: kind_real=REAL64

end module fv3jedi_lm_kinds_mod
