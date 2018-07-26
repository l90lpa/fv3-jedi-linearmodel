! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on pressure variables for fv3 
!> Daniel Holdaway, NASA/JCSDA

module fv_pressure_mod

use MAPL_ConstantsMod, only: kappa => MAPL_KAPPA

implicit none
public

contains

!----------------------------------------------------------------------------
! Compute all pressures needed as input to fv3 ------------------------------
!----------------------------------------------------------------------------

subroutine compute_fv3_pressures( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                  kappa, ptop, delp, pe, pk, pkz, peln )

 implicit none
 integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
 real, intent(in) :: kappa, ptop
 real, intent(in) :: delp(isd:ied, jsd:jed, npz)
 real, intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
 real, intent(out) :: pk(is:ie, js:je, npz+1)
 real, intent(out) :: peln(is:ie, npz+1, js:je)
 real, intent(out) :: pkz(is:ie, js:je, npz)
 integer :: i, j, k

 pe(:, :, :) = 0.0
 pe(:, 1, :) = ptop
 do k=2,npz+1
   do j=js,je
     do i=is,ie
       pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
     end do
   end do
 end do

 do k=1,npz+1
   do j=js,je
     do i=is,ie
       peln(i, k, j) = log(pe(i, k, j))
     end do
   end do
 end do

 do k=1,npz+1
   do j=js,je
     do i=is,ie
       pk(i, j, k) = exp(kappa*peln(i, k, j))
     end do
   end do
 end do

 do k=1,npz
   do j=js,je
     do i=is,ie
       pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1, j)-peln(i, k, j)))
     end do
   end do
 end do

end subroutine compute_fv3_pressures

subroutine compute_fv3_pressures_tlm( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                      kappa, ptop, delp, delp_tl, &
                                      pe, pe_tl, pk, pk_tl, pkz, pkz_tl, peln, peln_tl )

 implicit none
 integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
 real, intent(in) :: kappa, ptop
 real, intent(in) :: delp(isd:ied, jsd:jed, npz)
 real, intent(in) :: delp_tl(isd:ied, jsd:jed, npz)
 real, intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
 real, intent(out) :: pe_tl(is-1:ie+1, npz+1, js-1:je+1)
 real, intent(out) :: pk(is:ie, js:je, npz+1)
 real, intent(out) :: pk_tl(is:ie, js:je, npz+1)
 real, intent(out) :: peln(is:ie, npz+1, js:je)
 real, intent(out) :: peln_tl(is:ie, npz+1, js:je)
 real, intent(out) :: pkz(is:ie, js:je, npz)
 real, intent(out) :: pkz_tl(is:ie, js:je, npz)
 integer :: i, j, k

 pe(:, :, :) = 0.0
 pe(:, 1, :) = ptop
 pe_tl = 0.0
 do k=2,npz+1
   do j=js,je
     do i=is,ie
       pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
       pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
     end do
   end do
 end do

 peln_tl = 0.0
 do k=1,npz+1
   do j=js,je
     do i=is,ie
       peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
       peln(i, k, j) = log(pe(i, k, j))
     end do
   end do
 end do

 pk_tl = 0.0
 do k=1,npz+1
   do j=js,je
     do i=is,ie
       pk_tl(i, j, k) = kappa*peln_tl(i, k, j)*exp(kappa*peln(i, k, j))
       pk(i, j, k) = exp(kappa*peln(i, k, j))
     end do
   end do
 end do

 pkz_tl = 0.0
 do k=1,npz
   do j=js,je
     do i=is,ie
       pkz_tl(i, j, k) = ((pk_tl(i, j, k+1)-pk_tl(i, j, k))*kappa*(peln   (i, k+1, j)-peln   (i, k, j)) &
                         -(pk   (i, j, k+1)-pk   (i, j, k))*kappa*(peln_tl(i, k+1, j)-peln_tl(i, k, j))) &
                         /(kappa*(peln(i, k+1, j)-peln(i, k, j)))**2
       pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1, j)-peln(i, k, j)))
     end do
   end do
 end do

end subroutine compute_fv3_pressures_tlm

subroutine compute_fv3_pressures_bwd( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                      kappa, ptop, delp, delp_ad, &
                                      pe, pe_ad, pk, pk_ad, pkz, pkz_ad, peln, peln_ad)

  implicit none
  integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
  real, intent(in) :: kappa, ptop
  real, intent(in)    :: delp(isd:ied, jsd:jed, npz)
  real, intent(inout) :: delp_ad(isd:ied, jsd:jed, npz)
  real, intent(in)    :: pe(is-1:ie+1, npz+1, js-1:je+1)
  real, intent(inout) :: pe_ad(is-1:ie+1, npz+1, js-1:je+1)
  real, intent(in)    :: pk(is:ie, js:je, npz+1)
  real, intent(inout) :: pk_ad(is:ie, js:je, npz+1)
  real, intent(in)    :: peln(is:ie, npz+1, js:je)
  real, intent(inout) :: peln_ad(is:ie, npz+1, js:je)
  real, intent(in)    :: pkz(is:ie, js:je, npz)
  real, intent(inout) :: pkz_ad(is:ie, js:je, npz)
  integer :: i, j, k
  real :: temp_tj1, temp_ad1, temp_ad2

  do k=npz,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        temp_tj1 = kappa*(peln(i, k+1, j)-peln(i, k, j))
        temp_ad1 = pkz_ad(i, j, k)/temp_tj1
        temp_ad2 = -((pk(i, j, k+1)-pk(i, j, k))*kappa*temp_ad1/temp_tj1)
        pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad1
        pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad1
        peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad2
        peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad2
        pkz_ad(i, j, k) = 0.0
      end do
    end do
  end do

  do k=npz+1,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        peln_ad(i, k, j) = peln_ad(i, k, j) + exp(kappa*peln(i, k, j))*kappa*pk_ad(i, j, k)
        pk_ad(i, j, k) = 0.0
      end do
    end do
  end do

  do k=npz+1,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k, j)
        peln_ad(i, k, j) = 0.0
      end do
    end do
  end do

  do k=npz+1,2,-1
    do j=je,js,-1
      do i=ie,is,-1
        pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
        delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
        pe_ad(i, k, j) = 0.0
      end do
    end do
  end do

end subroutine compute_fv3_pressures_bwd

!----------------------------------------------------------------------------

end module pressure_vt_mod
