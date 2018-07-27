module fv3jedi_lm_turbulence_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

! BL Schemes
use bldriver

! Saturation table
use qsat_util

use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter, finalize_cp_iter
use tapenade_iter, only: cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray

!> Turbulence module

implicit none
private
public :: fv3jedi_lm_turbulence_type

!> Local trajectory object
type local_traj_turbulence
  real(kind_real), allocatable, dimension(:,:,:) :: akq, bkq, ckq
  real(kind_real), allocatable, dimension(:,:,:) :: aks, bks, cks
  real(kind_real), allocatable, dimension(:,:,:) :: akv, bkv, ckv
  real(kind_real), allocatable, dimension(:,:,:) :: pk
endtype local_traj_turbulence

!> Local constants object
type local_cnst_turbulence
  integer :: im,jm,lm
  real(kind_real), dimension(22)   :: TURBPARAMS
  integer, dimension(4) :: TURBPARAMSI
endtype local_cnst_turbulence

!> turbulence class (self)
type fv3jedi_lm_turbulence_type
 type(local_traj_turbulence) :: ltraj
 type(local_cnst_turbulence) :: lcnst
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type fv3jedi_lm_turbulence_type

contains

! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 real(kind_real), allocatable, dimension(:) :: pref
 integer :: l

 self%lcnst%im = conf%npx
 self%lcnst%jm = conf%npy
 self%lcnst%lm = conf%npz

 call allocate_ltraj(self%lcnst%im,self%lcnst%jm,self%lcnst%lm,self%ltraj)

 allocate(pref(0:self%lcnst%lm))
 DO l = 0,self%lcnst%lm
   pref(l) = conf%AK(l+1) + conf%BK(l+1)*p00
 enddo

 !Turbulence Parameters
 self%lcnst%TurbParams(1)  = 5.0_kind_real                    !LOUIS
 self%lcnst%TurbParams(2)  = 160.0_kind_real                  !LAMBDAM
 self%lcnst%TurbParams(3)  = 1.0_kind_real                    !LAMBDAM2
 self%lcnst%TurbParams(4)  = 160.0_kind_real                  !LAMBDAH
 self%lcnst%TurbParams(5)  = 1.0_kind_real                    !LAMBDAH2
 self%lcnst%TurbParams(6)  = 3000.0_kind_real                 !ZKMENV
 self%lcnst%TurbParams(7)  = 3000.0_kind_real                 !ZKHENV
 self%lcnst%TurbParams(8)  = 0.1_kind_real                    !MINTHICK
 self%lcnst%TurbParams(9)  = 0.0030_kind_real                 !MINSHEAR
 self%lcnst%TurbParams(10) = 2.5101471e-8_kind_real           !C_B
 self%lcnst%TurbParams(11) = 1500.0_kind_real                 !LAMBDA_B 
 self%lcnst%TurbParams(12) = 500.0_kind_real                  !AKHMMAX
 self%lcnst%TurbParams(13) = 1.0_kind_real                    !PRANDTLSFC
 self%lcnst%TurbParams(14) = 0.75_kind_real                   !PRANDTLRAD
 self%lcnst%TurbParams(15) = 0.50_kind_real                   !BETA_RAD
 self%lcnst%TurbParams(16) = 0.25_kind_real                   !BETA_SURF
 self%lcnst%TurbParams(17) = 0.85_kind_real                   !KHRADFAC
 self%lcnst%TurbParams(18) = 0.45_kind_real                   !KHSFCFAC
 self%lcnst%TurbParams(19) = 20.0_kind_real                   !TPFAC_SURF
 self%lcnst%TurbParams(20) = 1.5e-3_kind_real                 !ENTRATE_SURF
 self%lcnst%TurbParams(21) = 0.5_kind_real                    !PCEFF_SURF
 self%lcnst%TurbParams(22) = -999.0_kind_real                 !LOUIS_MEMORY
 self%lcnst%TurbParamsI(1) = count(pref < 50000.0_kind_real)  !KPBLMIN
 self%lcnst%TurbParamsI(2) = 1                                !LOCK_ON
 self%lcnst%TurbParamsI(3) = 1                                !PBLHT_OPTION
 self%lcnst%TurbParamsI(4) = 0                                !RADLW_DEP

 deallocate(pref)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout), target :: self
 type(fv3jedi_lm_traj), intent(inout) :: traj
 type(fv3jedi_lm_conf), intent(in) :: conf

 integer :: i,j,k
 type(local_traj_turbulence), pointer :: ltraj
 type(local_cnst_turbulence), pointer :: lcnst

 integer :: im,jm,lm

 !Set up the local trajectory
 call set_ltraj(conf,self%lcnst,traj,self%ltraj)

 !Convenience pointers
 ltraj => self%ltraj
 lcnst => self%lcnst

 im = lcnst%im
 jm = lcnst%jm
 lm = lcnst%lm

 traj%t = traj%t / ltraj%pk

 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,traj%u ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,traj%v ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%aks,ltraj%bks,ltraj%cks,traj%t ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,traj%qv,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,traj%qi,1,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,traj%ql,1,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,traj%o3,1,0)

 traj%t = traj%t * ltraj%pk 

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 integer :: i,j,k
 type(local_traj_turbulence), pointer :: ltraj
 type(local_cnst_turbulence), pointer :: lcnst

 integer :: im,jm,lm

 !Set up the local trajectory
 call set_ltraj(conf,self%lcnst,traj,self%ltraj)

 !Convenience pointers
 ltraj => self%ltraj
 lcnst => self%lcnst

 im = lcnst%im
 jm = lcnst%jm
 lm = lcnst%lm

 pert%t = pert%t / ltraj%pk 

 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,pert%u ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,pert%v ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%aks,ltraj%bks,ltraj%cks,pert%t ,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%qv,1,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%qi,1,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%ql,1,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%o3,1,0)

 pert%t = pert%t * ltraj%pk 

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 integer :: i,j,k
 type(local_traj_turbulence), pointer :: ltraj
 type(local_cnst_turbulence), pointer :: lcnst

 integer :: im,jm,lm

 !Set up the local trajectory
 call set_ltraj(conf,self%lcnst,traj,self%ltraj)

 !Convenience pointers
 ltraj => self%ltraj
 lcnst => self%lcnst

 im = lcnst%im
 jm = lcnst%jm
 lm = lcnst%lm

 pert%t = pert%t * ltraj%pk

 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,pert%u ,2,1)
 call vtrisolvepert(im,jm,lm,ltraj%akv,ltraj%bkv,ltraj%ckv,pert%v ,2,1)
 call vtrisolvepert(im,jm,lm,ltraj%aks,ltraj%bks,ltraj%cks,pert%t ,2,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%qv,2,1)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%qi,2,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%ql,2,0)
 call vtrisolvepert(im,jm,lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%o3,2,0)

 pert%t = pert%t / ltraj%pk

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self)

 implicit none
 class(fv3jedi_lm_turbulence_type), intent(inout) :: self

 call deallocate_ltraj(self%ltraj)

endsubroutine delete

! ------------------------------------------------------------------------------

subroutine set_ltraj(conf,lcnst,traj,ltraj)

 type(fv3jedi_lm_conf),       intent(in)    :: conf
 type(local_cnst_turbulence), intent(in)    :: lcnst
 type(fv3jedi_lm_traj),       intent(in)    :: traj
 type(local_traj_turbulence), intent(inout) :: ltraj

 integer :: i,j,l,im,jm,lm,isc,iec,jsc,jec 

 real(kind_real), allocatable, dimension(:,:,:) :: PTT1
 real(kind_real), allocatable, dimension(:,:)   :: ZPBL1, CT1
 real(kind_real), allocatable, dimension(:,:,:) :: EKV, FKV
 real(kind_real), allocatable, dimension(:,:,:) :: pet, pmt
 real(kind_real), allocatable, dimension(:,:,:) :: QIT1, QLT1
 real(kind_real), allocatable, dimension(:,:,:) :: fQi


 im = lcnst%im
 jm = lcnst%jm
 lm = lcnst%lm

 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 allocate(PTT1 (im,jm,lm))
 allocate(ZPBL1(im,jm))
 allocate(CT1  (im,jm))
 allocate(EKV  (im,jm,lm))
 allocate(FKV  (im,jm,lm))
 allocate(pet  (im,jm,0:lm))
 allocate(pmt  (im,jm,lm))
 allocate(QIT1 (im,jm,lm))
 allocate(QLT1 (im,jm,lm))
 allocate(fQi  (im,jm,lm))

 !Compute pressures from delp
 call compute_pressures(im,jm,lm,conf%ptop,traj%delp(isc:iec,jsc:jec,:),&
                        pet,pmt,ltraj%pk)

 !Use local copied to avoid overwrite
 ZPBL1(1:im,1:jm) = traj%ZPBL(isc:iec,jsc:jec)
 CT1(1:im,1:jm) = traj%CT(isc:iec,jsc:jec)

 ptt1(1:im,1:jm,:) = traj%t(isc:iec,jsc:jec,:) * ltraj%pk(1:im,1:jm,:)

 !Calculate total cloud ice and liquid trajectory
 if (conf%do_phy_mst == 0) then
 
    QIT1 = traj%QI(isc:iec,jsc:jec,:)
    QLT1 = traj%QL(isc:iec,jsc:jec,:)
 
 else
 
   DO l = 1,LM
     DO j = 1,JM
       DO i = 1,IM
         call IceFraction( traj%t(isc+i-1,jsc+j-1,l), fQi(i,j,l) )
       enddo
     enddo
  enddo
 
  QIT1(1:im,1:jm,:) = (traj%QLS(isc:iec,jsc:jec,:) + traj%QCN(isc:iec,jsc:jec,:)) * fQi(1:im,1:jm,:)
  QLT1(1:im,1:jm,:) = (traj%QLS(isc:iec,jsc:jec,:) + traj%QCN(isc:iec,jsc:jec,:)) * (1-fQi(1:im,1:jm,:))
 
 endif

 !Initialize tri-diagonal arrays
 ltraj%AKV = 0.0
 ltraj%BKV = 0.0
 ltraj%CKV = 0.0
 ltraj%AKS = 0.0
 ltraj%BKS = 0.0
 ltraj%CKS = 0.0
 ltraj%AKQ = 0.0
 ltraj%BKQ = 0.0
 ltraj%CKQ = 0.0
 
 !Call boundary layer routine. These routines will return the lower (AK*), main (BK*),
 !and upper (CK*) diagonals.
 
 call BL_DRIVER( IM           , &
                 JM           , &
                 LM           , &
                 conf%DT           , &
                 traj%U          , &
                 traj%V          , &
                 PTT1         , &
                 traj%QV          , &
                 PET          , &
                 QIT1         , &
                 QLT1         , &
                 traj%FRLAND       , &
                 traj%FROCEAN      , &
                 traj%VARFLT       , &
                 ZPBL1        , &
                 traj%CM           , &
                 CT1          , &
                 traj%CQ           , &
                 lcnst%TURBPARAMS   , &
                 lcnst%TURBPARAMSI  , &
                 traj%USTAR        , &
                 traj%BSTAR        , &
                 ltraj%AKS, ltraj%BKS, ltraj%CKS, &
                 ltraj%AKQ, ltraj%BKQ, ltraj%CKQ, &
                 ltraj%AKV, ltraj%BKV, ltraj%CKV, &
                 EKV, FKV       &
               )

 !Solver part 1, perform LU decomposition
 call VTRILUPERT(IM,JM,LM,ltraj%AKV,ltraj%BKV,ltraj%CKV)
 call VTRILUPERT(IM,JM,LM,ltraj%AKS,ltraj%BKS,ltraj%CKS)
 call VTRILUPERT(IM,JM,LM,ltraj%AKQ,ltraj%BKQ,ltraj%CKQ)

 deallocate(PTT1 )
 deallocate(ZPBL1)
 deallocate(CT1  )
 deallocate(EKV  )
 deallocate(FKV  )
 deallocate(pet  )
 deallocate(pmt  )
 deallocate(QIT1 )
 deallocate(QLT1 )
 deallocate(fQi  )

endsubroutine set_ltraj

! ------------------------------------------------------------------------------

subroutine allocate_ltraj(im,jm,lm,ltraj)

 integer, intent(in) :: im,jm,lm
 type(local_traj_turbulence), intent(inout) :: ltraj

 allocate(ltraj%pk(im,jm,lm))
 allocate(ltraj%akv(im,jm,lm))
 allocate(ltraj%bkv(im,jm,lm))
 allocate(ltraj%ckv(im,jm,lm))
 allocate(ltraj%aks(im,jm,lm))
 allocate(ltraj%bks(im,jm,lm))
 allocate(ltraj%cks(im,jm,lm))
 allocate(ltraj%akq(im,jm,lm))
 allocate(ltraj%bkq(im,jm,lm))
 allocate(ltraj%ckq(im,jm,lm))

endsubroutine allocate_ltraj

! ------------------------------------------------------------------------------

subroutine deallocate_ltraj(ltraj)

 type(local_traj_turbulence), intent(inout) :: ltraj

 deallocate(ltraj%pk )
 deallocate(ltraj%akv)
 deallocate(ltraj%bkv)
 deallocate(ltraj%ckv)
 deallocate(ltraj%aks)
 deallocate(ltraj%bks)
 deallocate(ltraj%cks)
 deallocate(ltraj%akq)
 deallocate(ltraj%bkq)
 deallocate(ltraj%ckq)

endsubroutine deallocate_ltraj

! ------------------------------------------------------------------------------

subroutine vtrilupert(im,jm,lm,a,b,c)

 !perform lu decomposition of tridiagonal system

 implicit none

 integer,                   intent(in   ) :: im, jm, lm
 real(kind_real), dimension(im,jm,lm), intent(in   ) :: c
 real(kind_real), dimension(im,jm,lm), intent(inout) :: a, b

 integer :: l

 b(:,:,1) = 1. / b(:,:,1)
 do l = 2,lm
    a(:,:,l) = a(:,:,l) * b(:,:,l-1)
    b(:,:,l) = 1. / ( b(:,:,l) - c(:,:,l-1) * a(:,:,l) )
 end do

end subroutine vtrilupert

! ------------------------------------------------------------------------------

subroutine vtrisolvepert(im,jm,lm,a,b,c,y,phase,ygswitch)

 !solve lu decomposed tridiagonal system

 implicit none

 !arguments
 integer,                   intent(in   ) :: im, jm, lm, phase, ygswitch
 real(kind_real), dimension(im,jm,lm), intent(in   ) :: a, b, c
 real(kind_real), dimension(im,jm,lm), intent(inout) :: y

 !locals
 integer :: l

 if (phase == 1) then !TLM

    !solve (lu)ynew = yold

    !sweep down, modifying rhs with multiplier a
    do l=2,lm
       y(:,:,l) = y(:,:,l) - a(:,:,l)*y(:,:,l-1)
    enddo

    !sweep up, solving for updated value. note b has the inverse of the main diagonal
    if (ygswitch == 1) then !winds, temperature and q
       y(:,:,lm) = y(:,:,lm)*b(:,:,lm) !ygprime = 0
    else !tracers
       y(:,:,lm) = y(:,:,lm)*b(:,:,lm-1)/(b(:,:,lm-1) - a(:,:,lm)*(1.0+c(:,:,lm-1)*b(:,:,lm-1) ))
    endif

    do l=lm-1,1,-1
       y(:,:,l) = b(:,:,l)*(y(:,:,l)-c(:,:,l)*y(:,:,l+1) )
    enddo

 elseif (phase == 2) then !ADM

    if (ygswitch == 1) then !can solve (lu)'ynew = u'l'ynew = yold

       !sweep down but with u' instead of l
       y(:,:,1) = y(:,:,1)*b(:,:,1)
       do l=2,lm
          y(:,:,l) = b(:,:,l) * (y(:,:,l) - c(:,:,l-1)*y(:,:,l-1))
       enddo

       !sweep up but with l' instead of u
       do l=lm-1,1,-1
          y(:,:,l) = y(:,:,l) - a(:,:,l+1)*y(:,:,l+1)
       enddo

    else !change for surface means transpose doesnt work so use line-by-line adjoint

       !adjoint of sweep up, solving for updated value. note b has the inverse of the main diagonal
       do l=1,lm-1,1
          y(:,:,l+1) = y(:,:,l+1) - c(:,:,l)*b(:,:,l)*y(:,:,l)
          y(:,:,l) = b(:,:,l)*y(:,:,l)
       end do

       !adjoint of surface fix
       y(:,:,lm) = b(:,:,lm-1)*y(:,:,lm)/(b(:,:,lm-1)-a(:,:,lm)*(c(:,:,lm-1)*b(:,:,lm-1)+1.0))

       !adjoint of sweep down, modifying rhs with multiplier a
       do l=lm,2,-1
          y(:,:,l-1) = y(:,:,l-1) - a(:,:,l)*y(:,:,l)
       end do

    endif

 endif

end subroutine vtrisolvepert

! ------------------------------------------------------------------------------

end module fv3jedi_lm_turbulence_mod
