module fv3jedi_lm_dynutils_mod

use fv3jedi_lm_kinds_mod
use fv_arrays_nlm_mod, only: fv_diag_type

implicit none

private
public init_fv_diag_type

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine init_fv_diag_type(diag)

implicit none
type(fv_diag_type), intent(inout) :: diag

diag%id_ps = 0
diag%id_slp = 0
diag%id_ua = 0
diag%id_va = 0
diag%id_pt = 0
diag%id_omga = 0
diag%id_vort = 0
diag%id_tm = 0
diag%id_pv = 0
diag%id_zsurf = 0
diag%id_oro = 0
diag%id_sgh = 0
diag%id_divg = 0
diag%id_w = 0
diag%id_ke = 0
diag%id_te = 0
diag%id_zs = 0
diag%id_ze = 0
diag%id_mq = 0
diag%id_vorts = 0
diag%id_us = 0
diag%id_vs = 0
diag%id_tq = 0
diag%id_rh = 0
diag%id_c15 = 0
diag%id_c25 = 0
diag%id_c35 = 0
diag%id_c45 = 0
diag%id_f15 = 0
diag%id_f25 = 0
diag%id_f35 = 0
diag%id_f45 = 0
diag%id_ctp = 0
diag%id_ppt = 0
diag%id_ts = 0
diag%id_tb = 0
diag%id_ctt = 0
diag%id_pmask = 0
diag%id_pmaskv2 = 0
diag%id_delp = 0
diag%id_delz = 0
diag%id_zratio = 0
diag%id_ws = 0
diag%id_iw = 0
diag%id_lw = 0
diag%id_pfhy = 0
diag%id_pfnh = 0
diag%id_qn = 0
diag%id_qn200 = 0
diag%id_qn500 = 0
diag%id_qn850 = 0
diag%id_qp = 0
diag%id_mdt = 0
diag%id_qdt = 0
diag%id_aam = 0
diag%id_amdt = 0
diag%id_acly = 0
diag%id_acl = 0
diag%id_acl2 = 0
diag%id_dbz = 0
diag%id_maxdbz = 0
diag%id_basedbz = 0
diag%id_dbz4km = 0
diag%id_vort200 = 0
diag%id_vort500 = 0
diag%id_w500 = 0
diag%id_w700 = 0
diag%id_vort850 = 0
diag%id_w850 = 0
diag%id_x850 = 0
diag%id_srh = 0
diag%id_srh25 = 0
diag%id_srh01 = 0
diag%id_uh03 = 0
diag%id_uh25 = 0
diag%id_theta_e = 0
diag%id_w200 = 0
diag%id_s200 = 0
diag%id_sl12 = 0
diag%id_sl13 = 0
diag%id_w5km = 0
diag%id_rain5km = 0
diag%id_w2500m = 0
diag%id_u_plev = 0
diag%id_v_plev = 0
diag%id_t_plev = 0
diag%id_h_plev = 0
diag%id_q_plev = 0
diag%id_omg_plev = 0
diag%id_rh10 = 0
diag%id_rh50 = 0
diag%id_rh100 = 0
diag%id_rh200 = 0
diag%id_rh250 = 0
diag%id_rh300 = 0
diag%id_rh500 = 0
diag%id_rh700 = 0
diag%id_rh850 = 0
diag%id_rh925 = 0
diag%id_rh1000 = 0
diag%id_rh1000_cmip = 0
diag%id_rh925_cmip = 0
diag%id_rh850_cmip = 0
diag%id_rh700_cmip = 0
diag%id_rh500_cmip = 0
diag%id_rh300_cmip = 0
diag%id_rh250_cmip = 0
diag%id_rh100_cmip = 0
diag%id_rh50_cmip = 0
diag%id_rh10_cmip = 0
diag%id_hght = 0
diag%id_u100m = 0
diag%id_v100m = 0
diag%id_w100m = 0
diag%ic_ps = 0
diag%ic_ua = 0
diag%ic_va = 0
diag%ic_ppt = 0
diag%ic_sphum = 0

diag%sphum = 0.0_kind_real
diag%liq_wat = 0.0_kind_real
diag%ice_wat = 0.0_kind_real
diag%rainwat = 0.0_kind_real
diag%snowwat = 0.0_kind_real
diag%graupel = 0.0_kind_real
diag%efx = 0.0_kind_real
diag%efx_sum = 0.0_kind_real
diag%efx_nest = 0.0_kind_real
diag%efx_sum_nest = 0.0_kind_real
diag%mtq = 0.0_kind_real
diag%mtq_sum = 0.0_kind_real
diag%steps = 0.0_kind_real

end subroutine init_fv_diag_type

! --------------------------------------------------------------------------------------------------

end module fv3jedi_lm_dynutils_mod
