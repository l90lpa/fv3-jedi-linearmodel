list( APPEND physics_files
fv3jedi_lm_physics_mod.F90
moist/cloud_ad.F90
moist/cloud.F90
moist/cloud_tl.F90
moist/convection_ad.F90
moist/convection.F90
moist/convection_tl.F90
moist/fv3jedi_lm_moist_mod.F90
moist/qsat_util.F90
turbulence/bldriver.F90
turbulence/blsimp.F90
turbulence/fv3jedi_lm_turbulence_mod.F90
)

PREPEND( _p_physics_files       "physics"       ${physics_files} )

set ( physics_src_files
      ${_p_physics_files}
      PARENT_SCOPE
)
