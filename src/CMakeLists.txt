list( APPEND fv3jedilm_src_files
fv3jedi_lm_mod.F90
physics/fv3jedi_lm_physics_mod.F90
physics/moist/cloud_ad.F90
physics/moist/cloud.F90
physics/moist/cloud_tl.F90
physics/moist/convection_ad.F90
physics/moist/convection.F90
physics/moist/convection_tl.F90
physics/moist/fv3jedi_lm_moist_mod.F90
physics/moist/qsat_util.F90
physics/turbulence/bldriver.F90
physics/turbulence/blsimp.F90
physics/turbulence/fv3jedi_lm_turbulence_mod.F90
utils/fv3jedi_lm_const_mod.F90
utils/fv3jedi_lm_kinds_mod.F90
utils/fv3jedi_lm_utils_mod.F90
utils/MAPL_Constants.F90
utils/tapenade/adBuffer.f
utils/tapenade/adStack.c
utils/tapenade/tapenade_iter.F90
)

function(PREPEND var prefix )
    set ( listVar "" )
    foreach (f ${ARGN})
        list (APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set ( ${var} "${listVar}" PARENT_SCOPE )
endfunction(PREPEND)

add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(utils)


list( APPEND fv3jedilm_src_files
    ${dynamics_src_files}
    ${physics_src_files}
    ${utils_src_files}
)

ecbuild_add_library( TARGET   fv3jedilm
                     SOURCES  ${fv3jedilm_src_files}
                     LIBS     ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${FMS_LIBRARIES}
                     INCLUDES ${FMS_INCLUDES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
