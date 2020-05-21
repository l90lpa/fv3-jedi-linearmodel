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
    fv3jedi_lm_mod.F90
    ${dynamics_src_files}
    ${physics_src_files}
    ${utils_src_files}
)

ecbuild_add_library( TARGET   fv3jedilm
                     SOURCES  ${fv3jedilm_src_files}
                     LIBS     ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                     INCLUDES ${FMS_INCLUDES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${FV3JEDILM_LINKER_LANGUAGE}
                    )
