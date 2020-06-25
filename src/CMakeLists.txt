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

#Requirement libraries
set(FV3JEDILM_LIB_DEP ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES})

if( NOT ${BUILD_WITH_MODEL} MATCHES "GEOS" AND NOT ${BUILD_WITH_MODEL} MATCHES "UFS" )
  list( APPEND FV3JEDILM_LIB_DEP fms )
endif()

ecbuild_add_library( TARGET   fv3jedilm
                     SOURCES  ${fv3jedilm_src_files}
                     LIBS     ${FV3JEDILM_LIB_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${FV3JEDILM_LINKER_LANGUAGE}
                    )

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    target_compile_options(fv3jedilm PRIVATE -fallow-argument-mismatch)
endif()
