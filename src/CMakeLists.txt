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

if( NOT FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT FV3_FORECAST_MODEL MATCHES "UFS" )
  list( APPEND FV3JEDILM_LIB_DEP fms )
endif()

ecbuild_add_library( TARGET   fv3jedilm
                     SOURCES  ${fv3jedilm_src_files}
                     PUBLIC_LIBS NetCDF::NetCDF_Fortran
                     PRIVATE_LIBS     ${FV3JEDILM_LIB_DEP}
                     INSTALL_HEADERS LISTED )
include_directories( fv3jedilm ${DEPEND_LIB_ROOT}/include_r8 )
target_include_directories( fv3jedilm INTERFACE "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
                                                "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>")

if(FV3_FORECAST_MODEL MATCHES "UFS")
  target_link_libraries(fv3jedilm PUBLIC FMS::fms_r8)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    target_compile_options(fv3jedilm PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
endif()

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( fv3jedilm INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )
