# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDILM
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedilm Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH} )

set( CMAKE_DIRECTORY_LABELS "fv3-jedi-lm" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( FV3JEDILM_LINKER_LANGUAGE Fortran )

################################################################################
# Dependencies
################################################################################

# MKL / LAPACK
# Default: use MKL as LAPACK implementation.
# Set -DENABLE_MKL=Off to disable use of MKL and prefer LAPACK
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

#openMP
find_package( OpenMP COMPONENTS Fortran )

# MPI
ecbuild_add_option( FEATURE MPI DEFAULT ON
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )
ecbuild_find_mpi( COMPONENTS Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_Fortran_LIBRARIES})

# NetCDF
set( NETCDF_F90 ON CACHE BOOL "Compile with Fortran NetCDF" )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# Set default model
if (NOT FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT FV3_FORECAST_MODEL MATCHES "UFS")
  set(FV3_FORECAST_MODEL "FV3CORE")
endif()

# fms
if (FV3_FORECAST_MODEL MATCHES "FV3CORE")
  ecbuild_use_package( PROJECT fms REQUIRED )
  include_directories( ${FMS_INCLUDE_DIRS} )
  list( APPEND FV3JEDILM_TPLS fms )
endif()

# geos
if(FV3_FORECAST_MODEL MATCHES "GEOS")
  find_package(MAPL REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake/MAPL )
  find_package(GEOSgcm REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake/GEOSgcm )
  if (FV3_PRECISION MATCHES DOUBLE OR NOT FV3_PRECISION)
    link_libraries(fms_r8)
  else()
    link_libraries(fms_r4)
  endif()
endif()

# ufs
if(FV3_FORECAST_MODEL MATCHES "UFS")
  # Should be: find_package(fms REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake )
  find_package(fms REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake )
  find_package(fv3atm REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake )
  find_package(ufs-weather-model-lib REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake )
  link_libraries(fms)
  link_libraries(NetCDF::NetCDF_Fortran)
  # Install includes
  include_directories( ${FV3_FORECAST_MODEL_ROOT}/mod )
  # Source incudes
  include_directories( ${FV3_FORECAST_MODEL_SRC}/FMS/include )
  include_directories( ${FV3_FORECAST_MODEL_SRC}/FMS/mpp/include )
  # Build includes
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FMS/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/cpl/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/gfsphysics/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/io/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/ipd/mod )
  include_directories( ${FV3_FORECAST_MODEL_BUILD}/FV3/stochastic_physics )
endif()

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDILM_TPLS MPI )

set( FV3JEDILM_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                            ${CMAKE_CURRENT_BINARY_DIR}/src
                            ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( FV3JEDILM_LIBRARIES fv3jedilm )

get_directory_property( FV3JEDILM_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${FV3JEDILM_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND FV3JEDILM_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND FV3JEDILM_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND FV3JEDILM_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include( fv3jedilm_compiler_flags )
include_directories( ${FV3JEDILM_INCLUDE_DIRS} ${FV3JEDILM_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME fv3jedilm )

# print the summary of the configuration
ecbuild_print_summary()
