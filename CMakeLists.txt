# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDILM
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedilm C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH} )

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
ecbuild_add_cxx11_flags()

set( FV3JEDILM_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################

# LAPACK
ecbuild_add_option( FEATURE LAPACK
                    DESCRIPTION "Linear Algebra PACKage"
                    REQUIRED_PACKAGES "LAPACK REQUIRED" )
if (NOT HAVE_LAPACK)
  ecbuild_error("LAPACK_LIBRARIES not defined")
endif()

# MPI
ecbuild_add_option( FEATURE MPI DEFAULT ON
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS CXX Fortran" )
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# NetCDF
set( NETCDF_F90 ON CACHE BOOL "Compile with Fortran NetCDF" )
find_package( NetCDF REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIR} )

# fms
include_directories( ${FMS_INCLUDE_DIRS} )

################################################################################
# Definitions
################################################################################

# Additional compiler definitions from user input
set( COMPDEFS_ "${COMPDEFS}" )
add_definitions( ${COMPDEFS_} )

message("Compiler definitions for fv3-jedi-lm are: ${COMPDEFS_}")

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDILM_TPLS LAPACK MPI fms )

set( FV3JEDILM_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
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
