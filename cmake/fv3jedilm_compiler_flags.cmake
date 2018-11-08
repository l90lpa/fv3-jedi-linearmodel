# (C) Copyright 2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Must have compiler definitions for FV3JEDI Linear Model
add_definitions( -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE -DOLDMPP )

if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

#######################################################################################
# Fortran
#######################################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_Fortran )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()

#######################################################################################
# C
#######################################################################################

if( CMAKE_C_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_C )
elseif( CMAKE_C_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_C )
elseif( CMAKE_C_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_C )
elseif( CMAKE_C_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_C )
elseif( CMAKE_C_COMPILER_ID MATCHES "Clang" )
  include( compiler_flags_Clang_C )
else()
  message( STATUS "C compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()
