# (C) Copyright 2018-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

####################################################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################################################

# OpenMP
if( HAVE_OMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
else( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp-stubs")
endif( )

# Dynamical core precision
if (FV3_PRECISION MATCHES DOUBLE OR NOT FV3_PRECISION)
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
endif()

# Traceback
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")

####################################################################################################
# RELEASE FLAGS
####################################################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -unroll -inline -no-heap-arrays" )

####################################################################################################
# DEBUG FLAGS
####################################################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check bounds -warn -heap-arrays -fpe-all=0 -fpe:0")

####################################################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -no-heap-arrays" )

####################################################################################################
# LINK FLAGS
####################################################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################################################
