# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( HAVE_OMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
else( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp-stubs")
endif( )

####################################################################
# FROM GMAO FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSPMD -DTIMING -DOLDMPP")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopt-report0 -ftz -align all -fno-alias -qno-offload -traceback -r8 -stack_temps -safe_cray_ptr -assume byterecl -fp-model source -ftz -align all -fno-alias -align dcommons -convert big_endian  -fPIC -fpe0 -fp-model source -heap-arrays 32 -assume noold_maxminloc  -align dcommons" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -ftz -align all -fno-alias -traceback -r8 -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -warn unused -init=snan,arrays" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -no-heap-arrays" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# todo
  
