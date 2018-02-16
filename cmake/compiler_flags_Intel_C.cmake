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
  set( CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -qopenmp")
else( )
  set( CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -qopenmp-stubs")
endif( )

####################################################################
# FROM GMAO FLAGS
####################################################################

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOVERLOAD_R4 -DMAPL_MODE -DEIGHT_BYTE -DSPMD -DTIMING -Duse_libMPI -Duse_netCDF")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_C_FLAGS_RELEASE     "-O3" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_C_FLAGS_DEBUG       "-O0 -g -traceback -fp-trap=common" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_C_FLAGS_BIT         "-O2" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_C_LINK_FLAGS        "" )

####################################################################

# Meaning of flags
# ----------------
# todo
  
