#!/bin/sh --login

# Get JEDI_ROOT from user input, defaults to level above fv3
read -p "Enter path to JEDI_ROOT: " JEDI_ROOT
JEDI_ROOT=${JEDI_ROOT:-$(dirname $(dirname $PWD))}

# Source / build environemnt
JEDI_SRC_FV3="$JEDI_ROOT/fv3"
JEDI_BUILD_FV3="$JEDI_SRC_FV3/build-fv3"
JEDI_BUILD_FMS="$JEDI_ROOT/fms/build-fms"

###########################################################

# Clean modules
source $MODULESHOME/init/sh

module purge

# Load Intel compilers, NetCDF and HDF5 libraries
module load intel/18.1.163
module load impi/5.1.2.150
module load hdf5/1.8.14
module load netcdf/4.3.0

# Load cmake
module use -a /contrib/modulefiles
module load cmake

module list

# Include ecbuild in path
export PATH="$PATH:$JEDI_ROOT/ecbuild/bin"

# Need correct MPIEXEC on Theia
MPIEXEC=$(which mpirun)

# Compiler definitions for FMS and FV3
#COMPDEFS="-DOLDMPP"
COMPDEFS=""

# Check built FMS exists
if [ ! -d $JEDI_BUILD_FMS ]; then
    echo "FMS must exist prior to building FV3 under:"
    echo "$JEDI_BUILD_FMS"
    echo "ABORT!"
    exit 99
fi

# Build FV3
rm -rf $JEDI_BUILD_FV3 ; mkdir -p $JEDI_BUILD_FV3 ; cd $JEDI_BUILD_FV3

ecbuild \
    --build=release \
    -DCMAKE_CXX_COMPILER=mpiicpc \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON \
    -DMPIEXEC=$MPIEXEC \
    -DNETCDF_ROOT=$NETCDF \
    -DCOMPDEFS=$COMPDEFS \
    -DFMS_PATH=$JEDI_BUILD_FMS \
    $JEDI_SRC_FV3

make VERBOSE=YES -j12
rc=$?
[[ $rc -ne 0 ]] && exit $rc
