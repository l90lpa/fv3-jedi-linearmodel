#!/bin/csh -f

#Script for compiling fv3 on the NASA Center for Climate Similation (NCCS) Discover cluster

#Clear path
set path = /bin:/usr/bin

#Shell
source /usr/share/modules/init/csh

#Jedi install config
setenv JEDI_ROOT /discover/nobackup/drholdaw/Jedi/
setenv JEDI_SRC $JEDI_ROOT/jedi-bundle/

module purge
module load other/cmake-3.8.2

if ($1 == "INT" || $1 == "Int" || $1 == "Intel"  || $1 == "intel") then

   module load other/comp/gcc-6.2
   module load comp/intel-18.0.1.163
   module load mpi/impi-18.0.1.163
   module load lib/mkl-18.0.1.163
   setenv BASEDIR /discover/swdev/mathomp4/Baselibs/ESMA-Baselibs-5.0.9/x86_64-unknown-linux-gnu/ifort_18.0.1.163-intelmpi_18.0.1.163/Linux/ 

   setenv MPIEXEC `which mpirun`
   setenv CPCcomp mpiicpc
   setenv CCcomp mpiicc
   setenv F90comp mpiifort
   
   setenv JEDI_BUILD $JEDI_ROOT/build_int/
   setenv GFDL_BUILD $JEDI_ROOT/build_gfdl_int/
  
else if ($1 == "GCC" || $1 == "gcc" || $1 == "GNU" || $1 == "gnu") then
   
   #Modules GCC
   module purge
   module load other/cmake-3.8.2
   module load other/comp/gcc-7.2
   module load other/mpi/openmpi/3.0.0-gcc-7.2 
   setenv BASEDIR /discover/swdev/mathomp4/Baselibs/ESMA-Baselibs-4.0.10/x86_64-unknown-linux-gnu/gfortran_7.2.0-openmpi_3.0.0/Linux

   setenv MPIEXEC `which mpirun`
   setenv CPCcomp mpicxx
   setenv CCcomp mpicc
   setenv F90comp mpifort

   setenv JEDI_BUILD $JEDI_ROOT/build_gcc/
   setenv GFDL_BUILD $JEDI_ROOT/build_gfdl_gcc/

else
   
   echo "No compiler specified"
   exit()

endif

#Check Jedi has been built first
if (! -d $JEDI_BUILD) then
   echo "Build Jedi first for NetCdf directories"
   exit() 
endif

#Source to be built
setenv SRC $JEDI_ROOT/fv3

#NetCDF lib search path
setenv NETCDF $JEDI_BUILD/netcdf
setenv NETCDF_INCLUDE_DIRS $JEDI_BUILD/netcdf/include #Need in order not to auto redefine NETCDF_LIBRARIES
setenv NETCDF_LIBRARIES "$JEDI_BUILD/netcdf/lib/libnetcdf.a;$JEDI_BUILD/netcdf/lib/libnetcdff.a;${BASEDIR}/lib/libhdf5_hl.a;${BASEDIR}/lib/libhdf5.a;${BASEDIR}/lib/libcurl.a;/usr/lib64/libcrypto.so;/usr/lib64/libssl.so;${BASEDIR}/lib/libmfhdf.a;${BASEDIR}/lib/libdf.a;${BASEDIR}/lib/libjpeg.a"

#FMS Path
setenv FMS_PATH "$GFDL_BUILD/fms/"

#Set defs from GMAO builds
setenv COMPDEFS "-DsysLinux;-DESMA64;-DHAS_NETCDF4;-DHAS_NETCDF3;-DH5_HAVE_PARALLEL;-DNETCDF_NEED_NF_MPIIO;-DMAPL_MODE;-DSPMD;-DTIMING;-DOLDMPP;-DHAVE_SHMEM"

#Add ecbuild to path
set path = (${path} ${JEDI_ROOT}/ecbuild/bin)

#If Build nonexistent the create
cd ${JEDI_ROOT}
if (! -d $GFDL_BUILD) then
  echo "Must build fms first"
  exit()
endif

cd $GFDL_BUILD

if (! -d fms/) then
  echo "Must build fms first"
  exit()
endif

if ($2 == "clean" || ! -d fv3/) then

   #Prepare build location
   rm -rf fv3/
   mkdir -p fv3
   cd fv3
   
   #Prepare to make
   ecbuild \
       --build=debug \
       -DCMAKE_CXX_COMPILER=$CPCcomp \
       -DCMAKE_C_COMPILER=$CCcomp \
       -DCMAKE_Fortran_COMPILER=$F90comp \
       -DNETCDF_INCLUDE_DIRS=$NETCDF_INCLUDE_DIRS \
       -DNETCDF_LIBRARIES=$NETCDF_LIBRARIES \
       -DNETCDF_PATH=$NETCDF \
       -DMPIEXEC=$MPIEXEC \
       -DCOMPDEFS=$COMPDEFS \
       -DFMSPATH=$FMS_PATH \
       $SRC

endif

#Build the model
cd $GFDL_BUILD/fv3
make VERBOSE=YES -j1
