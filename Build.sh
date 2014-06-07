#!/bin/bash
set -e

# This is not a makefile; it's just a simple script to build all programs.
# All of the classes are header-only, so no libraries are generated.

HDF5_INCDIR=/nfs/slac/g/exo_data4/users/cgd8d/HDF5/hdf5-1.8.13/hdf5/include
HDF5_LIBDIR=/nfs/slac/g/exo_data4/users/cgd8d/HDF5/hdf5-1.8.13/hdf5/lib

# I want to force lots of libraries to be linked statically; make that easier.
function static_link() {
echo "-Wl,-Bstatic -l$1 -Wl,-Bdynamic"
}

# Building noise correlation generation program.
g++ -O3 -o ComputeNoiseCorrelations -I. -I$HDF5_INCDIR `exo-config --cflags` -L`exo-config --libdir` -L$HDF5_LIBDIR ComputeNoiseCorrelations.C `static_link hdf5` -lEXOUtilities

# Building lightmap construction program.
gcc -O3 -o External/sqlite3.o -DSQLITE_THREADSAFE=0 -c External/sqlite3.c
g++ -O3 -o BuildLightMap -I. -I$HDF5_INCDIR -L$HDF5_LIBDIR BuildLightMap.C External/sqlite3.o `static_link hdf5` -lz -ldl
