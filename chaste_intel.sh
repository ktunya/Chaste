# Script to build Chaste executable with Intel build
#
# Before running this script you should create the Intel build with libraries (scons build=Intel chaste_libs=1)
#
# run from /home/eclipse/workspace/Chaste
#
mpicxx -CC=icpc -isystem /home/chaste/petsc-2.3.2-p4/bmake/linux-intel-opt-mkl -isystem /home/chaste/petsc-2.3.2-p4/include -isystem ../../../xsd-2.3.1-i686-linux-gnu/libxsd -Werror -DNDEBUG -O3 -I. -Icxxtest -Iglobal/src -Icancer/src -Icancer/src/common -Icancer/src/odes -Icancer/src/mesh -Icancer/src/tissue -Icancer/src/tissue/killers -Icancer/src/tissue/cell -Icancer/src/tissue/cell/cycle -Idealii/src -Idealii/src/common -Idealii/src/problem -Idealii/src/solver -Idealii/src/problem/cancer -Idealii/src/problem/elasticity -Idealii/src/problem/cardiac -Idealii/src/solver/common -Idealii/src/solver/elasticity -Ilinalg/src -Ilinalg/src/common -Iheart/src -Iheart/src/pdes -Iheart/src/odes -Iheart/src/stimulus -Iheart/src/convergence -Iheart/src/problem -Iheart/src/io -Iheart/src/solver -Ipde/src -Ipde/src/common -Ipde/src/problem -Ipde/src/solver -Ipde/src/problem/common -Ipde/src/solver/common -Imesh/src -Imesh/src/common -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/reader -Iode/src -Iode/src/common -Iode/src/problem -Iode/src/solver -Iio/src -Iio/src/writer -Iio/src/reader -c -o Chaste.o Chaste.cpp

mpicxx -CC=icpc -static-libcxa -o Chaste Chaste.o -Llinklib -Lheart/build/intel -Lheart -L/home/chaste/petsc-2.3.2-p4/lib/linux-intel-opt-mkl -L/opt/intel/cc/9.1.039/lib -L/home/chaste/petsc-2.3.2-p4/externalpackages/f2cblaslapack/linux-gnu -L/opt/intel/mkl/9.1.023/lib/32 -Llib -ltestheart -lheart -lode -lmesh -llinalg -lio -lglobal -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmkl_lapack -lmkl -lboost_serialization -lxerces-c
