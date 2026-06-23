## Compilation de la branche experimental sur math10

```
module load compiler-rt/2024.0.2 ifort/2024.0.2 mpi/latest tbb/latest mkl/latest vtune/latest itac/latest ScaleMP 

export ProjectDir=$PWD
export DepDir=$ProjectDir/ThirdParty  # choix non obligatoire pour DepDir !!
###
cd $ProjectDir/ThirdParty
mkdir hdf5-1.14.5-intel; cd hdf5-1.14.5-intel
tar xzf ../hdf5-1.14.5.tar.gz --strip=2
mkdir build;cd build
CC=icx FC=ifx cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/hdf5 -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_FORTRAN=ON
make -j8; make install  # 11min31s sur math10
###
cd $ProjectDir/ThirdParty
mkdir fftw-3.3.10-intel; cd fftw-3.3.10-intel
tar xzf ../fftw-3.3.10.tar.gz --strip=1
mkdir build;cd build
CC=icx FC=ifx cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/fftw3
make -j8; make install  # 1min10s sur math10
###
cd $ProjectDir/ThirdParty
mkdir 2decomp-fft-2.0.3-intel
cd 2decomp-fft-2.0.3-intel
tar xzf ../2decomp-fft-2.0.3.tar.gz --strip=1
mkdir build;cd build
CC=icx FC=ifx cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/decomp2d
make -j2;make install
###
cd $ProjectDir/slip_error
mkdir build;cd build
CC=icx FC=ifx cmake .. -DDECOMP2D_ROOT=$DepDir/lib/intel/decomp2d -DFFTW3_ROOT=$DepDir/lib/intel/fftw3 -DHDF5_ROOT=$DepDir/lib/intel/hdf5 
-DCMAKE_BUILD_TYPE=Release
make
####
rm -fr $ProjectDir/ThirdParty/hdf5-1.14.5-intel $ProjectDir/ThirdParty/2decomp-fft-2.0.3-intel $ProjectDir/ThirdParty/fftw-3.3.10-intel  #
 purge des install
```


## Compilation de la branche experimental sur AMU
```
module purge
module load userspace/all
module load intel-compiler/64/2020.2.254 openmpi/icc20/psm2/4.0.5 cmake/3.30.5
export ProjectDir=$PWD
export DepDir=$ProjectDir/ThirdParty  # choix non obligatoire pour DepDir !!

###
cd $ProjectDir/ThirdParty
mkdir hdf5-1.14.5-intel; cd hdf5-1.14.5-intel
tar xzf ../hdf5-1.14.5.tar.gz --strip=2
mkdir build;cd build
CC=icc FC=ifort cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/hdf5 -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_FORTRAN=ON
make -j8; make install
###
cd $ProjectDir/ThirdParty
mkdir fftw-3.3.10-intel; cd fftw-3.3.10-intel
tar xzf ../fftw-3.3.10.tar.gz --strip=1
mkdir build;cd build
CC=icc FC=ifort cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/fftw3
make -j8; make install
###
cd $ProjectDir/ThirdParty
mkdir 2decomp-fft-2.0.3-intel
cd 2decomp-fft-2.0.3-intel
tar xzf ../2decomp-fft-2.0.3.tar.gz --strip=1
mkdir build;cd build
CC=icc FC=ifort cmake .. -DCMAKE_INSTALL_PREFIX=$DepDir/lib/intel/decomp2d
make -j2;make install
###
cd $ProjectDir/slip_error
mkdir build;cd build
CC=icc FC=ifort cmake .. -DDECOMP2D_ROOT=$DepDir/lib/intel/decomp2d -DFFTW3_ROOT=$DepDir/lib/intel/fftw3 -DHDF5_ROOT=$DepDir/lib/intel/hdf
5 -DCMAKE_BUILD_TYPE=Release
make
###
rm -fr $ProjectDir/ThirdParty/hdf5-1.14.5-intel $ProjectDir/ThirdParty/2decomp-fft-2.0.3-intel $ProjectDir/ThirdParty/fftw-3.3.10-intel  #
 purge des install
###
srun -p skylake --ntasks-per-node=8  --pty bash -i
mpirun -np 8 ./rrb_rk3.x --input ../run-rrb-simple.in 
mpirun -np 8 ./rrb_hhi.x --input ../run-rrb-simple.in 
exit
```
