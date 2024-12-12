#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err
# install ABACUS with NTPoly support

# Build ABACUS by AOCC

ABACUS_DIR=$(pwd)/..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
#source $INSTALL_DIR/setup
cd $ABACUS_DIR

BUILD_DIR=build_abacus_aocc_ntpoly
rm -rf $BUILD_DIR/*

AOCL_DIR=/opt/AMD/aocl-linux-aocc-5.0.0/5.0.0/aocc
PREFIX=$BUILD_DIR
BLAS_LIBRARIES=$AOCL_DIR/lib/libblis-mt.a
LAPACK=$AOCL_DIR/lib/libflame.a
SCALAPACK_DIR=$AOCL_DIR/lib
ELPA=/opt/elpa/2024.05.001/openmpi-aocc-aocl
FFTW3=$AOCL_DIR/lib
CEREAL=$HOME/project/cereal/include/cereal
NTPoly_DIR=/opt/NTPoly/3.1.1/openmpi-aocc-aocl

#RAPIDJSON=$HOME/project/rapidjson/include/rapidjson
#LIBXC=$INSTALL_DIR/libxc-6.2.2
# LIBRI=$INSTALL_DIR/LibRI-0.2.1.0
# LIBCOMM=$INSTALL_DIR/LibComm-0.1.1
# LIBTORCH=$INSTALL_DIR/libtorch-2.1.2/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-1.0.1/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DMPI_CXX_COMPILER=mpicxx \
        -DLAPACK_DIR=$LAPACK \
        -DSCALAPACK_DIR=$SCALAPACK_DIR \
        -DELPA_DIR=$ELPA \
        -DFFTW3_DIR=$FFTW3 \
	-DFFTW3_INCLUDE_DIR=$AOCL_DIR/include \
        -DCEREAL_INCLUDE_DIR=$CEREAL \
        -DENABLE_LCAO=ON \
        -DUSE_OPENMP=ON \
        -DENABLE_NTPOLY=ON \
	-DNTPoly_DIR=$NTPoly_DIR \
        -DMPI_FORTRAN_LIBRARIES="-lmpi_mpifh" \
#        -DENABLE_RAPIDJSON=ON \
#        -DRapdidJSON_DIR=$RAPIDJSON \
#         -DENABLE_DEEPKS=1 \
#         -DTorch_DIR=$LIBTORCH \
#         -Dlibnpy_INCLUDE_DIR=$LIBNPY \
#         -DENABLE_LIBRI=ON \
#         -DLIBRI_DIR=$LIBRI \
#         -DLIBCOMM_DIR=$LIBCOMM \
# 	      -DDeePMD_DIR=$DEEPMD \
# 	      -DTensorFlow_DIR=$DEEPMD \
        -DCMAKE_VERBOSE_MAKEFILE=ON 

# # add mkl env for libtorch to link
# if one want to install libtorch, mkl should be load in build process
# for -lmkl when load libtorch
# module load mkl

# if one want's to include deepmd, your system gcc version should be >= 11.3.0 for glibc requirements

#cmake --build $BUILD_DIR -j `nproc` 
#cmake --install $BUILD_DIR 2>/dev/null
cd $BUILD_DIR && make VERBOSE=1 -j 16 && make install
# generate abacus_env.sh
cat << EOF > "${TOOL}/abacus_env.sh"
#!/bin/bash
source $INSTALL_DIR/setup
export PATH="${PREFIX}/bin":\${PATH}
EOF

# generate information
cat << EOF
========================== usage =========================
Done!
To use the installed ABACUS version
You need to source ${TOOL}/abacus_env.sh first !
"""
EOF
