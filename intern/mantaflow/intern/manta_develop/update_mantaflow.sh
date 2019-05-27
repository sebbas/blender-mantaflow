#!/bin/bash
#
# ========================================================================================
# UPDATING MANTAFLOW INSIDE BLENDER
# ========================================================================================

# ====================  1) ENVIRONMENT SETUP =============================================

# YOUR INSTALLATION PATHS GO HERE:
MANTA_INSTALLATION=/Users/sebbas/Developer/Mantaflow/mantaflowDevelop
BLENDER_INSTALLATION=/Users/sebbas/Developer/Blender/fluid-mantaflow

# ==================== 2) BUILD MANTAFLOW (OPENMP AND TBB) ===============================

# Need non-default (OpenMP enabled) compiler to build Mantaflow on OSX
if [[ "$OSTYPE" == "darwin"* ]]; then
  export CC=/usr/local/Cellar/gcc/9.1.0/bin/gcc-9
  export CXX=/usr/local/Cellar/gcc/9.1.0/bin/g++-9
fi

cd $MANTA_INSTALLATION
if cd mantaflowgit/; then git pull; else git clone git@bitbucket.org:thunil/mantaflowgit.git; cd mantaflowgit; fi
git checkout develop
MANTA_OMP_PATH=$MANTA_INSTALLATION/mantaflowgit/build_omp/
MANTA_TBB_PATH=$MANTA_INSTALLATION/mantaflowgit/build_tbb/
mkdir -p $MANTA_OMP_PATH $MANTA_TBB_PATH
cd $MANTA_OMP_PATH
cmake .. -DGUI=OFF -DOPENMP=ON -DBLENDER=ON -DPREPDEBUG=ON && make -j4
cd $MANTA_TBB_PATH
cmake .. -DGUI=OFF -DTBB=ON -DBLENDER=ON -DPREPDEBUG=ON && make -j4

# ==================== 3) COPY MANTAFLOW FILES TO BLENDER DIRECTORIES ====================

BLENDER_DEPENDENCIES_PATH=$BLENDER_INSTALLATION/blender/intern/mantaflow/intern/manta_develop/dependencies/
mkdir -p $BLENDER_DEPENDENCIES_PATH && cp -Rf $MANTA_INSTALLATION/mantaflowgit/dependencies/cnpy "$_"
echo "Copied Mantaflow dependencies to" $BLENDER_DEPENDENCIES_PATH

BLENDER_HELPER_PATH=$BLENDER_INSTALLATION/blender/intern/mantaflow/intern/manta_develop/helper/
mkdir -p $BLENDER_HELPER_PATH && cp -Rf $MANTA_INSTALLATION/mantaflowgit/source/util "$_"
mkdir -p $BLENDER_HELPER_PATH && cp -Rf $MANTA_INSTALLATION/mantaflowgit/source/pwrapper "$_"
echo "Copied Mantaflow helper files to" $BLENDER_HELPER_PATH

BLENDER_OMP_PATH=$BLENDER_INSTALLATION/blender/intern/mantaflow/intern/manta_develop/preprocessed/omp/
mkdir -p $BLENDER_OMP_PATH && cp -Rf $MANTA_INSTALLATION/mantaflowgit/build_omp/pp/source/. "$_"
echo "Copied Mantaflow OpenMP preprocessed files to" $BLENDER_OMP_PATH

BLENDER_TBB_PATH=$BLENDER_INSTALLATION/blender/intern/mantaflow/intern/manta_develop/preprocessed/tbb/
mkdir -p $BLENDER_TBB_PATH && cp -Rf $MANTA_INSTALLATION/mantaflowgit/build_tbb/pp/source/. "$_"
echo "Copied Mantaflow TBB preprocessed files to" $BLENDER_TBB_PATH

# ==================== 4) CHECK CMAKE SETUP ==============================================

# Make sure that all files copied from Mantaflow are listed in intern/mantaflow/CMakeLists.txt
# Especially if new source files / plugins were added to Mantaflow.
