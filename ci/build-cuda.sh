#!/bin/bash

set -xeuo pipefail


# make sure we keep the stage direcorty
spack config --scope=user add config:build_stage:$PWD


spack env create -d ./spack-env
# add local repository with current sirius recipe
spack -e ./spack-env repo add ./ci/sirius-uenv-recipe/repo

spack -e ./spack-env add $SPEC #cuda_arch==$CUDA_ARCH

cat ./spack-env/spack.yaml

# build sirius from source
spack -e ./spack-env develop -p . sirius@develop
spack -e ./spack-env concretize
spack -e ./spack-env install
builddir=$(spack -e ./spack-env location -b sirius)
# create a symlink to spack build directory (keep in artifacts)
ln -s $builddir builddir
