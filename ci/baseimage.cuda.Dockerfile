FROM ubuntu:24.04 as builder

ARG CUDA_ARCH=90

ENV DEBIAN_FRONTEND noninteractive

ENV FORCE_UNSAFE_CONFIGURE 1

ENV PATH="/spack/bin:${PATH}"

ENV MPICH_VERSION=3.4.3

ENV CMAKE_VERSION=3.31.2

RUN apt-get -y update

RUN apt-get install -y apt-utils

# install basic tools
RUN apt-get install -y --no-install-recommends gcc g++ gfortran git make unzip file \
  vim wget pkg-config python3-pip python3-dev cython3 python3-pythran curl tcl m4 cpio automake meson \
  xz-utils patch patchelf apt-transport-https ca-certificates gnupg software-properties-common perl tar bzip2 \
  liblzma-dev libbz2-dev

# install CMake
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-aarch64.tar.gz -O cmake.tar.gz && \
    tar zxvf cmake.tar.gz --strip-components=1 -C /usr

# get latest version of spack
RUN git clone -b releases/v0.23 https://github.com/spack/spack.git

# add local repo to spack
COPY ./spack /opt/spack
RUN spack repo add --scope system /opt/spack

# set the location of packages built by spack
RUN spack config add config:install_tree:root:/opt/local
# set cuda_arch for all packages
RUN spack config add packages:all:variants:cuda_arch=${CUDA_ARCH}

# find all external packages
RUN spack external find --all --exclude=openssl

# find compilers
RUN spack compiler find

# install yq (utility to manipulate the yaml files)
RUN wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_arm64 && chmod a+x /usr/local/bin/yq


# # change the fortran compilers: for gcc the gfortran is already properly set and the change has no effect; add it for clang
# RUN yq -i '.compilers[0].compiler.paths.f77 = "/usr/bin/gfortran"' /root/.spack/linux/compilers.yaml && \
#     yq -i '.compilers[0].compiler.paths.fc = "/usr/bin/gfortran"' /root/.spack/linux/compilers.yaml  && \
#     yq -i '.compilers[1].compiler.paths.f77 = "/usr/bin/gfortran"' /root/.spack/linux/compilers.yaml && \
#     yq -i '.compilers[1].compiler.paths.fc = "/usr/bin/gfortran"' /root/.spack/linux/compilers.yaml

# install MPICH
RUN spack install mpich@${MPICH_VERSION} %gcc

# install libvdwxc
RUN spack install libvdwxc %gcc +mpi ^mpich@${MPICH_VERSION}

# install openblas
RUN spack install openblas %gcc +fortran

RUN spack install magma %gcc +cuda +fortran ^openblas

RUN spack install nlcglib@master %gcc +cuda

# for the MPI hook
RUN echo $(spack find --format='{prefix.lib}' mpich) > /etc/ld.so.conf.d/mpich.conf
RUN ldconfig

# # create environments for several configurations and install dependencies
# RUN spack env create -d /sirius-env-clang && \
#     spack -e /sirius-env-clang add "sirius@develop %clang build_type=RelWithDebInfo ~cuda ~fortran +tests +pugixml ^openblas%gcc ^libxc%gcc ^mpich%gcc " && \
#     spack -e /sirius-env-clang develop -p /sirius-src sirius@develop && \
#     spack -e /sirius-env-clang install --only=dependencies --fail-fast

RUN spack env create -d /sirius-env-cuda && \
    spack -e /sirius-env-cuda add "sirius@develop %gcc build_type=RelWithDebInfo +scalapack +tests +pugixml +apps +cuda +magma +python ^netlib-scalapack ^mpich ^openblas threads=openmp" && \
    spack -e /sirius-env-cuda develop -p /sirius-src sirius@develop && \
    spack -e /sirius-env-cuda install --only=dependencies --fail-fast

RUN spack env create -d /sirius-env-elpa && \
    spack -e /sirius-env-elpa add "sirius@develop %gcc build_type=RelWithDebInfo +tests +pugixml +apps +cuda +scalapack +elpa ^netlib-scalapack ^mpich ^openblas ^elpa+cuda" && \
    spack -e /sirius-env-elpa develop -p /sirius-src sirius@develop && \
    spack -e /sirius-env-elpa install --only=dependencies --fail-fast

RUN spack env create -d /sirius-env-fp32 && \
    spack -e /sirius-env-fp32 add "sirius@develop %gcc build_type=RelWithDebInfo +tests +pugixml +apps +cuda ^mpich ^openblas ^elpa+cuda ^spfft+single_precision+cuda" && \
    spack -e /sirius-env-fp32 develop -p /sirius-src sirius@develop && \
    spack -e /sirius-env-fp32 install --only=dependencies --fail-fast

RUN spack env create -d /sirius-env-nlcg && \
    spack -e /sirius-env-nlcg add "sirius@develop %gcc build_type=RelWithDebInfo +fortran +tests +pugixml +apps +cuda +nlcglib ^openblas ^mpich" && \
    spack -e /sirius-env-nlcg develop -p /sirius-src sirius@develop && \
    spack -e /sirius-env-nlcg install --only=dependencies --fail-fast

# RUN spack env create -d /sirius-env-openmpi && \
#     spack -e /sirius-env-openmpi add "sirius@develop %gcc +tests +pugixml +apps +scalapack +fortran build_type=RelWithDebInfo ^netlib-scalapack ^openblas ^openmpi" && \
#     spack -e /sirius-env-openmpi develop -p /sirius-src sirius@develop && \
#     spack -e /sirius-env-openmpi install --only=dependencies --fail-fast

# RUN spack env create -d /sirius-env-cuda-sequential && \
#     spack -e /sirius-env-cuda-sequential add "sirius@develop %gcc +cuda +tests +pugixml +apps +fortran build_type=RelWithDebInfo ^openblas ^openmpi" && \
#     spack -e /sirius-env-cuda-sequential develop -p /sirius-src sirius@develop && \
#     spack -e /sirius-env-cuda-sequential install --only=dependencies --fail-fast

RUN spack env create -d /sirius-env-vdwxc-cuda && \
    spack -e /sirius-env-vdwxc-cuda add "sirius@develop %gcc build_type=RelWithDebInfo +fortran +tests +pugixml +apps +vdwxc +cuda +nlcglib ^openblas ^mpich +cuda" && \
    spack -e /sirius-env-vdwxc-cuda develop -p /sirius-src sirius@develop && \
    spack -e /sirius-env-vdwxc-cuda install --only=dependencies --fail-fast
