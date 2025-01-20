#!/bin/bash
set -e
set -u

_spack_version=${SPACK_VERSION:-develop}

_spack_folder=$1

if [ ! -d "${_spack_folder}" ]; then
    echo "Cloning spack"
    git clone --branch ${_spack_version} -c feature.manyFiles=true https://github.com/acts-project/spack.git ${_spack_folder}
    pushd ${_spack_folder} > /dev/null
    git config user.name 'CI'
    git config user.email '<>'
    popd > /dev/null
fi
