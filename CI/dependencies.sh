#!/bin/bash
set -e
set -u

function run() {
    set -x
    "$@"
    { set +x;  } 2> /dev/null
}



url=${1:-${DEPENDENCY_URL:-}}

if [ -n "${GITHUB_ACTIONS:-}" ]; then
    destination="${GITHUB_WORKSPACE}/dependencies"
    echo "DEPENDENCY_DIR=${destination}" >> $GITHUB_ENV
else
    destination=${2}
fi


if [ -z "${url}" ]; then
    echo "url is not set"
    exit 1
fi

echo "URL: $url"
echo "DESTINATION: $destination"

# check curl location
CURL=$(command -v curl)
if [ -z "$CURL" ]; then
    echo "curl is not available"
    exit 1
fi

UNZSTD=$(command -v unzstd)
if [ -z "$UNZSTD" ]; then
    echo "unzstd is not available"
    exit 1
fi

TAR=$(command -v tar)
if [ -z "$TAR" ]; then
    echo "tar is not available"
    exit 1
fi

run mkdir -p "${destination}"

run $CURL \
  --retry 5 \
  --connect-timeout 2 \
  --location $url \
  | unzstd \
  | tar \
    -x \
    --strip-components=1 \
    --directory "${destination}"

# Patch up geant4-config data install script
orig_share=$(dependencies/bin/geant4-config --datasets|head -n1|perl -pe 's|.*?(\/.*)\/share.*|\1|')
orig_share_escaped=$(echo $orig_share|perl -pe 's|/|\\/|g')
destination_escaped=$(echo "$destination"|perl -pe 's|/|\\/|g')
perl -pi.bak -e "s/$orig_share_escaped/$destination_escaped/g" dependencies/bin/geant4-config

if [ -n "${GITHUB_ACTIONS:-}" ]; then
    venv="${GITHUB_WORKSPACE}/venv"
    run "${destination}/bin/python3" -m venv "${venv}"
    run "${venv}/bin/python3" -m pip install pyyaml jinja2
    echo "PATH=${venv}/bin:${destination}/bin/:${PATH}" >> $GITHUB_ENV
    echo "CMAKE_PREFIX_PATH=${destination}" >> $GITHUB_ENV
    echo "LD_LIBRARY_PATH=${destination}/lib" >> $GITHUB_ENV
    echo "ROOT_INCLUDE_PATH=${destination}/include" >> $GITHUB_ENV
    # Geant4 puts CLHEP in a subdirectory
    echo "ROOT_INCLUDE_PATH=${destination}/include/Geant4" >> $GITHUB_ENV
    # Pythia8 looks for settings in this directory
    echo "PYTHIA8DATA=${destination}/share/Pythia8/xmldoc" >> $GITHUB_ENV
fi