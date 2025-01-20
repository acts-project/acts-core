#!/bin/bash
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  # Only set these properties if it's called as a subprocess
  set -e
  set -u
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export SPACK_COLOR=always

function set_env {
  key="$1"
  value="$2"

  echo "=> ${key}=${value}"

  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    echo "${key}=${value}" >> "$GITHUB_ENV"
  else
    export "${key}"="${value}"
  fi
}

function start_section() {
    local section_name="$1"
    if [ -n "${GITHUB_ACTIONS:-}" ]; then
        echo "::group::${section_name}"
    else
        echo "${section_name}"
    fi
}

function end_section() {
    if [ -n "${GITHUB_ACTIONS:-}" ]; then
        echo "::endgroup::"
    fi
}

# Parse command line arguments
while getopts "c:t:d:h" opt; do
  case ${opt} in
    c )
      compiler=$OPTARG
      ;;
    t )
      tag=$OPTARG
      ;;
    d )
      destination=$OPTARG
      ;;
    h )
      echo "Usage: $0 [-c compiler] [-t tag] [-d destination]"
      echo "Options:"
      echo "  -c <compiler>    Specify compiler (defaults to CXX env var)"
      echo "  -t <tag>         Specify dependency tag (defaults to DEPENDENCY_TAG env var)"
      echo "  -d <destination> Specify install destination (defaults based on CI environment)"
      echo "  -h              Show this help message"
      exit 0
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done

# Set defaults if not specified
if [ -z "${compiler:-}" ]; then
  compiler="${CXX:-default}"
fi

if [ -z "${tag:-}" ]; then
  tag="${DEPENDENCY_TAG:-}"
  if [ -z "${tag:-}" ]; then
    echo "No tag specified via -t or DEPENDENCY_TAG environment variable"
    exit 1
  fi
fi

if [ -z "${destination:-}" ]; then
  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    destination="${GITHUB_WORKSPACE}/dependencies"
  elif [ -n "${GITLAB_CI:-}" ]; then
    destination="${CI_PROJECT_DIR}/dependencies"
  else
    echo "No destination specified via -d and not running in CI"
    exit 1
  fi
fi



echo "Install tag: $tag"
echo "Install destination: $destination"

mkdir -p ${destination}

if [ -n "${GITLAB_CI:-}" ]; then
    _spack_folder=${CI_PROJECT_DIR}/spack
else
    _spack_folder=${PWD}/spack
fi

start_section "Install spack if not already installed"
if ! command -v spack &> /dev/null; then
  "${SCRIPT_DIR}/setup_spack.sh" "${_spack_folder}"
  source "${_spack_folder}/share/spack/setup-env.sh"
fi
end_section

if [ -n "${GITLAB_CI:-}" ]; then
  # Use the project spack config for GitLab CI so we can cache it
  mkdir -p ${CI_PROJECT_DIR}/.spack
  ln -s ${CI_PROJECT_DIR}/.spack ${HOME}/.spack
fi



if [ -n "${CI:-}" ]; then
  start_section "Add buildcache mirror"
  mirror_name="acts-spack-buildcache"
  mirror_url="oci://ghcr.io/acts-project/spack-buildcache"
  if [ -n "${GITLAB_CI:-}" ]; then
  # Use CERN mirror for non-Github Actions
    mirror_url="oci://registry.cern.ch/ghcr.io/acts-project/spack-buildcache"
  fi

  # Check if this buildcache is already configured
  if ! spack mirror list | grep -q ${mirror_name}; then
    echo "Adding buildcache ${mirror_name}"
    spack mirror add ${mirror_name} ${mirror_url} --unsigned
  fi
  end_section

  start_section "Locate OpenGL"
  "${SCRIPT_DIR}/opengl.sh"
  end_section
fi

start_section "Get spack lock file"
arch=$(spack arch --family)

env_dir="${destination}/env"
view_dir="${destination}/view"
mkdir -p ${env_dir}

lock_file_path="${destination}/spack.lock"
cmd=(
    "${SCRIPT_DIR}/select_lockfile.py"
    "--tag" "${tag}"
    "--arch" "${arch}"
    "--output" "${lock_file_path}"
)

if [ "${compiler}" != "default" ]; then
    cmd+=("--compiler-binary" "${compiler}")
fi

"${cmd[@]}"

end_section



start_section "Create spack environment"
time spack env create -d "${env_dir}" "${lock_file_path}" --with-view "$view_dir"
time spack -e "${env_dir}" spec -l
time spack -e "${env_dir}" find
end_section

start_section "Install spack packages"
if [ "$(uname)" = "Darwin" ]; then
  NCPUS=$(sysctl -n hw.ncpu)
else
  NCPUS=$(nproc)
fi
time "${SCRIPT_DIR}"/parallel.sh "$NCPUS" spack -e "${env_dir}" install --use-buildcache only \
  | tee install.log \
  | grep -v "^Waiting\|^\[+\]"
end_section

start_section "Patch up Geant4 data directory"
# ${SCRIPT_DIR}/with_spack_env.sh ${env_dir} geant4-config --install-datasets
geant4_dir=$(spack -e "${env_dir}" location -i geant4)
# Prepare the folder for G4 data, and symlink it to where G4 will look for it
mkdir -p "${geant4_dir}/share/Geant4"
ln -s "${geant4_dir}/share/Geant4/data" "${view_dir}/share/Geant4/data"
end_section


start_section "Prepare python environment"
ls -al
venv_dir="${view_dir}/venv"
"${view_dir}"/bin/python3 -m venv \
  --system-site-packages \
  "$venv_dir"

"${venv_dir}/bin/python3" -m pip install pyyaml jinja2

end_section

start_section "Set environment variables"
if [ -n "${GITHUB_ACTIONS:-}" ]; then
  echo "${view_dir}/bin" >> "$GITHUB_PATH"
  echo "${venv_dir}/bin" >> "$GITHUB_PATH"
fi
set_env PATH "${venv_dir}/bin:${view_dir}/bin/:${PATH}"
set_env CMAKE_PREFIX_PATH "${venv_dir}:${view_dir}"
set_env LD_LIBRARY_PATH "${view_dir}/lib"
set_env ROOT_INCLUDE_PATH "${view_dir}/include"
# Geant4 puts CLHEP in a subdirectory
set_env ROOT_INCLUDE_PATH "${view_dir}/include/Geant4"
end_section

# Pythia8 looks for settings in this directory
# set_env PYTHIA8DATA "${destination}/share/Pythia8/xmldoc"

set +u
