cmake_minimum_required(VERSION 3.14)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# determine project version; sets _acts_version and _acts_commit_hash
include(ActsRetrieveVersion)

project(Acts VERSION ${_acts_version} LANGUAGES CXX)

# build options

# all options and compile-time parameters must be defined here for clear
# visibility and to make them available everywhere
#
# NOTE if you are adding a new option make sure that is defined in such a way
#   that it is off/empty by default. if you think that is not possible, then
#   it probably is not an optional component.
option(ACTS_BUILD_EVERYTHING "Build with most options enabled (except HepMC3 and documentation)" OFF)
# core related options
set(ACTS_PARAMETER_DEFINITIONS_HEADER "" CACHE FILEPATH "Use a different (track) parameter definitions header")
set(ACTS_LOG_FAILURE_THRESHOLD "" CACHE STRING "Log level above which an exception should be automatically thrown")
# plugins related options
option(ACTS_BUILD_PLUGIN_AUTODIFF "Build the autodiff plugin" OFF)
option(ACTS_USE_SYSTEM_AUTODIFF "Use autodiff provided by the system instead of the bundled version" OFF)
option(ACTS_BUILD_PLUGIN_CUDA "Build CUDA plugin" OFF)
option(ACTS_BUILD_PLUGIN_DD4HEP "Build DD4hep plugin" OFF)
option(ACTS_BUILD_PLUGIN_DIGITIZATION "Build Digitization plugin" OFF)
option(ACTS_BUILD_PLUGIN_IDENTIFICATION "Build Identification plugin" OFF)
option(ACTS_BUILD_PLUGIN_JSON "Build json plugin" OFF)
option(ACTS_USE_SYSTEM_NLOHMANN_JSON "Use nlohmann::json provided by the system instead of the bundled version" OFF)
option(ACTS_BUILD_PLUGIN_LEGACY "Build legacy plugin" OFF)
option(ACTS_BUILD_PLUGIN_ONNX "Build ONNX plugin" OFF)
option(ACTS_BUILD_PLUGIN_SYCL "Build SYCL plugin" OFF)
option(ACTS_BUILD_PLUGIN_TGEO "Build TGeo plugin" OFF)
# fatras related options
option(ACTS_BUILD_FATRAS "Build FAst TRAcking Simulation package" OFF)
option(ACTS_BUILD_FATRAS_GEANT4 "Build Geant4 Fatras package" OFF)
# alignment related options
option(ACTS_BUILD_ALIGNMENT "Build Alignment package" OFF)
# examples related options
option(ACTS_BUILD_EXAMPLES "Build standalone examples" OFF)
option(ACTS_BUILD_EXAMPLES_DD4HEP "Build DD4hep-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_GEANT4 "Build Geant4-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_HEPMC3 "Build HepMC3-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_PYTHIA8 "Build Pythia8-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_SCRIPTS "Build Analysis applications in the examples" OFF)
# test related options
option(ACTS_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(ACTS_BUILD_INTEGRATIONTESTS "Build integration tests" OFF)
option(ACTS_BUILD_UNITTESTS "Build unit tests" OFF)
# other options
option(ACTS_BUILD_DOCS "Build documentation" OFF)
option(ACTS_USE_SYSTEM_BOOST "Use a system-provided boost" ON)
option(ACTS_USE_SYSTEM_EIGEN3 "Use a system-provided eigen3" ON)

# handle option inter-dependencies and the everything flag
# NOTE: ordering is important here. dependencies must come before dependees
include(ActsOptionHelpers)
# options without inter-dependencies
set_option_if(ACTS_BUILD_BENCHMARKS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_UNITTESTS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_INTEGRATIONTESTS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_DD4HEP ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_GEANT4 ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_HEPMC3 ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_PYTHIA8 ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_FATRAS_GEANT4 ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_FATRAS ACTS_BUILD_FATRAS_GEANT4)
set_option_if(ACTS_BUILD_ALIGNMENT ACTS_BUILD_EVERYTHING)
# any examples component activates the general examples option
set_option_if(
  ACTS_BUILD_EXAMPLES
  ACTS_BUILD_EXAMPLES_DD4HEP
  OR ACTS_BUILD_EXAMPLES_GEANT4
  OR ACTS_BUILD_EXAMPLES_HEPMC3
  OR ACTS_BUILD_EXAMPLES_PYTHIA8
  OR ACTS_BUILD_EVERYTHING)
# core plugins might be required by examples or depend on each other
set_option_if(
  ACTS_BUILD_PLUGIN_DD4HEP
  ACTS_BUILD_EXAMPLES_DD4HEP OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_PLUGIN_TGEO
  ACTS_BUILD_PLUGIN_DD4HEP OR ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_PLUGIN_IDENTIFICATION
  ACTS_BUILD_PLUGIN_TGEO OR ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_PLUGIN_DIGITIZATION
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_PLUGIN_JSON
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_FATRAS
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_PLUGIN_LEGACY ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_PLUGIN_AUTODIFF ACTS_BUILD_EVERYTHING)

# feature tests
include(CheckCXXSourceCompiles)

# function that tests if the root installation is compatible (i.e., compiled with C++17)
function(check_root_compatibility)
  get_target_property(ROOT_INCLUDE_DIR ROOT::Core INTERFACE_INCLUDE_DIRECTORIES)
  set(CMAKE_REQUIRED_FLAGS -std=c++17)
  set(CMAKE_REQUIRED_INCLUDES ${ROOT_INCLUDE_DIR})
  check_cxx_source_compiles("#include <string>\n#include <TString.h>\nint main(){}" ROOT_COMPATIBILITY_CHECK)
  if(NOT ROOT_COMPATIBILITY_CHECK)
    message(FATAL_ERROR "Root installation is misconfigured. Ensure that your Root installation was compiled with C++17.")
  endif()
endfunction()

# additional configuration and tools
include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include(ActsCompilerOptions) # default compile options
include(ActsComponentsHelpers) # handle components via add_..._if commands
# place build products in `<build_dir>/bin` and `<build_dir>/lib` for simple use
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# minimal dependency versions. they are defined here in a single place so
# they can be easily upgraded, although they might not be used if the
# dependency is included via `add_subdirectory(...)`.
set(_acts_autodiff_version 0.5.11)
set(_acts_boost_version 1.71.0)
set(_acts_dd4hep_version 1.11)
set(_acts_doxygen_version 1.8.15)
set(_acts_eigen3_version 3.3.7)
set(_acts_hepmc3_version 3.2.1)
set(_acts_nlohmanjson_version 3.2.0)
set(_acts_root_version 6.20)
set(_acts_tbb_version 2020.1)

# required packages

if (ACTS_USE_SYSTEM_BOOST)
  # NOTE we use boost::filesystem instead of std::filesystem since the later
  # is not uniformly supported even on compilers that nominally support C++17
  # NOTE FindBoost.cmake looks for BoostConfig.cmake first, before running it's own logic.
  find_package(Boost ${_acts_boost_version} REQUIRED COMPONENTS filesystem program_options unit_test_framework)
else()
  add_subdirectory(thirdparty/boost)
endif()

if (ACTS_USE_SYSTEM_EIGEN3)
  find_package(Eigen3 ${_acts_eigen3_version} REQUIRED CONFIG)
else()
  add_subdirectory(thirdparty/eigen3)
endif()

# the `<project_name>_VERSION` variables set by `setup(... VERSION ...)` have
# only local scope, i.e. they are not accessible her for dependencies added
# via `add_subdirectory`. this overrides the `project(...)` funcion for
# sub-projects such that the resulting `<project_name>_VERSION` has
# global scope and is accessible within the main project later on.
cmake_policy(SET CMP0048 NEW)
macro(project)
  _project(${ARGN})
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}" CACHE INTERNAL "")
endmacro()

# optional packages
#
# find packages explicitly for each component even if this means searching for
# the same package twice. This avoids having complex if/else trees to sort out
# when a particular package is actually needed.
# plugin dependencies
if(ACTS_BUILD_PLUGIN_AUTODIFF)
  if(ACTS_USE_SYSTEM_AUTODIFF)
    find_package(autodiff ${_acts_autodiff_version} CONFIG REQUIRED)
    message(STATUS "Using system installation of autodiff")
  else()
    add_subdirectory(thirdparty/autodiff)
    add_library(autodiff::autodiff ALIAS autodiff)
    message(STATUS "Using bundled autodiff")
  endif()
endif()
if(ACTS_BUILD_PLUGIN_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14 CACHE STRING "CUDA C++ standard to use")
  set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL
    "Force the C++ standard requirement")
  set(CMAKE_CUDA_ARCHITECTURES "35;52;75" CACHE STRING
    "CUDA architectures to generate code for")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G")
endif()
if(ACTS_BUILD_PLUGIN_DD4HEP)
  find_package(DD4hep ${_acts_dd4hep_version} REQUIRED CONFIG COMPONENTS DDCore DDDetectors)
endif()
if(ACTS_BUILD_PLUGIN_JSON)
  if(ACTS_USE_SYSTEM_NLOHMANN_JSON)
    message(STATUS "Using system installation of nlohmann::json")
    find_package(nlohmann_json ${_acts_nlohmanjson_version} REQUIRED CONFIG)
  else()
    message(STATUS "Using bundled nlohmann::json")
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(thirdparty/nlohmann_json)
  endif()
endif()
if(ACTS_BUILD_PLUGIN_ONNX)
  find_package(OnnxRuntime REQUIRED)
endif()
if(ACTS_BUILD_PLUGIN_SYCL)
  find_package(SYCL REQUIRED)
endif()
if(ACTS_BUILD_PLUGIN_TGEO)
  find_package(ROOT ${_acts_root_version} REQUIRED CONFIG COMPONENTS Geom)
  check_root_compatibility()
endif()
# examples dependencies
if(ACTS_BUILD_EXAMPLES)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  # for simplicity always request all potentially required components.
  find_package(ROOT ${_acts_root_version} REQUIRED CONFIG COMPONENTS Core Geom GenVector Hist Tree TreePlayer)
  check_root_compatibility()
  # newer DD4hep version require TBB and search interally for TBB in
  # config-only mode. to avoid mismatches we explicitely search using
  # config-only mode first to be sure that we find the same version.
  find_package(TBB ${_acts_tbb_version} CONFIG)
  if(NOT TBB_FOUND)
    # no version check possible when using the find module
    find_package(TBB ${_acts_tbb_version} MODULE REQUIRED)
  endif()
  add_subdirectory(thirdparty/dfelibs)
endif()
if(ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_GEANT4)
  find_package(DD4hep ${_acts_dd4hep_version} REQUIRED CONFIG COMPONENTS DDCore DDG4 DDDetectors)
elseif(ACTS_BUILD_EXAMPLES_DD4HEP)
  find_package(DD4hep ${_acts_dd4hep_version} REQUIRED CONFIG COMPONENTS DDCore DDDetectors)
endif()
if(ACTS_BUILD_EXAMPLES_GEANT4)
  find_package(Geant4 REQUIRED CONFIG COMPONENTS gdml)
endif()
if(ACTS_BUILD_EXAMPLES_HEPMC3)
  find_package(HepMC3 ${_acts_hepmc3_version} REQUIRED CONFIG)
endif()
if(ACTS_BUILD_EXAMPLES_PYTHIA8)
  find_package(Pythia8 REQUIRED)
endif()
# other dependencies
if(ACTS_BUILD_DOCS)
  find_package(Doxygen ${_acts_doxygen_version} REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

if(ACTS_CUSTOM_SCALARTYPE)
    message(STATUS "Building Acts with custom scalar type: ${ACTS_CUSTOM_SCALARTYPE}")
endif()

# core library, core plugins, and other components
add_component(Core Core)
add_subdirectory(Plugins)
add_component_if(Fatras Fatras ACTS_BUILD_FATRAS)
add_component_if(Alignment Alignment ACTS_BUILD_ALIGNMENT)

# examples
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)

# automated tests and benchmarks
if(ACTS_BUILD_BENCHMARKS OR ACTS_BUILD_INTEGRATIONTESTS OR ACTS_BUILD_UNITTESTS)
  enable_testing() # must be set in the main CMakeLists.txt
  add_subdirectory(Tests)
endif()

# documentation
add_subdirectory_if(docs ACTS_BUILD_DOCS)


# create cmake configuration files and environment setup script
include(ActsCreatePackageConfig)
include(ActsCreateSetup)
