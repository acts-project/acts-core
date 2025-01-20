cmake_minimum_required(VERSION 3.25)

message(STATUS "Building Eigen ${_acts_eigen3_version}")

include(ExternalProject)

set(Eigen3_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen3-prefix/src/eigen3)

ExternalProject_Add(
    eigen3
    SYSTEM
    "${ACTS_EIGEN3_SOURCE};CONFIGURE_COMMAND;;BUILD_COMMAND;;INSTALL_COMMAND;"
)

add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
add_dependencies(Eigen3::Eigen eigen3)
target_include_directories(Eigen3::Eigen INTERFACE ${Eigen3_INCLUDE_DIR})
install(
    DIRECTORY ${Eigen3_INCLUDE_DIR}/Eigen
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
