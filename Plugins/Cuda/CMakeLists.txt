# Set up the "CUDA 1" plugin.
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

add_library(
  ActsPluginCuda SHARED
  src/Seeding/Kernels.cu
)
target_include_directories(
  ActsPluginCuda
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(
  ActsPluginCuda
  PUBLIC ActsCore ${CUDART_LIBRARY})

# Set up the "CUDA 2" plugin.
add_library(
  ActsPluginCuda2 STATIC
  include/Acts/Plugins/Cuda/Seeding2/Details/CountDublets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/FindDublets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/FindTriplets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/Types.hpp
  include/Acts/Plugins/Cuda/Seeding2/SeedFinder.hpp
  include/Acts/Plugins/Cuda/Seeding2/SeedFinder.ipp
  include/Acts/Plugins/Cuda/Seeding2/TripletFilterConfig.hpp
  include/Acts/Plugins/Cuda/Utilities/Arrays.hpp
  include/Acts/Plugins/Cuda/Utilities/Info.hpp
  include/Acts/Plugins/Cuda/Utilities/MemoryManager.hpp
  include/Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp
  include/Acts/Plugins/Cuda/Utilities/StreamWrapper.hpp
  src/Seeding2/CountDublets.cu
  src/Seeding2/FindDublets.cu
  src/Seeding2/FindTriplets.cu
  src/Utilities/Arrays.cu
  src/Utilities/ErrorCheck.cuh
  src/Utilities/ErrorCheck.cu
  src/Utilities/Info.cu
  src/Utilities/MatrixMacros.hpp
  src/Utilities/MemoryManager.cu
  src/Utilities/NewMemoryManager.cu
  src/Utilities/StreamHandlers.cuh
  src/Utilities/StreamWrapper.cu
)
target_include_directories(
  ActsPluginCuda2
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginCuda2
  PUBLIC ActsCore)
set_target_properties(ActsPluginCuda2 PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON)

# Install all CUDA plugins.
install(
  TARGETS ActsPluginCuda ActsPluginCuda2
  EXPORT ActsPluginCudaTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
