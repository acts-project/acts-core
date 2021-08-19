add_library(
  ActsExamplesTrackFinding SHARED
  src/SeedingAlgorithm.cpp
  src/SeedingOrthogonalAlgorithm.cpp
  src/SpacePointMaker.cpp
  src/SpacePointMakerOptions.cpp
  src/TrackFindingAlgorithm.cpp
  src/TrackFindingAlgorithmFunction.cpp
  src/TrackFindingOptions.cpp
  src/TrackParamsEstimationAlgorithm.cpp)
target_include_directories(
  ActsExamplesTrackFinding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTrackFinding
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesIoJson
    ActsExamplesMagneticField
    Boost::program_options)

install(
  TARGETS ActsExamplesTrackFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
