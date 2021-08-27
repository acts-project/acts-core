add_library(
  ActsExamplesCommon SHARED
  src/CommonGeometry.cpp
  src/CommonOptions.cpp
  src/GeometryExampleBase.cpp
  src/MaterialMappingBase.cpp
  src/MaterialValidationBase.cpp
  src/PropagationExampleBase.cpp
  src/ReconstructionBase.cpp)
target_include_directories(
  ActsExamplesCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesCommon
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    ActsExamplesDetectorsCommon ActsExamplesPropagation
    ActsExamplesMaterialMapping ActsExamplesIoCsv ActsExamplesIoJson
    ActsExamplesIoRoot ActsExamplesIoObj ActsExamplesIoPerformance
    ActsExamplesDigitization ActsExamplesTrackFitting
    ActsExamplesTrackFinding ActsExamplesTruthTracking
    ActsExamplesIoCsv)

install(
  TARGETS ActsExamplesCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
