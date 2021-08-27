add_library(
  ActsExamplesCommon SHARED
  src/CommonGeometry.cpp
  src/CommonOptions.cpp
  src/CommonSimulation.cpp
  src/GeometryExampleBase.cpp
  src/CommonMaterialMapping.cpp
  src/MaterialValidationBase.cpp
  src/PropagationExampleBase.cpp)
target_include_directories(
  ActsExamplesCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesCommon
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    ActsExamplesDetectorsCommon ActsExamplesPropagation ActsExamplesGenerators
    ActsExamplesMaterialMapping ActsExamplesIoCsv ActsExamplesIoJson
    ActsExamplesIoRoot ActsExamplesIoObj ActsExamplesIoPerformance)

install(
  TARGETS ActsExamplesCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
