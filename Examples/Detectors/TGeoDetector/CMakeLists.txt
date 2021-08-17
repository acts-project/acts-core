add_library(
  ActsExamplesDetectorTGeo SHARED
  src/BuildTGeoDetector.cpp
  src/TGeoDetector.cpp)
  
target_include_directories(
  ActsExamplesDetectorTGeo
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorTGeo
  PUBLIC
    ActsCore ActsPluginIdentification ActsPluginDigitization ActsPluginTGeo
    ActsExamplesFramework ActsExamplesDetectorsCommon
    ActsExamplesDetectorGeneric ActsExamplesIoJson)

install(
  TARGETS ActsExamplesDetectorTGeo
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
