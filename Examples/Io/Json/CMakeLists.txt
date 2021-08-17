add_library(
  ActsExamplesIoJson SHARED
  src/JsonGeometryList.cpp
  src/JsonMaterialWriter.cpp
  src/JsonSurfacesWriter.cpp
  src/JsonOptionsWriter.cpp
  src/JsonDigitizationConfig.cpp
  src/JsonTGeoDetectorConfig.cpp)
target_include_directories(
  ActsExamplesIoJson
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoJson
  PUBLIC ActsCore ActsPluginJson ActsExamplesDigitization ActsExamplesFramework ActsExamplesMaterialMapping Boost::program_options)
# Needed to include tgeo plugin specific configs
if(ACTS_BUILD_PLUGIN_TGEO)
  target_link_libraries(ActsExamplesIoJson PUBLIC ActsPluginTGeo)
endif()

install(
  TARGETS ActsExamplesIoJson
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
