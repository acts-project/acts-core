add_library(
  ActsPluginTGeo SHARED
  src/TGeoCylinderDiscSplitter.cpp
  src/TGeoDetectorElement.cpp
  src/TGeoITkModuleSplitter.cpp
  src/TGeoLayerBuilder.cpp
  src/TGeoParser.cpp
  src/TGeoPrimitivesHelper.cpp
  src/TGeoSurfaceConverter.cpp)
target_include_directories(
  ActsPluginTGeo
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginTGeo
  PUBLIC ActsCore ActsPluginIdentification ROOT::Geom)

install(
  TARGETS ActsPluginTGeo
  EXPORT ActsPluginTGeoTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
