add_library(
  ActsExamplesIoCsv SHARED
  src/CsvMeasurementReader.cpp
  src/CsvMeasurementWriter.cpp
  src/CsvOptionsReader.cpp
  src/CsvOptionsWriter.cpp
  src/CsvParticleReader.cpp
  src/CsvParticleWriter.cpp
  src/CsvPlanarClusterReader.cpp
  src/CsvPlanarClusterWriter.cpp
  src/CsvSimHitReader.cpp
  src/CsvSimHitWriter.cpp
  src/CsvTrackingGeometryWriter.cpp
  src/CsvPropagationStepsWriter.cpp
  src/CsvMultiTrajectoryWriter.cpp)
  
target_include_directories(
  ActsExamplesIoCsv
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoCsv
  PRIVATE
    ActsCore ActsPluginDigitization ActsPluginIdentification
    ActsExamplesFramework ActsExamplesDigitization
    Threads::Threads Boost::program_options dfelibs)

install(
  TARGETS ActsExamplesIoCsv
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
