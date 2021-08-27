add_library(
  ActsExamplesGeant4Common SHARED
  Geant4.cpp)
target_include_directories(
  ActsExamplesGeant4Common
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesGeant4Common
  PUBLIC
    ActsCore
    ActsExamplesFramework 
    ActsExamplesCommon
    ActsExamplesGenerators
    ActsExamplesMagneticField
    ActsExamplesGeant4    
    ActsExamplesIoCsv 
    ActsExamplesIoRoot
    ${Geant4_LIBRARIES}
    Boost::program_options)

install(
  TARGETS ActsExamplesGeant4Common
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
