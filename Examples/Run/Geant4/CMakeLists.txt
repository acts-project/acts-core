add_subdirectory(Common)

add_executable(
  ActsExampleGeantinoRecordingGdml
  GdmlGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleGeantinoRecordingGdml
  PRIVATE
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4Common
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS 
    ActsExampleGeantinoRecordingGdml 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(HepMC ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_HEPMC3)