add_library(
  ActsExamplesGeant4 SHARED
  src/EventStoreRegistry.cpp
  src/EventAction.cpp
  src/GdmlDetectorConstruction.cpp
  src/GeantinoRecording.cpp
  src/Geant4Options.cpp
  src/Geant4Simulation.cpp
  src/PrimaryGeneratorAction.cpp
  src/RunAction.cpp
  src/SteppingAction.cpp)
target_compile_definitions(
  ActsExamplesGeant4
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesGeant4
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesGeant4
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesGeant4
  PUBLIC ActsCore ActsExamplesFramework ${Geant4_LIBRARIES}
  PRIVATE Boost::program_options)

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
