add_library(
  ActsExamplesFatras SHARED
  src/FatrasSimulation.cpp
  src/FatrasOptions.cpp)
target_include_directories(
  ActsExamplesFatras
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesFatras
  PUBLIC
    ActsCore ActsFatras
    ActsExamplesFramework ActsExamplesMagneticField
  PRIVATE
    Boost::program_options)

install(
  TARGETS ActsExamplesFatras
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
