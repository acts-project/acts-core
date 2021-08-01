add_library(
  ActsExamplesAlignmentCommon SHARED
  DetectorAlignment.cpp)
target_include_directories(
  ActsExamplesAlignmentCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesAlignmentCommon
  PUBLIC
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesMagneticField
    ActsExamplesIoPerformance
    ActsExamplesAlignment)

install(
  TARGETS ActsExamplesAlignmentCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
