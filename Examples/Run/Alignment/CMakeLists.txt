# shared code 
add_subdirectory(Common)

# Alignment for Generic Contextual detector 
add_executable(
  ActsExampleDetectorAlignContextual
  ContextualAlignment.cpp)
target_link_libraries(
    ActsExampleDetectorAlignContextual 
  PRIVATE
    ActsExamplesAlignmentCommon
    ActsExamplesDetectorContextual
)


install(
  TARGETS
  ActsExampleDetectorAlignContextual
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
