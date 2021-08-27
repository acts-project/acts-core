add_library(ActsAlignment SHARED
  src/Kernel/detail/AlignmentEngine.cpp)

target_compile_features(
  ActsAlignment
  PUBLIC cxx_std_17)

target_include_directories(
  ActsAlignment
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(
  ActsAlignment
  PUBLIC ActsCore)

install(
  TARGETS ActsAlignment
  EXPORT ActsAlignmentTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/ActsAlignment
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

