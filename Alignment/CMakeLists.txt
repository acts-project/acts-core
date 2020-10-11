add_library(ActsAlignment INTERFACE)

target_include_directories(ActsAlignment
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(ActsAlignment INTERFACE
  ActsCore
  )

install(
  TARGETS  ActsAlignment 
  EXPORT ActsAlignmentTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

install(
  DIRECTORY include/ActsAlignment
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
