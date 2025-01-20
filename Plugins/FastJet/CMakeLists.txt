add_library(ActsPluginFastJet INTERFACE)

target_include_directories(
    ActsPluginFastJet
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPluginFastJet INTERFACE FastJet)

install(
    TARGETS ActsPluginFastJet
    EXPORT ActsPluginFastJetTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
