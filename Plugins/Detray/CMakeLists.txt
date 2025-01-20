add_library(
    ActsPluginDetray
    SHARED
    src/DetrayConversionUtils.cpp
    src/DetrayConverter.cpp
    src/DetrayGeometryConverter.cpp
    src/DetrayMaterialConverter.cpp
    src/DetraySurfaceGridsConverter.cpp
)

add_dependencies(ActsPluginDetray detray::core covfie::core vecmem::core)

target_include_directories(
    ActsPluginDetray
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    ActsPluginDetray
    PUBLIC
        ActsCore
        ActsPluginJson
        detray::core
        detray::core_array
        detray::io
        detray::detectors
        detray::test_utils
        detray::validation_utils
        vecmem::core
)

install(
    TARGETS ActsPluginDetray
    EXPORT ActsPluginDetrayTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
