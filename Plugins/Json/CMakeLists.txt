include(ActsTargetLinkLibrariesSystem)

add_library(
    ActsPluginJson
    SHARED
    src/AlgebraJsonConverter.cpp
    src/DetectorJsonConverter.cpp
    src/DetectorVolumeJsonConverter.cpp
    src/ExtentJsonConverter.cpp
    src/GridJsonConverter.cpp
    src/DetectorVolumeFinderJsonConverter.cpp
    src/IndexedSurfacesJsonConverter.cpp
    src/JsonMaterialDecorator.cpp
    src/MaterialMapJsonConverter.cpp
    src/MaterialJsonConverter.cpp
    src/PortalJsonConverter.cpp
    src/ProtoDetectorJsonConverter.cpp
    src/ProtoAxisJsonConverter.cpp
    src/SurfaceBoundsJsonConverter.cpp
    src/SurfaceJsonConverter.cpp
    src/UtilitiesJsonConverter.cpp
    src/VolumeBoundsJsonConverter.cpp
    src/VolumeJsonConverter.cpp
    src/AmbiguityConfigJsonConverter.cpp
    src/DetrayJsonHelper.cpp
    src/JsonDetectorElement.cpp
    src/JsonSurfacesReader.cpp
)
target_include_directories(
    ActsPluginJson
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsPluginJson
    PUBLIC ActsCore nlohmann_json::nlohmann_json
)
target_compile_definitions(
    ActsPluginJson
    PRIVATE JSON_DISABLE_ENUM_SERIALIZATION=1
)

install(
    TARGETS ActsPluginJson
    EXPORT ActsPluginJsonTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
