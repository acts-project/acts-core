add_library(ActsExamplesDetectorGeoModel SHARED src/GeoModelDetector.cpp)

target_include_directories(
    ActsExamplesDetectorGeoModel
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorGeoModel
    PUBLIC
        ActsCore
        ActsExamplesFramework
        ActsPluginGeoModel
        ActsExamplesDetectorsCommon
)

if(ACTS_BUILD_EXAMPLES_GEANT4)
    find_library(GeoModel2G4_LIBRARY GeoModel2G4 REQUIRED)

    target_link_libraries(
        ActsExamplesDetectorGeoModel
        PUBLIC ActsExamplesGeant4 ${GeoModel2G4_LIBRARY}
    )

    target_sources(
        ActsExamplesDetectorGeoModel
        PUBLIC
            src/GeoModelGeant4DetectorConstruction.cpp
            src/GeoModelDetectorGeant4.cpp
    )
else()
    target_sources(
        ActsExamplesDetectorGeoModel
        PUBLIC src/GeoModelDetectorGeant4Stub.cpp
    )
endif()

install(
    TARGETS ActsExamplesDetectorGeoModel
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
