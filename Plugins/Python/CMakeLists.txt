add_library(ActsPluginPython INTERFACE)

install(
  TARGETS ActsPluginPython
  EXPORT ActsPluginPythonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(ActsPythonUtilities INTERFACE)
target_include_directories(ActsPythonUtilities INTERFACE include)

set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(ActsPythonBindings 
  src/ModuleEntry.cpp
  src/Base.cpp
  src/Detector.cpp
  src/Material.cpp
  src/Geometry.cpp
  src/ExampleAlgorithms.cpp
  src/MagneticField.cpp
  src/Output.cpp
  src/Input.cpp
  src/Propagation.cpp
  src/Generators.cpp
  src/TruthTracking.cpp
  src/TrackFitting.cpp
  src/TrackFinding.cpp
  src/Vertexing.cpp
)
install(TARGETS ActsPythonBindings DESTINATION .)

set_target_properties(ActsPythonBindings PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
set_target_properties(ActsPythonBindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)

target_link_libraries(ActsPythonBindings PUBLIC 
  ActsCore 
  ActsPythonUtilities
  ActsExamplesFramework
  ActsExamplesDetectorsCommon
  ActsExamplesDetectorGeneric
  ActsExamplesDetectorContextual
  ActsExamplesDetectorTGeo
  ActsExamplesMagneticField
  ActsExamplesIoRoot
  ActsExamplesIoNuclearInteractions
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoJson
  ActsExamplesIoPerformance
  ActsExamplesGenerators
  ActsExamplesTrackFinding
  ActsExamplesTrackFitting
  ActsExamplesVertexing
  ActsExamplesFatras
  ActsExamplesPrinters
)

set(py_files 
  __init__.py 
  examples/__init__.py
  _adapter.py
)

if(ACTS_BUILD_PLUGIN_JSON)
  target_link_libraries(ActsPythonBindings PUBLIC ActsPluginJson)
  target_sources(ActsPythonBindings PRIVATE src/Json.cpp)
else()
  target_sources(ActsPythonBindings PRIVATE src/JsonStub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
  pybind11_add_module(ActsPythonBindingsDD4hep src/DD4hepComponent.cpp)
  target_link_libraries(ActsPythonBindingsDD4hep PUBLIC 
    ActsPythonUtilities
    ActsExamplesDetectorDD4hep)
  add_dependencies(ActsPythonBindings ActsPythonBindingsDD4hep)

  install(TARGETS ActsPythonBindingsDD4hep DESTINATION .)
  set_target_properties(ActsPythonBindingsDD4hep PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
  set_target_properties(ActsPythonBindingsDD4hep PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
  list(APPEND py_files examples/dd4hep.py)
endif()

if(ACTS_BUILD_EXAMPLES_PYTHIA8)
  target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesGeneratorsPythia8)
  target_sources(ActsPythonBindings PRIVATE src/Pythia8.cpp)
else()
  target_sources(ActsPythonBindings PRIVATE src/Pythia8Stub.cpp)
endif()

if(ACTS_BUILD_PLUGIN_DIGITIZATION AND ACTS_BUILD_PLUGIN_IDENTIFICATION)
  target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesDigitization)
  target_sources(ActsPythonBindings PRIVATE src/Digitization.cpp)
else()
  target_sources(ActsPythonBindings PRIVATE src/DigitizationStub.cpp)
endif()

if(ACTS_BUILD_EXAMPLES_GEANT4)
  pybind11_add_module(ActsPythonBindingsGeant4 src/Geant4Component.cpp)
  target_link_libraries(ActsPythonBindingsGeant4 PUBLIC 
    ActsExamplesGeant4 
    ActsPythonUtilities)
  add_dependencies(ActsPythonBindings ActsPythonBindingsGeant4)

  install(TARGETS ActsPythonBindingsGeant4 DESTINATION .)
  set_target_properties(ActsPythonBindingsGeant4 PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
  set_target_properties(ActsPythonBindingsGeant4 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
  list(APPEND py_files examples/geant4/__init__.py)

  if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
    pybind11_add_module(ActsPythonBindingsGeant4DD4hep src/Geant4DD4hepComponent.cpp)
    target_link_libraries(ActsPythonBindingsGeant4DD4hep PUBLIC 
      ActsPythonUtilities
      ActsExamplesGeant4
      ActsExamplesGeant4DD4hep 
      ActsExamplesDetectorDD4hep)
    add_dependencies(ActsPythonBindings ActsPythonBindingsGeant4DD4hep)

    install(TARGETS ActsPythonBindingsGeant4DD4hep DESTINATION .)
    set_target_properties(ActsPythonBindingsGeant4DD4hep PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
    set_target_properties(ActsPythonBindingsGeant4DD4hep PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
    list(APPEND py_files examples/geant4/dd4hep.py)
  endif()

  if(ACTS_BUILD_EXAMPLES_HEPMC3)
    target_link_libraries(ActsPythonBindingsGeant4 PUBLIC ActsExamplesEventRecording)
    target_sources(ActsPythonBindingsGeant4 PRIVATE src/Geant4HepMC3.cpp)
    list(APPEND py_files examples/geant4/hepmc3.py)
  else()
    target_sources(ActsPythonBindingsGeant4 PRIVATE src/Geant4HepMC3Stub.cpp)
  endif()

endif()

if(ACTS_BUILD_EXAMPLES_HEPMC3)
  target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesHepMC3)
  target_sources(ActsPythonBindings PRIVATE src/HepMC3.cpp)
  list(APPEND py_files examples/hepmc3.py)
else()
  target_sources(ActsPythonBindings PRIVATE src/HepMC3Stub.cpp)
endif()


add_custom_target(ActsPythonGlueCode)
configure_file(setup.sh.in ${_python_dir}/setup.sh)


foreach(f ${py_files})
set(_target ${_python_dir}/acts/${f})
get_filename_component(_dir ${_target} DIRECTORY)
file(MAKE_DIRECTORY ${_dir})

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/${f} ${_target})
endforeach()

add_dependencies(ActsPythonBindings ActsPythonGlueCode)

  