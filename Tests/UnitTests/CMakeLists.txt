# add a unittest executable w/ default dependencies and register it.

# the common libraries which are linked to every unittest can be
# extended by setting the `unittest_extra_libraries` variables before
# calling the macro.
macro(add_unittest _name)
  # automatically prefix the target name
  set(_target "ActsUnitTest${_name}")
  add_executable(${_target} ${ARGN})
  # define required BOOST_TEST_... macros here to ensure consistent names
  target_compile_definitions(
    ${_target}
    PRIVATE "-DBOOST_TEST_DYN_LINK" "-DBOOST_TEST_MODULE=${_target}")
  target_include_directories(
    ${_target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(
    ${_target}
    PRIVATE
      ActsCore
      ActsTestsCommonHelpers
      Boost::unit_test_framework
      ${unittest_extra_libraries})
  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

add_subdirectory(Core)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory_if(Benchmarks ACTS_BUILD_BENCHMARKS)
add_subdirectory_if(Fatras ACTS_BUILD_FATRAS)
add_subdirectory(Plugins)
add_subdirectory_if(Alignment ACTS_BUILD_ALIGNMENT)
