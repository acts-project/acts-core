# Third party software

This contains software that is usually not available through the system
package manager.

**Note** Only include software here that is either a header-only library or
can be built as a static library. Do not use any of this software as part of
the public interface and only include it in pure implementation files to avoid
issues with missing files after installation.

`nlohmann_json` is exempted from this rule as it is handled specially.

## dfelibs

A copy of [dfelibs](https://github.com/msmk0/dfelibs) v20200416, with the
following folders removed:

-   examples
-   unittests

## nlohmann_json

A copy of [nlohmann::json](https://github.com/nlohmann/json) revision 84f19d6, with the
following folders removed:

-   benchmarks
-   doc
-   test
-   include
-   third_party

## autodiff

A copy of [autodiff](https://github.com/autodiff/autodiff), v0.5.12 (however, the 
version given in the `CMakeLists.txt` of autodiff has not been changed by the maintainers
and is still v0.5.11). In the `CMakeLists.txt` the commands `add_subdirectory(tests)` and 
`add_subdirectory(examples)` have been commented out. All folders/files have been 
removed except the following ones:

-   autodiff (contains the header files)
-   cmake
-   CMakeLists.txt

## boost 

For convenience, it's possible to use the ACTS build system to build the minimum
required version of [boost](https://www.boost.org/) (currently 1.71.0).  No source is
bundled here, and if requested via "-DACTS_USE_SYSTEM_BOOST=OFF", only the filesystem,
program_options, and test libraries will be built.

Warning: during installation, the built boost libraries will be installed alongside the
ACTS libraries, with a version suffix. This location may be known to the system linker.

## eigen3

For convenience, it's possible to use the ACTS build system to build
the minimum required version of [Eigen](https://eigen.tuxfamily.org)
(currently 3.3.7), with "-DACTS_USE_SYSTEM_EIGEN3=OFF".

## pybind11

[Pybind11](https://github.com/pybind/pybind11) is used to create python bindings for 
the examples. The following files/directories are removed:
-   docs
-   tests
-   .github