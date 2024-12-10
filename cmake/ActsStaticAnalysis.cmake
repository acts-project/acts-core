if(ACTS_RUN_CLANG_TIDY)
    find_program(CLANG_TIDY_COMMAND NAMES clang-tidy)
    if(NOT CLANG_TIDY_COMMAND)
        message(
            WARNING
            "ACTS_RUN_CLANG_TIDY is ON but clang-tidy is not found!"
        )
        set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE)
    else()
        message(STATUS "Setting up clang-tidy run")

        set(_chks "")
        list(APPEND _chks "-*")
        list(APPEND _chks "clang-analyzer-optin.cplusplus.UninitializedObject")
        list(APPEND _chks "cppcoreguidelines-init-variables")
        list(APPEND _chks "cppcoreguidelines-pro-type-member-init")
        list(APPEND _chks "google-readability-casting")
        list(APPEND _chks "google-explicit-constructor")
        list(APPEND _chks "modernize-concat-nested-namespaces")
        list(APPEND _chks "modernize-use-equals-default")
        list(APPEND _chks "modernize-use-default-member-init")
        list(APPEND _chks "modernize-use-nullptr")
        list(APPEND _chks "modernize-use-override")
        list(APPEND _chks "modernize-use-using")
        list(APPEND _chks "performance-for-range-copy")
        list(APPEND _chks "performance-move-const-arg")
        list(APPEND _chks "performance-unnecessary-value-param")
        list(APPEND _chks "readability-braces-around-statements")
        list(APPEND _chks "readability-container-size-empty")
        list(APPEND _chks "readability-implicit-bool-cast")
        list(APPEND _chks "readability-implicit-bool-conversion")
        list(APPEND _chks "readability-inconsistent-declaration-parameter-name")
        list(APPEND _chks "readability-named-parameter")
        list(APPEND _chks "readability-operators-representation")
        list(JOIN _chks "," CLANG_TIDY_CHECKS)

        message(STATUS "Configured checks")
        foreach(_chk ${_chks})
            message(STATUS "|-> ${_chk}")
        endforeach()

        set(_errs "")
        list(JOIN _errs "," CLANG_TIDY_ERRORS)

        message(STATUS "Enabled errors:")
        foreach(_err ${_errs})
            message(STATUS "|-> ${_err}")
        endforeach()

        set(CLANG_TIDY_HEADER_FILTER ".*")

        set(CMAKE_CXX_CLANG_TIDY
            "${CLANG_TIDY_COMMAND};-checks=${CLANG_TIDY_CHECKS};-header-filter=${CLANG_TIDY_HEADER_FILTER};-warnings-as-errors=${CLANG_TIDY_ERRORS}"
        )
    endif()
endif()
