# Compiler flag definitions common to Chaste & third party library builds

set(CMAKE_USER_MAKE_RULES_OVERRIDE
   ${CMAKE_CURRENT_LIST_DIR}/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_CURRENT_LIST_DIR}/cxx_flag_overrides.cmake)

#For GUI configs. Change C, and CXX compiler flags dynamically to static, debug build.
#The cmake CMAKE_USER_MAKE_RULES_OVERRIDEs above takes care of non GUI builds.
foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
   if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
   endif(${flag_var} MATCHES "/MD")
endforeach(flag_var)

# A bunch of compiler flags. Since we now propagate the compiler and linker flags to all external projects
# this is one place to set all the flags.
add_definitions(
    -Z7     # => embed debugging info in library as opposed to using an external .pdb database
    -wd4996 # => disable insecure api warnings
    -wd4267 # => disable "possible loss of data due to 'narrowing' conversion, e.g. size_t to int"
    -wd4290 # => disable warning "C++ exception specification ignored except to indicate a function is not __declspec(nothrow)"
    -wd4005 # => 'identifier' : macro redefinition
    -wd4018 # => 'expression' : signed/unsigned mismatch
    -wd4244 # => 'argument' : conversion from 'type1' to 'type2', possible loss of data
    -wd4101 # => 'identifier' : unreferenced local variable
    -wd4661 # => 'identifier' : no suitable definition provided for explicit template instantiation request
)

if(NOT "${PROJECT_NAME}" MATCHES "Chaste")
    # Some third party libaries need further suppressions
    add_definitions(
        -MTd    # => static debug build
        -Yu     # => use precompiled headers: this prevents a linker warning: "H5FDdirect.obj : warning LNK4221: This object
                # file does not define any previously undefined public symbols, so it will not be used by any link operation
                # that consumes this library".
        -wd4554 # => possible operator precedence error warning
        -wd4305 # => truncation from type1 to type2
        -wd4133 # => 'type' : incompatible types - from 'type1' to 'type2'
    )
endif()
