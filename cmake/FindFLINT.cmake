find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FLINT QUIET flint)
endif()

find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    HINTS ${PC_FLINT_INCLUDEDIR}
    PATHS
    /opt/homebrew/include
    /usr/local/include
    /usr/include
)

find_library(FLINT_LIBRARY
    NAMES flint
    HINTS ${PC_FLINT_LIBDIR}
    PATHS
    /opt/homebrew/lib
    /usr/local/lib
    /usr/local/lib64
    /usr/lib
    /usr/lib64
)

set(FLINT_LIBRARIES ${FLINT_LIBRARY})
set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG FLINT_LIBRARIES
    FLINT_INCLUDE_DIRS)

mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY)
