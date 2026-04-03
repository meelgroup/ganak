find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FLINT QUIET flint)
endif()

# On macOS with Homebrew, query the exact install prefix for flint.
# brew --prefix returns e.g. /opt/homebrew/opt/flint which has lib/ and include/.
if(APPLE)
    execute_process(
        COMMAND brew --prefix flint
        OUTPUT_VARIABLE HOMEBREW_FLINT_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    HINTS ${PC_FLINT_INCLUDEDIR} ${HOMEBREW_FLINT_PREFIX}/include
    PATHS
    /opt/homebrew/include
    /usr/local/include
    /usr/include
)

find_library(FLINT_LIBRARY
    NAMES flint
    HINTS ${PC_FLINT_LIBDIR} ${HOMEBREW_FLINT_PREFIX}/lib
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
