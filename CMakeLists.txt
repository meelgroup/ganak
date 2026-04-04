# Copyright (c) 2017, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# CMP0104 (CUDA architectures default) was introduced in 3.18; guard needed
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 NEW)
endif()

# -----------------------------------------------------------------------------
# Make RelWithDebInfo the default build type if otherwise not set
# -----------------------------------------------------------------------------
set(build_types Debug Release RelWithDebInfo MinSizeRel)
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "You can choose the type of build, options are:${build_types}")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Options are ${build_types}"
        FORCE
    )

    # Provide drop down menu options in cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${build_types})
endif()
message(STATUS "Doing a ${CMAKE_BUILD_TYPE} build")

# -----------------------------------------------------------------------------
# Option to enable/disable assertions
# -----------------------------------------------------------------------------
# Filter out definition of NDEBUG from the default build configuration flags.
# We will add this ourselves if we want to disable assertions
foreach(build_config ${build_types})
    string(TOUPPER ${build_config} upper_case_build_config)
    foreach(language CXX C)
        set(VAR_TO_MODIFY "CMAKE_${language}_FLAGS_${upper_case_build_config}")
        string(REGEX REPLACE "(^| )[/-]D *NDEBUG($| )"
                             " "
                             replacement
                             "${${VAR_TO_MODIFY}}")
        set(${VAR_TO_MODIFY} "${replacement}" CACHE STRING "Default flags for ${build_config} configuration" FORCE)
    endforeach()
endforeach()

set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)
set(CMAKE_EXPORT_PACKAGE_REGISTRY OFF)
project(ganak)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(GenerateExportHeader)
include(GNUInstallDirs)

message(STATUS "LIB directory is '${CMAKE_INSTALL_LIBDIR}'")
message(STATUS "BIN directory is '${CMAKE_INSTALL_BINDIR}'")

# contains some library search cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# generate JSON file of compile commands -- useful for code extension
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Python extension — defined early so the GMP/MPFR/flint discovery below can use it.
option(BUILD_PYTHON_EXTENSION "Build the pyganak Python extension module" OFF)

# static compilation
option(BUILD_SHARED_LIBS "Build the shared library" ON)
# STATIC_BINARY controls whether the ganak executable is linked fully statically.
# Defaults to ON when BUILD_SHARED_LIBS=OFF, but can be set to OFF independently
# (e.g. when building the Python wheel: BUILD_SHARED_LIBS=OFF embeds static .a
# deps in the .so, but we don't want a fully static ganak binary).
if(NOT BUILD_SHARED_LIBS)
    option(STATIC_BINARY "Link the ganak binary fully statically" ON)
else()
    option(STATIC_BINARY "Link the ganak binary fully statically" OFF)
endif()

if((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
    if(STATIC_BINARY AND ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        message(STATUS "Linking ganak binary statically")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
include(CheckLinkerFlag)
macro(add_cxx_flag_if_supported flagname)
  check_cxx_compiler_flag("${flagname}" HAVE_FLAG_${flagname})

  if(HAVE_FLAG_${flagname})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flagname}")
  endif()
endmacro()
macro(add_c_flag_if_supported flagname)
  check_c_compiler_flag("${flagname}" HAVE_FLAG_${flagname})

  if(HAVE_FLAG_${flagname})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flagname}")
  endif()
endmacro()
macro(add_link_flag_if_supported flagname)
  check_linker_flag(CXX "${flagname}" HAVE_FLAG_${flagname})

  if(HAVE_FLAG_${flagname})
    add_link_options("${flagname}")
  endif()
endmacro()

# Note: O3 gives slight speed increase, 1 more solved from SAT Comp'14 @ 3600s
if(NOT MSVC)
    add_compile_options(-g)

    add_compile_options("$<$<CONFIG:RelWithDebInfo>:-O3>")

    add_compile_options("$<$<CONFIG:Release>:-O3>")
    add_compile_options("$<$<CONFIG:Release>:-g0>")

    add_compile_options("$<$<CONFIG:Debug>:-O0>")

    foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG)
        if(${flag_var} MATCHES "-O3")
            string(REGEX REPLACE "-O3" "-O0" ${flag_var} "${${flag_var}}")
        endif()
        if(${flag_var} MATCHES "-O2")
            string(REGEX REPLACE "-O2" "-O0" ${flag_var} "${${flag_var}}")
        endif()
    endforeach()

else()
    # see https://msdn.microsoft.com/en-us/library/fwkeyyhe.aspx for details
    # /ZI = include debug info
    # /Wall = all warnings

    add_compile_options("$<$<CONFIG:RelWithDebInfo>:/Ox>")
    add_compile_options("$<$<CONFIG:RelWithDebInfo>:/ZI>")

    add_compile_options("$<$<CONFIG:Release>:/Ox>")
    add_compile_options("$<$<CONFIG:Release>:/DNDEBUG>")
    add_compile_options("$<$<CONFIG:Release>:/ZI>")

    add_compile_options("$<$<CONFIG:Debug>:/Od>")

    # buffers security check
    add_compile_options(/GS)

    # Proper warning level
    add_compile_options(/W1)

    # Disable STL used in DLL-boundary warning
    add_compile_options(/wd4251)
    add_compile_options(/D_CRT_SECURE_NO_WARNINGS)

    # Wall is MSVC's Weverything, so annoying unless used from the start
    # and with judiciously used warning disables
    # add_compile_options(/Wall)

    # /Za = only ansi C98 & C++11
    # /Za is not recommended for use, not tested, etc.
    # see: http://stackoverflow.com/questions/5489326/za-compiler-directive-does-not-compile-system-headers-in-vs2010
    # add_compile_options(/Za)

    add_compile_options(/fp:precise)

    # exception handling. s = The exception-handling model that catches C++ exceptions only and tells the compiler to assume that functions declared as extern "C" may throw an exception.
    # exception handling. c = If used with s (/EHsc), catches C++ exceptions only and tells the compiler to assume that functions declared as extern "C" never throw a C++ exception.
    add_compile_options(/EHsc)

    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /PDBCOMPRESS")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:1572864")

    set(DEF_INSTALL_CMAKE_DIR CMake)
endif()


option(ENABLE_ASSERTIONS "Build with assertions enabled" ON)
message(STATUS "build type is ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(ENABLE_ASSERTIONS OFF)
endif()

if(ENABLE_ASSERTIONS)
    # NDEBUG was already removed.
else()
    # Note this definition doesn't appear in the cache variables.
    add_compile_definitions(NDEBUG)
    add_cxx_flag_if_supported("-fno-stack-protector")
    add_compile_definitions(_FORTIFY_SOURCE=0)
endif()

# -----------------------------------------------------------------------------
# Flags
# -----------------------------------------------------------------------------
if(NOT WIN32)
    if(NOT ENABLE_TESTING AND ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        #add_cxx_flag_if_supported("-fvisibility=hidden")
    endif()

    add_cxx_flag_if_supported("-mpopcnt")
    add_cxx_flag_if_supported("-msse4.2")
    add_cxx_flag_if_supported("-mpclmul")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
      # add_cxx_flag_if_supported("-flto")
      # add_c_flag_if_supported("-flto")
      # add_link_flag_if_supported("-flto")
    else()
      add_cxx_flag_if_supported("-Wall")
      add_cxx_flag_if_supported("-Wextra")
      add_cxx_flag_if_supported("-Wunused")
      add_cxx_flag_if_supported("-Wsign-compare")
      add_cxx_flag_if_supported("-fno-omit-frame-pointer")
      add_cxx_flag_if_supported("-ggdb3")
      add_cxx_flag_if_supported("-g")
      add_cxx_flag_if_supported("-Wtype-limits")
      add_cxx_flag_if_supported("-Wuninitialized")
      add_cxx_flag_if_supported("-Wstrict-aliasing")
      add_cxx_flag_if_supported("-Wpointer-arith")
      add_cxx_flag_if_supported("-Wheader-guard")
      add_cxx_flag_if_supported("-Wformat-nonliteral")
      add_cxx_flag_if_supported("-Winit-self")
      add_cxx_flag_if_supported("-Wparentheses")
      add_cxx_flag_if_supported("-Wunreachable-code")
      add_cxx_flag_if_supported("-Wno-class-memaccess")
      add_cxx_flag_if_supported("-Wextra-semi-stmt")
      add_cxx_flag_if_supported("-Wnoweak-vtables")
      add_cxx_flag_if_supported("-Wdouble-promotion")
      add_cxx_flag_if_supported("-Warith-conversion")
    endif()
endif()

if(NOT WIN32)
    message(STATUS "${CMAKE_C_FLAGS}")
    add_c_flag_if_supported("-Wall")
    add_c_flag_if_supported("-msse4.2")
    add_c_flag_if_supported("-mpclmul")
    add_c_flag_if_supported("-mpopcnt")
endif()

option(SANITIZE "Use Clang sanitizers. You MUST use clang++ as the compiler for this to work" OFF)
if(SANITIZE)
    message(WARNING " --Using clang sanitizers -- you MUST use clang++ or the compile WILL fail")
    add_compile_options("-fsanitize=address")
    add_link_options("-fsanitize=address")
    add_link_options("-fsanitize=undefined")
    # add_compile_options("-fsanitize=integer")
    add_compile_options("-fsanitize=undefined")

    add_compile_options("-fsanitize=null")
    add_compile_options("-fsanitize=alignment")
    add_compile_options("-fno-sanitize-recover")

    add_compile_options("-fsanitize=return")
    add_compile_options("-fsanitize=bounds")
    # add_compile_options("-fsanitize=float-divide-by-zero")
    # add_compile_options("-fsanitize=integer-divide-by-zero")
    # add_compile_options("-fsanitize=unsigned-integer-overflow")
    # add_compile_options("-fsanitize=signed-integer-overflow")
    add_compile_options("-fsanitize=bool")
    add_compile_options("-fsanitize=enum")
    add_compile_options("-fsanitize=float-cast-overflow")
    add_compile_options("$<$<CONFIG:RelWithDebInfo>:-D_GLIBCXX_ASSERTIONS>")

    #add_compile_options("-Weverything")
    #add_compile_options("-Wshorten-64-to-32")
    #add_compile_options("-Wweak-vtables")
    #add_compile_options("-Wsign-conversion")
    #add_compile_options("-Wconversion")
endif()

option(ANALYZE "Run Clang static analyzer. You MUST use clang++ as the compiler for this to work" OFF)
if(ANALYZE)
    message(WARNING " --Using clang static analyzer -- you MUST use clang++ or the compile WILL fail")
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=clang-analyzer-*")
    set(CMAKE_C_CLANG_TIDY "clang-tidy;-checks=clang-analyzer-*")
endif()

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
# When building statically, prefer .a over .so for all find_library/find_package calls
if(NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
endif()
find_package(PkgConfig REQUIRED)
if(APPLE)
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE _homebrew_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_homebrew_prefix)
        set(ENV{PKG_CONFIG_PATH} "${_homebrew_prefix}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
endif()
# Ensure /usr/local/lib/pkgconfig is searched (for locally-installed libraries).
if(EXISTS "/usr/local/lib/pkgconfig")
    set(ENV{PKG_CONFIG_PATH} "/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
# GMP, MPFR, and flint are external pre-built system libraries.  On Linux, a
# source-built static .a typically lacks -fPIC and cannot be linked into a shared
# object (.so).  Restore the platform-default suffix order (.so first) when building
# the Python extension on Linux so the shared versions are found instead.
# On macOS, Homebrew always builds static libs with PIC so the .a works fine, and
# using the .dylib would drag in a macOS-version-specific deployment target that
# breaks the wheel's declared minimum OS version.
if(BUILD_PYTHON_EXTENSION AND UNIX AND NOT APPLE)
    set(_saved_find_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
    unset(CMAKE_FIND_LIBRARY_SUFFIXES)
endif()
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx HINTS ${GMP_LIBRARY_DIRS})
if(GMPXX_LIBRARY)
    get_property(_gmp_iface_libs TARGET PkgConfig::GMP PROPERTY INTERFACE_LINK_LIBRARIES)
    set_property(TARGET PkgConfig::GMP PROPERTY INTERFACE_LINK_LIBRARIES "${GMPXX_LIBRARY}" ${_gmp_iface_libs})
endif()
pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
pkg_check_modules(flint REQUIRED IMPORTED_TARGET flint)
if(BUILD_PYTHON_EXTENSION AND UNIX AND NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES "${_saved_find_suffixes}")
endif()

option(NOZLIB "Don't use zlib" OFF)
if(NOT NOZLIB AND NOT (NOT BUILD_SHARED_LIBS AND WIN32))
    find_package(ZLIB)
    if(ZLIB_FOUND)
        message(STATUS "OK, Found ZLIB!")
        add_compile_definitions(USE_ZLIB)
    else()
        message(STATUS "WARNING: Did not find ZLIB, gzipped file support will be disabled")
    endif()
endif()

include(FetchContent)
# FetchContent_Populate with declared details is deprecated in CMake 3.30+ (CMP0169).
# We call it intentionally to clone-only without add_subdirectory, so silence the warning.
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()


# ---- CryptoMiniSat5 ----
set(cryptominisat5_DIR "" CACHE PATH "cryptominisat5 install/build prefix (contains lib/cmake/cryptominisat5/). Auto-resolved if empty.")
if(NOT TARGET cryptominisat5)
    if(cryptominisat5_DIR)
        find_package(cryptominisat5 CONFIG REQUIRED HINTS "${cryptominisat5_DIR}")
        message(STATUS "CryptoMiniSat5: using pre-built at ${cryptominisat5_DIR}")
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../cryptominisat/CMakeLists.txt")
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../cryptominisat"
                         "${CMAKE_CURRENT_BINARY_DIR}/cryptominisat-build")
        message(STATUS "CryptoMiniSat5: using sibling repo")
    else()
        FetchContent_Declare(cryptominisat5
            GIT_REPOSITORY https://github.com/msoos/cryptominisat
            GIT_TAG        master
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(cryptominisat5)
        message(STATUS "CryptoMiniSat5: fetched from GitHub")
    endif()
endif()

# ---- Arjun ----
set(arjun_DIR "" CACHE PATH "arjun install/build prefix (contains lib/cmake/arjun/). Auto-resolved if empty.")
if(NOT TARGET arjun)
    if(arjun_DIR)
        find_package(arjun CONFIG REQUIRED HINTS "${arjun_DIR}")
        message(STATUS "Arjun: using pre-built at ${arjun_DIR}")
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../arjun/CMakeLists.txt")
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../arjun"
                         "${CMAKE_CURRENT_BINARY_DIR}/arjun-build")
        message(STATUS "Arjun: using sibling repo")
    else()
        FetchContent_Declare(arjun
            GIT_REPOSITORY https://github.com/meelgroup/arjun.git
            GIT_TAG        master
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(arjun)
        message(STATUS "Arjun: fetched from GitHub")
    endif()
endif()

# ---- ApproxMC ----
set(approxmc_DIR "" CACHE PATH "approxmc install/build prefix (contains lib/cmake/approxmc/). Auto-resolved if empty.")
if(NOT TARGET approxmc)
    if(approxmc_DIR)
        find_package(approxmc CONFIG REQUIRED HINTS "${approxmc_DIR}")
        message(STATUS "ApproxMC: using pre-built at ${approxmc_DIR}")
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../approxmc/CMakeLists.txt")
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../approxmc"
                         "${CMAKE_CURRENT_BINARY_DIR}/approxmc-build")
        message(STATUS "ApproxMC: using sibling repo")
    else()
        FetchContent_Declare(approxmc
            GIT_REPOSITORY https://github.com/meelgroup/ApproxMC.git
            GIT_TAG        master
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(approxmc)
        message(STATUS "ApproxMC: fetched from GitHub")
    endif()
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Collect all compile definitions set so far into a string for embedding
# via @COMPILE_DEFINES@ in GitSHA1.cpp.in
get_directory_property(DirDefs COMPILE_DEFINITIONS)
set(COMPILE_DEFINES "")
foreach(d ${DirDefs})
    set(COMPILE_DEFINES "${COMPILE_DEFINES} -D${d}")
endforeach()
message(STATUS "All defines at startup: ${COMPILE_DEFINES}")

macro(ganak_add_public_header LIBTARGET HEADER)
    get_target_property(EXISTING_PUBLIC_HEADERS ${LIBTARGET} PUBLIC_HEADER)
    if(EXISTING_PUBLIC_HEADERS)
        list(APPEND EXISTING_PUBLIC_HEADERS "${HEADER}")
    else()
        # Do not append to empty list
        set(EXISTING_PUBLIC_HEADERS "${HEADER}")
    endif()
    set_target_properties(
        ${LIBTARGET}
        PROPERTIES
        PUBLIC_HEADER "${EXISTING_PUBLIC_HEADERS}"
     )
endmacro()


set(DEF_INSTALL_CMAKE_DIR lib/cmake/ganak)
set(GANAK_INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
    "Installation directory for ganak CMake files")

# -----------------------------------------------------------------------------
# Add GIT version
# -----------------------------------------------------------------------------
function(SetVersionNumber PREFIX VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
  set(${PREFIX}_VERSION_MAJOR ${VERSION_MAJOR} PARENT_SCOPE)
  set(${PREFIX}_VERSION_MINOR ${VERSION_MINOR} PARENT_SCOPE)
  set(${PREFIX}_VERSION_PATCH ${VERSION_PATCH} PARENT_SCOPE)
  set(${PREFIX}_VERSION
        "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
        PARENT_SCOPE)
endfunction()

find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  message(STATUS "GIT hash found: ${GIT_SHA1}")
else()
  set(GIT_SHA1 "GIT-hash-notfound")
endif()
set(GANAK_FULL_VERSION "1.1.0")

string(REPLACE "." ";" GANAK_FULL_VERSION_LIST ${GANAK_FULL_VERSION})
SetVersionNumber("PROJECT" ${GANAK_FULL_VERSION_LIST})
message(STATUS "PROJECT_VERSION: ${PROJECT_VERSION}")
message(STATUS "PROJECT_VERSION_MAJOR: ${PROJECT_VERSION_MAJOR}")
message(STATUS "PROJECT_VERSION_MINOR: ${PROJECT_VERSION_MINOR}")
message(STATUS "PROJECT_VERSION_PATCH: ${PROJECT_VERSION_PATCH}")


# -----------------------------------------------------------------------------
# Provide an export name to be used by targets that wish to export themselves.
# -----------------------------------------------------------------------------
set(GANAK_EXPORT_NAME "ganakTargets")

add_subdirectory(src)

# -----------------------------------------------------------------------------
# Add uninstall target for makefiles
# -----------------------------------------------------------------------------
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)

add_custom_target(uninstall_ganak
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

option(ENABLE_TESTING "Enables testing" ON)
if(ENABLE_TESTING)
  find_package(Python3 COMPONENTS NumPy Interpreter Development REQUIRED)
  message(STATUS "Python 3 -- Python3_EXECUTABLE=${Python3_EXECUTABLE}")
  message(STATUS "Python 3 -- Python3_LIBRARIES=${Python3_LIBRARIES}")
  message(STATUS "Python 3 -- Python3_INCLUDE_DIRS=${Python3_INCLUDE_DIRS}")
  message(STATUS "Python 3 -- Python3_VERSION=${Python3_VERSION}")
  message(STATUS "Python 3 -- Python3_NumPy_INCLUDE_DIRS=${Python3_NumPy_INCLUDE_DIRS}")
  message(STATUS "Python 3 -- Python3_NumPy_VERSION=${Python3_NumPy_VERSION}")
  enable_testing()
  add_subdirectory(tests/cnf-files)
endif()

# -----------------------------------------------------------------------------
# Export our targets so that other CMake based projects can interface with
# the build of ganak in the build-tree
# -----------------------------------------------------------------------------
set(GANAK_TARGETS_FILENAME "ganakTargets.cmake")
set(GANAK_CONFIG_FILENAME "ganakConfig.cmake")

# Export targets
set(MY_TARGETS ganak)
export(
    TARGETS ${MY_TARGETS}
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${GANAK_TARGETS_FILENAME}"
)

# Build-tree config (absolute paths, for use without installation)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ganakConfig-buildtree.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${GANAK_CONFIG_FILENAME}"
    @ONLY
)

# Export this package to the CMake user package registry
# Now the user can just use find_package(ganak) on their system
export(PACKAGE ganak)

# Install-tree config and version files (relocatable via CMakePackageConfigHelpers)
set(GANAK_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(GANAK_VERSION_FILENAME "ganakConfigVersion.cmake")

include(CMakePackageConfigHelpers)
configure_package_config_file(ganakConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${GANAK_CONFIG_FILENAME}"
    INSTALL_DESTINATION "${GANAK_INSTALL_CMAKE_DIR}"
    PATH_VARS GANAK_INSTALL_INCLUDEDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${GANAK_VERSION_FILENAME}"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${GANAK_CONFIG_FILENAME}"
    "${CMAKE_CURRENT_BINARY_DIR}/${GANAK_VERSION_FILENAME}"
    DESTINATION "${GANAK_INSTALL_CMAKE_DIR}"
)

# Install the export set for use with the install-tree
install(
    EXPORT ${GANAK_EXPORT_NAME}
    DESTINATION "${GANAK_INSTALL_CMAKE_DIR}"
)
