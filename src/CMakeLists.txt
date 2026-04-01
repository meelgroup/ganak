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

if(NOT WIN32)
    #add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    add_cxx_flag_if_supported("-Wduplicated-cond")
    add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

set(SOURCES
    comp_analyzer.cpp
    comp_types/comp_archetype.cpp
    comp_management.cpp
    counter.cpp
    outer_counter.cpp
    statistics.cpp
    clauseallocator.cpp
    buddy-2.4/src/bddio.c
    buddy-2.4/src/bddop.c
    buddy-2.4/src/bvec.c
    buddy-2.4/src/cache.c
    buddy-2.4/src/fdd.c
    buddy-2.4/src/imatrix.c
    buddy-2.4/src/kernel.c
    buddy-2.4/src/pairs.c
    buddy-2.4/src/prime.c
    buddy-2.4/src/reorder.c
    buddy-2.4/src/tree.c
    buddy-2.4/src/cppext.cxx
    ganak.cpp
    comp_cache_if.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

add_library(ganak ${SOURCES})

target_include_directories(ganak
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/buddy-2.4/
    ${CMAKE_CURRENT_SOURCE_DIR}/buddy-2.4/src/
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ganak
  PRIVATE
    approxmc
    treedecomp
    arjun
  PUBLIC
    Threads::Threads
    cryptominisat5
)

set_target_properties(ganak PROPERTIES
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

add_executable(ganak-bin
  mpqi.c
  main.cpp
)
add_executable(example-bin example.cpp)

target_link_libraries(ganak-bin
    PRIVATE
      ganak
      arjun
      PkgConfig::flint
      PkgConfig::MPFR
      PkgConfig::GMP
      ${MPFI_LIBRARIES}
      approxmc
)
target_link_libraries(example-bin
    PRIVATE
      ganak
      arjun
      PkgConfig::GMP
)

target_include_directories(ganak-bin PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPFI_INCLUDES}
)
target_include_directories(example-bin PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ZLIB_FOUND)
    target_include_directories(ganak-bin PRIVATE ${ZLIB_INCLUDE_DIR})
    target_link_libraries(ganak-bin PUBLIC ${ZLIB_LIBRARY})
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set_target_properties(ganak-bin PROPERTIES
    OUTPUT_NAME ganak
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
  target_link_options(ganak-bin PRIVATE
    "SHELL:-s WASM=1"
    "SHELL:-s ALLOW_MEMORY_GROWTH=1"
    "SHELL:-s EXPORTED_RUNTIME_METHODS='[\"callMain\", \"ccall\", \"cwrap\", \"FS\", \"print\"]'"
    "SHELL:-s FORCE_FILESYSTEM=1"
  )
else()
  set_target_properties(ganak-bin PROPERTIES
    OUTPUT_NAME ganak
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
endif()

set_target_properties(example-bin PROPERTIES
  OUTPUT_NAME example
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/ganak.hpp)
ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/counter_config.hpp)
ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/lit.hpp)

# Ensure headers from sub-builds are copied before we compile against them.
foreach(_dep_copy_target
        CopyPublicHeaders_cryptominisat5
        CopyPublicHeaders_arjun
        CopyPublicHeaders_approxmc
        CopyPublicHeaders_treedecomp
        CopyPublicHeaders_sbva
        CopyPublicHeaders)
    if(TARGET ${_dep_copy_target})
        foreach(_tgt ganak ganak-bin example-bin)
            add_dependencies(${_tgt} ${_dep_copy_target})
        endforeach()
    endif()
endforeach()

# Public headers are fully known now; the CopyPublicHeaders target below reads
# this property, so it must be set before that target is created.
get_target_property(ganak_public_headers ganak PUBLIC_HEADER)

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The ganakConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/ganak")
add_custom_target(CopyPublicHeaders_ganak ALL)
foreach(public_header ${ganak_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders_ganak PRE_BUILD
                       COMMAND ${CMAKE_COMMAND} -E make_directory
                               "${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E echo
                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E
                           copy_if_different
                           ${public_header}
                           "${HEADER_DEST}"
                      )
endforeach()

foreach(_tgt ganak-bin example-bin)
    add_dependencies(${_tgt} CopyPublicHeaders_ganak)
endforeach()

install(TARGETS ganak
    EXPORT ${GANAK_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ganak"
)

install(TARGETS ganak-bin
    EXPORT ${GANAK_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ---------------------------------------------------------------
# Python extension module (pyganak)
# Activated by scikit-build-core via -DBUILD_PYTHON_EXTENSION=ON.
# Requires CMake >= 3.18 (python_add_library).
# (option defined early in root CMakeLists.txt so GMP/MPFR/flint discovery can use it)
# ---------------------------------------------------------------
if(BUILD_PYTHON_EXTENSION)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    python_add_library(pyganak MODULE WITH_SOABI
        ${PROJECT_SOURCE_DIR}/python/src/pyganak.cpp
    )

    target_compile_definitions(pyganak PRIVATE
        PYGANAK_VERSION="${PROJECT_VERSION}"
    )

    target_include_directories(pyganak PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
    )

    target_link_libraries(pyganak PRIVATE
        ganak
        arjun
        PkgConfig::GMP
    )

    target_compile_features(pyganak PRIVATE cxx_std_20)

    install(TARGETS pyganak DESTINATION . COMPONENT python)
endif()
