add_executable(test_mpqi
    ${PROJECT_SOURCE_DIR}/src/mpqi.c
    ${PROJECT_SOURCE_DIR}/tests/test_mpqi.cpp
)
target_include_directories(test_mpqi PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
    ${GMP_INCLUDE_DIR}
    ${MPFI_INCLUDES}
)
target_link_libraries(test_mpqi PRIVATE
    ${MPFI_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
)
set_target_properties(test_mpqi PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)
add_test(NAME mpqi_unit_tests
    COMMAND test_mpqi
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)
