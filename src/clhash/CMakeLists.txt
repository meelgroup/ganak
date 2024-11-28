SET(clhash_SOURCES clhash.cpp)
ADD_LIBRARY(clhash ${clhash_SOURCES})

# set(clhash_public_headers clhash.h)

# set_target_properties(clshash PROPERTIES
#   #     PUBLIC_HEADER "${clhash_public_headers}"
#     VERSION 1.0
#     SOVERSION 1.0
# )

install(TARGETS clhash
    EXPORT ${GANAK_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ganak"
)

