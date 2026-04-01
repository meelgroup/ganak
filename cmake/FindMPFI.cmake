# Try to find the MPFI library
# See https://gitlab.inria.fr/mpfi/mpfi
#
# Once done this will define
#
#  MPFI_FOUND - system has MPFI lib
#  MPFI_INCLUDES - the MPFI include directory
#  MPFI_LIBRARIES - the MPFI library

# Set MPFI_INCLUDES

find_path(MPFI_INCLUDES
  NAMES
  mpfi.h
  PATHS
  $ENV{GMPDIR}
  ${INCLUDE_INSTALL_DIR}
)

# Set MPFI_LIBRARIES

find_library(MPFI_LIBRARIES mpfi PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

# Epilogue

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFI DEFAULT_MSG
                                  MPFI_INCLUDES MPFI_LIBRARIES)
mark_as_advanced(MPFI_INCLUDES MPFI_LIBRARIES)
