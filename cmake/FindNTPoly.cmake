###############################################################################
# - Find NTPoly
# Find NTPoly solver
#
#  NTPoly_FOUND          - True if NTPoly is found.
#  NTPoly_INCLUDE_DIR    - Where to find NTPoly headers.
#  NTPoly_LIBRARY        - NTPoly library.

find_path(NTPoly_INCLUDE_DIR
    NAMES PSMatrix.h
    HINTS ${NTPoly_DIR}
    PATH_SUFFIXES "include"
)

find_library(NTPoly_LIBRARY
    NAMES NTPoly
    HINTS ${NTPoly_DIR}
    PATH_SUFFIXES "lib"
)


# Handle the QUIET and REQUIRED arguments and
# set NTPoly_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NTPoly DEFAULT_MSG NTPoly_LIBRARY NTPoly_INCLUDE_DIR)


# Copy the results to the output variables and target.
mark_as_advanced(NTPoly_LIBRARY NTPoly_INCLUDE_DIR)

