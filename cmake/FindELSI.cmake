###############################################################################
# - Find ELSI
# Find the native ELSI headers and libraries.
#
#  ELSI_FOUND        - True if libelsi is found.
#  ELSI_LIBRARIES    - List of libraries when using libyaml
#  ELSI_INCLUDE_DIR - Where to find ELSI headers.
#

find_path(ELSI_INCLUDE_DIR
    elsi.h
    HINTS ${ELSI_DIR}
    PATH_SUFFIXES "include" 
    )
find_library(ELSI_LIBRARY
	NAMES elsi
	HINTS ${ELSI_DIR}
	PATH_SUFFIXES "lib"
	)

string(REGEX REPLACE "libelsi.a" "" ELSI_LIBRARY_DIR ${ELSI_LIBRARY})
# Handle the QUIET and REQUIRED arguments and
# set ELSI_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ELSI DEFAULT_MSG ELSI_LIBRARY_DIR ELSI_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(ELSI_FOUND)
    set(ELSI_LIBRARIES "-Wl,--start-group -L${ELSI_LIBRARY_DIR} -lelsi -lelpa -lfortjson -lMatrixSwitch -lNTPoly -lOMM -lpexsi -lptscotch -lptscotcherr -lptscotchparmetis -lscotch -lscotcherr -lscotchmetis -lsuperlu_dist -Wl,--end-group")
    # set(ELSI_LIBRARIES ${ELSI_LIBRARY})
    message(STATUS "ELSI_INCLUDE_DIR is ${ELSI_INCLUDE_DIR}")
    message(STATUS "ELSI_LIBRARY is ${ELSI_LIBRARY}")
    message(STATUS "ELSI_LIBRARY_DIR is ${ELSI_LIBRARY_DIR}")
    message(STATUS "ELSI_LIBRARIES is ${ELSI_LIBRARIES}")
    #set(ELSI_INCLUDE_DIR ${ELSI_INCLUDE_DIR})

    if(NOT TARGET ELSI::ELSI)
        add_library(ELSI::ELSI UNKNOWN IMPORTED)
        set_target_properties(ELSI::ELSI PROPERTIES
           IMPORTED_LINK_INTERFACE_LANGUAGES "C"
           IMPORTED_LOCATION "${ELSI_LIBRARY}"
           INTERFACE_INCLUDE_DIRECTORIES "${ELSI_INCLUDE_DIR}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ELSI_INCLUDE_DIR})

mark_as_advanced(ELSI_INCLUDE_DIR ELSI_LIBRARY_DIR)
