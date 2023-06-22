##
#  CMake script for the step-<XYZ> tutorial program:
##

# Set the name of the project and target:
SET(TARGET "cook_membrane")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc
  # You can specify additional files here!
  )

SET(CLEAN_UP_FILES
  # a custom list of globs, e.g. *.log *.vtk
  *.vtk
)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)

# seek eigen
find_package(Eigen3 REQUIRED)

# include eigen in dealii
include_directories(${EIGEN3_INCLUDE_DIR})

FIND_PACKAGE(deal.II 9.3.0 QUIET 
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.ucd *.d2 *.vtu *.pvtu)
MESSAGE(STATUS "deal.II_DIR: ${deal.II_DIR}")
MESSAGE(STATUS "DEAL_II_DIR: ${DEAL_II_DIR}")

