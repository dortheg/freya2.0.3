#----------------------------------------------------------------------------
#
# Jerome Verbeke
#
# Created 08/26/2014
#

# Setup the project
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

PROJECT(events CXX)

# Locate source and headers for this project
SET(CPPSRC ${PROJECT_SOURCE_DIR}/events.cpp)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/../../../include)

# Locate libFission
FIND_LIBRARY(LFISSION Fission
             PATHS ${PROJECT_SOURCE_DIR}/../../../lib
             DOC "LLNL Fission Library"
             NO_DEFAULT_PATH
            )
IF(NOT LFISSION)
  MESSAGE(FATAL_ERROR "Could not find LLNL Fission Library libFission")
ENDIF()

# Build executable
ADD_EXECUTABLE(events ${CPPSRC})
TARGET_LINK_LIBRARIES(events ${LFISSION})
ADD_DEFINITIONS(-DFREYA)
