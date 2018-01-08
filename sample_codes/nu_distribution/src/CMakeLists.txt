#----------------------------------------------------------------------------
#
# Jerome Verbeke
#
# Created 07/28/2014
#

# Setup the project
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

PROJECT(nudist CXX)

# Performance profile build
set(PROF false CACHE BOOL "performance profile of sample code")

# Locate source and headers for this project
SET(CPPSRC ${PROJECT_SOURCE_DIR}/nu_dist.cpp)
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

# Add custom performance profiler build
#
if(PROF)
  set (CMAKE_CXX_FLAGS "-pg")
  if(FORT AND CMAKE_Fortran_COMPILER)
    get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    if (${Fortran_COMPILER_NAME} MATCHES "gfortran.*")
      set (CMAKE_Fortran_FLAGS "-pg")
      message(STATUS "FORTRAN compiler is gfortran, adding full '-pg' flags")
    endif()
  endif()
endif()

# Build executable
ADD_EXECUTABLE(nudist ${CPPSRC})
TARGET_LINK_LIBRARIES(nudist ${LFISSION})
ADD_DEFINITIONS(-DFREYA)
