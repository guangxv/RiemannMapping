cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Set and enforce C++-11 flags
set( CMAKE_CXX_STANDARD_REQUIRED TRUE )
set( CMAKE_CXX_STANDARD 11 )

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.0" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.9" OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "4.9")
                message(WARNING "Your version of GCC contains an optimizer bug. Please verify that you do not use -O3!")
                string(REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE_NEW "${CMAKE_CXX_FLAGS_RELEASE}")
                set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE_NEW}" CACHE STRING "" FORCE)
        endif()
endif()


# SET( CMAKE_BUILD_TYPE Debug CACHE STRING "set build type to debug" ) #Release
SET(ProjectName RiemannMapping)

PROJECT(${ProjectName})

SET(OpenMesh_DIR /usr/local/build/OpenMesh)

FIND_PACKAGE(OpenMesh REQUIRED)
FIND_PACKAGE(GLUT)
FIND_PACKAGE(OpenGL)
#set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#INCLUDE(ACGCommon)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Where to put the executables") 
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Where to put the libraries") 

INCLUDE_DIRECTORIES(
  ${OpenMesh_INCLUDE_DIRS}
  /usr/local/include/eigen3/
  ./include/
  )

ADD_SUBDIRECTORY( src )
