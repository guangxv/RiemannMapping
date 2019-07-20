CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

SET(CMAKE_MODULE_PATH /usr/local/src/OpenMesh/cmake/ )
include (ACGCommon)

SET(targetName SphericalMapping)
PROJECT(${targetName})

# collect all header and source files
acg_append_files (headers "*.h" .)
acg_append_files (sources "*.cpp" .)

INCLUDE_DIRECTORIES(
 "/usr/local/include/OpenMesh/"
 "/home/liguangxu/Projects/Include/"
 )
acg_add_executable (${targetName} ${headers} ${sources})

LINK_DIRECTORIES("/usr/local/lib/")
TARGET_LINK_LIBRARIES (${targetName}
  -lOpenMeshCore 
  -lOpenMeshTools
)
