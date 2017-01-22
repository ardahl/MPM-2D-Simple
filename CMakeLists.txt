cmake_minimum_required (VERSION 2.6)
project(MPM)

###
# Add build type that uses O3 but includes debugging flag
###
SET( CMAKE_CXX_FLAGS_OPTDEBUG "-Wall -Wno-deprecated-declarations -O3" CACHE STRING
    "Flags used by the C++ compiler during optimized debugging builds."
    FORCE )
SET( CMAKE_C_FLAGS_OPTDEBUG "-Wall -Wno-deprecated-declarations -pedantic -O3" CACHE STRING
    "Flags used by the C compiler during optimized debugging builds."
    FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_OPTDEBUG
    "-Wl" CACHE STRING
    "Flags used for linking binaries during maintainer builds."
    FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_OPTDEBUG
    "-Wl" CACHE STRING
    "Flags used by the shared libraries linker during maintainer builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_OPTDEBUG
    CMAKE_C_FLAGS_OPTDEBUG
    CMAKE_EXE_LINKER_FLAGS_OPTDEBUG
    CMAKE_SHARED_LINKER_FLAGS_OPTDEBUG )

###
# User Options
###
# Options are Release, Debug, and Optdebug
set( CMAKE_BUILD_TYPE Release )
#Source directory
set( SRC "src")
#mpm
option(MPM_GRAPHICS "Compile the frame-grabbing mpm viewer" ON)
#mpm_sim
option(MPM_SIM "Compile the simulation that just dumps data to file for later viewing" ON)
#Haven't written these parts yet
#mpm_sim viewer
option(MPM_VIEW "Compile the OpenGL viewer for the data files" ON)
#mpm_sim viewer use opencv
option(MPM_VIEW_CV "Viewer compiles with OpenCV to save frames to video" ON)
#sampler
option(SAMPLER "Simple template code for creating input particle files" ON)

###
# build details
###
set( CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR} )
#set( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall" )
#set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-sign-compare")
#set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-sign-compare")

######
# For the simulation
######
#Eigen library
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(PARTIO_LIBRARY "partio -L/usr/local/lib")
#set(PARTIO_LIBRARY "partio")
#zlib library
set(ZLIB "z")

######
# For the viewer
######
#OpenGL and GLUT (Eventually make this optional)
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OpenGL not found!")
endif(NOT OPENGL_FOUND)

#OpenCV
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
if(NOT OpenCV_FOUND)
    message(ERROR " OpenCV not found!")
endif(NOT OpenCV_FOUND)

###
# Sources
###
set (SOURCES
	${SRC}/grid.cpp
	${SRC}/world.cpp
	${SRC}/jsoncpp.cpp
)

#Simulation and Viewer
add_library(mpmlib ${SOURCES})

if(MPM_GRAPHICS)
    add_executable(mpm_graphics ${SRC}/main_graphics.cpp)
    target_link_libraries(mpm_graphics ${PARTIO_LIBRARY} ${ZLIB} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OpenCV_LIBS} mpmlib)
endif(MPM_GRAPHICS)

#Standalone sim
if(MPM_SIM)
    add_executable(mpm ${SRC}/main.cpp)
    target_link_libraries(mpm ${PARTIO_LIBRARY} ${ZLIB} mpmlib)
endif(MPM_SIM)

#Data viewer
if(MPM_VIEW)
    add_executable(viewer ${SRC}/viewer.cpp)
    target_link_libraries(viewer ${PARTIO_LIBRARY} ${ZLIB} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OpenCV_LIBS} mpmlib)
endif(MPM_VIEW)

#Simple Sampling Tool
if(SAMPLER)
    add_executable(sampler ${SRC}/sampler.cpp)
    target_link_libraries(sampler ${PARTIO_LIBRARY} ${ZLIB} mpmlib)
endif(SAMPLER)

