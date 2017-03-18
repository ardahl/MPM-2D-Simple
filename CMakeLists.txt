cmake_minimum_required (VERSION 2.6)
project(MPM)

###
# User Options
###
# Options are Release, Debug
set( CMAKE_BUILD_TYPE Release )
#Uncomment line below to still use O3 optimization with debugging
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations" )

#Source directory
set( SRC "src")
#mpm graphics
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
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")	#Adds the cmake directory for the FindEigen3.cmake file

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

#set(PARTIO_LIBRARY "partio -L/usr/local/lib")
find_library(PARTIO_LIBRARY partio)
#zlib library
set(ZLIB "z")

######
# For the viewer
######
#OpenGL and GLUT (Eventually make this optional)
find_package(GLUT QUIET)
if(NOT GLUT_FOUND)
    message(STATUS "GLUT not found! Turning off mpm_graphics and viewer")
    set(MPM_GRAPHICS OFF)
    set(MPM_VIEW OFF)
else(NOT GLUT_FOUND)
    include_directories(${GLUT_INCLUDE_DIRS})
    link_directories(${GLUT_LIBRARY_DIRS})
    add_definitions(${GLUT_DEFINITIONS})
endif(NOT GLUT_FOUND)

find_package(OpenGL QUIET)
if(NOT OPENGL_FOUND)
    message(STATUS "OpenGL not found! Turning off mpm_graphics and viewer")
    set(MPM_GRAPHICS OFF)
    set(MPM_VIEW OFF)
else(NOT OPENGL_FOUND)
    include_directories(${OpenGL_INCLUDE_DIRS})
    link_directories(${OpenGL_LIBRARY_DIRS})
    add_definitions(${OpenGL_DEFINITIONS})
endif(NOT OPENGL_FOUND)

#OpenCV
find_package( OpenCV QUIET )
include_directories( ${OpenCV_INCLUDE_DIRS} )
if(NOT OpenCV_FOUND)
    #message(ERROR " OpenCV not found!")
    message(STATUS "OpenCV not found! Turning off mpm_graphics and viewer")
    set(MPM_GRAPHICS OFF)
    set(MPM_VIEW OFF)
endif(NOT OpenCV_FOUND)

###
# Sources
###
set (SOURCES
	${SRC}/world.cpp
	${SRC}/jsoncpp.cpp
)

#Simulation and Viewer
add_library(mpmlib ${SOURCES})

if(MPM_GRAPHICS)
    add_executable(mpm_graphics ${SRC}/main_graphics.cpp)
    target_link_libraries(mpm_graphics mpmlib ${PARTIO_LIBRARY} ${ZLIB} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OpenCV_LIBS})
endif(MPM_GRAPHICS)

#Standalone sim
if(MPM_SIM)
    add_executable(mpm ${SRC}/main.cpp)
    target_link_libraries(mpm mpmlib ${PARTIO_LIBRARY} ${ZLIB})
endif(MPM_SIM)

#Data viewer
if(MPM_VIEW)
    add_executable(viewer ${SRC}/viewer.cpp)
    target_link_libraries(viewer mpmlib ${PARTIO_LIBRARY} ${ZLIB} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OpenCV_LIBS})
endif(MPM_VIEW)

#Simple Sampling Tool
if(SAMPLER)
    add_executable(sampler ${SRC}/sampler.cpp)
    target_link_libraries(sampler mpmlib ${PARTIO_LIBRARY} ${ZLIB})
endif(SAMPLER)

