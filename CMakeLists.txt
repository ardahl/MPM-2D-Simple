cmake_minimum_required (VERSION 2.6)
project(MPM)

###
# User Options
###
set( CMAKE_BUILD_TYPE Debug )
#Source directory
set( SRC "src")
#mpm
option(MPM_GRAPHICS "Compile the frame-grabbing mpm viewer" ON)
#mpm_dump 
option(MPM "Compile the simulation that just dumps data to file for later viewing" ON)
#Haven't written these parts yet
#mpm_dump viewer
option(MPM_DUMP_VIEW "Compile the OpenGL viewer for the data files" ON)
#mpm_dump viewer use opencv
option(MPM_DUMP_VIEW_CV "Viewer compiles with OpenCV to save frames to video" ON)

###
# build details
###
set( CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR} )
#set( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -DNDEBUG")

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
	${SRC}/configparser.cpp
	${SRC}/grid.cpp
	${SRC}/world.cpp
	${SRC}/jsoncpp.cpp
)

#Simulation and Viewer
if(MPM_GRAPHICS)
    add_executable(mpm ${SOURCES} ${SRC}/main_graphics.cpp)
    target_link_libraries(mpm ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OpenCV_LIBS} )
endif(MPM_GRAPHICS)
#Data Dumper
if(MPM)
    add_executable(mpm ${SOURCES} ${SRC}/main.cpp)
endif(MPM_DUMP)
#Data viewer
if(MPM_DUMP_VIEW)
    add_executable(viewer ${SOURCES} ${SRC}/viewer.cpp)
    target_link_libraries(viewer ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
endif(MPM_DUMP_VIEW)
