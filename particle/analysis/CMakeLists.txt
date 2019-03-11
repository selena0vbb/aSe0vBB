# CMakeLists.txt for aSe0vBB analysis

#----------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(aseanalysis)
set(CMAKE_CXX_FLAGS "-std=c++11")
message(${CMAKE_CXX_FLAGS})
#----------------------------------------------------------
# Location of root fies
list(APPEND CMAKE_PREFIX_PATH $ROOTSYS/etc/cmake/FindROOT.cmake})
find_package(ROOT REQUIRED)

#----------------------------------------------------------
# Include root and project directories
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------
# Locate source and header files for project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


#----------------------------------------------------------
# Create an executable and link to root libraries
add_executable(aseanalysis main.cxx ${sources} ${headers})
target_link_libraries(aseanalysis ${ROOT_LIBRARIES})
