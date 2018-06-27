# CMakeList.txt for Selenium Detctor

#------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(selenium)

#------------------------------------------------------------------
# Find Geant4 package
# Default is to build with all available UI and Visualization drivers
# By setting WITH_GEANT4_UIVIS to OFF in the command line, this will build a batch mode only executable

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis Drivers" ON)
if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

#------------------------------------------------------------------
# Setup Geant4 include directories
# Setup Selenium include directories

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#------------------------------------------------------------------
# Locate sources and headers for selenium project

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#------------------------------------------------------------------
# Add executable and link it to the Geant4 libraries

add_executable(selenium selenium.cc ${sources} ${headers})
target_link_libraries(selenium ${Geant4_LIBRARIES})

#------------------------------------------------------------------
# Copy useful files from source to build directories

set(SELENIUM_SCRIPTS
  122.mac
  136.mac
  plotHist.c
  plotCo57.c
)

foreach(_script ${SELENIUM_SCRIPTS})
	configure_file(
		${PROJECT_SOURCE_DIR}/${_script}
		${PROJECT_BINARY_DIR}/${_script}
		COPYONLY
	)
endforeach()


