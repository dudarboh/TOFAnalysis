if(NOT WIN32)
  string(ASCII 27 Esc)
  set(Reset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12)

project(TOFAnalysis VERSION 1.0)

# Include include folder so source files can find headers
# Write all source files in the DIRECTORY to the VARIABLE
# Add shared library with the NAME based on those source files
include_directories(${PROJECT_SOURCE_DIR}/include)

aux_source_directory(${PROJECT_SOURCE_DIR}/src source_files)
add_library(${PROJECT_NAME} SHARED ${source_files})


### DEPENDENCIES ###
include($ENV{ILCSOFT}/ILCSoft.cmake)
find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)
# Add compler optimization flag to the default ones
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-effc++")

find_package(Marlin 1.15 REQUIRED)
find_package(MarlinUtil 1.4 REQUIRED)
FOREACH(pkg Marlin MarlinUtil)
    include_directories(${${pkg}_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${${pkg}_LIBRARIES})
ENDFOREACH()


# Print compiler flags for compiling C++ sources with a nice coloring
message(STATUS "${Yellow}CMAKE_CXX_FLAGS:  ${CMAKE_CXX_FLAGS} ${Reset}")

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)
