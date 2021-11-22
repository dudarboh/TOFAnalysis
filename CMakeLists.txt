cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(TOFAnalysis)

find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)

# Remove this ugly warning
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-effc++")
message(STATUS "C++ flags are set to: " ${CMAKE_CXX_FLAGS})


include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCE_DIR}/src/TOFAnaUtils.cpp ${PROJECT_SOURCE_DIR}/src/TOFAnalysis.cpp)


### DEPENDENCIES ###
find_package(Marlin REQUIRED)
find_package(MarlinUtil REQUIRED)
find_package(MarlinTrk REQUIRED)

# this is broken as fck
# find_package (MarlinTrk PATHS /afs/desy.de/user/d/dudarboh/iLCSoft/MarlinTrk/ NO_DEFAULT_PATH)

FOREACH(pkg Marlin MarlinUtil MarlinTrk CED DD4hep)
   include_directories(${${pkg}_INCLUDE_DIRS})
   target_link_libraries(${PROJECT_NAME} ${${pkg}_LIBRARIES})
ENDFOREACH()

# include_directories(/afs/desy.de/user/d/dudarboh/iLCSoft/MarlinTrk/include)
# target_link_libraries(${PROJECT_NAME} /afs/desy.de/user/d/dudarboh/iLCSoft/MarlinTrk/lib/libMarlinTrk.so)

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)
