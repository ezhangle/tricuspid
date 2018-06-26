project(tricuspid)
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Gui REQUIRED) # It's required because one of the include files includes QPen.
find_package(Bezitopo REQUIRED)

add_executable(tricuspid tricuspid.cpp)
target_link_libraries(tricuspid ${Bezitopo_LIBRARY} Qt5::Widgets Qt5::Core)

set(TRICUSPID_MAJOR_VERSION 0)
set(TRICUSPID_MINOR_VERSION 1)
set(TRICUSPID_PATCH_VERSION 3)
set(TRICUSPID_VERSION ${TRICUSPID_MAJOR_VERSION}.${TRICUSPID_MINOR_VERSION}.${TRICUSPID_PATCH_VERSION})
set(COPY_YEAR 2018)

include_directories(${PROJECT_BINARY_DIR})
include_directories(${Bezitopo_INCLUDE_DIR})
configure_file (config.h.in config.h)

set(CPACK_PACKAGE_VERSION_MAJOR ${TRICUSPID_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${TRICUSPID_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${TRICUSPID_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)
