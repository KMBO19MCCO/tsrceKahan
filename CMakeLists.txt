cmake_minimum_required(VERSION 3.23)
project(tsrceKahan)

set(CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/excerpt)
link_directories(${PROJECT_SOURCE_DIR}/excerpt)
add_executable(${CMAKE_PROJECT_NAME} main.cpp ${PROJECT_SOURCE_DIR}/excerpt/excerpt.cpp ${PROJECT_SOURCE_DIR}/tsrceKahan.cpp ${PROJECT_SOURCE_DIR}/tsrceKahan.h)
target_link_libraries(${CMAKE_PROJECT_NAME} OpenMP::OpenMP_CXX)