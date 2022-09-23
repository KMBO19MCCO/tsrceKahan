cmake_minimum_required(VERSION 3.23)
project(tsrceKahan)

set(CMAKE_CXX_STANDARD 20)
add_subdirectory(excerpt)
include_directories(${tsrceKahan_SOURCE_DIR}/excerpt/inc)
add_executable(tsrceKahan main.cpp)
target_link_libraries(tsrceKahan excerpt)
