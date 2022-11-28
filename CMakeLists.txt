cmake_minimum_required(VERSION 3.23)
project(tsrceKahan)

set(CMAKE_CXX_STANDARD 23)


find_package(OpenMP REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
add_executable(${CMAKE_PROJECT_NAME} main.cpp)
if(WIN32)
    set(LIB ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll)
else()
    set(LIB ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.so)
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} OpenMP::OpenMP_CXX ${LIB})
file(COPY ${LIB} DESTINATION ${PROJECT_BINARY_DIR})
