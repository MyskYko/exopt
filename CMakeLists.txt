cmake_minimum_required(VERSION 3.8)
project(exopt CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_DEBUG "Build debug executable for the project" OFF)

add_compile_options(-g -O3 -Wall -Wextra)

add_subdirectory(lib)

file(GLOB LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
add_library(exopt_lib STATIC ${LIB_SOURCES})
target_include_directories(exopt_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(exopt_lib PUBLIC aig kissat argparse cadical)

add_executable(exopt ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_include_directories(exopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(exopt PUBLIC exopt_lib)

if(BUILD_DEBUG)
	add_executable(exopt_debug ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${LIB_SOURCES})
	target_include_directories(exopt_debug PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
	target_link_libraries(exopt_debug PUBLIC aig kissat argparse cadical)
	target_compile_options(exopt_debug PRIVATE -DDEBUG)
endif()
