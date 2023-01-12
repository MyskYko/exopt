cmake_minimum_required(VERSION 3.8)
project(exopt CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-g -O3 -Wall -Wextra)

add_subdirectory(lib)

file(GLOB FILENAMES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(exopt ${FILENAMES})
target_include_directories(exopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(exopt PUBLIC aig kissat argparse cadical)

add_executable(exopt_debug ${FILENAMES})
target_include_directories(exopt_debug PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(exopt_debug PUBLIC aig kissat argparse cadical)
target_compile_options(exopt_debug PRIVATE -DDEBUG)
