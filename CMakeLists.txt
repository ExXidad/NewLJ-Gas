cmake_minimum_required(VERSION 3.17)
project(newLJ_gas)

set(CMAKE_CXX_STANDARD 17)

add_executable(newLJ_gas Solver.cpp Solver.h main.cpp)
target_link_libraries(newLJ_gas -lgsl)