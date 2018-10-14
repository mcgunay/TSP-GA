cmake_minimum_required(VERSION 3.10)
project(TSP-GA C)

set(CMAKE_C_STANDARD 99)

set(test1_SOURCES main.c)
set(test2_SOURCES run_experiments.c)
add_executable(test1 ${test1_SOURCES})
add_executable(test2 ${test2_SOURCES})

