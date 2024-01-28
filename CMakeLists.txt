cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_CX_STANDARD 17)
set(CMAKE_CX_STANDARD_REQUIRED ON)
set(CMAKE_CX_EXTENSIONS OFF)

project("cacheperf" CXX)

add_executable(cacheperf src/cacheperf.cpp)
