cmake_minimum_required(VERSION 2.6.1)
project(MetisTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lmetis")


set(SOURCE_FILES
    src/metis.h
    src/MetisNode.cpp
    src/MetisNode.h
    src/MetisTest.cpp)

add_executable(MetisTest ${SOURCE_FILES})