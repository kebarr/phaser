cmake_minimum_required(VERSION 3.3)
project(phaser)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.cpp src/graph.h src/graph.cpp src/haplotype_scorer.cpp src/haplotype_scorer.h)
add_executable(phaser ${SOURCE_FILES})
