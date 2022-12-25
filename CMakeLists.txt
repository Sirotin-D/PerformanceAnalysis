cmake_minimum_required(VERSION 3.22)
project(PerfomanceAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(PerfomanceAnalysis main.cpp SlowSolution.h OptimizedSolution.h StartConditions.h)
