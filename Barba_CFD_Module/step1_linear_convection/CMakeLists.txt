cmake_minimum_required(VERSION 3.21.1)
project(step1_linear_convection)


set(CMAKE_CXX_STANDARD 17)


include_directories(src)

add_executable(step1_linear_convection src/main.cpp)
