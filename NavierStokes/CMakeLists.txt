cmake_minimum_required(VERSION 3.21.1)
project(NavierStokes)

set(CMAKE_CXX_STANDARD 17)

enable_testing()
add_subdirectory(test)

find_package (Eigen3 3.3 REQUIRED)

add_executable(NavierStokes "src/main.cpp")

add_subdirectory("src/simulation")


target_link_libraries(NavierStokes PUBLIC Eigen3::Eigen)
target_link_libraries(NavierStokes PUBLIC simulation)

