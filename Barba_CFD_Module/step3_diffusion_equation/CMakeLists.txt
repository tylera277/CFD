cmake_minimum_required(VERSION 3.21.1)
project(step3_diffusion_equation)


set(CMAKE_CXX_STANDARD 17)


find_package (Eigen3 3.3 REQUIRED)



#add_executable(step1_linear_convection src/main.cpp)

#target_link_libraries(step1_linear_convection PUBLIC supplemental_files)
#target_link_libraries(step1_linear_convection Eigen3::Eigen)


add_executable(step3_diffusion_equation "src/main.cpp")

target_link_libraries(step3_diffusion_equation PUBLIC Eigen3::Eigen)
