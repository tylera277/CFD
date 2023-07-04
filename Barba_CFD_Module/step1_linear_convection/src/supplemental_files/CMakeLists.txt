

add_library(supplemental_files STATIC supplemental_files.cpp supplemental_files.hpp)

target_link_libraries(supplemental_files PRIVATE Eigen3::Eigen)
