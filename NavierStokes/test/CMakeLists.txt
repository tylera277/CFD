set(CATCH_INCLUDE_DIR </catch.hpp>)

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})



#Adding test subdirectories

#add_subdirectory(test_main)
add_subdirectory(test_simulation)
