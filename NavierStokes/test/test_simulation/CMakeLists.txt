add_executable(test_simulation)

target_sources(test_simulation PRIVATE 
test_simulation.cpp 
../../src/simulation/simulation.cpp
../../src/simulation/simulation.hpp 
)


target_link_libraries(test_simulation Catch)