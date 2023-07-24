
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "simulation/simulation.hpp"
#include "simulation_parameters.hpp"


using namespace navier_stokes;

int main()
{

    std::string outputFileName = "../outputData/";

    double current_time_step = 1;

    // I need to think of a better way to pass all of these parameters.
    // My current motivation for doing this is to be able to more easily test each unit of a class
    // by being able to create an instance of the class.
    simulation::Simulation sim(
        sim_params::imin, sim_params::imax, sim_params::jmin, sim_params::jmax,
        sim_params::dxi, sim_params::dyi,
        sim_params::viscosity, sim_params::density, sim_params::time_step_size
    );


    sim.CreateLaplacianV2();

    sim.PrintLaplacian("../outputData/22jul_laplacianv2.csv");

    while(current_time_step < sim_params::total_simulation_time)
    {   


        // Applying boundary conditions to the problem domain.
        // Values are entered directly into function simulation.cpp definition,
        // I will change this in the future
        sim.ApplyBoundaryConditions("generic");

        

        // Compute an intermediate velocity by solving the momentum eqn.
        // but omitting the effect of pressure
        sim.predictor_step();

        sim.PrintVelocities(outputFileName + "22jul_star_velocities.csv");

        sim.corrector_step();



        current_time_step += sim_params::time_step_size;

    }
    sim.PrintVelocities(outputFileName + "22jul_velocities.csv");
    sim.PrintPressure(outputFileName + "22jul_pressure.csv");


    return 0;
}