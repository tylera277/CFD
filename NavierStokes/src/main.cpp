
#include <iostream>
#include <Eigen/Dense>

#include "simulation/simulation.hpp"


using namespace navier_stokes;

int main()
{



    //Eigen::Matrix<double, ny, nx> x;
    //Eigen::VectorXd x = Eigen::VectorXd::LinSpaced ( 0, simulation_parameters::length_x, simulation_parameters::x_grid_points + 1);
    //Eigen::VectorXd y = Eigen::VectorXd::LinSpaced ( 0, simulation_parameters::length_y, simulation_parameters::y_grid_points + 1);

    /* Eigen::VectorXd xm = 0.5 * (x( Eigen::seq(imin, imax) ) + x( Eigen::seq(imin+1, imax+1) ));
    Eigen::VectorXd ym = 0.5 * (y( Eigen::seq(jmin, jmax) ) + y( Eigen::seq(jmin+1, jmax+1) )); */
    double current_time_step = 0;
    double total_time = 100;
    double time_step = 1;

    simulation::Simulation sim;

    //sim.applyBoundaryConditions("PUT SOMETHING HERE FOR THE FUTURE");

    while(current_time_step < total_time)
    {
        sim.predictor_step();

        current_time_step += time_step;

    }



   
    





    return 0;
}