
#include <iostream>
#include <cmath>

#include "simulation.hpp"
#include "../simulation_parameters.hpp"

using namespace navier_stokes;


// Basic constructor
simulation::Simulation::Simulation():
us(sim_params::imax+1, sim_params::jmax+1), vs(sim_params::imax+1, sim_params::jmax+1),
u(sim_params::imax+1, sim_params::jmax+1), v(sim_params::imax+1, sim_params::jmax+1),
p(sim_params::imax+1, sim_params::jmax+1)
{ }


void simulation::Simulation::predictor_step()
{
    for(int i = sim_params::imin; i < sim_params::imax; i++)
    {   
        for(int j = sim_params::jmin; j < sim_params::jmax; j++)
        {
            us(i, j) = this->u_star(i, j);
            vs(i, j) = this->v_star(i, j);


        }
    }
}

double simulation::Simulation::u_star(int xPos, int yPos)
{
    int i = xPos;
    int j = yPos;


    double d2u_dx2 = (u(i-1, j) - 2 * u(i, j) + u(i+1, j)) * pow(sim_params::dxi, 2);
    double d2u_dy2 = (u(i, j-1) - 2 * u(i, j) + u(i, j+1)) * pow(sim_params::dyi, 2);

    double u_du_dx = u(i,j) * ( (u(i+1, j) - u(i-1, j)) * (2 * sim_params::dxi));
    double v_du_dy = (1.0/4.0) * (v(i-1, j) + v(i, j) + v(i-1, j+1) + v(i, j+1)) * (u(i, j+1) - u(i,j-1))*(2*sim_params::dyi);

    double value =  u(i,j) + sim_params::time_step_size * ((sim_params::viscosity)*(d2u_dx2 + d2u_dy2) - (u_du_dx + v_du_dy));

    return value;
}

double simulation::Simulation::v_star(int xPos, int yPos)
{
    int i = xPos;
    int j = yPos;

    double d2v_dx2 = (v(i-1,j) - 2 * v(i,j) + v(i+1,j)) * pow(sim_params::dxi, 2);
    double d2v_dy2 = (v(i,j-1) - 2*v(i,j) + v(i,j+1)) * pow(sim_params::dyi, 2);

    double udv_dx = (1.0/4.0) * (u(i, j-1) + u(i, j) + u(i+1, j-1) + u(i+1, j)) * (v(i+1,j)-v(i-1,j))*(2*sim_params::dxi);
    double vdv_dy = v(i,j) * (v(i,j+1) - v(i,j-1)) * (2 * sim_params::dyi);

    double value = v(i,j) + sim_params::time_step_size * (sim_params::viscosity*(d2v_dx2 + d2v_dy2) - (udv_dx + vdv_dy));

    return value;
}