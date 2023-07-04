
// Program: Simulating a 1d Linear Convection equation
// Date: 02-Jul-2023


#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "supplemental_files/supplemental_files.hpp"


int main()
{
    // File where results are going to be output to
    std::string outputFilename = "../src/output_data/04Jul2023_data.csv";

    const double nx = 40;
    const double dx = 2 / (nx-1);
    double counter = 0;

    const int nt = 25;
    const int dt = 1;

    const int c = 1;

    Eigen::VectorXf u( (int(nx)) );
    Eigen::VectorXf un( (int(nx)) );

    // Initializing the Initial conditions
    while(counter<nx)
    {
        if(counter>=0.5 && counter <=1)
        {
            u.coeffRef(counter) = 2.0;

        }
        else
        {
            u.coeffRef(counter) = 1.0;

        }

        counter += 1;
    }


    // Main loop 
    for(int i = 0; i < nt; i += dt)
    {
        un = u;

        for(int j  = 1; j < nx; j++)
        {
            u.coeffRef(j) = un.coeffRef(j) - c * ( dt / dx ) * ( un.coeffRef(j) - un.coeffRef(j-1));

        }


    }

    // Handling the output of the calculated values in order to be plotted with 
    // python in a separate program
    SupplementalFiles sf;

    sf.printFilesToCSV(outputFilename, u);

    

    return 0;
}