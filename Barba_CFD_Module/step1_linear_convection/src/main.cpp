
// Program: Simulating a 1d Linear Convection equation
// Date: 02-Jul-2023


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

#include "supplemental_files/supplemental_files.hpp"


int main()
{
    // File where results are going to be output to
    std::string outputFilename = "../src/output_data/04Jul2023_data.csv";

    const double nx = 40;
    const double dx = 2 / (nx-1);
    double counter = 0;
    

    const int nt = 50;
    const double dt = 0.025;

    const int c = 1;

    // Vectors for storing the calculated and reference values
    Eigen::VectorXf u( (int(nx)) );
    Eigen::VectorXf un( (int(nx)) );

    // Initializing the Initial conditions
    while(counter<nx)
    {
        double positionOfPoint = counter * dx;
        if(positionOfPoint>=0.5 && positionOfPoint <=1)
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
    for(int i = 0; i < nt; i += 1)
    {
        un = u;

        for(int j  = 1; j < nx; j++)
        {
            u.coeffRef(j) = un.coeffRef(j) - (c * ( dt / dx ) ) * ( un.coeffRef(j) - un.coeffRef(j-1));

            std::cout << "New Value: " << u.coeffRef(j) << "\n";

        }
        std::cout << "------------\n"; 


    }

    /*
    ** Struggling on getting this to compile for some reason

    // Handling the output of the calculated values in order to be plotted with 
    // python in a separate program
    SupplementalFiles sf;

    sf.printFilesToCSV(outputFilename, u);
    */

    std::ofstream output;
    output.open(outputFilename);

    for(int iterator = 0; iterator < u.size(); iterator++)
    {
        if(iterator==(u.size()-1))
        {
            output<<0;
        }
        else
        {
            output << 0 << ",";

        }


    }
    output<<"\n";
    for(int iterator = 0; iterator < u.size(); iterator++)
    {
        if(iterator==(u.size()-1))
        {
             output << u.coeffRef(iterator);
        }
        else
        {
            output << u.coeffRef(iterator) << ",";

        }

    }
    // Closing file after having finished printed out to it
    output.close();
    

    return 0;
}