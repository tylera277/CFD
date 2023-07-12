
// Program: Simulating a 1d NON-Linear Convection equation
// Date: 06-Jul-2023


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>


int main()
{
    // File where results are going to be output to
    std::string outputFilename = "../src/output_data/06Jul2023_data.csv";

    const double nx = 40;
    const double dx = 2 / (nx-1);
    double counter = 0;
    

    const int nt = 25;
    const double dt = 0.025;

    // User defined parameter.
    // VISCOSITY
    const double nu = 0.3;



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


    // File handling stuff
    std::ofstream output;
    output.open(outputFilename);

    // +This loops sole purpose is to make it quicker for me 
    // to plot the results in python with also a  little less work.
    // +Outputs a row of zeros as a buffer, so python pandas takes these
    // as the column names instead of my actual data.
    for(int iterator = 1; iterator < (u.size()-1); iterator++)
    {
        if(iterator==(u.size()-2))
        {
            output<<0;
        }
        else
        {
            output << 0 << ",";
        }
    }
    output<<"\n";


    // Main loop 
    for(int i = 0; i < nt; i += 1)
    {
        un = u;

        for(int j  = 1; j < (nx-1); j++)
        {
            u.coeffRef(j) = un.coeffRef(j) + nu * (dt/pow(dx,2)) * ( un.coeffRef(j+1) + un.coeffRef(j-1) - 2*un.coeffRef(j));
        }

        // Print out results each certain number of loops
        if(i == 24){
            for(int iterator = 1; iterator < (u.size()-1); iterator++)
            {
                if(iterator==(u.size()-2))
                {
                    output << u.coeffRef(iterator);
                }
                else
                {
                    output << u.coeffRef(iterator) << ",";
                }
        
                }
            }
    //output << "\n";
    }

    

    
    
    
    // Closing file after having finished printing out to it
    output.close();
    

    return 0;
}