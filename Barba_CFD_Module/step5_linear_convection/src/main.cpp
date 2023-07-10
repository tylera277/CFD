
// Program: Simulating a 1d NON-Linear Convection equation
// Date: 08-Jul-2023


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>


int main()
{
    // File where results are going to be output to
    std::string outputFilename = "../src/output_data/09Jul2023_data.csv";

    const int nx = 40;
    const int ny = 40;

    const double dx = 2 / ((double)nx - 1);
    const double dy = 2 / ((double)ny - 1);
    double counter_x = 0;
    double counter_y = 0;
    

    const int nt = 25;
    const double dt = 0.025;

    // User defined parameter.
    const double c = 1;



    // Vectors for storing the calculated and reference values
    //Eigen::Vector2d u( (int(nx)), (int(ny)) );
    //Eigen::Vector2d un( (int(nx)), (int(ny)) );

    Eigen::Matrix<double, ny, nx> u_2d;
    Eigen::Matrix<double, ny, nx> un_2d;

    
    // Initializing the Initial conditions
    while(counter_x < nx)
    {
        double positionOfPoint_x = (counter_x * dx);
        counter_y = 0;

        while(counter_y < ny)
        {
             double positionOfPoint_y = (counter_y * dy);
            
            if(( positionOfPoint_x>=0.5 && positionOfPoint_x <=1) && \
                (positionOfPoint_y>=0.5 && positionOfPoint_y <=1))
            {
                u_2d((int)counter_y, (int)counter_x) = 2.0;

            }
            else
            {
                u_2d((int)counter_y, (int)counter_x) = 1.0;

            } 
            counter_y += 1;
        }

        counter_x += 1;
    }


    // File handling stuff
    std::ofstream output;
    output.open(outputFilename);

    // +This loops sole purpose is to make it quicker for me 
    // to plot the results in python with also a  little less work.
    // +Outputs a row of zeros as a buffer, so python pandas library takes these
    // as the column names instead of my actual data.
    
    for(int iterator = 0; iterator < nx; iterator++)
    {
        if(iterator==(nx-1))
        {
            output<<0;
        }
        else
        {
            output << 0 << ",";
        }
    }
    output<<"\n";


    //Main loop 
    for(int currentTime = 0; currentTime < nt; currentTime += 1)
    {
        un_2d = u_2d;

        for(int y_position  = 1; y_position < ny; y_position++)
        {
            for(int x_position = 1; x_position < nx; x_position++)
            {
                double first_part = un_2d(x_position, y_position);
                double second_part = -c * (dt/dx) * (un_2d(y_position, x_position) - un_2d(y_position-1, x_position));
                double third_part = -c * (dt/dy) * (un_2d(y_position, x_position) - un_2d(y_position, x_position-1));
                double result = first_part + second_part + third_part;

                u_2d(y_position, x_position) = result;
            }
        }

        // Print out resulting calculated vector into a designated csv file
        if(currentTime == 0){
            for(int y_position_printing  = 0; y_position_printing < ny; y_position_printing++)
            {
                for(int x_position_printing = 0; x_position_printing < nx; x_position_printing++)
                {

                    if( x_position_printing == (nx-1))
                    {
                        output << u_2d(y_position_printing, x_position_printing);
                    }
                    else
                    {
                        output << u_2d(y_position_printing, x_position_printing) << ",";
                    }
                }
                output << "\n";
            }
        }
    }

    

    
    
    
    // Closing file after having finished printing out to it
    output.close();
    

    return 0;
}