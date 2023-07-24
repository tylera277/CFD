

#include "simulation.hpp"

using namespace navier_stokes;


// Basic constructor
// (which looks attrocious at the moment, I know...)
simulation::Simulation::Simulation(int imin, int imax, int jmin, int jmax, 
double dxi, double dyi, double viscosity, double density,
double time_step_size):
imin(imin), imax(imax), jmin(jmin), jmax(jmax),
dxi(dxi), dyi(dyi),
viscosity(viscosity), density(density), time_step_size(time_step_size),

us(imax+2, jmax+2), vs(imax+2, jmax+2),
u(imax+2, jmax+2), v(imax+2, jmax+2),
p(imax+2, jmax+2),
//laplacian( (imax)*(jmax), (imax)*(jmax)),
laplacian(21*21, 21*21),
rhs( (imax)*(jmax),1 )
{ }

/////////////////////////////////////////
//////////  BOUNDARY CONDITIONS  /////////
/////////////////////////////////////////

void simulation::Simulation::ApplyBoundaryConditions(std::string bc_name)
{

    if(bc_name == "generic")
    {
        this->ApplyGeneric_BC();
    }
    else if(bc_name == "lid_driven_cavity")
    {
        this->ApplyLidDrivenCavity_BC();
    }       
    else
    {
        throw std::invalid_argument("You did not enter one of the boundary condition options!");
    }
            
}

void simulation::Simulation::ApplyGeneric_BC()
{   

    // HACK: I will enter these specific boundary conditions in a better way
    // later on
    double u_bot = 0;
    double u_top = 1;
    double v_left = 0;
    double v_right = 0;


    for(int it = 1; it < imax; it++)
    {
        u(imin-1, it) = (2 * u_top) - u(imin, it);
        u(imax+1, it) = (2 * u_bot) - u(imax, it);

        v(it, jmin-1) = (2 * v_left) - v(it, jmin);
        v(it, jmax+1) = (2 * v_right) - v(it, jmax); 
    }
}

void simulation::Simulation::ApplyLidDrivenCavity_BC()
{


}

/////////////////////////////////////////
//////////  MAIN CALCULATIONS  //////////
/////////////////////////////////////////

void simulation::Simulation::predictor_step()
{
    //std::cout << "BEEPERS \n ";
    for(int i = imin; i <= imax; i++)
    {   
        for(int j = jmin; j <= jmax; j++)
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


    double d2u_dx2 = (u(i-1, j) - 2 * u(i, j) + u(i+1, j)) * pow(dxi, 2);
    double d2u_dy2 = (u(i, j-1) - 2 * u(i, j) + u(i, j+1)) * pow(dyi, 2);

    double u_du_dx = u(i,j) * ( (u(i+1, j) - u(i-1, j)) * (2 * dxi));
    double v_du_dy = (1.0/4.0) * (v(i-1, j) + v(i, j) + v(i-1, j+1) + v(i, j+1)) * (u(i, j+1) - u(i,j-1))*(2*dyi);

    double value =  u(i,j) + time_step_size * ((viscosity)*(d2u_dx2 + d2u_dy2) - (u_du_dx + v_du_dy));

    return value;
}

double simulation::Simulation::v_star(int xPos, int yPos)
{
    int i = xPos;
    int j = yPos;

    double d2v_dx2 = (v(i-1,j) - 2 * v(i,j) + v(i+1,j)) * pow(dxi, 2);
    double d2v_dy2 = (v(i,j-1) - 2*v(i,j) + v(i,j+1)) * pow(dyi, 2);

    double udv_dx = (1.0/4.0) * (u(i, j-1) + u(i, j) + u(i+1, j-1) + u(i+1, j)) * (v(i+1,j)-v(i-1,j))*(2*dxi);
    double vdv_dy = v(i,j) * (v(i,j+1) - v(i,j-1)) * (2 * dyi);

    double value = v(i,j) + time_step_size * (viscosity*(d2v_dx2 + d2v_dy2) - (udv_dx + vdv_dy));

    return value;
}



void simulation::Simulation::corrector_step()
{
    CalculateRHS();
    
    PrintPressureVector(rhs, "../outputData/22jul_rhs.csv" );
    // PrintLaplacian("../outputData/22jul_laplacianV2.csv");
    std::cout << "ROWS, COLS: " << rhs.rows() << "," << rhs.cols() << std::endl;
    std::cout << "ROWS, COLS: " << laplacian.rows() << "," << laplacian.cols() << std::endl;

    std::cout << "DET:" << laplacian.determinant() << std::endl;


    Eigen::PartialPivLU<Eigen::MatrixXd> lu(laplacian);
    Eigen::MatrixXd pv = lu.solve(rhs);
    

    std::cout << "ROWS, COLS: " << pv.rows() << "," << pv.cols() << std::endl;
    std::cout << "----------------" << std::endl;

    PrintPressureVector(pv, "../outputData/22jul_pv.csv");
    //this->ConvertPressureVectorIntoMatrix(pv);


    //calculateUpdatedVelocities();



}

void simulation::Simulation::CreateLaplacian()
{
    int nx = imax - 1;
    int ny = jmax - 1;

    

    for(int j = 1; j < ny; j++)
    {
        for(int i = 1; i < nx; i++)
        {
            laplacian(i+(j-1)*nx, i+(j-1)*nx) = 2 * pow(dxi, 2) + 2 * pow(dyi,2);

            for(int ii=(i-1); ii<(i+2); ii++)
            {
                if(ii>0 && ii<=nx)
                {
                    laplacian(i+(j-1)*nx, ii+(j-1)*nx) = -pow(dxi,2);
                }
                else
                {
                    laplacian(i+(j-1)*nx, i+(j-1)*nx) = laplacian(i+(j-1)*nx, i+(j-1)*nx)-pow(dxi,2);
                }
            }
            for(int jj=(j-1); jj<(j+2); jj++)
            {
                if(jj>0 && jj<=ny)
                {
                    laplacian(i+(j-1)*nx, i+(jj-1)*nx) = -pow(dyi,2);
                }
                else
                {
                    laplacian(i+(j-1)*nx, i+(j-1)*nx) = laplacian(i+(j-1)*nx, i+(j-1)*nx) - pow(dyi,2);
                }

            }


        }
    }
    laplacian(1, 1) = 1;

}

void simulation::Simulation::CreateLaplacianV2()
{
    
    int nx = imax - 1;
    int ny = jmax - 1;

    double dx = 1 / dxi;
    double dy = 1 / dyi;

    double invDeltaX2 = 1.0 / (dx * dx);
    double invDeltaY2 = 1.0 / (dy * dy);

    for(int i = 0; i < nx; ++i)
    {
        for(int j = 0; j < ny; ++j)
        {
            int currentPoint = i * ny + j;

            laplacian(currentPoint, currentPoint) = -2.0 * (invDeltaX2 + invDeltaY2);
             if(i > 0)
            {
                laplacian(currentPoint, currentPoint - ny) = invDeltaX2;
            }
            if (i < nx - 1)
            {
                laplacian(currentPoint, currentPoint + ny) = invDeltaX2;
            }
            if (j > 0)
            {
                laplacian(currentPoint, currentPoint - 1) = invDeltaY2;
            }
            if(j < ny - 1)
            {
                laplacian(currentPoint, currentPoint + 1) = invDeltaY2;
            } 
        }
    }
}

void simulation::Simulation::CalculateRHS()
{
    int n = 0;

    for(int j = jmin; j <= jmax; j++)
    {
        for(int i = imin; i <= imax; i++)
        {
            rhs(n,0) = density/time_step_size * ((us(i+1, j) - us(i, j)) * dxi +\
            (vs(i, j+1) - vs(i,j))*dyi);
            
            std::cout << rhs(n,0) << std::endl;
            n++;
        }

    }

}

void simulation::Simulation::ConvertPressureVectorIntoMatrix(Eigen::MatrixXd pv)
{
    int n = 0;
    
    for(int i = imin; i <= imax; i++)
    {
        for(int j = jmin; j <= jmax; j++)
        {
            p(i,j) = pv(n);

            n++;
        }
    }
}

void simulation::Simulation::calculateUpdatedVelocities()
{
    for(int i = imin; i < imax; i++)
    {
        for(int j = jmin; j < jmax; j++)
        {
            u(i,j) = us(i,j) - (time_step_size / density) * (p(i,j) - p(i-1,j)) * dxi;
            v(i,j) = vs(i,j) - (time_step_size / density) * (p(i,j) - p(i,j-1)) * dyi;
        }
    }


}

/////////////////////////////////////////
//////////  UTILITY FUNCTIONS  //////////
/////////////////////////////////////////

void simulation::Simulation::PrintPressure(std::string outputFile)
{
    this->PrintMatrixToCSVFile(p, outputFile);
}

void simulation::Simulation::PrintVelocities(std::string outputFile)
{
    this->PrintMatrixToCSVFile(us, outputFile);
}

void simulation::Simulation::PrintLaplacian(std::string outputFile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputFile);

    for(int i = imin; i < ((imax)*(jmax)); i++)
    {
        for(int j = jmin; j < ((imax)*(jmax)); j++)
        {
            file << laplacian(i, j) << ",";
        }
        file << "\n";
    }

    file.close();
}

void simulation::Simulation::PrintPressureVector(Eigen::MatrixXd pv, std::string outputFile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputFile);

    for(int i = imin; i < ((imax-1)*(jmax-1)); i++)
    {
        file << pv(i, 0) << "\n";
    }

    file.close();
}

void simulation::Simulation::PrintMatrixToCSVFile(Eigen::MatrixXd matrix, std::string outputfile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputfile);

    // Used to make it easier to use pandas in python for plotting the results
    for(int filler = imin; filler < imax; filler++)
    {
        file << 0 << ",";
    }
    file << "\n";


    for(int i = imin; i < imax; i++)
    {
        for(int j = jmin; j < jmax; j++)
        {
            file << matrix(i, j) << ",";
        }
        file << "\n";
    }

    file.close();


}