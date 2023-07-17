// 16Jul2023
// Houses the main simulation body class, which will
// be interfaced via the main.cpp file

#include <Eigen/Dense>
#include <string>

namespace navier_stokes{

    namespace simulation{


        class Simulation
        {
        
            // X and Y velocity components. 
            //(Values which are based at the faces of each cell.
            // s is for predictor step values, no s is corrector step solution)
            Eigen::MatrixXd us;
            Eigen::MatrixXd vs;

            Eigen::MatrixXd u;
            Eigen::MatrixXd v;
            // Pressure. 
            // (A value which is stored at the center of each cell)
            Eigen::MatrixXd p;


            public:
                // Basic constructor which initializes the matrices that are going to be used
                Simulation();


                /**
                 * @brief Sets the boundary conditions of the problem that is trying to be solved
                 * 
                 * @param boundary_conditions_name  The name of the particular B.C.'s one wants to utilize
                 *
                 */
                void applyBoundaryConditions(std::string boundary_conditions_name);

                /**
                 * @brief Calculates the velocity components for the initial predictor step
                 * (u and v star).
                 * REMINDER: the pressure component is ignored for this step
                 *
                 */
                void predictor_step();

                double u_star(int xPos, int yPos);

                double v_star(int xPos, int yPos);



                void corrector_step();



        };



    } // namespace: simulation
} //namespace: navier_stokes