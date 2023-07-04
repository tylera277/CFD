
#include <fstream>

#include "supplemental_files.hpp"



void SupplementalFiles::printFilesToCSV(std::string outputFileName, Eigen::VectorXf vector)
{

    std::ofstream output;
    output.open(outputFileName);

    for(int iterator = 0; iterator < vector.size(); iterator++)
    {
        output << vector.coeffRef(iterator) << ", ";


    }
    output<<"\n";

    // Closing file after having finished printed out to it
    output.close();

}