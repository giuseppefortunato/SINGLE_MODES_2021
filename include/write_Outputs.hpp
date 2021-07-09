#ifndef WRITE_OUTPUTS_HPP
#define WRITE_OUTPUTS_HPP

#include "Reconstruction.hpp"
#include "read_Inputs.hpp"

const int CGNS_STRING_SIZE = 33;

void Config_stream ( prob_settings settings );



string write_reconstruction_file ( const Eigen::MatrixXd &Rec,
                                   const Eigen::MatrixXd &Coords,
                                   int nt,
                                   int nC,
                                   const int Nm,
                                   prob_settings settings,
                                   Eigen::VectorXi &ID,
                                   bool flag);




#endif // WRITE_OUTPUTS_HPP