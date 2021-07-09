#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "Extract_Basis.hpp"
#include "manifold_learning.hpp"
#include "Extract_Basis.hpp"
#include <cmath>
#include "Surrogates/rbf.h"

using namespace smartuq::surrogate;


smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ); 


Eigen::MatrixXd Reconstruction_POD ( const Eigen::MatrixXd &param,
                                    const Eigen::RowVectorXd &norme,
                                    const Eigen::VectorXd &K_pc,
                                    const Eigen::VectorXd &lam,
                                    const Eigen::MatrixXd &Coeffs,
                                    const Eigen::MatrixXd &phi,
                                    double param_rec_1,
                                    double param_rec_2,
                                    int Nrec,
                                    int Np,
                                    std::string flag_interp,
                                    std::string analysis,
                                    const Eigen::RowVectorXd &supercritical);



std::vector<double> ISOMAP_surr_coeff( const Eigen::MatrixXd &param,
                                       const Eigen::RowVectorXd &norme,
                                       const Eigen::MatrixXd &Coeffs,
                                       double param_rec_1,
                                       double param_rec_2,
                                       const int Np,
                                       std::string flag_interp,
                                       std::string analysis,
                                       const Eigen::RowVectorXd supercritical);






std::vector<rbf> getSurrCoefs ( const std::vector<double> &t_vec,
                            const Eigen::MatrixXd &Coeffs,
                            std::string flag_interp = "LINEAR" );






#endif // RECONSTRUCTION_HPP

