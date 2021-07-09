#include "Reconstruction.hpp"



#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ) {

    if (key_string == "LINEAR")
        return smartuq::surrogate::LINEAR;
    else if (key_string == "CUBIC")
        return smartuq::surrogate::CUBIC;
    else if (key_string == "GAUSSIAN")
        return smartuq::surrogate::GAUSSIAN;
    else if (key_string == "THIN_PLATE")
        return smartuq::surrogate::THIN_PLATE;
    else if (key_string == "MULTIQUADRATICS")
        return smartuq::surrogate::MULTIQUADRATICS;

    return smartuq::surrogate::LINEAR;

}


std::vector<double> ISOMAP_surr_coeff( const Eigen::MatrixXd &param,
                                       const Eigen::RowVectorXd &norme,
                                       const Eigen::MatrixXd &Coeffs,
                                       const double param_rec_1,
                                       const double param_rec_2,
                                       const int Np,
                                       std::string flag_interp,
                                       std::string analysis,
                                       const Eigen::RowVectorXd supercritical) {

    std::vector <std::vector<double>> T(param.rows(), std::vector<double>(param.cols()));
    //for supercritical
    //std::vector<double> all_rec(Np);
    //for (int i=0; i<Np; i++)
    //    all_rec[i]=supercritical(i);
    std::vector<double> all_rec(Np);
    if(analysis== "UNSTEADY"){
        all_rec[0]=param_rec_1;
    }else{
        all_rec[0] = param_rec_1;
        all_rec[1] = param_rec_2;}


    for (int i = 0; i < param.rows(); i++) {
        for(int j=0; j< param.cols(); j++)
            T[i][j] = param(i, j)/ norme(j);

    }
    //std::cout<<"inside RBF function, T="<<std::endl;
    //std::cout<<T<<std::endl;

    std::vector<double> t(Np);
        for (int i=0; i<Np; i++)
        t[i] = all_rec[i]/ norme(i);


    //Vector of surrogate coefficients
    std::vector<double> coefs_intrp(Coeffs.rows());

    double avgDt=0.0;

    // Create surrogates for coefficients
    std::vector <rbf> surr_coefs = {};
    RBF_CONSTANTS rbf_const{avgDt, 0.0};

    for (int i = 0; i < Coeffs.rows(); i++) {

        std::vector<double> coefs;
        for (int j = 0; j < param.rows(); j++)
            coefs.push_back(Coeffs(i, j));

        surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
        surr_coefs[i].build();
        surr_coefs[i].evaluate(t, coefs_intrp[i]);
    }

    return coefs_intrp;

}




Eigen::MatrixXd Reconstruction_POD ( const Eigen::MatrixXd &param,
                                     const Eigen::RowVectorXd &norme,
                                     const Eigen::VectorXd &K_pc,
                                     const Eigen::VectorXd &lam,
                                     const Eigen::MatrixXd &Coeffs,
                                     const Eigen::MatrixXd &phi,
                                     const double param_rec_1,
                                     const double param_rec_2,
                                     int Nrec,
                                     int Np,
                                     std::string flag_interp,
                                     std::string analysis,
                                     const Eigen::RowVectorXd &supercritical) {

    std::vector<std::vector<double> > T(param.rows(), std::vector<double>(param.cols()));
    //for supercritical
    //std::vector<double> all_rec(Np);
    //for (int i=0; i<Np; i++)
    //    all_rec[i]=supercritical(i);
    std::vector<double> all_rec(Np);
    if(analysis== "UNSTEADY"){
        all_rec[0]=param_rec_1;
    }else{
        all_rec[0] = param_rec_1;
        all_rec[1] = param_rec_2;}


    double avgDt = 0.0;
    for (int i = 0; i < param.rows(); i++) {
        for (int j = 0; j < param.cols(); j++) {
            T[i][j] = param(i, j) / norme(j);
        }
    }


    std::vector<double> t(Np);
        for (int i = 0; i < Np; i++)
            t[i] = all_rec[i] / norme(i);

        //for ( int i = 1; i < t_vec.size(); i++ ) {
        //    avgDt += t_vec[i] - t_vec[i-1];
        //}

        //avgDt = avgDt/(double)(t_vec.size()-1);

        //Vector of surrogate coefficients
        std::vector<double> coefs_intrp(Nrec);

        // Create surrogates for coefficients
        std::vector<rbf> surr_coefs = {};
        RBF_CONSTANTS rbf_const{avgDt, 0.0};


        for (int i = 0; i < Nrec; i++) {

            std::vector<double> coefs;
            for (int j = 0; j < param.rows(); j++)
                coefs.push_back(Coeffs(i, j));

            surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
            surr_coefs[i].build();
            surr_coefs[i].evaluate(t, coefs_intrp[i]);
        }

        Eigen::VectorXd coefs_t(Nrec);

        for (int i = 0; i < Nrec; i++)
            coefs_t(i) = coefs_intrp[i];

        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);

        for (int i = 0; i < Nrec; i++)
            Sig(i, i) = std::sqrt(lam(i));

        Eigen::MatrixXd Rec_field = phi.leftCols(Nrec) * Sig * coefs_t.head(Nrec);

        return Rec_field;

    }



    std::vector<rbf> getSurrCoefs(const std::vector<double> &t_vec,
                                  const Eigen::MatrixXd &Coeffs,
                                  std::string flag_interp) {
        int Ns = Coeffs.rows();
        int Nrec = Coeffs.cols();

        std::vector<std::vector<double> > T(t_vec.size(), std::vector<double>(1));
        double avgDt = 0.0;

        for (int i = 0; i < t_vec.size(); i++) {
            T[i][0] = t_vec[i];
        }

        for (int i = 1; i < t_vec.size(); i++) {
            avgDt += t_vec[i] - t_vec[i - 1];
        }

        avgDt = avgDt / (double) (t_vec.size() - 1);

        // Create surrogates for coefficients
        std::vector<rbf> surr_coefs{};
        RBF_CONSTANTS rbf_const{avgDt, 0.0};

        for (int i = 0; i < Nrec; i++) {

            std::vector<double> coefs{};
            for (int j = 0; j < t_vec.size(); j++)
                coefs.push_back(Coeffs(j, i));

            surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
            surr_coefs[i].build();

        }

        return surr_coefs;

    }










#pragma clang diagnostic pop
