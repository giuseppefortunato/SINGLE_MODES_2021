//
// Created by giuseppe on 08/06/21.
//


#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"
#include "manifold_learning.hpp"

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------Single-MODES start-------------" << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];



    //Reading configuration file
    Read_cfg( filecfg, settings );
    Config_stream ( settings );

    // Calculate number of grid points


    std::string file_1 = settings.in_file[0];


    //DEFINE PARAMETERS
    Eigen::MatrixXd param_t(settings.Ns, settings.Np);
    std::string file_parameters = argv[2];
    define_parameters(settings, param_t, file_parameters);


    //DEFINE ALL POINTS DESIGN OF DESIGN SPACE WHERE TO PERFORM
    bool design_space_var = false;
    Eigen::MatrixXd axis = Eigen::MatrixXd::Zero(15, settings.Np);
    Eigen::MatrixXd design_space= Eigen::MatrixXd::Zero(axis.rows()*axis.rows(), settings.Np);
    if (design_space_var) {

        axis.col(0)= Eigen::VectorXd::LinSpaced(15, 50, 100); //Pe
        axis.col(1)= Eigen::VectorXd::LinSpaced(15, 0.5, 3.5); //Ta

        int interval = 0;
        for (int i=0; i < axis.rows(); i++) {
            for (int j=0; j < axis.rows(); j++) {
                design_space(interval+j, 0) = axis(i,0);
                design_space(interval+j, 1) = axis(j,1);
            }
            interval = interval + axis.rows();
        }

        for (int i=0;i<design_space.rows(); i++){
            settings.param_rec_1.push_back(design_space(i,0));
            settings.param_rec_2.push_back(design_space(i,1));}

    }


    Eigen::RowVectorXd supercritical(settings.Np);
    for (int i=0; i<settings.Np; i++)
        supercritical(i)= settings.supercritical_rec[i];


    Eigen::RowVectorXd normes_parameters_t(settings.Np);
    //COMPUTING PARAMETERS NORM
    for (int i=0; i<settings.Np; i++)
        normes_parameters_t(i)= param_t.col(i).norm();

    int Nr;
    Eigen::VectorXi ID = N_gridpoints ( file_1, Nr );
    std::cout<<"ID point"<<std::endl;
    std::cout<<ID<<std::endl;
    std::cout << "Number of grid points : " << Nr << std::endl;
    int nC = settings.Cols.size();
    int s_Nf = 1;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    int nfj = 0;

    // Reading coordinates
    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );          //file_1 in origine, file_cord se supercritical
    std::cout << "Done " << std::endl;

    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set_t = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                                     settings.Cols,
                                                     settings.in_file,
                                                     settings.flag_prob);

    Eigen::VectorXd mean = sn_set_t.rowwise().mean();
    std::cout<<"done"<<std::endl;

    std::cout << "Storing test Matrix ... \n ";
    Eigen::MatrixXd sn_set_test = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                                     settings.Cols,
                                                     settings.test_file,
                                                     settings.flag_prob);

    std::cout<<"done"<<std::endl;





    Eigen::VectorXd K_pc(settings.Ns);


    if (settings.flag_method[0] == "ISOMAP") {
        std::cout<<"ISOMAP OPTIMAL ANALYSIS"<<std::endl;



        for (int number_r = 0; number_r < settings.r_isomap.size(); number_r++) {
            std::cout << " ISOMAP for \t" << settings.r_isomap[number_r] << " \t dimension of the manifold"
                      << std::endl;
            //Eigen::MatrixXd snapshot_star_total = Eigen::MatrixXd::Zero(Nr, nC);
            Eigen::MatrixXd sn_set_partial(Nr, sn_set_t.cols());
            Eigen::VectorXd norma(sn_set_t.cols());
            Eigen::MatrixXd y(settings.r_isomap[number_r] * nC, sn_set_t.cols());
            Eigen::VectorXi k(nC);
            Eigen::VectorXi Rank = Eigen::VectorXi::Ones(nC);
            Eigen::VectorXd kruskal_min = Eigen::VectorXd::Ones(nC);




            //Loop for each conservative variable---> ISOMAP
            for (int q = 0; q < nC; q++) {
                std::cout << "ISOMAP FOR CONSERVATIVE VARIABLE NUMBER " << "\t" << q << std::endl;
                ISOMAP(sn_set_t, settings, Nr, q, k, kruskal_min, y, number_r, Rank);
            }

            std::cout << "DONE " << std::endl;
            std::cout << "Kruskal stress minimo \t " << kruskal_min << "\t per k=" << k << std::endl;
            std::cout << y << std::endl << std::endl << std::endl;
            std::cout << "END OF ISOMAP..." << std::endl << std::endl;


            std::cout << "------------START BACK-MAPPING SECTION--------------" << std::endl << std::endl;

            Eigen::MatrixXd final_decision= Eigen::MatrixXd::Zero(nC,settings.in_file.size()-2);
            for (int rec = 2; rec < settings.in_file.size(); rec++) {
                Eigen::MatrixXd rec_test= Eigen::MatrixXd::Zero(nC*Nr, settings.test_file.size());
                Eigen::MatrixXd diff_test= Eigen::MatrixXd::Zero(nC*Nr, settings.test_file.size());
                settings.k_rec[0] = rec;
                if(settings.test_file.size() != settings.param_rec_1.size()){
                    std::cout<<"your test files doesn't match the number of parameters"<<std::endl;
                    exit (EXIT_FAILURE);}
                for (int nt = 0; nt < settings.param_rec_1.size(); nt++) {
                    std::cout << "RECONSTRUCTING FIELD FOR:\t" << settings.param_rec_1[nt]
                              << "\t" << settings.param_rec_2[nt] << std::endl;

                    for (int q = 0; q < nC; q++) {

                        sn_set_partial = sn_set_t.block(Nr * q, 0, Nr, sn_set_t.cols());
                        Eigen::MatrixXd y_part(settings.r_isomap[number_r], sn_set_t.cols());
                        for (int j = 0; j < settings.r_isomap[number_r]; j++) {
                            y_part.row(j) = y.row((q * settings.r_isomap[number_r]) + j);
                        }


                        std::vector<double> y_star = ISOMAP_surr_coeff(param_t, normes_parameters_t, y_part,
                                                                       settings.param_rec_1[nt],
                                                                       settings.param_rec_2[nt],
                                                                       settings.Np, settings.flag_interp,
                                                                       settings.analysis,
                                                                       supercritical);

                        std::cout << "values RBF" << std::endl;
                        for (int i = 0; i < y_star.size(); i++)
                            std::cout << y_star[i] << std::endl;

                        int f=0;
                        Eigen::MatrixXd snapshot_star = create_rec_fields(settings, y_star, y_part, sn_set_partial, number_r, f);

                        std::cout << "size  snapstar [" << snapshot_star.rows() << ", " << snapshot_star.cols() << "]" << std::endl;
                        rec_test.block( q*Nr,nt,Nr,1) = snapshot_star;

                    }

                }

                //qua dovrei aver finito di riempire la matrice delle ricostruzioni per un valore unico di punti vicini
                //std::cout << "size rec test [" << rec_test.rows() << ", " << rec_test.cols() << "]" << std::endl;
                //std::cout << "size sn_set_test [" << sn_set_test.rows() << ", " << sn_set_test.cols() << "]" << std::endl;
                //std::cout << "size  diff_test [" << diff_test.rows() << ", " << diff_test.cols() << "]" << std::endl;

//                for (int i=0; i<diff_test.rows(); i++){
//                    for(int j=0; j<diff_test.cols(); j++){
//                        diff_test(i,j)= (sn_set_test(i,j) - rec_test(i,j))/sn_set_test(i,j);
//                    }
//                }
                diff_test= (sn_set_test - rec_test);


                // trovare un modo per calcolare l'errore
                Eigen::MatrixXd decision(nC,diff_test.cols());

                for (int q=0; q<nC; q++){
                    Eigen::MatrixXd tmp = diff_test.block(q*Nr,0,Nr,diff_test.cols());
                    //std::cout<<tmp<<std::endl;
                    for(int j=0; j<diff_test.cols(); j++){
                        double sum=0;
                        for (int i=0; i<Nr; i++) {
                            sum = sum + (tmp(i,j)*tmp(i,j));
                        }
                        //std::cout<<sum<<std::endl;
//                        decision(q,j)= sqrt((sum));
                        decision(q,j)=sqrt((sum))/sn_set_test.block(q*Nr,j,Nr,1).norm();

                        //std::cout<<"decision"<<std::endl;
                        //std::cout<<decision(q,j)<<std::endl;
                    }
                        std::cout<<"error matrix for k_rec= "<< rec<<std::endl;
                        std::cout<<decision<<std::endl;
                        final_decision(q,rec-2)= decision.row(q).mean();

                }

                //salvare in modo furbo il numero di nearest neighbors da utilizzare a seconda della variabile conservata
            }

            std::cout<<"FINAL DECISION"<<std::endl;
            for (int i=0; i<nC; i++){
                int j;
                std::cout<<" for variable number "<< i <<" min value= "<< final_decision.row(i).minCoeff(&j)<<std::endl;
                std::cout<<" for number of neighbors equal to: "<< j+2 <<std::endl;
            }
            std::cout<<final_decision<<std::endl;

        }

    }

    std::cout << std::endl;
    std::cout<<"------------------------------------------" <<std::endl;
    std::cout<<"------------Optimal ISO end--------------" <<std::endl;
    std::cout<<"------------------------------------------" <<std::endl;

    return 0;


}
