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
    std::cout << "-----------Single-MODES-leave-one-out-start-------------" << std::endl << std::endl;

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

    if (settings.in_file.size() != param_t.rows() ){
        std::cout<<"Number of files as input is different from the list in param "<<std::endl;
        exit(EXIT_FAILURE);
    }



    Eigen::RowVectorXd supercritical(settings.Np);
    for (int i=0; i<settings.Np; i++)
        supercritical(i)= settings.supercritical_rec[i];


    bool flag= true;
    int Nr ;
    Eigen::VectorXi ID = N_gridpoints ( file_1, Nr );
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
//  Eigen::VectorXd Ic = IC(settings, ,Nr)
    std::cout<<"done"<<std::endl;

    std::cout << "Storing test Matrix ... \n ";
    Eigen::MatrixXd sn_set_test = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                                        settings.Cols,
                                                        settings.test_file,
                                                        settings.flag_prob);

    std::cout<<"done"<<std::endl;



    Eigen::VectorXd K_pc(settings.Ns-1);
   if (settings.flag_method[0] == "POD"){
        std::cout<<"POD ANALYSIS"<<std::endl;

        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);


        //std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::MatrixXd eig_vec(settings.Ns-1, settings.Ns-1);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now

       Eigen::MatrixXd rec_test= Eigen::MatrixXd::Zero(nC*Nr, settings.in_file.size());
       Eigen::MatrixXd diff_test= Eigen::MatrixXd::Zero(nC*Nr, settings.in_file.size());
     for (int leave=0; leave<settings.Ns; leave++) {
        std::cout << "LEAVE ONE OUT FOR SNAPSHOT " << leave << std::endl;
        //removing the column that correspond at the snapshot
        Eigen::MatrixXd sn_set(sn_set_t.rows(), sn_set_t.cols() - 1);
        std::vector<Eigen::MatrixXd> coefs = {};
        int count = 0;
        for (int j = 0; j < sn_set_t.cols(); j++){
            if (leave != j) {
                sn_set.col(count) = sn_set_t.col(j);
                count++;
            }
        }

        for (int number_r = 0; number_r < settings.r.size(); number_r++) {
            std::cout << "PROCESSING BASIS CONSIDERING RANK = " << settings.r[number_r] << std::endl;
            for (int i = 0; i < nC; i++) {
                Phi[i] = Eigen::MatrixXd::Zero(Nr, settings.Ns-1 );
                lambda[i] = Eigen::VectorXd::Zero(settings.Ns-1);
            }

            std::cout << "Extracting basis ... " << std::endl << std::endl;
            for (int ncons = 0; ncons < nC; ncons++) {
                std::cout << "Processing conservative variable " << ncons <<"...";
                Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                        lambda[ncons], K_pc, eig_vec,
                                        Nf[nfj],
                                        settings.flag_bc,
                                        settings.flag_filter,
                                        settings.sigma);


                N_notZero = Phi[ncons].cols();
                if (settings.r[number_r] == 0) Nm[ncons] = Nmod(settings.En, K_pc);
                else Nm[ncons] = std::min(settings.r[number_r], N_notZero);
                coefs.push_back(eig_vec.transpose());
                std::cout<<"done"<<std::endl;
            }

            Eigen::MatrixXd Field(Nr, nC);

            settings.param_rec_1[0]=param_t(leave,0);
            settings.param_rec_2[0]=param_t(leave,1);

            // devo escludere il parametro dal file per il leave one out
            Eigen::MatrixXd param(param_t.rows()-1, settings.Np);
            Eigen::RowVectorXd normes_parameters(settings.Np);
            count=0;
            for (int w=0; w<param_t.rows(); w++) {
                if (leave != w) {
                    param(count,0)=param_t(w,0);
                    param(count,1)= param_t(w,1);
                    count ++;
                }
            }
            for (int i=0; i<settings.Np; i++)
                normes_parameters(i)= param.col(i).norm();


            std::cout << "Reconstructing fields for:" << settings.param_rec_1[0]<< " " << settings.param_rec_2[0] <<"...";
            for (int ncons = 0; ncons < nC; ncons++) {
                Eigen::MatrixXd Rec = Reconstruction_POD(param, normes_parameters, K_pc, lambda[ncons],
                                                             coefs[ncons],
                                                             Phi[ncons],
                                                             settings.param_rec_1[0],
                                                             settings.param_rec_2[0],
                                                             Nm[ncons],
                                                             settings.Np,
                                                             settings.flag_interp,
                                                             settings.analysis,
                                                             supercritical);

                Field.col(ncons) = Rec;
                rec_test.block( ncons*Nr,leave,Nr,1) = Rec;
            }
            std::cout << "Done" << std::endl;

            std::string filename;
            //filename = write_reconstruction_file(Field, Coords, leave, nC, Nm[0], settings, ID, flag);


            //SU2_DTR(settings, filename, su2_conf, nt, number_r);
            //RMS_residuals(settings, number_r, nt, Nr);
            //Write_History_ResError_global(settings, number_r);
        }
     }
     diff_test= (sn_set_t - rec_test);
     Eigen::MatrixXd decision(nC,diff_test.cols());
     for (int q=0; q<nC; q++){
           Eigen::MatrixXd tmp = diff_test.block(q*Nr,0,Nr,diff_test.cols());
           for(int j=0; j<diff_test.cols(); j++){
               double sum=0;
               for (int i=0; i<Nr; i++) {
                   sum = sum + (tmp(i,j)*tmp(i,j));
               }
               decision(q,j)=sqrt((sum))/sn_set_t.block(q*Nr,j,Nr,1).norm();
           }
     }
     std::cout<<"POD error reconstruction matrix"<<std::endl;
     std::cout<<decision<<std::endl;
     std::string filename;
     filename="ERROR_POD_leave_newset2D.csv";
     std::ofstream  flow_data;
     flow_data.open(filename);
     if (!flow_data.is_open()){
         std::cout<<"I can't open the Error file and write into it"<<std::endl;
         exit (EXIT_FAILURE);}
     flow_data<<"\"alpha\"" << ",";
     flow_data<<"\"Reynolds\"" << ",";
     flow_data<<"\"error_rec_Pressure\"" << ",";
     flow_data<<"\"error_rec_Cp\"" << ",";
     flow_data<<"\"error_rec_skin_friction_coeff_x\"" << ",";
     flow_data<<"\"error_rec_skin_friction_coeff_y\"" ;
     flow_data << std::endl;
     for (int i=0; i<param_t.rows(); i++) {
           for (int j = 0; j < param_t.cols(); j++)
               flow_data << std::setprecision(12) << std::scientific << param_t(i, j) << ", ";
           for( int j=0; j<nC; j++)
               flow_data << std::setprecision(12) << std::scientific << decision(j,i) << ", ";

           flow_data<< std::endl;
     }

     flow_data.close();




    std::cout<<"------------------END OF POD--------------------" <<std::endl;

   }

   if (settings.flag_method[0] == "ISOMAP") {
       std::cout<<"ISOMAP ANALYSIS"<<std::endl;
       Eigen::MatrixXd rec_test= Eigen::MatrixXd::Zero(nC*Nr, settings.in_file.size());
       Eigen::MatrixXd diff_test= Eigen::MatrixXd::Zero(nC*Nr, settings.in_file.size());
       for( int leave =0; leave<settings.Ns; leave++) {
         std::cout << "LEAVE ONE OUT FOR SNAPSHOT " << leave << std::endl;
         //removing the column that correspond at the snapshot
         Eigen::MatrixXd sn_set(sn_set_t.rows(), sn_set_t.cols() - 1);
         int count = 0;
         for (int j = 0; j < sn_set_t.cols(); j++){
             if (leave != j) {
                 sn_set.col(count) = sn_set_t.col(j);
                 count++;
             }
         }


         for (int number_r = 0; number_r < settings.r_isomap.size(); number_r++) {
             std::cout << " ISOMAP for " << settings.r_isomap[number_r] <<" dimension of the manifold" << std::endl;
             Eigen::MatrixXd snapshot_star_total = Eigen::MatrixXd::Zero(Nr, nC);
             Eigen::MatrixXd sn_set_partial(Nr, sn_set.cols());
             Eigen::VectorXd norma(sn_set.cols());
             Eigen::MatrixXd y(settings.r_isomap[number_r] * nC, sn_set.cols());
             Eigen::VectorXi k(nC);
             Eigen::VectorXi Rank = Eigen::VectorXi::Ones(nC);
             Eigen::VectorXd kruskal_min = Eigen::VectorXd::Ones(nC);


             //Loop for each conservative variable---> ISOMAP
             for (int q = 0; q < nC; q++) {
                 std::cout << "ISOMAP FOR CONSERVATIVE VARIABLE NUMBER   "<< q <<"...";
                 ISOMAP(sn_set, settings, Nr, q, k, kruskal_min, y, number_r, Rank);
             }

             std::cout << "DONE " << std::endl;
             std::cout << "Kruskal stress minimo \t " << kruskal_min << "\t per k=" << k << std::endl;
             std::cout << y << std::endl << std::endl << std::endl;
             std::cout << "END OF ISOMAP..." << std::endl << std::endl;


             std::cout << "------------START BACK-MAPPING SECTION--------------" << std::endl << std::endl;


             settings.param_rec_1[0]=param_t(leave,0);
             settings.param_rec_2[0]=param_t(leave,1);
             //devo escludere il parametro dal file per il leave one out
             Eigen::MatrixXd param(param_t.rows()-1, settings.Np);
             Eigen::RowVectorXd normes_parameters(settings.Np);
             count=0;
             for (int j=0; j<param_t.rows(); j++) {
                 if (leave != j) {
                     param(count,0)=param_t(j,0);
                     param(count,1)=param_t(j,1);
                     count ++;
                 }
             }
             std::cout<<param<<std::endl;
             for (int i=0; i<settings.Np; i++)
                 normes_parameters(i)= param.col(i).norm();


             std::cout<<"RECONSTRUCTING FIELD FOR:\t"<<settings.param_rec_1[0]<<" "<<settings.param_rec_2[0]<<std::endl;

                 for (int q = 0; q < nC; q++) {

                     sn_set_partial = sn_set.block(Nr * q, 0, Nr, sn_set.cols());
                     Eigen::MatrixXd y_part(settings.r_isomap[number_r], sn_set.cols());
                     for (int j = 0; j < settings.r_isomap[number_r]; j++) {
                         y_part.row(j) = y.row((q * settings.r_isomap[number_r]) + j);
                     }
                     //std::cout<<y_part<<std::endl;

                     std::vector<double> y_star = ISOMAP_surr_coeff(param, normes_parameters, y_part,
                                                                    settings.param_rec_1[0],
                                                                    settings.param_rec_2[0],
                                                                    settings.Np, settings.flag_interp,
                                                                    settings.analysis,
                                                                    supercritical);

                     std::cout << "values RBF" << std::endl;
                     for (int i = 0; i < y_star.size(); i++)
                         std::cout << y_star[i] << std::endl;

                     Eigen::MatrixXd snapshot_star = create_rec_fields(settings, y_star, y_part, sn_set_partial,number_r,q);
                     snapshot_star_total.col(q) = snapshot_star;
                     rec_test.block( q*Nr,leave,Nr,1) = snapshot_star;
                 }

                 //std::cout << "snap*" << std::endl << snapshot_star_total.topRows(5) << std::endl;
                 std::string filename;

                 //filename = write_reconstruction_file(snapshot_star_total, Coords, leave, nC, settings.r_isomap[number_r],
                       //                               settings, ID, flag);

                 //SU2_DTR(settings, filename, su2_conf, nt, number_r);
                 //RMS_residuals(settings, number_r, nt, Nr);
                 //Write_History_ResError_global(settings, number_r);


             std::cout << "-----END OF ISOMAP+BACK-MAPPING------" << std::endl << std::endl;
         }
       }
       diff_test= (sn_set_t - rec_test);
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
               //decision(q,j)= sqrt((sum));
               decision(q,j)=sqrt((sum))/sn_set_t.block(q*Nr,j,Nr,1).norm();
           }
       }
       std::cout<<"ISOMAP error reconstruction matrix"<<std::endl;
       std::cout<<decision<<std::endl;
       std::string filename;
       filename="ERROR_ISOMAP_leave_newset2D.csv";
       std::ofstream  flow_data;
       flow_data.open(filename);
       if (!flow_data.is_open()){
           std::cout<<"I can't open the Error file and write into it"<<std::endl;
           exit (EXIT_FAILURE);}
       flow_data<<"\"alpha\"" << ",";
       flow_data<<"\"Reynolds\"" << ",";
       flow_data<<"\"error_rec_Pressure\"" << ",";
       flow_data<<"\"error_rec_Cp\"" << ",";
       flow_data<<"\"error_rec_skin_friction_coeff_x\"" << ",";
       flow_data<<"\"error_rec_skin_friction_coeff_y\"" << ",";
       flow_data << std::endl;
       for (int i=0; i<param_t.rows(); i++) {
           for (int j = 0; j < param_t.cols(); j++)
               flow_data << std::setprecision(12) << std::scientific << param_t(i, j) << ", ";
           for( int j=0; j<nC; j++)
               flow_data << std::setprecision(12) << std::scientific << decision(j,i) << ", ";

           flow_data<< std::endl;
       }

       flow_data.close();


   }

    std::cout << std::endl;
    std::cout<<"------------------------------------------" <<std::endl;
    std::cout<<"------------Single-MODES end--------------" <<std::endl;
    std::cout<<"------------------------------------------" <<std::endl;

   return 0;


}
