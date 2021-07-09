//
// Created by giuseppe on 01/10/2020.
//
#include "Post-Process.hpp"



void Write_History_ResError_global(prob_settings settings, int number_r){

    //Defining Common Variables
    std::cout<<"Write global file starting from the each history file"<<std::endl;
    std::string header_history, value_history; //to write headers in global files
    std::ifstream history_su2;
    std::ofstream history_global;


    //create the name for the history global file
    std::string filename_history_global={};
    if(settings.flag_method[0]== "ISOMAP")
        filename_history_global= "history_global_rbm_"+settings.flag_method[0]+"_"+std::to_string(settings.r_isomap[number_r])+".csv";
    else if (settings.flag_method[0]=="POD")
        filename_history_global= "history_global_rbm_"+settings.flag_method[0]+"_"+std::to_string(settings.r[number_r])+".csv";


    //Getting filename
    std::string filename_history_su2 ={};
    history_global.open(filename_history_global);
    for (int nt=0; nt < settings.param_rec_1.size(); nt++) {
       if (settings.flag_method[0] == "ISOMAP") {
           filename_history_su2 = "history_" + settings.flag_method[0] + "_" + std::to_string(settings.r_isomap[number_r]) + "_" +
                                  std::to_string(nt + 1) + ".csv";
       }else if (settings.flag_method[0] == "POD"){
           filename_history_su2 = "history_" + settings.flag_method[0] + "_" + std::to_string(settings.r[number_r]) + "_" +
                                  std::to_string(nt + 1) + ".csv";
       }


        //se nt=0, apri il file globale e scrivi l'header
        history_su2.open(filename_history_su2);
        if(history_su2.is_open()){
            std::string linedata1, linedata2;
            //Getting row of headers
            getline(history_su2,linedata1);
            header_history = linedata1;
            //Getting values
            getline(history_su2,linedata2);
            value_history = linedata2;
        } else {
            std::cout << "Unable to open SU2 single history file. Exiting ..." << std::endl;
            exit(EXIT_FAILURE);
        }
        history_su2.close();
        //history_global.open(filename_history_global);
        if(history_global.is_open()){
            if( nt==0 ) {
                history_global << header_history << std::endl;
                history_global << std::to_string(nt+1) << "_" << std::to_string(settings.param_rec_1[nt]) << "_" << std::to_string(settings.param_rec_2[nt])
                               << "," << value_history << std::endl;
            }else{
                history_global << std::to_string(nt + 1) << "_" << std::to_string(settings.param_rec_1[nt]) << "_"
                               << std::to_string(settings.param_rec_2[nt])
                               << "," << value_history << std::endl;
            }
        }else{
            std::cout << "Unable to open global file for writing. Exiting... " << std::endl;
            exit(EXIT_FAILURE);
        }

        //String for removing useless solutions (single history files that now are included in the global one)
        std::string rmf_string = "rm -f " + filename_history_su2;
        int len_s = rmf_string.length();
        char rmf_sys_call[len_s + 20];
        strcpy(rmf_sys_call, rmf_string.c_str());
        std::system(rmf_sys_call);

        //string for removing the output.log from SU2_CFD / SU2_SOL---->if I want to see these files, comment these lines.
        std::string rmf_string_su2log= "rm -f SU2_"+settings.flag_method[0]+"_"+std::to_string(nt+1)+".log";
        len_s= rmf_string_su2log.length();
        rmf_sys_call[len_s + 20];
        strcpy(rmf_sys_call, rmf_string_su2log.c_str());
        std::system(rmf_sys_call);


        //string for removing the rec_flow file that comes out before SU2.
        std::string rmf_string_rec_flow= "rm -f rec_flow*.csv";
        len_s= rmf_string_rec_flow.length();
        rmf_sys_call[len_s + 20];
        strcpy(rmf_sys_call, rmf_string_rec_flow.c_str());
        std::system(rmf_sys_call);

    }


    history_global.close();
    std::cout<<"DONE"<<std::endl;
}


void RMS_residuals(prob_settings settings, int number_r, int nt, int Nr){

    //Defining Common Variables
    std::cout<<"Write the residual global file with the RMS of residual of each point of the mesh...";
    std::ifstream rec_flow_from_su2;
    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5)  << std::to_string(nt+1);

    //LEGGI I FILE DI P_SU2_rec_flow
    std::string filename_rec_flow_from_su2 ={};
    if (settings.flag_method[0] == "ISOMAP") {
            filename_rec_flow_from_su2 =
                    "P_SU2_rec_flow_" + std::to_string(settings.r_isomap[number_r]) + "_" + settings.flag_method[0] +
                    "_" + buffer.str() + ".csv";
    } else if (settings.flag_method[0] == "POD") {
            filename_rec_flow_from_su2 =
                    "P_SU2_rec_flow_" + std::to_string(settings.r[number_r]) + "_" + settings.flag_method[0] + "_" +
                    buffer.str() + ".csv";
    }


    // LEGGO IL FILE DI RESTART RIGA PER RIGA E LO SALVO IN UNA MATRICE EIGEN
    rec_flow_from_su2.open(filename_rec_flow_from_su2);
    if(!rec_flow_from_su2.is_open()){
        std::cout<<"unable to open the P_SU2_rec_flow"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line_flow_data;
    // Read row of headers--->understand the number of colums
    getline(rec_flow_from_su2, line_flow_data);
    int Ncol=1;
    char c;
    for(int i=0; i<line_flow_data.size(); i++) {
        c= line_flow_data[i];
        if ( c == ',')
        Ncol++;
    }

    Eigen::MatrixXd P_SU2_rec_flow(Nr,Ncol);
    int n_row=0;
    while(getline(rec_flow_from_su2,line_flow_data))
    {
        Eigen::RowVectorXd ID(Ncol);
        std::istringstream ss(line_flow_data);
        std::string token;
        long double value;
        int col =0;
        while (getline(ss,token, ','))
            {
                value = std::stold(token);
                ID(col)=value;
                col++;
            }

        P_SU2_rec_flow.row(n_row)= ID;
        n_row++;
    }

    rec_flow_from_su2.close();

    //CALCOLO DEL RMS DEI RESIDUI PER UN BOX SCELTO
    Eigen::RowVectorXd RMS = Eigen::RowVectorXd::Zero(settings.number_residuals.size());
    if( settings.ndim == 2 ){
        for (int j=0; j<settings.number_residuals.size(); j++){
            int id=1;
            for( int i=0; i<Nr; i++){
                if (settings.xmin <= P_SU2_rec_flow(i,1)  && P_SU2_rec_flow(i,1)<=settings.xmax) {
                    if (settings.ymin <= P_SU2_rec_flow(i,2) && P_SU2_rec_flow(i,2) <= settings.ymax) {
                        RMS(j) = RMS(j) + pow(P_SU2_rec_flow(i, settings.number_residuals[j]), 2);
                        id = id + 1;
                    }
                }
            }
            RMS(j)= sqrt(RMS(j)/id);
        }
    }else if (settings.ndim == 3){
        for (int j=0; j<settings.number_residuals.size(); j++){
            int id=1;
            for( int i=0; i<Nr; i++){
                if(settings.xmin <= P_SU2_rec_flow(i,1) && P_SU2_rec_flow(i,1)<=settings.xmax){
                    if (settings.ymin <= P_SU2_rec_flow(i,2) && P_SU2_rec_flow(i,2)<=settings.ymax) {
                        if (settings.zmin <= P_SU2_rec_flow(i,3) && P_SU2_rec_flow(i,3) <= settings.zmin) {
                            RMS(j) = RMS(j) + pow(P_SU2_rec_flow(i,settings.number_residuals[j]), 2);                                  //check the column that you really want
                            id = id + 1;
                        }
                    }
                }
            }
            RMS(j)= sqrt(RMS(j)/id);
        }

    }

    //SAVE THE RMS IN A FILE
    std::ofstream box_residualerror;

    //create the name for the history global file
    std::string filename_box_residualerror={};
    if(settings.flag_method[0]== "ISOMAP")
        filename_box_residualerror= "box_residual_error"+settings.flag_method[0]+"_"+std::to_string(settings.r_isomap[number_r])+".csv";
    else if (settings.flag_method[0]=="POD")
        filename_box_residualerror= "box_residual_error"+settings.flag_method[0]+"_"+std::to_string(settings.r[number_r])+".csv";

    if (nt==0)
        box_residualerror.open(filename_box_residualerror);
    else
        box_residualerror.open(filename_box_residualerror, std::ofstream::out | std::ofstream::app);

    if(!box_residualerror.is_open())
        std::cout<<"unable to open the file of box_residuals"<<std::endl;
    else{
        if (settings.flag_method[0] == "POD")
                box_residualerror << settings.r[number_r]<<"  "<< settings.param_rec_1[nt]<<"  "<< settings.param_rec_2[nt]<<"  "
                                  << RMS<<std::endl;
        else if (settings.flag_method[0] == "ISOMAP")
                box_residualerror << settings.r_isomap[number_r]<< "  " << settings.param_rec_1[nt]<<"  "<< settings.param_rec_2[nt]<< "  "
                                  << RMS<<std::endl;
    }



    //DELETE THE P_SU2_rec_flow----> to much memory required otherwise
    //std::string rmf_string = "rm -f P_SU2_rec_flow*";
    //int len_s = rmf_string.length();
    //char rmf_sys_call[len_s + 20];
    //strcpy(rmf_sys_call, rmf_string.c_str());
    //std::system(rmf_sys_call);


    std::cout<<"done"<<std::endl;
}



