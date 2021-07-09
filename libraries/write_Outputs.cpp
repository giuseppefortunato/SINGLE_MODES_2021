#include "write_Outputs.hpp"


void Config_stream ( prob_settings settings )
{

    std::cout << " Number of snapshots : " << settings.Ns << std::endl;
    std::cout << " Delta T between snapshots ( only EQUISPACED for now) : " << settings.Dt_cfd*settings.Ds << std::endl;
    std::cout << " Starting snapshot index : " << settings.nstart << std::endl;
    std::cout << " Data-type to be processed : " << settings.flag_prob << std::endl << std::endl;

    std::cout << "----------- Performing " << settings.flag_method[0][0] << " --------------" << std::endl;
    std::cout << " Subtracting mean from snapshot = " << settings.flag_mean << std::endl << std::endl;

    if ( settings.flag_method[0] == "SPOD" )
    {
        std::cout << " Filter size : " << settings.Nf << std::endl;
        std::cout << " Filter type : " << settings.flag_filter << std::endl;
        std::cout << " Energy level desired : " << settings.En*100 << "%" << std::endl;
    }   

    if ( settings.flag_method[0] == "DMD" || settings.flag_method[0] == "fbDMD" || settings.flag_method[0] == "HODMD" || settings.flag_method[0] == "mrDMD" || settings.flag_method[0] == "RDMD")
    {
        std::cout << " Rank for the reduced dynamic (if -1 pick all modes, if 0 do SVHT) : " << settings.r[0] << std::endl;
        std::cout << " Method to compute coefficients : " << settings.dmd_coef_flag << std::endl;
    }

    if ( settings.flag_method[0] == "RDMD" )
    {
        std::cout << " Rank for the RDMD (if 0 do the recursion until the desired energetic content ( max = 3*Ns)) : " << settings.r_RDMD << std::endl;
        std::cout << " Energy level desired (to be considered only if rank = 0): " << settings.En*100 << "%" << std::endl;
    }

    if ( settings.flag_method[0] == "mrDMD" )
    {
        std::cout << " Max levels of the multi resolution : " << settings.max_levels << std::endl;
        std::cout << " Number of samples per time bin : " << settings.max_cycles << std::endl;
    }

    if ( settings.flag_method[0] == "HODMD" )
    {
        std::cout << " Number of levels for high order : " << settings.d << std::endl;
    }

    std::cout << std::endl;

}



//SingleModes
string write_reconstruction_file (const Eigen::MatrixXd &Rec,
                                  const Eigen::MatrixXd &Coords,
                                  int nt,
                                  int nC,
                                  const int Nm,
                                  prob_settings settings,
                                  Eigen::VectorXi &ID,
                                  bool flag) {

    std::cout << "creating reconstructed fields file...";
    int Nr = Coords.rows();
    std::string root_outfile;
    root_outfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string out_format;
    //out_format.assign ( settings.out_file, settings.out_file.size() - 3, 3);
    if (!flag)
        out_format.assign( settings.test_file[nt], settings.test_file[nt].size() - 9, 9 );
    else
        out_format.assign( settings.in_file[nt], settings.in_file[nt].size() - 9, 9 );

    //std::stringstream buffer;
    //buffer << std::setfill('0') << std::setw(5)  << std::to_string(nt+1);
    std::string filename;
    if(settings.flag_method[0]== "POD") {
        filename = root_outfile + "_" + std::to_string(Nm)+ "_" + settings.flag_method[0] + "_"  + out_format;
    }else if (settings.flag_method[0] =="ISOMAP"){
        //filename = root_outfile + "_"+std::to_string(Nm)+"_" + settings.flag_method[0] + "_" + buffer.str() + "." + out_format;
        filename = root_outfile + "_"+std::to_string(Nm)+"_" + settings.flag_method[0] + "_" + out_format;
    }
    std::cout<<filename<<std::endl;
    std::ofstream  flow_data;
    flow_data.open(filename);
    if (!flow_data.is_open()){
        std::cout<<"I can't open the recostruction file and write into it"<<std::endl;
        exit (EXIT_FAILURE);}
    if (Coords.cols()==1){                          //case for advection diffusion equation
        // Write row of Headers
        flow_data << "\"PointID\"" << ",";
        flow_data << "\"x\"" << ",";
        for(int i=0; i < settings.Cols.size(); i++){
            switch (settings.Cols[i]){
                case 1:
                    flow_data<< "\"Rec_density\"" << ",";
                    break;
                case 2:
                    flow_data<< "\"Rec_momentum_x\"" << ",";
                    break;
                case 3:
                    flow_data<< "\"Rec_energy_t\"" << ",";
                    break;
                case 4:
                    flow_data<< "\"Rec_velocity\"" << ",";
                    break;
                case 5:
                    flow_data<< "\"Rec_Pressure\"" << ",";
                    break;
                case 6:
                    flow_data<< "\"Rec_entropy\"" << ",";
                    break;
                case 7:
                    flow_data<<"\"Rec_internal_energy\"" << ",";
                    break;

            }
        }
        flow_data << std::endl;

        //Write fields
        for(int i=0; i<Nr; i++){
            flow_data << ID(i) << ", ";
            for(int j=0; j<Coords.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Coords(i,j)  << ", ";
            for (int j=0; j<Rec.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Rec(i,j)  << ", " ;

            flow_data<< std::endl;
        }
        flow_data.close();

    }
    if (Coords.cols()==2){
        // Write row of Headers
        flow_data << "\"PointID\"" << ",";
        flow_data << "\"x\"" << ",";
        flow_data << "\"y\"" << ",";
        for(int i=0; i < settings.Cols.size(); i++){
            switch (settings.Cols[i]){
                case 3:
                    flow_data<< "\"Rec_density\"" << ",";
                break;
                case 4:
                    flow_data<< "\"Rec_MomentumX\"" << ",";
                break;
                case 5:
                    flow_data<< "\"Rec_MomentumY\"" << ",";
                break;
                case 6:
                    flow_data<< "\"Rec_Energy\"" << ",";
                break;
                case 7:
                    flow_data<< "\"Rec_Turb_var1\"" << ",";
                break;
                case 8:
                    flow_data<< "\"Rec_Pressure\"" << ",";
                break;
                case 9:
                    flow_data<< "\"Rec_Temperature\"" << ",";
                break;
                case 10:
                    flow_data<< "\"Rec_Mach\"" << ",";
                break;
                case 11:
                    flow_data<< "\"Rec_Pressure_coeff\"" << ",";
                break;
                case 12:
                    flow_data<< "\"Rec_laminar_visc\"" << ",";
                break;
                case 13:
                    flow_data<< "\"rec_skin_friction_x\"" << ",";
                break;
                case 14:
                    flow_data<< "\"rec_skin_friction_y\"" << ",";
                break;
            }
        }
        flow_data << std::endl;

        //Write fields
        for(int i=0; i<Nr; i++){
            flow_data << ID(i) << ", ";
            for(int j=0; j<Coords.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Coords(i,j)  << ", ";
            for (int j=0; j<Rec.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Rec(i,j)  << ", " ;

            flow_data<< std::endl;
        }
        flow_data.close();

    }else if (Coords.cols()==3){
        // Write row of Headers
        flow_data << "\"PointID\"" << ",";
        flow_data << "\"x\"" << ",";
        flow_data << "\"y\"" << ",";
        flow_data << "\"z\"" << ",";
        for(int i=0; i < settings.Cols.size(); i++){
            switch (settings.Cols[i]){
                case 4:
                    flow_data<< "\"Rec_density\"" << ",";
                    break;
                case 5:
                    flow_data<< "\"Rec_MomentumX\"" << ",";
                    break;
                case 6:
                    flow_data<< "\"Rec_MomentumY\"" << ",";
                    break;
                case 7:
                    flow_data<< "\"Rec_MomentumZ\"" << ",";
                    break;
                case 8:
                    flow_data<< "\"Rec_Energy\"" << ",";
                    break;
                case 9:
                    flow_data<< "\"Rec_Turb_var1\"" << ",";
                    break;
                case 10:
                    flow_data <<"\"Rec_Turb_var2\"" << ",";
                break;
                case 13:
                    flow_data <<"\"Rec_pressure_coeff\"" << ",";
                break;
                case 15:
                    flow_data <<"\"Rec_skin_fric_x\"" << ",";
                break;
                case 16:
                    flow_data <<"\"Rec_skin_fric_y\"" << ",";
                break;
                case 17:
                    flow_data <<"\"Rec_skin_fric_z\"" << ",";
                break;
            }
        }
        flow_data << std::endl;

        //Write fields
        for(int i=0; i<Nr; i++){
            flow_data << ID(i) << ", ";
            for(int j=0; j<Coords.cols(); j++){
                flow_data << std::setprecision(12) << std::scientific << Coords(i,j)  << ", ";
            }
            for (int j=0; j<Rec.cols(); j++)
                flow_data << std::setprecision(12) << std::scientific << Rec(i,j)  << ", " ;

            flow_data<< std::endl;
        }
        flow_data.close();
    }

    std::cout<<"done"<<std::endl;

return filename;
}





