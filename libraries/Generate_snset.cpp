#include "read_Inputs.hpp"
#include "Generate_snset.hpp"

Eigen::MatrixXd generate_snap_matrix( const int Nr, const int Ns, const int ds, const int init,
                                        std::vector<int> Cols,
                                        std::vector<std::string> inputfile,
                                        std::string flag_prob,
                                        std::string solver )
{

    std::string file_temp;
    int k = 0;

    if ( solver == "SU2")
    {   
        Eigen::MatrixXd field(Nr, Cols.size());                

        if (flag_prob == "VECTOR-2D")
        {
        
            Eigen::MatrixXd snap(2*Nr, inputfile.size());

            for( int i=0; i< inputfile.size(); i++ )
            {


                file_temp= inputfile[i];
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols);
                std::cout << "Complete!" << std::endl;


                Eigen::VectorXd gx = field.col(0);
                Eigen::VectorXd gy = field.col(1);

                snap.col(k) << gx,
                                gy;
                
                k++;
            }

            return snap;

        } else if ( flag_prob == "SCALAR" )
        {

            Eigen::MatrixXd snap(Cols.size()*Nr, inputfile.size());

            for( int i=0; i<inputfile.size(); i++)
            {


                file_temp=inputfile[i];
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols);
                std::cout << "Complete!" << std::endl;
                if (Cols.size()== 1)
                {
                  snap.col(k)<<field.col(0);

                }else if (Cols.size()== 2)
                {
                    snap.col(k) << field.col(0),
                                    field.col(1);

                }else if ( Cols.size()==3)
                {
                    snap.col(k) << field.col(0),
                                    field.col(1),
                                    field.col(2);

                }else if ( Cols.size() == 4 )
                {
                    snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3);

                } else if ( Cols.size() == 5 )
                {
                    snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3),
                            field.col(4);
                } else if ( Cols.size() == 6 )
                {
                    snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3),
                            field.col(4),
                            field.col(5);
                } else if ( Cols.size() == 7 )
                {
                    snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3),
                            field.col(4),
                            field.col(5),
                            field.col(6);
                } else if (Cols.size()== 8)
                {   snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3),
                            field.col(4),
                            field.col(5),
                            field.col(6),
                            field.col(7);
                }else
                {
                    std::cout << "Check the number of conservative variables in use " << std::endl;
                }

                k++;

            }

            return snap;


        } else if ( flag_prob == "CONSERVATIVE" )
        {

            Eigen::MatrixXd snap(Cols.size()*Nr, inputfile.size());

            for( int i=0; i<inputfile.size(); i++)
            {

                file_temp = inputfile[i];
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols);
                std::cout << "Complete!" << std::endl;

                if ( Cols.size() == 4 )
                {
                    snap.col(k) << field.col(0),
                                field.col(1),
                                field.col(2),
                                field.col(3);
		        } else if ( Cols.size() == 5 )
		        {
		            snap.col(k) << field.col(0),
                                field.col(1),
                                field.col(2),
                                field.col(3),
				                field.col(4);
		        } else if ( Cols.size() == 6 )
                {
                    snap.col(k) << field.col(0),
                                field.col(1),
                                field.col(2),
                                field.col(3),
                                field.col(4),
                                field.col(5);
                } else if ( Cols.size() == 7 )
                {
                    snap.col(k) << field.col(0),
                                field.col(1),
                                field.col(2),
                                field.col(3),
                                field.col(4),
                                field.col(5),
                                field.col(6);
                } else if (Cols.size()== 8)
                {   snap.col(k) << field.col(0),
                            field.col(1),
                            field.col(2),
                            field.col(3),
                            field.col(4),
                            field.col(5),
                            field.col(6),
                            field.col(7);
                }else
                {
                    std::cout << "Check the number of conservative variables in use " << std::endl;
                }
                            
                k++;

            }

            return snap;


        } else {

            std::cout << "Set well flag_prob! Now Exiting ..." << std::endl;
            exit (EXIT_FAILURE);

        }

    } else if ( solver == "CS3D")
    {/*
        int dum = 0;
        file_temp = std::to_string(init) + ".q";
        plot3d_info Info = read_plot3d_info (file_temp);

        for ( int iblock = 0; iblock < Info.nblocks; iblock++ )
        {
            dum += Info.ni[iblock]*Info.nj[iblock]*Info.nk[iblock];
        }

        if ( flag_prob == "VELOCITY-3D" )
        {
            
            Eigen::MatrixXd snap(3*dum, Ns);
            for( int i = init; i < (Ns*ds + init); i += ds )
            {
                // int k = 0;
                file_temp = std::to_string(i) + ".q"; 
                std::cout << "Reading fields from : " << file_temp << "\t";
                std::vector<Eigen::VectorXd> data_fields = read_plot3d (file_temp, Info);
                std::cout << "Complete!" << std::endl;

                int dum1 = 0;
                for ( int iblock = 0; iblock < Info.nblocks; iblock++ )
                {
                    int np_block = Info.ni[iblock]*Info.nj[iblock]*Info.nk[iblock];
                    snap.middleRows(dum1, np_block).col(k) = data_fields[iblock].middleRows(np_block, np_block).cwiseQuotient(data_fields[iblock].head(np_block));
                    snap.middleRows(dum + dum1, np_block).col(k) = data_fields[iblock].middleRows(2*np_block, np_block).cwiseQuotient(data_fields[iblock].head(np_block));
                    snap.middleRows(2*dum + dum1, np_block).col(k) = data_fields[iblock].middleRows(3*np_block, np_block).cwiseQuotient(data_fields[iblock].head(np_block));
                    dum1 += Info.ni[iblock]*Info.nj[iblock]*Info.nk[iblock];
                }

                k++;

            }

            return snap;
        } else {
            
            std::cout << "Only VELOCITY-3D implemented for CS3D so far. Set well FLAG_PROB! Now Exiting ..." << std::endl;
            exit (EXIT_FAILURE);
        }
    */} else {
        std::cout << "Set well problem solver " << std::endl;
        exit (EXIT_FAILURE);
    }
}




