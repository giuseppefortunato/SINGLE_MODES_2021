
#ifndef MODES_POST_PROCESS_HPP
#define MODES_POST_PROCESS_HPP

#include "read_Inputs.hpp"



void Write_History_ResError_global(prob_settings settings, int number_r);

void RMS_residuals(prob_settings settings, int number_r, int nt, int Nr);


#endif //MODES_POST_PROCESS_HPP
