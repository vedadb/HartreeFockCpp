#include <iostream>
#include <cmath>
#include "headers/Eigen/Dense"
#include "headers/SCF_CS.h"
#include "headers/Param_settings.h"

int main(){

    read_param parameters;
    Eigen::VectorXd alpha(parameters.alphamat[0].size());

    for(int i=0;i<parameters.alphamat[0].size();i++)
        alpha(i)=parameters.alphamat[0][i];

    SCF_CS calc;

    calc.set_elnum(parameters.N_el);
    calc.set_atoms(&parameters.atomvec);
    // calc.set_alpha(&alpha);
    calc.set_alphamat(&parameters.alphamat);

    // auto t1=std::chrono::high_resolution_clock::now();
    // calc.SCF_SinglePoint();
    // auto t2=std::chrono::high_resolution_clock::now();
    
    // std::chrono::duration<double,std::milli> ms_double= t2-t1;
    // std::cout<<"Time elapsed: "<<ms_double.count()<<"ms"<<std::endl;
    return 0;
}