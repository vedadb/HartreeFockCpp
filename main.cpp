#include <iostream>
#include <cmath>
#include "headers/Eigen/Dense"
#include "headers/SCF_CS.h"
#include "headers/Param_settings.h"

int main(){

    read_param parameters;
    std::cout<<parameters.alphavec.size()<<std::endl;

    int N_el=parameters.N_el;
    Eigen::MatrixXd nucl(parameters.atomvec.size()/4,3);
    Eigen::VectorXd alpha(parameters.alphavec.size());
    Eigen::VectorXi Z_el(parameters.atomvec.size()/4);

    for(int i=0;i<parameters.alphavec.size();i++)
        alpha(i)=parameters.alphavec[i];
    for(int i=0;i<parameters.atomvec.size()/4;i++){
        Z_el(i)=parameters.atomvec[i*4];
        for(int j=0;j<3;j++)
            nucl(i,j)=parameters.atomvec[i*4+j];
    }

    SCF_CS calc;

    calc.set_elnum(N_el);
    calc.set_nucl(&nucl);
    calc.set_alpha(&alpha);
    calc.set_Z(&Z_el);

    auto t1=std::chrono::high_resolution_clock::now();
    calc.SCF_SinglePoint();
    auto t2=std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double,std::milli> ms_double= t2-t1;
    std::cout<<"Time elapsed: "<<ms_double.count()<<"ms"<<std::endl;
    return 0;
}