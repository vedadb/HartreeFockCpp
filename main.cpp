#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include "SCF_CS.h"


int main(){

    Eigen::MatrixXd nucl {{1,0,0},{0,1,0},{-1,0,0},{0,0,0}};
    Eigen::VectorXd alpha {{0.233136, 1.309757, 2, 3, 4}};
    Eigen::VectorXi Z_el {{1,1,1,1,1}};
    int N_el=2;

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