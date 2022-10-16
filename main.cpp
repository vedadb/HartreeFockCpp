#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include "SCF_CS.h"


int main(){

    Eigen::MatrixXd nucl {{1,0,0},{2,3,4},{5,6,7}};
    Eigen::VectorXd alpha {{0.233, 1.33}};
    int N_el=2;

    SCF_CS calc;
    calc.set_elnum(N_el);
    calc.set_nucl(&nucl);
    calc.set_alpha(&alpha);

    calc.SCF_SinglePoint();


    return 0;
}