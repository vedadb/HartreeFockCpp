#include <iostream>
#include "SCF.h"
#include "Eigen/Dense"
 
using Eigen::MatrixXd;

int main(){

    double alpha[]={0.233136, 1.309757};
    double nucl[3][3] = { 
                          
                          {-1, 0, 0} , 
                          {0 , 1, 0},
                          {1,0,0}
    };

    int N_el = 2;
    int N_exp= sizeof(alpha)/sizeof(double);
    int N_nucl=sizeof(nucl)/sizeof(nucl[0]);

    SCF_cycle calc(alpha,N_exp,&(nucl[0][0]),N_nucl,N_el);
    calc.SCF_StartCalc(alpha,N_exp,&(nucl)[3],N_nucl);

    return 0;
}