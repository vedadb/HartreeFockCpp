#include <iostream>
#include "Eigen/Dense"
using Eigen::MatrixXd;

class SCF_cycle
{
private:

    

    void calcBasisMat(double*,double*,int,double*,int);
    void calcXHMMat(double*,double*,double*,double*,int,int);

public:
    SCF_cycle(double*,int,double*[3],int, int);
    void SCF_StartCalc(double*,int,double*[3],int);
    
};