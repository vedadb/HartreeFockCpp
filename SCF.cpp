#include "SCF.h"
#include <iostream>
#include "Eigen/Dense"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;



SCF_cycle::SCF_cycle(double *alpha, int n_exp, double *(nucl)[3], int n_nucl, int n_el){

    std::cout << "Initialized SCF"<<std::endl;
    std::cout<<"-------------------"<<std::endl;
    std::cout<<"The given exponentials are: "<<n_exp<<std::endl;

    for(int i=0;i<n_exp;i++){
        
        std::cout<<i+1<<": "<<*(alpha+i)<<std::endl;
        
    };
    std::cout<<"-------------------"<<std::endl;
    std::cout<<"The given nuclei are: "<<n_nucl<<std::endl;
    std::cout<<"(x),(y),(z) "<<std::endl;
    for(int i=0;i<n_nucl;i++){

        std::cout<<i+1<<": "<<*(nucl+i)<<","<<*(nucl+1+i*3)<<","<<*(nucl+2+3*i)<<std::endl;
    }; 
    std::cout<<"-------------------"<<std::endl;
    std::cout<<"The given number of electrons are: "<<n_el<<std::endl;
    std::cout<<"-------------------"<<std::endl;
};





void SCF_cycle::SCF_StartCalc(double *alpha, int n_exp, double *nucl[3], int n_nucl){
    int n_basis=n_exp*n_nucl;
    double BasisMat[n_basis][4];
    double XMat[n_basis][n_basis];
    double HMat[n_basis][n_basis];
    double MMat[n_basis*n_basis][n_basis*n_basis];

    std::cout << "Initialized B,C,H,M matrices"<<std::endl;
    std::cout<<"calculating basis matrix"<<std::endl;
    calcBasisMat(&(BasisMat[0][0]), alpha, n_exp,nucl,n_nucl);
    std::cout<<"calculating X-H-M matrices"<<std::endl;
    calcXHMMat(&(BasisMat[0][0]),&(XMat[0][0]),&(HMat[0][0]),&(MMat[0][0]),n_exp,n_nucl);
};






void SCF_cycle::calcBasisMat(double* BasisMat, double *alpha, int n_exp, double *nucl[3], int n_nucl)
{

   // for(int i=0;i<n_exp;i++){
   //     for(int j=0;j<n_nucl;j++){
   //         *(BasisMat+i*n_nucl+j)=*(alpha+i);
   //         *(BasisMat+i*n_nucl+j+1)=*(nucl+j*3);
   //         *(BasisMat+i*n_nucl+j+2)=*(nucl+1+j*3);
   //         *(BasisMat+i*n_nucl+j+3)=*(nucl+2+j*3);

   //         std::cout<<*(BasisMat+i*n_nucl+j)<<","<<*(BasisMat+i*n_nucl+j+1)<<","<<*(BasisMat+i*n_nucl+j+2)<<","<<*(BasisMat+i*n_nucl+j+3)<<","<<std::endl;
   //     }
   // }
    
};

void SCF_cycle::calcXHMMat(double* BasisMat,double* XMat,double* HMat,double* MMat,int n_exp,int n_nucl)
{   
    int n_basis=n_exp*n_nucl;
    MatrixXd SMat(n_basis,n_basis);
    MatrixXd TMat(n_basis,n_basis);
    MatrixXd VMat(n_basis,n_basis);
    std::cout<<"Number of basis functions: "<<n_basis<<std::endl;

    double a_i, a_j;
    VectorXd R_i(3), R_j(3);

    for(int i=0;i<n_basis;i++){
        a_i=*(BasisMat+i);
        for(int k=0;k<3;k++){
            R_i(k)=*(BasisMat+k+1+i*4);
        }
        
        std::cout<<R_i<<std::endl;

        for(int j=0;j<n_basis;j++){
            a_j=*(BasisMat+j);

            SMat(i,j)=pow(4*a_i*a_j/pow(a_i+a_j,2),0.75)*exp(-a_i*a_j/(a_i+a_j));
        
        }   
    }

    //std::cout<<SMat<<std::endl;
};