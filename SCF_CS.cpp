#include <iostream>
#include "Eigen/Dense"
#include <cmath>
#include "SCF_CS.h"


SCF_CS::SCF_CS(){

    std::cout<<"Closed Shell Hartree Fock calculation initialized"<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::set_elnum(int elnum){
    N_el=elnum;
    std::cout<<"Number of electrons: "<<N_el<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::set_nucl(Eigen::MatrixXd* nuclref){

    nucl=*nuclref;

    std::cout<<"Nuclear coordinates: "<<nucl<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::set_alpha(Eigen::VectorXd* alpharef){

    alpha=*alpharef;

    std::cout<<"Gaussian coefficients: "<<alpha<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::SCF_SinglePoint(){

    std::cout<<"Closed Shell Hartree Fock Single Point Calculation Starting"<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
    SCF_CS::calcXHMMat();


};

void SCF_CS::calcXHMMat(){

    int n_exp=alpha.rows();
    int n_nucl=nucl.rows();
    int n_basis=n_exp*n_nucl;
    Eigen::MatrixXd BasisM(n_basis,4);

    for(int i=0;i<n_exp;i++){
        for(int j=0;j<n_nucl;j++){

            BasisM(i*n_nucl+j,0)=alpha(i);
            BasisM(i*n_nucl+j,1)=nucl(j,0);
            BasisM(i*n_nucl+j,2)=nucl(j,1);
            BasisM(i*n_nucl+j,3)=nucl(j,2);
        }
    }
   
    std::cout<<"Calculated Basis matrix"<<std::endl;


    Eigen::MatrixXd Sij(n_basis,n_basis);
    double a_i,a_j;

    for(int i=0;i<n_basis;i++){
        a_i=BasisM(i,0);
        for(int j=0;j<n_basis;j++){
            a_j=BasisM(j,0);
            Sij(i,j)=4*(a_i*a_j)/pow(a_i+a_j,2);
        }
    }

    std::cout<<Sij<<std::endl;
};