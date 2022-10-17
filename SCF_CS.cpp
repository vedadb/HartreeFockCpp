#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <cmath>
#include "SCF_CS.h"
#define _USE_MATH_DEFINES
#include "unsupported/Eigen/MatrixFunctions"

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

    std::cout<<"Nuclear coordinates:\n "<<nucl<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::set_alpha(Eigen::VectorXd* alpharef){

    alpha=*alpharef;

    std::cout<<"Gaussian coefficients:\n "<<alpha.transpose()<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
};

void SCF_CS::SCF_SinglePoint(){

    std::cout<<"Closed Shell Hartree Fock Single Point Calculation Starting"<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
    
    int n_exp=alpha.rows();
    int n_nucl=nucl.rows();
    int n_basis=n_exp*n_nucl;

    Eigen::MatrixXd Hij(n_basis,n_basis);
    Eigen::MatrixXd Xij(n_basis,n_basis);
    Eigen::MatrixXd Mijkl(n_basis*n_basis,n_basis*n_basis);

 

    SCF_CS::calcXHMMat(&Xij, &Hij, &Mijkl);

    Eigen::MatrixXd Pn=Eigen::MatrixXd::Identity(n_basis,n_basis);
    double Ee, Etot;
    double error=1;
    double tol=1e-12;
    int iter=0;
    while(error>tol)
    {
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        Eigen::MatrixXd P(n_basis,n_basis);
        Eigen::MatrixXd G(n_basis,n_basis);
        Eigen::MatrixXd Ft(n_basis,n_basis);
        
        P=Pn;
        for(int i=0;i<n_basis;i++){
            for(int j=0;j<n_basis;j++){
                G(i,j)=0;
                for(int k=0;k<n_basis;k++){
                    for(int l=0;l<n_basis;l++){
                        G(i,j)=G(i,j)+P(k,l)*Mijkl(i*n_basis+j,k*n_basis+l);
                    }
                }
            }
        }
        
        Ft=Xij*(Hij+G)*Xij;
        es.compute(Ft);
        Eigen::VectorXd eval=es.eigenvalues();
        Eigen::MatrixXd evec=es.eigenvectors();

        Eigen::MatrixXd Cij=Xij*evec;
        Eigen::MatrixXd Cijt(n_basis,N_el/2);

        
        for(int i=0;i<n_basis;i++){
            for(int j=0;j<N_el/2;j++){
                Cijt(i,j)=Cij(i,j);
            }
        }
        
        Pn=2*Cijt*(Cijt.transpose());
        error=0;
        Ee=0;
        for(int i=0;i<n_basis;i++){
            for(int j=0;j<n_basis;j++){   
                error=pow(P(i,j)-Pn(i,j),2);
                Ee=Ee+0.5*(Pn(i,j)*(2*Hij(i,j)+G(i,j)));
            }
        }
        error=sqrt(error);
        Etot=Ee;
        for(int i=0;i<n_nucl;i++){
            for(int j=i+1;j<n_nucl;j++){
                Etot=Etot+pow(pow(nucl(i,0)-nucl(j,0),2)+pow(nucl(i,1)-nucl(j,1),2)+pow(nucl(i,2)-nucl(j,2),2),-0.5);
            }
        }
        std::cout<<"Iteration: "<<++iter<<": E_el = "<<Ee<<": E_tot = "<<Etot<<std::endl;
         
    }



};

void SCF_CS::calcXHMMat(Eigen::MatrixXd *Xij,Eigen::MatrixXd *Hij,Eigen::MatrixXd *Mijkl){

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

    Eigen::MatrixXd Sij(n_basis,n_basis);
   
    {
        Eigen::MatrixXd Tij(n_basis,n_basis);
        Eigen::MatrixXd Vij(n_basis,n_basis);
        Eigen::Vector3d Rp;

        double a_i,a_j,RiRj,RiRc,a;

        for(int i=0;i<n_basis;i++){
            a_i=BasisM(i,0);

            Sij(i,i)=1;
            Tij(i,i)=1.5*a_i;
            Vij(i,i)=0;
            
            for(int k=0;k<n_nucl;k++){
                
                RiRc=sqrt(pow(BasisM(i,1)-nucl(k,0),2)+pow(BasisM(i,2)-nucl(k,1),2)+pow(BasisM(i,3)-nucl(k,2),2));
                
                if(RiRc==0)
                    Vij(i,i)=Vij(i,i)-2*sqrt(2*a_i/M_PI);
                else
                    Vij(i,i)=Vij(i,i)-1/(RiRc)*erf(sqrt(2*a_i)*RiRc);

            }
            

            for(int j=i+1;j<n_basis;j++){
                a_j=BasisM(j,0);
                a=a_i*a_j/(a_i+a_j);
                
                RiRj=pow(BasisM(i,1)-BasisM(j,1),2)+pow(BasisM(i,2)-BasisM(j,2),2)+pow(BasisM(i,3)-BasisM(j,3),2);
                
                Sij(j,i)=pow(4*a/(a_i+a_j),0.75)*exp(-a*RiRj);
                Sij(i,j)=Sij(j,i);

                Tij(i,j)=Sij(i,j)*a*(3-2*a*RiRj);
                Tij(j,i)=Tij(i,j);
                
                for(int l=1;l<4;l++){
                    Rp(l-1)=(a_i*BasisM(i,l)+a_j*BasisM(j,l))/(a_i+a_j);
                }
                
                for(int k=0;k<n_nucl;k++){
                    
                    RiRc=sqrt(pow(Rp(0)-nucl(k,0),2)+pow(Rp(1)-nucl(k,1),2)+pow(Rp(2)-nucl(k,2),2));
                    
                    if(RiRc==0)
                        Vij(i,j)=Vij(i,j)-2*Sij(i,j)*sqrt((a_i+a_j)/M_PI);
                    else
                        Vij(i,j)=Vij(i,j)-Sij(i,j)/(RiRc)*erf(sqrt(a_i+a_j)*RiRc);
                } 

                Vij(j,i)=Vij(i,j);
                
            }  
        }

        *(Hij)=Tij+Vij;      
    }
    
    {
        Eigen::MatrixXd Mijkl1(n_basis*n_basis,n_basis*n_basis); 
        Eigen::MatrixXd M0ijkl(n_basis*n_basis,n_basis*n_basis);
        Eigen::Vector3d Rp,Rq;
        double a_i,a_j,a_k,a_l,c,RpRq;

        for(int i=0;i<n_basis;i++){
            a_i=BasisM(i,0);
            for(int j=0;j<n_basis;j++){
                a_j=BasisM(j,0);

                
                Rp(0)=(a_i*BasisM(i,1)+a_j*BasisM(j,1))/(a_i+a_j);
                Rp(1)=(a_i*BasisM(i,2)+a_j*BasisM(j,2))/(a_i+a_j);
                Rp(2)=(a_i*BasisM(i,3)+a_j*BasisM(j,3))/(a_i+a_j);

                for(int k=0;k<n_basis;k++){
                    a_k=BasisM(k,0);
                    for(int l=0;l<n_basis;l++){
                        a_l=BasisM(l,0);
                        c=sqrt((a_i+a_j)*(a_k+a_l)/(a_i+a_j+a_k+a_l));

                        Rq(0)=(a_k*BasisM(k,1)+a_l*BasisM(l,1))/(a_k+a_l);
                        Rq(1)=(a_k*BasisM(k,2)+a_l*BasisM(l,2))/(a_k+a_l);
                        Rq(2)=(a_k*BasisM(k,3)+a_l*BasisM(l,3))/(a_k+a_l);

                        RpRq=sqrt(pow(Rp(0)-Rq(0),2)+pow(Rp(1)-Rq(1),2)+pow(Rp(2)-Rq(2),2));
                        if(RpRq==0)
                            M0ijkl(i*n_basis+j,k*n_basis+l)=2*Sij(i,j)*Sij(k,l)*c/sqrt(M_PI);
                        else
                            M0ijkl(i*n_basis+j,k*n_basis+l)=Sij(i,j)*Sij(k,l)/RpRq*erf(c*RpRq);

                            
                    }
                }
            }
        }

        for(int i=0;i<n_basis;i++){
            for(int j=0;j<n_basis;j++){
                for(int k=0;k<n_basis;k++){
                    for(int l=0;l<n_basis;l++){
                        Mijkl1(i*n_basis+j,k*n_basis+l)=M0ijkl(i*n_basis+j,k*n_basis+l)-0.5*M0ijkl(i*n_basis+k,j*n_basis+l);
                    }
                }
            }
        }
        *(Mijkl)=Mijkl1;
    }

    *(Xij)=Sij.pow(-0.5);
};