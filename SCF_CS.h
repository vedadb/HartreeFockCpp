#include "Eigen/Dense"

class SCF_CS
{
private:

    

    
    void calcXHMMat(Eigen::MatrixXd*, Eigen::MatrixXd*,Eigen::MatrixXd*);
    
    int N_el;
    Eigen::MatrixXd nucl;
    Eigen::VectorXd alpha;
    Eigen::VectorXi Z_at;



public:

    SCF_CS();


    void set_nucl(Eigen::MatrixXd*);
    void set_alpha(Eigen::VectorXd*);
    void set_elnum(int);
    void set_Z(Eigen::VectorXi*);
    void SCF_SinglePoint();
};