#include "Eigen/Dense"

class SCF_CS
{
private:

    

    
    void calcXHMMat();

    int N_el;
    Eigen::MatrixXd nucl;
    Eigen::VectorXd alpha;


public:

    SCF_CS();


    void set_nucl(Eigen::MatrixXd*);
    void set_alpha(Eigen::VectorXd*);
    void set_elnum(int);

    void SCF_SinglePoint();
};