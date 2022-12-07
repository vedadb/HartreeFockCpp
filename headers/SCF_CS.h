#include "Eigen/Dense"
#include <vector>
class SCF_CS
{
private:

    

    
    void calcXHMMat(Eigen::MatrixXd*, Eigen::MatrixXd*,Eigen::MatrixXd*);
    
    int N_el;
    Eigen::MatrixXd nucl;
    Eigen::MatrixXd alphamat;
    Eigen::VectorXd alpha;
    Eigen::VectorXi Z_at;



public:

    SCF_CS();
    void set_alphamat(std::vector<std::vector<double> >*);
    void set_atoms(std::vector<double>*);

    void set_alpha(Eigen::VectorXd*);

    void set_elnum(int);
    void SCF_SinglePoint();
};