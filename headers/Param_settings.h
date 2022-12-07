#include <iostream>
#include "Eigen/Dense"
#include <vector>

class read_param
{
private:
    void read_file(std::string);
public:
    read_param();
    read_param(std::string);
    int N_el;
    std::vector<std::vector<double> > alphamat;
    std::vector<double> atomvec;
};


