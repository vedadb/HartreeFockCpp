#include "headers/Param_settings.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "headers/Eigen/Dense"


const std::string WHITESPACE = " \n\r\t\f\v";
 
std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

//https://stackoverflow.com/questions/10058606/splitting-a-string-by-a-character
#include <vector>
#include <sstream>




using namespace std;

read_param::read_param(){
    read_file("param.txt");
}

read_param::read_param(std::string filename){
    read_file(filename);

}

void read_param::read_file(std::string filename){
    string line;

   // cout<<"inside constructor"<<endl;

    ifstream myfile(filename);

    const char comment='/';
    const char blocksgn='#';

    string optPar="";

    bool Nel=false;
    bool Nexp=false;
    bool Natoms=false;
    vector<double> seglist;

    if (myfile.is_open()){
        while (getline(myfile,line))
        {
            
            if(line[0]==comment){
               // cout<<"Found a comment"<<endl;
            }
            else if(line[0]==blocksgn){
                    string newline (line.substr(1));
                    //cout<<newline<<endl;

                    if(newline.compare("end")){
                        optPar=newline;
                        if(!optPar.compare("num_el")){
                        
                        if (Nel==false){
                                Nel=true;
                            }
                            else{
                                 cout<<"Error, number of electrons already defined"<<endl;
                            }
                        }
                        else if(!optPar.compare("exponents")){
                            
                            if (Nexp==false){
                                    Nexp=true;
                                }
                                else{
                                    cout<<"Error, exponent already defined"<<endl;//Error, throw error msg
                                }
                        }
                        else if(!optPar.compare("atoms")){
                            
                            if (Natoms==false){
                                    Natoms=true;
                                }
                                else{
                                    cout<<"Error, atoms already defined"<<endl;
                                }
                        }
                        else{
                            cout<<"Error, unknown command"<<endl;
                        }
                    }else{
                    
                        if(!optPar.compare("num_el")){
                            if(seglist.size()==1){
                                N_el=seglist[0];
                            }
                        }
                        else if(!optPar.compare("exponents")){
                            if(seglist.size()>0){
                                alphavec=seglist;
                            }
                        }
                        else if(!optPar.compare("atoms")){
                            if(seglist.size()>0 && seglist.size()%4==0){
                                atomvec=seglist;
                            }
                        }   

                        seglist.clear();
                    }
            }
            else{
                if(!optPar.empty()){
                    if(!trim(line).empty()){
                        string newline=trim(line);
                        stringstream ss;
                        ss<<newline;
                       
                        string segment;
                        while(getline(ss,segment, ','))
                        {
                            //cout<<segment<<endl;
                            seglist.push_back(stod(segment));
                        }
                        if(!optPar.compare("atoms") && seglist.size()%4!=0){
                            cout<<"wrong number of elements in row"<<endl;
                            //throw error
                        }
                        // cout<<"size of seglist "<<seglist.size()<<endl;
                    }
                }


            }
        }

    }
    myfile.close();
}
