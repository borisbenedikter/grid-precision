#include "import.h"
#include <fstream>

void read_data(std::string filename, std::vector<double> &t, std::vector<std::vector<double>> &X_MAT, std::vector<std::vector<double>> &U_MAT, std::vector<double> &p)
{
    std::ifstream ist{filename};

    const int n_state = 6;
    const int n_ctrl = 0;
    const int n_param = 0;

    // Skip the first line
    const int n_lines_to_skip = 1;
    for(int i = 0; i < n_lines_to_skip; i++)
    {
        std::string line;
        std::getline(ist, line);
    }

    double dum;
    t = std::vector<double>();
    X_MAT = std::vector<std::vector<double>>();

    while(true){
        std::vector<double> vec_dum;
        if (ist >> dum)
        {
            t.push_back(dum);
            for (int k = 0; k < n_state; k++)
            {
                ist >> dum;
                vec_dum.push_back(dum);
            }
            X_MAT.push_back(vec_dum);
        }
        else
        {
            break;
        }
    }
}
