#include "odes.h"
#include <math.h>


std::vector<double> ode_fun(const double t, const std::vector<double> &x,
                            const std::vector<double> &u, const std::vector<double> &p)
{
    const double mu = 1.;
    const double r = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
    std::vector<double> dx(x.size(), 0.);
    dx[0] = x[3];
    dx[1] = x[4];
    dx[2] = x[5];
    dx[3] = -mu / (pow(r, 3)) * x[0];
    dx[4] = -mu / (pow(r, 3)) * x[1];
    dx[5] = -mu / (pow(r, 3)) * x[2];
    return dx;
}

// std::vector<double> ode_fun_tanh_vector(const double t)
// {
//     extern std::vector<double> t_in;
//     std::vector<double>::iterator it = std::lower_bound(t_in.begin(), t_in.end(), t);
//     int j = std::distance(t_in.begin(), it);


// }


double ode_fun_tanh(const double t)
{
    return 0.;
}
