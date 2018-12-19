#include "numerical_functions.h"
#include "math.h"


double linear_interpolation(double t, double x0, double xf, double t0, double tf)
{
    double x = x0 + (xf - x0) / (tf - t0) * (t - t0);
    return x;
}



Spline_t::Spline_t(const std::vector<double> f, const std::vector<double> df, const double t_0, const double tf)
{
    t0 = t_0;
    double S0 = f[0];
    double dS0 = df[0];
    double Sf = f[1];
    double dSf = df[1];
    a0 = S0;
    a1 = dS0;
    a2 = 3. / pow((tf - t0), 2) * (Sf - S0) - (2. * a1 + dSf) / (tf - t0);
    a3 = (dSf - a1) / (3. * pow((tf - t0), 2)) - (2. / 3.) * a2 / (tf - t0);
}

Spline_t::~Spline_t()
{
}

double Spline_t::eval(const double t) const
{
    return a0 + a1 * (t - t0) + a2 * pow((t - t0), 2) + a3 * pow((t - t0), 3);
}

double Spline_t::eval_prime(const double t) const
{
    return a1 + 2 * a2 * (t - t0) + 3 * a3 * pow((t - t0), 2);
}


