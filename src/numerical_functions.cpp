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
    a0 = *f.begin();
    a1 = *df.begin();
    a2 = 3. / pow((tf - t0), 2) * (*f.end() - *f.begin()) - (2. * a1 + *df.end()) / (tf - t0);
    a3 = (*df.end() - a1) / (3. * pow((tf - t0), 2)) - (2. / 3.) * a2 / (tf - t0);
    // a3 = 2. / pow((tf - t0), 3) * (*f.begin() - *f.end()) + 2. * a1 / pow((tf - t0), 2) + (*df.end() - *df.begin()) / pow((tf - t0), 2);
    // a2 = 1. / (2 * (tf - t0)) * (*df.end() - *df.begin()) - 3. / 2. * a3 * (tf - t0);
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


