#pragma once

#include <vector>

double linear_interpolation(double t, double x0, double xf, double t0, double tf);

class Spline_t
{
private:
    double a0;
    double a1;
    double a2;
    double a3;
    double t0;
public:
    Spline_t(const std::vector<double> f, const std::vector<double> df, const double t_0, const double tf);
    ~Spline_t();

    double eval(const double t) const;
    double eval_prime(const double t) const;
};

