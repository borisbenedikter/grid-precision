#pragma once

#include <vector>

std::vector<double> ode_fun(const double t, const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &p);
double ode_fun_tanh(const double t);