// #include <iostream>
// #include <string>
// #include <fstream>
// #include <cmath>
#include <algorithm>
#include <iterator>
// #include <vector>

#include "odes.h"
#include "import.h"
#include "numerical_functions.h"
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <fstream>

using boost::math::quadrature::tanh_sinh;
using namespace boost::math::quadrature;

int main(int argc, char const *argv[])
{
    // Print width and precision
    const int width = 20;
    const int prec = 9;

    // Import solution
    std::string filename{"/home/boris/Documenti/grid-precision/matlab/mat_results.dat"};
    std::vector<double> t;
    std::vector<std::vector<double>> X_MAT;
    std::vector<std::vector<double>> U_MAT;
    std::vector<double> p;
    read_data(filename, t, X_MAT, U_MAT, p);
    const int n_dis = X_MAT.size();
    const int n_state = X_MAT[0].size();
    int n_ctrl_dum;
    if (U_MAT.size() > 0)
    {
        n_ctrl_dum = U_MAT[0].size();
    }
    else
    {
        n_ctrl_dum = 0;
    }
    const int n_ctrl = n_ctrl_dum;

    // Compute derivatives
    std::vector<std::vector<double>> dX_MAT;
    for (int j = 0; j < n_dis; j++)
    {
        std::vector<double> x = X_MAT[j];
        std::vector<double> u;
        if (n_ctrl > 0)
        {
            u = U_MAT[j];
        }
        std::vector<double> dx = ode_fun(t[j], x, u, p);
        dX_MAT.push_back(dx);
    }

    // Interpolate solution
    // std::vector<std::vector<boost::math::cubic_b_spline<double>>> X_spline_MAT;
    std::vector<std::vector<Spline_t>> X_spline_MAT;
    for (int j = 0; j < n_dis - 1; j++)
    {
        if (j == 38)
        {
            double test = 0;
        }
        std::vector<Spline_t> X_spline_j;
        // std::vector<boost::math::cubic_b_spline<double>> X_spline_j;
        for (int k = 0; k < n_state; k++)
        {
            double x0 = X_MAT[j][k];
            double x1 = X_MAT[j + 1][k];
            std::vector<double> x_nodes{x0, x1};
            double t0 = t[j];
            double t1 = t[j + 1];
            // double h = t[j + 1] - t0;
            double dx0 = dX_MAT[j][k];
            double dx1 = dX_MAT[j + 1][k];
            std::vector<double> dx_nodes{dx0, dx1};
            Spline_t X_spline_jk{x_nodes, dx_nodes, t0, t1};
            // boost::math::cubic_b_spline<double> X_spline_jk(x_nodes.begin(), x_nodes.end(), t0, h, dx0, dx1);
            X_spline_j.push_back(X_spline_jk);
        }
        X_spline_MAT.push_back(X_spline_j);
    }

    std::string file_spline{"/home/boris/Documenti/grid-precision/results/spline00.dat"};
    std::ofstream ost_spline{file_spline};
    const int n_points_spline = 10001;
    for (int j = 0; j < n_points_spline; j++)
    {
        double time = t[38] + j * (t[38 + 1] - t[38]) / (n_points_spline - 1);
        ost_spline << std::setw(width) << std::scientific << std::setprecision(prec) << time
                   << std::setw(width) << std::scientific << std::setprecision(prec) << X_spline_MAT[38][4].eval(time) << std::endl;
    }

    // Assess precision
    std::string file_out{"/home/boris/Documenti/grid-precision/results/grid_precision.dat"};
    std::ofstream ost{file_out};

    ost << std::setw(width) << "Segment (t0)"
        << std::setw(width) << "Dx"
        << std::setw(width) << "Dy"
        << std::setw(width) << "Dz"
        << std::setw(width) << "Dvx"
        << std::setw(width) << "Dvy"
        << std::setw(width) << "Dvz" << std::endl;

    double tol = 1e-9;
    double err;
    for (int j = 0; j < n_dis - 1; j++)
    {
        if (j > 36)
        {
            double test = 0;
        }
        ost << std::setw(width) << std::scientific << std::setprecision(prec) << t[j];
        for (int k = 0; k < n_state; k++)
        {
            auto eps_tanh = [&, t, X_spline_MAT, U_MAT, p, j, k, n_state, n_ctrl](double time) {
                // Spline derivatives
                double dx_spline = X_spline_MAT[j][k].eval_prime(time);

                // State derivatives (computed on the spline)
                // Find closest nodes
                auto it = std::lower_bound(t.begin(), t.end(), time);
                int jj = std::distance(t.begin(), it) - 1;

                // State vector
                std::vector<double> x(n_state, 0.);
                for (int kk = 0; kk < n_state; kk++)
                {
                    // x[k] = X_spline_MAT[jj][k](time);
                    x[kk] = X_spline_MAT[jj][kk].eval(time);
                }

                // Control vector (linear)
                std::vector<double> u;
                if (n_ctrl > 0)
                {
                    double t0 = t[jj];
                    double t1 = t[jj + 1];
                    for (int kk = 0; kk < n_ctrl; kk++)
                    {
                        double u0 = U_MAT[jj][kk];
                        double u1 = U_MAT[jj + 1][kk];
                        u.push_back(linear_interpolation(time, u0, u1, t0, t1));
                    }
                }

                // Derivatives
                std::vector<double> dx = ode_fun(time, x, u, p);

                // Error
                double eps = dx_spline - dx[k];
                return fabs(eps);
            };
            // Numerical quadrature
            double eta_jk = gauss_kronrod<double, 15>::integrate(eps_tanh, t[j], t[j + 1], 15, tol, &err);
            //    tanh_sinh<double> tanh_integrator;
            //    double eta_jk = tanh_integrator.integrate(eps_tanh, t[j], t[j + 1]);
            //    double eta_jk = tanh_integrator.integrate(eps_tanh, t[j], t[j + 1], tol, &err);
            ost << std::setw(width) << std::scientific << std::setprecision(prec) << eta_jk;
        }
        ost << std::endl;
    }
    return 0;
}
