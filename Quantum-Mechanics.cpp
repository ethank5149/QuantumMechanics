///\file Quamtum-Mechanics.cpp
///\author Ethan Knox
///\date 8/2/2020.

#include <iostream>
#include <fstream>
#include <limits>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <boost/program_options.hpp>

namespace opt = boost::program_options;

struct Params {
    const double m = 1.0;
    const double k = 1.0;
    const double hbar = 1.0;
    double E{};
    double (Params::* V)(double) {};
    double (Params::* dV)(double) {};

    double omega{};

    Params(double input_E) {
        E = input_E;
        omega = sqrt(k / m);
    }

    double qho(double x) { return 0.5 * m * pow(omega * x, 2.0); }
    double d_qho(double x) { return m * pow(omega, 2.0) * x; }

    double isw(double x) { return 0.0; }
    double d_isw(double x) { return 0.0; }
};


int tise(double x, const double psi[], double dpsi_dx[], void* params)
{
    struct Params* p = (Params*)params;

    dpsi_dx[0] = psi[1];
    dpsi_dx[1] = (2.0 * p->m / pow(p->hbar, 2.0)) * ((p->*(p->V))(x) - p->E) * psi[0];

    return GSL_SUCCESS;
}


int jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
    struct Params* p = (Params*)params;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix* mat = &dfdy_mat.matrix;

    gsl_matrix_set(mat, 0, 0, 0.0);
    gsl_matrix_set(mat, 0, 1, 1.0);
    gsl_matrix_set(mat, 1, 0, (2.0 * p->m / pow(p->hbar, 2)) * ((p->*(p->V))(t) - p->E));
    gsl_matrix_set(mat, 1, 1, 0.0);

    dfdt[0] = 0.0;
    dfdt[1] = (p->*(p->dV))(t);

    return GSL_SUCCESS;
}


int main(int argc, const char* argv[])
{
    long nsteps;
    double L, psi0, dpsi0, E, h_0, eps_abs, eps_rel;
    char Vfunc;

    opt::options_description params("Simulation Parameters");

    params.add_options()
        ("help,h", "Show usage")
        ("nsteps,n", opt::value<long>(&nsteps)->default_value(1000), "Number of descretization steps")
        ("L,l", opt::value<double>(&L)->default_value(10.0), "Characteristic length")
        ("psi0", opt::value<double>(&psi0)->default_value(0.0), "Initial value")
        ("dpsi0", opt::value<double>(&dpsi0)->default_value(0.0), "Initial value")
        ("E,e", opt::value<double>(&E)->default_value(0.5), "Energy")
        ("V,v", opt::value<char>(&Vfunc)->default_value('Q'), "Potential function")
        ("h_0,h", opt::value<double>(&h_0)->default_value(1.0e-6), "Initial step size")
        ("eps_abs", opt::value<double>(&eps_abs)->default_value(1.0e-6), "Allowed absolute error")
        ("eps_rel", opt::value<double>(&eps_rel)->default_value(0.0), "Allowed relative error")
        ;

    opt::variables_map vm;
    opt::store(opt::parse_command_line(argc, argv, params), vm);

    if (vm.count("help")) {
        std::cout << params << std::endl;
        return 1;
    }
    else {
        opt::notify(vm);

        Params p(E);
        gsl_odeiv2_system sys = { tise, jac, 2, &p };
        gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h_0, eps_abs, eps_rel);

        p.V = (Vfunc == 'Q') ? &Params::qho : &Params::isw;
        p.dV = (Vfunc == 'Q') ? &Params::d_qho : &Params::d_isw;
        
        const double xi = -0.5 * L;
        const double xf = 0.5 * L;
        double x_current = xi;
        double x_next;
        int status;

        psi0 = (psi0 == 0.0) ? DBL_EPSILON : psi0;
        dpsi0 = (dpsi0 == 0.0) ? DBL_EPSILON : dpsi0;

        double psi[2] = { psi0, dpsi0 };
        double* psi_series = (double *)malloc(nsteps * sizeof(double));
        double* psi2_series = (double*)malloc(nsteps * sizeof(double));
        double* x_series = (double*)malloc(nsteps * sizeof(double));

        psi_series[0] = psi[0];
        psi2_series[0] = pow(psi[0], 2.0);
        x_series[0] = x_current;

        std::ofstream file("output.dat");
        file << "# E = " << E << "\n";
        file << "# x, V(x), E + psi, E + |psi|^2\n";

        for (int i = 1; i < nsteps; i++)
        {
            x_next = xi + i * (xf - xi) / (nsteps - 1.0);

            status = gsl_odeiv2_driver_apply(d, &x_current, x_next, psi);

            if (status != GSL_SUCCESS)
            {
                std::cout << "Error!, return value = " << status << std::endl;
                break;
            }
            
            psi_series[i] = psi[0];
            psi2_series[i] = pow(psi[0], 2.0);
            x_series[i] = x_current;
        }

        gsl_odeiv2_driver_free(d);


        gsl_interp_accel* acc = gsl_interp_accel_alloc();
        gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
        gsl_spline_init(spline, x_series, psi2_series, nsteps);

        double normalization = gsl_spline_eval_integ(spline, xi, xf, acc);

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        for (int i = 0; i < nsteps; i++)
        {
            file << std::scientific << x_series[i] << ", " << (p.*(p.V))(x_series[i]) << ", " << (&p)->E + psi_series[i] / sqrt(normalization) << ", " << (&p)->E + psi2_series[i] / normalization << "\n";
        }

        return 0;
    }
}


