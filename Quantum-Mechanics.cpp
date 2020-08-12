///\file Quamtum-Mechanics.cpp
///\author Ethan Knox
///\date 8/2/2020.

#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include <gsl/gsl_const_mksa.h>

#include "Quantum-Potentials.h"
#include "determine_domain_from_potential.h"

using namespace std::placeholders;
namespace opt = boost::program_options;

typedef std::vector<double> state_t;
typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_t> stepper_rkf78_t;
typedef boost::numeric::odeint::runge_kutta4<state_t> rk4;


void assign_potential(std::function<double(double)> &V, std::function<double(double)> &dV, std::string Vfunc) {
    if (Vfunc == "QHO") {
        V = std::bind(quantum_harmonic_oscillator, _1, m, omega);
        dV = std::bind(d_quantum_harmonic_oscillator, _1, m, omega);
    }
    else { // Vfunc == "ISW"
        V = infinite_square_well;
        dV = d_infinite_square_well;
    }
}


int main(int argc, const char* argv[])
{
    double L, psi0, dpsi0, E, h_0, eps_abs, eps_rel;
    std::string Vfunc, filename, units;

    opt::options_description params("Simulation Parameters");

    params.add_options()
        ("help,h", "Show usage")
        ("L,l", opt::value<double>(&L)->default_value(10.0), "Characteristic length")
        ("psi0", opt::value<double>(&psi0)->default_value(1.0e-8), "Initial value")
        ("dpsi0", opt::value<double>(&dpsi0)->default_value(1.0e-8), "Initial value")
        ("E,e", opt::value<double>(&E)->default_value(6.5), "Energy")
        ("V,v", opt::value<std::string>(&Vfunc)->default_value("QHO"), "Potential function")
        ("h_0,h", opt::value<double>(&h_0)->default_value(1.0e-2), "Initial step size")
        ("eps_abs", opt::value<double>(&eps_abs)->default_value(1.0e-12), "Allowed absolute error")
        ("eps_rel", opt::value<double>(&eps_rel)->default_value(0.0), "Allowed relative error")
        ("filename,f", opt::value<std::string>(&filename)->default_value("output.dat"), "Output file name")
        ("units,u", opt::value<std::string>(&units)->default_value("natural"), "Use 'natural' or 'mks' units.")
        ;

    opt::variables_map vm;
    opt::store(opt::parse_command_line(argc, argv, params), vm);

    if (vm.count("help")) {
        std::cout << params << std::endl;
        return 1;
    }
    else {
        opt::notify(vm);

        // Units
        if (units == "natural") {
            double m = 1.0;
            double k = 1.0;
            double hbar = 1.0;
            double omega = 1.0;
        }
        else { // units == 'mks'
            double m = GSL_CONST_MKSA_MASS_ELECTRON;
            // GSL_CONST_MKSA_MASS_PROTON
            // GSL_CONST_MKSA_MASS_NEUTRON

            double k = 158.2;
            /* k \approx \lambda \cdot 2 r_0
            | Material | Configuration | Young's Modulus [GPa] | Atomic Radii [pm] | k [N/m] |
            | -------- | ------------- | --------------------- | ----------------- | ------- |
            | Mg       |               | 45                    | 150               | 13.5    |
            | Al       |               | 69                    | 125               | 17.25   |
            | Ti       |               | 110.3                 | 140               | 30.884  |
            | Cu       |               | 117                   | 135               | 31.59   |
            | Si       | (Crystal)     | 130-185               | 110               | 34.65   |
            | Be       |               | 287                   | 105               | 60.27   |
            | Mo       |               | 329-330               | 145               | 95.555  |
            | W        |               | 400-410               | 135               | 109.35  |
            | Os       |               | 525-562               | 130               | 141.31  |
            | C        | (Graphene)    | 1050                  | 70                | 147     |
            | C        | (Diamond)     | 1050-1210             | 70                | 158.2   |
            | C        | (Carbyne)     | 32100                 | 70                | 294     |
            */
            double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
            double omega = sqrt(k / m);
        }

        // Potential
        std::function<double(double)> V;
        std::function<double(double)> dV;
        assign_potential(V, dV, Vfunc);

        // Domain
        std::pair<double, double> domain = assign_domain(Vfunc, L);
        const double xi = domain.first;
        const double xf = domain.second;
        double x = xi;

        // Initial Condition
        state_t psi{ psi0, dpsi0 };
        state_t dpsi_dx(2);

        // Data Series
        std::vector<double> psi_series;
        std::vector<double> dpsi_series;
        std::vector<double> x_series;
        // Data Recorder
        auto observer = [&](const state_t& psi, const double x) {
            psi_series.push_back(psi[0]);
            dpsi_series.push_back(psi[1]);
            x_series.push_back(x);
        };

        // ODE System
        std::function<void(const state_t&, state_t&, double)> sys = [&](const state_t& psi, state_t& dpsi_dx, const double x) {
            dpsi_dx[0] = psi[1];
            dpsi_dx[1] = (2.0 * m / pow(hbar, 2.0)) * (V(x) - E) * psi[0];
        };

        boost::numeric::odeint::integrate_const(rk4(), sys, psi, x, xf, h_0, observer);
        //boost::numeric::odeint::integrate_adaptive(make_controlled(eps_abs, eps_rel, stepper_rkf78_t()), sys, psi, x, xf, h_0, observer);

        std::vector<double> psi_series_copy(psi_series);
        std::vector<double> dpsi_series_copy(dpsi_series);
        std::vector<double> x_series_copy(x_series);

        auto spline = boost::math::interpolators::cubic_hermite(std::move(x_series_copy), std::move(psi_series_copy), std::move(dpsi_series_copy));
        double norm = sqrt(boost::math::quadrature::trapezoidal([&](double x) {return pow(spline(x), 2.0); }, xi, xf));

        std::ofstream file(filename);
        file << "# E = " << E << "\n";
        file << "# x, V(x), psi, |psi|^2\n";

        for (size_t i = 0; i < x_series.size(); i++) {
            psi_series[i] /= norm;
            file << std::scientific << x_series[i] << " " << V(x_series[i]) << " " << psi_series[i] + E << " " << pow(psi_series[i], 2.0) + E << "\n";
        }
    };
}

///////////////////////////////////////////////////////////////////////////////

///\file Stationary-Quantum-EigenE.cpp
///\author Ethan Knox
///\date 8/2/2020.

#include "pch.h"
#include "framework.h"

// TODO: This is an example of a library function
void fnStationaryQuantumEigenE()
{
}

#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <functional>

#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

#include "determine_domain_from_potential.h"

using namespace std::placeholders;

const double m = 1.0;
const double k = 1.0;
const double hbar = 1.0;
const double omega = 1.0;

typedef std::vector<double> state_t;
typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_t> stepper_rkf78_t;
typedef boost::numeric::odeint::runge_kutta4<state_t> rk4;

void TISE(state_t& psi, state_t& dpsi_dx, double x, std::function<double(double)> V, double E) {
    dpsi_dx[0] = psi[1];
    dpsi_dx[1] = (2.0 * m / pow(hbar, 2.0)) * (V(x) - E) * psi[0];
}

int eigenE_shooting(state_t& psi, double boundary_condition, std::string Vfunc, double L, double h_0, double eps_rel, double eps_abs, double E)
{
    // Potential
    std::function<double(double)> V;
    std::function<double(double)> dV;
    if (Vfunc == "QHO") {
        V = std::bind(quantum_harmonic_oscillator, _1, m, omega);
        dV = std::bind(d_quantum_harmonic_oscillator, _1, m, omega);
    }
    else { // Vfunc == "ISW"
        V = infinite_square_well;
        dV = d_infinite_square_well;
    }

    // Domain
    std::pair<double, double> domain = determine_domain_from_potential(Vfunc, L);
    const double xi = domain.first;
    const double xf = domain.second;
    double x = xi;

    // Initial Condition
    state_t psi_copy(psi);
    state_t dpsi_dx(2);
    double h_0 = 1.0e-4;

    std::function<double(double)> defect = [x, xf, h_0, psi_copy, &V](double E) {
        boost::numeric::odeint::integrate_const(rk4(), std::bind(TISE, _1, _2, _3, V, E), psi_copy, x, xf, h_0);
        return psi_copy[0];
    };
    boost::math::tools::eps_tolerance<double> tol;
    std::pair<double, double> r = boost::math::tools::bisect(defect, 0.0, 10.0, tol);
    return r.first + (r.second - r.first) / 2;
}
