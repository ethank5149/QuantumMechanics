///\file shooting-method.cpp
///\author Ethan Knox
///\date 8/2/2020.

#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <functional>

#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

#include "potentials.h"
#include "determine_domain_from_potential.h"

using namespace std::placeholders;

const double m = 1.0;
const double k = 1.0;
const double hbar = 1.0;
const double omega = 1.0;

typedef std::vector<double> state_t;
typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_t> stepper_rkf78_t;
typedef boost::numeric::odeint::runge_kutta4<state_t> rk4;

void TISE(state_t &psi, state_t &dpsi_dx, double x, std::function<double(double)> V, double E) {
    dpsi_dx[0] = psi[1];
    dpsi_dx[1] = (2.0 * m / pow(hbar, 2.0)) * (V(x) - E) * psi[0];
}

int eigenE_shooting(state_t& psi, double boundary_condition, std::string Vfunc, double L, double h_0, double eps_rel, double eps_abs, double E)
{
    // Potential
    std::function<double(double)> V;
    std::function<double(double)> dV;
    if (Vfunc == "QHO") {
        V = quantum_harmonic_oscillator;
        dV = d_quantum_harmonic_oscillator;
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
