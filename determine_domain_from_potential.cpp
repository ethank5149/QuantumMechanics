#include "determine_domain_from_potential.h"

std::pair<double, double> determine_domain_from_potential(std::string &Vfunc, double L) {
	std::vector<std::string> centered{"QHO", "ISW"};
	std::vector<std::string> shifted{};

	if (std::any_of(centered.begin(), centered.end(), [&](std::string s) {return Vfunc == s; })) {
		return std::make_pair(-0.5*L, 0.5*L);
	}
	else if (std::any_of(shifted.begin(), shifted.end(), [&](std::string s) {return Vfunc == s; })) {
		return std::make_pair(0.0, L);
	}
	else {
		std::cout << "Error! Potential Function Not Recognized!" << std::endl;
		return std::make_pair(-1.0, -1.0);
	}
}