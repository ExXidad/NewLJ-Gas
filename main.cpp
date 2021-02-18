//
// Created by xidad on 12.02.2021.
//

#include "Solver.h"
#include <fstream>
#include <string>
#include <ios>

int main(int argc, char *argv[])
{
	int N = 729;
	double rho = 0.8;
	double tMax = 30;
	double dt = 0.001;
	int NSteps = static_cast<int>(1. * tMax / dt);
	Solver solver(dt, 2, N, rho, -1);

	std::fstream file;
	std::fstream energyFile;
	std::fstream pressureFile;
	energyFile.open("../energy/energy", std::ios::out);

	for (int i = 0; i <= NSteps; ++i)
	{
		if (dt * i >= 0)
		{
			file.open("../pos&vel/" + std::to_string(i) + ".xyz", std::ios::out);
			solver.xyzOut(file, true, true);
			file.close();

			solver.saveEnergyToFile(energyFile);
		}

		if (i % (NSteps / 100) == 0)
		{
			std::cout << "Progress: " << 1. * i / NSteps * 100. << "%" << std::endl;
		}

		solver.step();
	}

	energyFile.close();
	pressureFile.close();
	return 0;
}