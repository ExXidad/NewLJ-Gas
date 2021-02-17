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
	double L = pow(N / rho,0.333333);
	double tMax = 10;
	double dt = 0.001;
	int NSteps = tMax / dt;
	Solver solver(N, L, dt, 2);

	std::fstream file;
	std::fstream energyFile;
	energyFile.open("../energy/energy", std::ios::out);

	for (int i = 0; i <= NSteps; ++i)
	{
		if (dt * i >= 0)
		{
			file.open("../pos&vel/" + std::to_string(i) + ".xyz", std::ios::out);
			solver.xyzOut(file, true, false);
			file.close();

			solver.saveEnergyToFile(energyFile);
		}

		if (i % (NSteps / 100) == 0)
		{
			std::cout << "Progress: " << 1.* i / NSteps*100. << "%" << std::endl;
		}

		solver.step();
	}

	energyFile.close();
	return 0;
}