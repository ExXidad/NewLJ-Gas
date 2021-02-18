//
// Created by xidad on 12.02.2021.
//

#include "Solver.h"
#include <fstream>
#include <string>
#include <ios>

int main(int argc, char *argv[])
{
	std::fstream file;
	std::fstream energyFile;
	std::fstream pressureFile;
	std::fstream PkPfpPkFile;
	energyFile.open("../energy/energy", std::ios::out);
	pressureFile.open("../pressure&density/rCut&pressure", std::ios::out);

	double rCutStep = 7. / 20;
	for (int rCutIndex = 0; rCutIndex <= 20; ++rCutIndex)
	{
		double rCut = rCutStep * rCutIndex + 1.5;
		int N = 512;
		double rho = 0.7;
		double tMax = 10;
		double dt = 0.001;
		int NSteps = static_cast<int>(1. * tMax / dt);

		Solver solver(dt, 2, N, rho, -1, rCut);

		double pressure = 0;

		for (int i = 0; i <= NSteps; ++i)
		{
			if (dt * i >= tMax * 0.2)
			{
//			file.open("../pos&vel/" + std::to_string(i) + ".xyz", std::ios::out);
//			solver.xyzOut(file, true, true);
//			file.close();

//			solver.saveEnergyToFile(energyFile);
//			solver.savePressureToFile(pressureFile);

				pressure += solver.getPressure();
			}

			if (i % (NSteps / 100) == 0)
			{
				std::cout << "Progress: " << 1. * i / NSteps * 100. << "%" << std::endl;
			}

			solver.step();
		}

		pressure /= NSteps / 2.;
		pressureFile << rCut << "\t" << pressure << std::endl;
	}
	energyFile.close();
	pressureFile.close();
	return 0;
}