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
	pressureFile.open("../pressure&density/pressure&density", std::ios::out);
	PkPfpPkFile.open("../pressure&density/PkPfpPk&density", std::ios::out);

	double rhoStep = 0.9 / 20;
	for (int rhoIndex = 1; rhoIndex <= 20; ++rhoIndex)
	{
		double rho = rhoStep * rhoIndex;
		int N = 512;
//	double rho = 0.8;
		double tMax = 10;
		double dt = 0.001;
		int NSteps = static_cast<int>(1. * tMax / dt);

		Solver solver(dt, 2, N, rho, -1);

		double pressure = 0;
		double PkPfpPk = 0;

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
				PkPfpPk += solver.getPkPfpPk();
			}

			if (i % (NSteps / 100) == 0)
			{
				std::cout << "Progress: " << 1. * i / NSteps * 100. << "%" << std::endl;
			}

			solver.step();
		}

		pressure /= NSteps / 2.;
		PkPfpPk /= NSteps / 2.;
		pressureFile << rho << "\t" << pressure << std::endl;
		PkPfpPkFile << rho << "\t" << PkPfpPk << std::endl;
	}
	energyFile.close();
	pressureFile.close();
	PkPfpPkFile.close();
	return 0;
}