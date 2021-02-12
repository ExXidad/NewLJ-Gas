//
// Created by xidad on 12.02.2021.
//

#include "Solver.h"
#include <fstream>
#include <string>
#include <ios>

int main(int argc, char *argv[])
{
	int N = 64;
	double L = 4;
	double dt = 0.01;
	int NSteps = 2000;
	Solver solver(N, L, dt);

	std::fstream file;

	for (int i = 0; i < NSteps; ++i)
	{
		file.open("../results/" + std::to_string(i)+".xyz", std::ios::out);
		solver.xyzOut(file, true, false);
		file.close();

		solver.step();
	}
}