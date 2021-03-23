//
// Created by xidad on 12.02.2021.
//

#include "Solver.h"
#include <fstream>
#include <string>
#include <ios>
#include <cstring>
#include <iomanip>

/* Prints usage information */
void usage()
{
	fprintf(stdout, "mdlj usage:\n");
	fprintf(stdout, "mdlj [options]\n\n");
	fprintf(stdout, "Options:\n");
	std::cout << "It's mandatory to set 2 quantities of the following:\n";
	fprintf(stdout, "\t -N [integer]\t\tNumber of particles\n");
	fprintf(stdout, "\t -rho [real]\t\tDensity\n");
	fprintf(stdout, "\t -L [real]\t\tBox side length\n");
	std::cout << std::endl;
	fprintf(stdout, "\t -rCut [real]\t\tCutoff radius, base = 2.5\n");
	fprintf(stdout, "\t -T0 [real]\t\tInitial temperature, base = 1\n");
	fprintf(stdout, "\t -dt [real]\t\tTime step, base = 0.001\n");
	fprintf(stdout, "\t -tMax [real]\t\tIntegrate to this time, base = 0\n");
	fprintf(stdout, "\t -saveStartTime [real]\t\tStart saving values from this time, base = 0\n");
	fprintf(stdout, "\t -saveStep [integer]\t\tSample frequency, base=1\n");
	fprintf(stdout, "\t -saveVelocities\t\tInclude velocities in output file, base = 1\n");
	fprintf(stdout, "\t -savePressure\t\tInclude pressure in output file, base = 0\n");
	fprintf(stdout, "\t -saveEnergy\t\tInclude energy in output file, base = 1\n");
//	fprintf(stdout, "\t -Tb [real]\t\tBerendsen therm. temperature\n");
//	fprintf(stdout, "\t -tau [real]\t\tBerendsen therm. freq. \n");
//	fprintf(stdout, "\t -sf [a|w]\t\tAppend or write config output file\n");
//	fprintf(stdout, "\t -icf [string]\t\tInitial configuration file\n");
//	fprintf(stdout, "\t -seed [integer]\tRandom number generator seed\n");
	fprintf(stdout, "\t -resetDirs\t\tDeletes and creates new dirs for coordinates, energy, pressure\n");
	fprintf(stdout, "\t -unfold\t\tPrint unfolded coordinates in output files\n");
	fprintf(stdout, "\t -h\t\tPrint this info\n");
}

int main(int argc, char *argv[])
{
	int N = -1;
	double rho = -1;
	double tMax = 0;
	double L = -1;
	double dt = 0.001;
	double T0 = 1;
	double rCut = 2.5;
	int unfold = 0;
	int saveStep = 1;
	int saveVelocities = 1;
	double saveStartTime = 0;
	int saveEnergy = 1;
	int savePressure = 0;

	/* Here we parse the command line arguments;  If
	 you add an option, document it in the usage() function! */
	for (int i = 1; i < argc; ++i)
	{
		if (!strcmp(argv[i], "-N")) N = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-rho")) rho = atof(argv[++i]);
		else if (!strcmp(argv[i], "-L")) L = atof(argv[++i]);
//		else if (!strcmp(argv[i], "-T")) T = atof(argv[++i]);
		else if (!strcmp(argv[i], "-dt")) dt = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rCut")) rCut = atof(argv[++i]);
//		else if (!strcmp(argv[i], "-ns")) nSteps = atoi(argv[++i]);
//		else if (!strcmp(argv[i], "-so")) short_out = 1;
		else if (!strcmp(argv[i], "-T0")) T0 = atof(argv[++i]);
		else if (!strcmp(argv[i], "-saveVelocities")) saveVelocities = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-savePressure")) savePressure = atoi(argv[++i]);
//		else if (!strcmp(argv[i], "-Tb")) Tb = atoi(argv[++i]);
//		else if (!strcmp(argv[i], "-tau")) tau = atof(argv[++i]);
//		else if (!strcmp(argv[i], "-fs")) fSamp = atoi(argv[++i]);
//		else if (!strcmp(argv[i], "-sf")) wrt_code_str = argv[++i];
//		else if (!strcmp(argv[i], "-icf")) init_cfg_file = argv[++i];
//		else if (!strcmp(argv[i], "-ecorr")) use_e_corr = 1;
//		else if (!strcmp(argv[i], "-seed")) Seed = (unsigned long) atoi(argv[++i]);
		else if (!strcmp(argv[i], "-unfold")) unfold = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-tMax")) tMax = atof(argv[++i]);
		else if (!strcmp(argv[i], "-saveStep")) saveStep = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-saveStartTime")) saveStartTime = atof(argv[++i]);
		else if (!strcmp(argv[i], "-saveEnergy")) saveEnergy = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-resetDirs"))
		{
			system("bash add_rm_pressure");
			system("bash add_rm_pv");
			system("bash add_rm_energy");
			exit(0);
		} else if (!strcmp(argv[i], "-h"))
		{
			usage();
			exit(0);
		} else
		{
			fprintf(stderr, "Error: Command-line argument '%s' not recognized.\n",
					argv[i]);
			exit(-1);
		}
	}

//	N = 512;
//	rho = 0.7;
//	tMax = 10;
//	dt = 0.001;
//	T0 = 2;
//	rCut = 3;
	int NSteps = static_cast<int>(1. * tMax / dt);
	if (saveStep < 0)
	{
		saveStep = NSteps / 100 == 0 ? 1 : NSteps / 100;
	}
	int saveStartTimeStep = static_cast<int>(1. * saveStartTime / dt);

	Solver solver(dt, T0, N, rho, L, rCut);

	std::cout << "Solver configuration:"
			  << "\nN: " << N
			  << "\nL: " << L
			  << "\nrho: " << rho
			  << "\nT0: " << T0
			  << "\nrCut: " << rCut
			  << "\n" << std::endl;

	std::cout << "Time & save configuration: "
			  << "\ndt: " << dt
			  << "\ntMax: " << tMax
			  << "\nNSteps: " << NSteps
			  << "\nSave starts at t = " << saveStartTime << " or step = " << saveStartTimeStep
			  << "\nSaving every " << saveStep << " steps"
			  << (saveVelocities ? "\nVelocities are included" : "\nVelocities are not included")
			  << (unfold ? "\nCoordinates are unfolded" : "\nCoordinates are folded")
			  << (saveEnergy ? "\nEnergy is included" : "\nEnergy is not included")
			  << (savePressure ? "\nPressure is included" : "\nPressure is not included")
			  << "\n" << std::endl;

	std::fstream file;
	std::fstream energyFile;
	std::fstream pressureFile;

	if (saveEnergy)
		energyFile.open("./energy/energy", std::ios::out);

	if (saveEnergy)
		pressureFile.open("./pressure&density/pressure&density", std::ios::out);


	for (int i = 0; i <= NSteps; ++i)
	{
		if (i >= saveStartTimeStep)
		{
			if (i % saveStep == 0)
			{
				std::cout << "Progress: " << 1. * i / NSteps * 100. << "%" << std::endl;

				file.open("./pos&vel/" + std::to_string(i) + ".xyz", std::ios::out);
				solver.xyzOut(file, saveVelocities, unfold);
				file.close();

				if (saveEnergy)
					solver.saveEnergyToFile(energyFile);

				if (savePressure)
					solver.savePressureToFile(pressureFile);
			}
		}


		solver.step();
	}

	if (saveEnergy)
		energyFile.close();

	if (savePressure)
		pressureFile.close();
	return 0;
}