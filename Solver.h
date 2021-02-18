//
// Created by xidad on 12.02.2021.
//

#ifndef NEWLJ_GAS_SOLVER_H
#define NEWLJ_GAS_SOLVER_H

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Solver
{
private:
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

	// Amount of particles
	int N;
	bool initFileVelocitiesIncluded;

	// Cell size
	double rho;
	double L, V;

	// Coordinates, velocities and forces arrays
	double *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz;

	// Boundary crossing counter arrays
	int *xBCCounter, *yBCCounter, *zBCCounter;

	// Time step, squared time step, internal time counter
	double dt, dt2, t;

	double rCut2, energyCut = 0;

	// Energies, virial
	double PE, KE, TotalE, virial = 0;

	// Initial temperature
	double T0;


private:
	void initializeArrays();

	void updateForces();

	void initializeSystem();

public:
	Solver(std::fstream &initFile, const double &L, const double &dt);

	Solver(const double &dt,
		   const double &T0 = 1,
		   const int &N = -1,
		   const double &rho = -1,
		   const double &L = -1,
		   const double &rCut = 10
	);

	~Solver();

	void step();

	void saveEnergyToFile(std::fstream &file);

	void savePressureToFile(std::fstream &file);

	void xyzOut(std::fstream &file, bool includeVelocities, bool unfold);

	void xyzIn(std::fstream &fstream);

	double getPressure();

	double getPkPfpPk();
};


#endif //NEWLJ_GAS_SOLVER_