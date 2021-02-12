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
	double L = 1;

	// Coordinates, velocities and forces arrays
	double *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz;

	// Boundary crossing counter arrays
	int *xBCCounter, *yBCCounter, *zBCCounter;

	// Time step, squared time step
	double dt, dt2;

	double rCut2 = 2.5, energyCut = 0;

	double PE, KE, TotalE, T0 = -1;


public:
	Solver(std::fstream &initFile, const double &L, const double &dt);
	Solver(const int &N, const double &L, const double &dt);

	~Solver();

	void step();

	void initializeArrays();

	void updateForces();

	void xyzOut(std::fstream &file, bool includeVelocities, bool unfold);

	void xyzIn(std::fstream &fstream);

};


#endif //NEWLJ_GAS_SOLVER_H
// what is e_corr & tail correction