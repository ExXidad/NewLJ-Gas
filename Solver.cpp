//
// Created by xidad on 12.02.2021.
//

#include "Solver.h"

void Solver::step()
{
//	std::cout << "before step" << std::endl;
//	for (int i = 0; i < 5; ++i)
//	{
//		std::cout << x[i] << "\t" <<
//				  y[i] << "\t" <<
//				  z[i] << "\t" <<
//				  vx[i] << "\t" <<
//				  vy[i] << "\t" <<
//				  vz[i] << "\t" <<
//				  fx[i] << "\t" <<
//				  fy[i] << "\t" <<
//				  fz[i] << "\t" << std::endl;
//	}
	t += dt;
	for (int i = 0; i < N; ++i)
	{
		// First integration step
		x[i] += vx[i] * dt + 0.5 * dt2 * fx[i];
		y[i] += vy[i] * dt + 0.5 * dt2 * fy[i];
		z[i] += vz[i] * dt + 0.5 * dt2 * fz[i];

		vx[i] += 0.5 * fx[i] * dt;
		vy[i] += 0.5 * fy[i] * dt;
		vz[i] += 0.5 * fz[i] * dt;

		/* Apply periodic boundary conditions */
		if (x[i] < 0.0)
		{
			x[i] += L;
			xBCCounter[i]--;
		}
		if (x[i] > L)
		{
			x[i] -= L;
			xBCCounter[i]++;
		}
		if (y[i] < 0.0)
		{
			y[i] += L;
			yBCCounter[i]--;
		}
		if (y[i] > L)
		{
			y[i] -= L;
			yBCCounter[i]++;
		}
		if (z[i] < 0.0)
		{
			z[i] += L;
			zBCCounter[i]--;
		}
		if (z[i] > L)
		{
			z[i] -= L;
			zBCCounter[i]++;
		}
	}
	updateForces();

	// Second integration step
	KE = 0.0;
	for (int i = 0; i < N; ++i)
	{
		vx[i] += 0.5 * dt * fx[i];
		vy[i] += 0.5 * dt * fy[i];
		vz[i] += 0.5 * dt * fz[i];
		KE += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	}
	KE *= 0.5;
	TotalE = KE + PE;
}

Solver::~Solver()
{
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] vx;
	delete[] vy;
	delete[] vz;
	delete[] fx;
	delete[] fy;
	delete[] fz;
	delete[] xBCCounter;
	delete[] yBCCounter;
	delete[] zBCCounter;
}

void Solver::updateForces()
{
	// Zero all the forces
	for (int i = 0; i < N; ++i)
	{
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
	}

	double dx, dy, dz;
	double r2, r6Inverse, hL = L / 2.0;
	double energy = 0, force;

	for (int i = 0; i < (N - 1); ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			dx = (x[i] - x[j]);
			dy = (y[i] - y[j]);
			dz = (z[i] - z[j]);

			// Minimum image convention
			if (dx > hL) dx -= L;
			else if (dx < -hL) dx += L;
			if (dy > hL) dy -= L;
			else if (dy < -hL) dy += L;
			if (dz > hL) dz -= L;
			else if (dz < -hL) dz += L;

			// Update forces & compute energy
			r2 = dx * dx + dy * dy + dz * dz;
			if (r2 < rCut2)
			{
				r6Inverse = 1.0 / (r2 * r2 * r2);
				energy += 4 * (r6Inverse * r6Inverse - r6Inverse) - energyCut;
				force = 48 * (r6Inverse * r6Inverse - 0.5 * r6Inverse);
				fx[i] += dx * force / r2;
				fx[j] -= dx * force / r2;
				fy[i] += dy * force / r2;
				fy[j] -= dy * force / r2;
				fz[i] += dz * force / r2;
				fz[j] -= dz * force / r2;
			}
		}
	}
	PE = energy;
}

void Solver::xyzOut(std::fstream &file, bool includeVelocities, bool unfold)
{
	if (!file.is_open())
	{
		std::cout << "Can't open file to write" << std::endl;
		exit(1);
	}

	file << N << std::endl << std::endl;

	for (int i = 0; i < N; ++i)
	{
		file << i << "\t" <<
			 x[i] + (unfold ? (xBCCounter[i] * L) : 0.0) << "\t" <<
			 y[i] + (unfold ? (yBCCounter[i] * L) : 0.0) << "\t" <<
			 z[i] + (unfold ? (zBCCounter[i] * L) : 0.0) << "\t";
		if (includeVelocities)
		{
			file <<
				 vx[i] << "\t" <<
				 vy[i] << "\t" <<
				 vz[i] << "\t" << std::endl;
		} else
			file << std::endl;
	}
}

void Solver::xyzIn(std::fstream &file)
{
	if (!file.is_open())
	{
		std::cout << "Can't open file to read" << std::endl;
		exit(1);
	}

	std::string line;
	file >> N >> initFileVelocitiesIncluded;
	std::getline(file, line);

	for (int i = 0; i < N; ++i)
	{
		std::getline(file, line);
		std::istringstream ss(line);
		ss >> x[i] >> y[i] >> z[i];
		if (ss >> vx[i])
		{
			if (ss >> vy[i])
			{
				if (ss >> vz[i])
				{

				}
			}
		}
	}
}

void Solver::initializeArrays()
{
	x = new double[N];
	y = new double[N];
	z = new double[N];
	vx = new double[N];
	vy = new double[N];
	vz = new double[N];
	fx = new double[N];
	fy = new double[N];
	fz = new double[N];

	xBCCounter = new int[N];
	yBCCounter = new int[N];
	zBCCounter = new int[N];
}

Solver::Solver(const int &N, const double &L, const double &dt, const double &T0)
{
	this->N = N;
	this->L = L;
	this->dt = dt;
	this->dt2 = dt * dt;
	this->T0 = T0;

	initializeArrays();

	// Generate uniform grid
	int n = std::ceil(std::cbrt(N));

	/* Assign particle positions */
	int xIndex = 0, yIndex = 0, zIndex = 0;
	for (int i = 0; i < N; ++i)
	{
		x[i] = (xIndex + 0.5) * L / n;
		y[i] = (yIndex + 0.5) * L / n;
		z[i] = (zIndex + 0.5) * L / n;
		vx[i] = gsl_ran_exponential(r, 1.0);
		vy[i] = gsl_ran_exponential(r, 1.0);
		vz[i] = gsl_ran_exponential(r, 1.0);
		++xIndex;
		if (xIndex == n)
		{
			xIndex = 0;
			++yIndex;
			if (yIndex == n)
			{
				yIndex = 0;
				++zIndex;
			}
		}
	}

	double cmvx = 0, cmvy = 0, cmvz = 0;
	/* Take away any center-of-mass drift; compute initial KE */
	for (int i = 0; i < N; ++i)
	{
		cmvx += vx[i];
		cmvy += vy[i];
		cmvz += vz[i];
	}
	KE = 0;
	for (int i = 0; i < N; ++i)
	{
		vx[i] -= cmvx / N;
		vy[i] -= cmvy / N;
		vz[i] -= cmvz / N;
		KE += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
	}
	KE *= 0.5;
	/* if T0 is specified, scale velocities */
	double T, fac;
	if (T0 > 0)
	{
		T = KE / N * 2. / 3.;
		fac = sqrt(T0 / T);
		KE = 0;
		for (int i = 0; i < N; ++i)
		{
			vx[i] *= fac;
			vy[i] *= fac;
			vz[i] *= fac;
			KE += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
		}
		KE *= 0.5;
	}
	updateForces();
}

void Solver::saveEnergyToFile(std::fstream &file)
{
	if (file.is_open())
	{
		file << t << "\t" << KE << "\t" << PE << "\t" << TotalE << std::endl;
	} else
	{
		std::cout << "File to write energy not found" << std::endl;
	}
}
