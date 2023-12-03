#ifndef DUMBVISCO_H
#define DUMBVISCO_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class dumbVisco
{
public:
	dumbVisco(elasticRod &m_rod, timeStepper &m_stepper, int idx);
	~dumbVisco();
	void computeFv();
	void computeJv();

	VectorXd ForceVec;
	bool isReleasing;

private:
	elasticRod *rod;
	timeStepper *stepper;

	int ind1, ind2, mappedInd1, mappedInd2;
	double f, jac;
	int rod_idx;
	double k_dumb_visco;
	double ind;
};

#endif
