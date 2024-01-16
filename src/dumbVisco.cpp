#include "dumbVisco.h"

dumbVisco::dumbVisco(elasticRod &m_rod, timeStepper &m_stepper, int idx, double k_dumb_visco)
{
	rod = &m_rod;
	stepper = &m_stepper;

	ForceVec = VectorXd::Zero(rod->ndof);
	rod_idx = idx;
	isReleasing = true;
	k_visco = k_dumb_visco; // defined in world.cpp
}

dumbVisco::~dumbVisco()
{
	;
}

void dumbVisco::computeFv()
{
	int ind;

	ForceVec = VectorXd::Zero(rod->ndof);

	for (int i = 0; i < rod->ndof; i++)
	{

		f = (k_visco * rod->refLen[0] * (rod->x[i] - rod->x0[i]) / rod->dt) * isReleasing;

		ForceVec(i) = f;

		stepper->addForce(i, f, rod_idx);

	}
}

void dumbVisco::computeJv()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		jac = k_visco * rod->refLen[0] / ((rod->dt)) * isReleasing;
		stepper->addJacobian(i, i, jac, rod_idx);
	}
}
