#include "dumbVisco.h"

dumbVisco::dumbVisco(elasticRod &m_rod, timeStepper &m_stepper, int idx)
{
	rod = &m_rod;
	stepper = &m_stepper;

	ForceVec = VectorXd::Zero(rod->ndof);
	rod_idx = idx;
	isReleasing = false;

	double k_dumb_visco;
	k_dumb_visco = 5;
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

		f = (k_dumb_visco * rod->refLen[0] * rod->u[i]) * isReleasing;

		ForceVec(i) = f;

		stepper->addForce(i, f, rod_idx);
	}
}

void dumbVisco::computeJv()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		jac = k_dumb_visco * rod->refLen[0] / ((rod->dt));
		stepper->addJacobian(i, i, jac, rod_idx);
	}
}
