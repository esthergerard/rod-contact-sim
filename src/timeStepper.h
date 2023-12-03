#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticRod.h"

#include "eigenIncludes.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

// extern "C" void dgbsv_( int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info );

/* PARDISO prototype. */
// extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
// extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
//                   double *, int    *,    int *, int *,   int *, int *,
//                      int *, double *, double *, int *, double *);
// extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
// extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
// extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
//                            double *, int *);

/**
 * @class timeStepper
 * @brief Class for managing the time integration of a set of elastic rods.
 *
 * This class provides methods for calculating the force and the jacobian,
 * adding forces and jacobians, and integrating the equations of motion.
 */
class timeStepper
{
public:
	/**
	 * @brief Constructor for the timeStepper class.
	 * @param m_rod_vec A vector of pointers to the elasticRod objects to integrate.
	 */
	timeStepper(std::vector<elasticRod *> m_rod_vec);

	/**
	 * @brief Destructor for the timeStepper class.
	 */
	~timeStepper();

	/**
	 * @brief Retrieves the total force on all rods.
	 * @return A pointer to the force array.
	 */
	double *getForce();

	/**
	 * @brief Retrieves the total jacobian for all rods.
	 * @return A pointer to the jacobian array.
	 */
	double *getJacobian();

	/**
	 * @brief Resets the force and jacobian to zero.
	 */
	void setZero();

	/**
	 * @brief Computes the number of free degrees of freedom.
	 */
	void computeFreeDOF();

	/**
	 * @brief Adds a force to a specific degree of freedom.
	 * @param ind The index of the degree of freedom.
	 * @param p The force to add.
	 * @param idx The index of the rod.
	 */
	void addForce(int ind, double p, int idx);

	/**
	 * @brief Adds a value to the jacobian at a specific location.
	 * @param ind1 The index of the first dimension.
	 * @param ind2 The index of the second dimension.
	 * @param p The value to add.
	 * @param idx The index of the rod.
	 */
	void addJacobian(int ind1, int ind2, double p, int idx);

	/**
	 * @brief Adds a value to the jacobian at a specific location.
	 * @param ind1 The index of the first dimension.
	 * @param ind2 The index of the second dimension.
	 * @param p The value to add.
	 * @param idx1 Column offset
	 * @param idx2 Row offset
	 */
	void addJacobian(int ind1, int ind2, double p, int idx1, int idx2);

	/**
	 * @brief Integrates the equations of motion for a time step.
	 */
	void integrator();
	double *dx;
	VectorXd DX;
	VectorXd Force;
	VectorXd force;
	MatrixXd Jacobian;

private:
	elasticRod *rod;
	elasticRod *rod1;

	std::vector<elasticRod *> rod_vec;
	int kl, ku, freeDOF;

	double *totalForce;
	double *jacobian;

	vector<int> freeDOF_vec;

	// utility variables
	int mappedInd, mappedInd1, mappedInd2;
	int row, col, offset;
	int NUMROWS;
	int jacobianLen;
	int nrhs;
	int *ipiv;
	int info;
	int ldb;

	void pardisoSolver();
};

#endif
