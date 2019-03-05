/*! SQOPT */
/**
 *  \brief CitoSQOPT class consists of functions for using SQOPT
 *
 *  This class defines functions that are used to convert the convex subproblems
 *  used in SCvx to the appropriate form for SQOPT.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_SQOPT_H
#define CITO_SQOPT_H

#include "cito_params.h"
#include "snoptProblem.hpp"

class CitoSQOPT
{
public:
    /// Constructor
    CitoSQOPT();
    /// Destructor
    ~CitoSQOPT() {}
    /// This function solves the convex subproblem
    void solveCvx(double *xTraj, double r, const stateVecThread X, const ctrlVecThread U,
                  const stateDerThread Fx, const ctrlDerThread Fu, int *isJFree, int *isAFree,
                  double *qposLB, double *qposUB, double *tauLB, double *tauUB);
private:
    /// This function sets linear and constant cost terms
    void setCost(const stateVecThread X, const ctrlVecThread U,
                 double *ru, double *cObj, double& ObjAdd);
    /// This function sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
    void setBounds(double r, const stateVecThread X, const ctrlVecThread U,
                   double *bl, double *bu, int *isJFree, int *isAFree,
                   double *qposLB, double *qposUB, double *tauLB, double *tauUB);
    /* This function creates the sparse A matrix for linearized dynamics, auxiliary
     *  variables, and trust region constraints */
    void setA(double *valA, int *indA, int *locA,
              const stateDerThread Fx, const ctrlDerThread Fu);
    /// This function modifies A and the bounds such that non-zero elements in H come first
    void sortToMatch(double *valA, int *indA, int *locA, int *indMove, double *bl, double *bu);
    /// This function moves column iMove in A to left
    void moveColA(double *valA, int *indA, int *locA, int iMove);
    /// This function moves row iMove in bounds to top
    void moveRowBounds(double *bl, double *bu, int iMove);
    /// This function sorts decision variables back to original order
    void sortX(double *x, int *indMove);
    /// SQOPT parameters
    int nnH     = 6 + 6 + NTS*NPAIR;    // number of non-zero elements of the Hessian
    int lencObj = nnH;                  // number of non-zero elements of the linear term
    int lenru   = 4;                    // number of weights
    int neA = NTS*N*N + (NTS+1)*N + NTS*N*M + ((NTS+1)*N+NTS*M)*5;
    int n   = ((NTS+1)*N + NTS*M)*2;    // *2 is for auxiliary variables for l1-norm
    int nc  = (NTS+1)*N + ((NTS+1)*N+NTS*M)*2 + 1;
    int iObj = -1;
    int nS, nInf;
    double *cObj   = new double[lencObj];
    double *ru     = new double[lenru];
    double ObjAdd = 0, sInf = 0, objective;
    double infBnd = 1.0e20;
    int Cold = 0, Basis = 1, Warm = 2;
    /// Cost parameters
    int controlJointDOF0;
    Eigen::Matrix<double, 6, 1> desiredPose, desiredVelo, dPose, dVelo;
    kConVecThread dKCon;
    double dKConSN;
    /// Sort parameters
    int nMove    = nnH;
    int *indMove = new int[nMove];
    double *xTemp;
    /// setBounds parameters
    int dUOffset  = (NTS+1)*N;
    int auxOffset = (NTS+1)*N+NTS*M;
    /// SQOPT problem object
    sqoptProblem cvxProb;
};

#endif //CITO_SQOPT_H
