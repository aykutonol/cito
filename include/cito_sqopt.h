// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SQOPT class consists of functions that setup the CITO problem for the
// convex programming solver SQOPT.

// ***** CLASS TYPE ************************************************************
// Solver specific

#ifndef CITO_SQOPT_H
#define CITO_SQOPT_H

#include "cito_params.h"
#include "snoptProblem.hpp"

class CitoSQOPT
{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
    CitoSQOPT();
    ~CitoSQOPT() {}
    // ***** FUNCTIONS *************************************************************
    void solveCvx(double *xTraj, double r, const stateVecThread X, const ctrlVecThread U,
                  const stateDerThread Fx, const ctrlDerThread Fu, int *isJFree, int *isAFree,
                  double *qposLB, double *qposUB, double *tauLB, double *tauUB);
private:
    // ***** FUNCTIONS *************************************************************
    void setCost(const stateVecThread X, const ctrlVecThread U,
                 double *ru, double *cObj, double& ObjAdd);
    void setBounds(double r, const stateVecThread X, const ctrlVecThread U,
                   double *bl, double *bu, int *isJFree, int *isAFree,
                   double *qposLB, double *qposUB, double *tauLB, double *tauUB);
    void setA(double *valA, int *indA, int *locA,
              const stateDerThread Fx, const ctrlDerThread Fu);
    void sortToMatch(double *valA, int *indA, int *locA, int *indMove, double *bl, double *bu);
    void moveColA(double *valA, int *indA, int *locA, int iMove);
    void moveRowBounds(double *bl, double *bu, int iMove);
    void sortX(double *x, int *indMove);
    // ***** PARAMETERS ************************************************************
    // solver parameters
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
    // setCost parameters
    int controlJointDOF0;
    Eigen::Matrix<double, 6, 1> dPose, dVelo, desiredPose, desiredVelo;
    kConVecThread dKCon;
    double dKConSN;
    // sort parameters
    int    nMove    = nnH;
    int    *indMove = new int[nMove];
    double *xTemp   = new double[n];
    // setBounds parameters
    int dUOffset  = (NTS+1)*N;
    int auxOffset = (NTS+1)*N+NTS*M;
    // ***** OBJECTS ***************************************************************
    sqoptProblem cvxProb;
};

#endif //CITO_SQOPT_H
