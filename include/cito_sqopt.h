/*! SQOPT Interface */
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
    CitoSQOPT(const mjModel* model);
    /// Destructor
    ~CitoSQOPT() {}
    /// This function solves the convex subproblem
    void solveCvx(double *xTraj, double r, const stateTraj X, const ctrlTraj U,
                  const stateDerTraj Fx, const ctrlDerTraj Fu, int *isJFree, int *isAFree,
                  double *qposLB, double *qposUB, double *tauLB, double *tauUB);
private:
    /// This callback function sets Hx to the H*x part of the quadratic cost to be multiplied by x'
    static void qpHx(int *nnH, double x[], double Hx[], int *nState,
                     char cu[], int *lencu, int iu[], int *leniu,
                     double ru[], int *lenru);
    /// This function sets linear and constant cost terms of the cost
    void setCObj(const stateTraj X, const ctrlTraj U,
                 double *ru, double *cObj, double &ObjAdd);
    /// This function sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
    void setBounds(double r, const stateTraj X, const ctrlTraj U,
                   double *bl, double *bu, int *isJFree, int *isAFree,
                   double *qposLB, double *qposUB, double *tauLB, double *tauUB);
    /* This function creates the sparse A matrix for linearized dynamics, auxiliary
     *  variables, and trust region constraints */
    void setA(double *valA, int *indA, int *locA,
              const stateDerTraj Fx, const ctrlDerTraj Fu);
    /// This function modifies A and the bounds such that non-zero elements in H come first
    void sortToMatch(double *valA, int *indA, int *locA, int *indMove, double *bl, double *bu);
    /// This function moves column iMove in A to left
    void moveColA(double *valA, int *indA, int *locA, int iMove);
    /// This function moves row iMove in bounds to top
    void moveRowBounds(double *bl, double *bu, int iMove);
    /// This function sorts decision variables back to original order
    void sortX(double *x, int *indMove);
    /// MuJoCo model
    const mjModel* m;
    /// SQOPT parameters
    int nnH, neA, n, nc, lencObj, lenru, leniu, *iu;
    double *cObj, *ru;
    int iObj = -1;
    int nS, nInf;
    double ObjAdd = 0, sInf = 0, objective;
    double infBnd = 1.0e20;
    int Cold = 0, Basis = 1, Warm = 2;
    /// Cost parameters
    int controlJointDOF0;
    eigDbl desiredPos, desiredVel, deltaPos, deltaVel, dKCon;
    /// Sort parameters
    int *indMove;
    double *xTemp;
    /// setBounds parameters
    int dUOffset, auxOffset;
    double kCon0;
    /// Objects
    CitoParams cp;
    sqoptProblem cvxProb;
};

#endif //CITO_SQOPT_H
