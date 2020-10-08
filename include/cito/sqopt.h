/*! SQOPT Interface */
/**
 *  \brief SQOPT consists of methods for using SQOPT in CITO.
 *
 *  This class defines functions that are used to convert the convex subproblems 
// used in the SCVX method to the appropriate form for SQOPT.
 *
 *  \author Aykut Onol
 */

#ifndef SQOPT_H
#define SQOPT_H

#include "cito/params.h"
#include "snoptProblem.hpp"

class SQOPT
{
public:
    /// Constructor
    SQOPT(const mjModel *m_, Params *cp_);
    /// Destructor
    ~SQOPT();
    /// This function solves the convex subproblem
    void solveCvx(double *xTraj, double r, const eigMd &X, const eigMd &U,
                  const eigTd &Fx, const eigTd &Fu, int *isJFree, int *isAFree,
                  double *qposLB, double *qposUB, double *tauLB, double *tauUB);

private:
    /// This callback function sets Hx to the H*x part of the quadratic cost to be multiplied by x'
    static void qpHx(int *nnH, double x[], double Hx[], int *nState,
                     char cu[], int *lencu, int iu[], int *leniu,
                     double ru[], int *lenru);
    /// This function sets linear and constant cost terms of the cost
    void setCObj(const eigMd &X, const eigMd &U,
                 double *ru, double *cObj, double &ObjAdd);
    /// This function sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
    void setBounds(double r, const eigMd &X, const eigMd &U,
                   double *bl, double *bu, int *isJFree, int *isAFree,
                   double *qposLB, double *qposUB, double *tauLB, double *tauUB);
    /* This function creates the sparse A matrix for linearized dynamics, auxiliary
     *  variables, and trust region constraints */
    void setA(double *valA, int *indA, int *locA, const eigTd &Fx, const eigTd &Fu);
    /// This function modifies A and the bounds such that non-zero elements in H come first
    void sortToMatch(double *valA, int *indA, int *locA, int *indMove, double *bl, double *bu);
    /// This function moves column iMove in A to left
    void moveColA(double *valA, int *indA, int *locA, int iMove);
    /// This function moves row iMove in bounds to top
    void moveRowBounds(double *bl, double *bu, int iMove);
    /// This function sorts decision variables back to original order
    void sortX(double *x, int *indMove);
    /// MuJoCo model
    const mjModel *m;
    /// SQOPT parameters
    int nnH, neA, n, nc, lencObj, lenru, leniu, *iu;
    double *cObj, *ru;
    int iObj = -1;
    int nS, nInf;
    double ObjAdd = 0, sInf = 0, objective;
    double infBnd = 1.0e20;
    int Cold = 0, Basis = 1, Warm = 2;
    /// Cost parameters
    eigVd deltaPos;
    /// Sort parameters
    int nMove, *indMove;
    double *xTemp;
    /// setBounds parameters
    int dUOffset, auxOffset;
    double kCon0;
    /// Objects
    Params *cp;
    sqoptProblem cvxProb;
};

#endif //SQOPT_H
