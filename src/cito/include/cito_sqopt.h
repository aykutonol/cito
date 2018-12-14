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

class CitoSQOPT
{
private:
    // solver parameters
    int nnH     = 6 + NTS;      // number of non-zero elements of the Hessian
    int lencObj = 6 + NTS;      // number of non-zero elements of the linear term
    int lenru   = 3;            // number of weights
    double *cObj = new double[lencObj]; double *ru   = new double[lenru];
    int neA = NTS*N*N + (NTS+1)*N + NTS*N*M + ((NTS+1)*N+NTS*M)*5;
    int n   = ((NTS+1)*N + NTS*M)*2;    // *2 is for auxiliary variables for l1-norm
    int nc  = (NTS+1)*N + ((NTS+1)*N+NTS*M)*2 + 1;
    int iObj = -1;
    int nS, nInf;
    double ObjAdd  = 0, sInf = 0, objective;
    double infBnd = 1.0e20;
    int Cold = 0, Basis = 1, Warm = 2;
    // setBounds parameters
    int dU_offset = (NTS+1)*N;
    int aux_offset  = (NTS+1)*N+NTS*M;
    // sortX parameters
    double xtemp[n];

public:
    CitoSQOPT()  {}
    ~CitoSQOPT() {}
};

#endif //CITO_SQOPT_H
