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

class CitoNumDiff
{
private:
    int nnH     = 6 + NTS;
    int lencObj = 6 + NTS;
    int lenru = 4;
    double *cObj = new double[lencObj]; double *ru   = new double[lenru];
    int neA = NTS*N*N + (NTS+1)*N + NTS*N*M + ((NTS+1)*N+NTS*M)*5;
    int n   = ((NTS+1)*N + NTS*M)*2;              // x2 is for auxiliary variables
    int nc  = (NTS+1)*N + ((NTS+1)*N+NTS*M)*2 + 1;
    int iObj = -1;
    int nS, nInf;
    double ObjAdd  = 0, sInf = 0, objective;
    double infBnd = 1.0e20;
    int Cold = 0, Basis = 1, Warm = 2;

public:
    CitoSqopt();
    ~CitoSqpopt();
};

#endif //CITO_SQOPT_H
