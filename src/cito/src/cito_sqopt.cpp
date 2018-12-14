// setBounds: sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
void CitoSqopt::setBounds(double *bl, double *bu, double *qpos_lb, double *qpos_ub,
                          double *tau_lb, double *tau_ub, int *isJFree, int *isAFree,
                          int n, int nc, const stateVecThread X, const ctrlVecThread U, double r)
{
    // decision variables
    // states
    int dU0i = (NTS+1)*N;
    int y0i  = (NTS+1)*N+NTS*M;
    for( int i=0; i<NTS+1; i++ )
    {
        for( int j=0; j<NV; j++ )
        {
            // change in free joint positions: unbounded
            if( isJFree[j] == 0 )
            {
                // change in joint positions: bounded by joint limits - current value
                bl[i*N+j] = qpos_lb[j] - X[i](j);
                bu[i*N+j] = qpos_ub[j] - X[i](j);
            }
            // change in joint velocities: unbounded (already set)
        }
        if( i < NTS )
        {
            // change in joint torques: bounded by torque limits - current value
            for( int j=0; j<ndof; j++ )
            {
                bl[dU0i+i*M+j] = tau_lb[j] - U[i](j);
                bu[dU0i+i*M+j] = tau_ub[j] - U[i](j);
            }
            bl[dU0i+i*M+ndof] = 0 - U[i](ndof);
            bu[dU0i+i*M+ndof] = kcon0 - U[i](ndof);
        }
    }
    // auxiliary variables > 0
    for( int i=y0i; i<n; i++ )
    {
        bl[i] = 0;
    }

    // constraints
    // dynamics == 0
    for( int i=n; i<n+(NTS+1)*N; i++ )
    {
        bl[i] = 0;
        bu[i] = 0;
    }
    // 0 <= absolute value constraints <= inf
    for( int i=n+(NTS+1)*N; i<n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2; i++ )
    {
        bl[i] = 0;
    }
    // 0 <= trust region <= r
    bl[n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2] = 0;
    bu[n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2] = r;
}

// setA: creates the sparse A matrix for linearized dynamics, auxiliary
// variables, and trust region constraints
void CitoSqopt::setA(double *valA, int *indA, int *locA,
                     const stateMatThread Fx, const ctrlMatThread Fu)
{
    int colNo = 0, indNo = 0, indTS = 0;//, auxNo1 = 0, auxNo2 = 0;

    locA[0] = 0;
    // columns associated with dx[1,...,NTS]
    for( int i=0; i<NTS*N; i++ )
    {
        indA[indNo] = i;
        valA[indNo] = -1;
        indNo += 1;
        // time step index
        indTS = floor(i/N);
        for( int j=0; j<N; j++ )
        {
            indA[indNo] = (indTS+1)*N+j;
            valA[indNo] = Fx[indTS](j,i%N);
            indNo += 1;
        }

        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo += 1; //auxNo1 += 1;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo += 1; //auxNo2 += 1;

        colNo += 1;
        locA[colNo] = indNo;
    }

    // columns associated with dx[NTS+1]
    for( int i=0; i<N; i++ )
    {
        indA[indNo] = NTS*N+i;
        valA[indNo] = -1;
        indNo += 1;

        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo += 1; //auxNo1 += 1;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo += 1; //auxNo2 += 1;

        colNo += 1;
        locA[colNo] = indNo;
    }

    // columns associated with du[1,...,NTS]
    for( int i=0; i<NTS*M; i++ )
    {
        // time step index
        indTS = floor(i/M);
        for( int j=0; j<N; j++ )
        {
            indA[indNo] = (indTS+1)*N+j;
            valA[indNo] = Fu[indTS](j,i%M);
            indNo += 1;
        }

        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo += 1; //auxNo1 += 1;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo += 1; //auxNo2 += 1;

        colNo += 1;
        locA[colNo] = indNo;
    }

    // columns associated with auxiliary variables
    int auxNo1 = 0, auxNo2 = 0;
    for( int i=0; i<(NTS+1)*N+NTS*M; i++ )
    {
        indA[indNo] = (NTS+1)*N + auxNo1;
        valA[indNo] = +1.0;
        indNo += 1; auxNo1 += 1;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + auxNo2;
        valA[indNo] = +1.0;
        indNo += 1; auxNo2 += 1;
        // for l1-norm
        indA[indNo] = (NTS+1)*N + 2*((NTS+1)*N + NTS*M);
        valA[indNo] = +1.0;
        indNo += 1;

        colNo += 1;
        locA[colNo] = indNo;
    }
}

// moveColA: moves iMove to beginning
void CitoSqopt::moveColA(double *valA, int *indA, int *locA, int iMove, int neA, int n)
{
    int indAtemp[neA], neCol[n]; double valAtemp[neA];

    for( int i=0; i<n; i++ )
    {
        neCol[i] = locA[i+1]-locA[i];
    }
    for( int i=0; i<neCol[iMove]; i++)
    {
        indAtemp[i] = indA[locA[iMove]+i];
        valAtemp[i] = valA[locA[iMove]+i];
    }
    for( int i=0; i<locA[iMove]; i++ )
    {
        indAtemp[neCol[iMove]+i] = indA[i];
        valAtemp[neCol[iMove]+i] = valA[i];
    }

    locA[0] = 0;  locA[1] = neCol[iMove];
    for( int i=1; i<iMove; i++ )
    {
        locA[i+1] = locA[i] + neCol[i-1];
    }
    for( int i=0; i<locA[iMove+1]; i++ )
    {
        indA[i] = indAtemp[i]; valA[i] = valAtemp[i];
    }
}

// moveRowBounds: moves iMove to top
void CitoSqopt::moveRowBounds(double *bl, double *bu, int iMove)
{
    double bl_temp[iMove], bu_temp[iMove];
    bl_temp[0] = bl[iMove];
    bu_temp[0] = bu[iMove];
    for( int i=0; i<iMove+1; i++ )
    {
        if( i<iMove )
        {
            bl_temp[i+1] = bl[i];
            bu_temp[i+1] = bu[i];
        }
        bl[i] = bl_temp[i];
        bu[i] = bu_temp[i];
    }
}

// sortX: sorts decision variables back to original
void CitoSqopt::sortX(double *x, int *rowsMove, int nMove, int n)
{
    double xtemp[n];
    int k = 0;
    for( int i=0; i<n; i++ )
    {
        xtemp[i] = x[nMove+k];
        for( int j=0; j<nMove; j++ )
        {
            if( i == rowsMove[j] )
            {
                xtemp[rowsMove[j]] = x[nMove-j-1];
                k = k - 1;
            }
        }
        k = k + 1;
    }
    for( int i=0; i<n; i++ )
    {
        x[i] = xtemp[i];
    }
}

// Test
// // print A
// int neCol[n];
// for( int i=0; i<n; i++ )
// {
//   neCol[i] = locA[i+1]-locA[i];
// }
// int indNo = 0, colNo = 0;
// for( int i=0; i<n; i++ )
// {
//   std::cout << "Column " << i << ":\n";
//   for( int j=0; j<neCol[i]; j++ )
//   {
//     std::cout << "\tindNo: " << indNo << ", indA: " << indA[indNo] << ", valA: " << valA[indNo] << '\n';
//     indNo+=1;
//   }
//   colNo += 1;
//   std::cout << "\tlocA: " << locA[colNo] << '\n';
// }
// std::cout << "n: " << n << ", neA: " << neA << '\n';

// // test shifts
// int ndene = 10;
// double *dene   = new double[ndene];
// double *deneDummy  = new double[ndene];
// double assignDene[10] = {1, 2, 3, 4, 5, 6 ,7, 8, 9, 10};
// for( int i=0; i<ndene; i++ )
// {
//   dene[i] = assignDene[i];
//   deneDummy[i] = assignDene[i];
// }
// for( int i=0; i<ndene; i++ )
//   std::cout << dene[i] << ' ';
// std::cout << "\n\n\n";
// const int nshift = 5;
// int *shift = new int[nshift];
// int assignShift[nshift] = {3, 9, 4, 7, 8};
// for( int i=0; i<nshift; i++ )
//   shift[i] = assignShift[i]-1;
// int less_counter = 0;
// for( int i=0; i<nshift; i++ )
// {
//   if( i>0 && shift[i]<shift[i-1] )
//   {
//     less_counter += 1;
//   }
//   int rMove = shift[i]+less_counter;
//   sc.moveRowBounds(dene, deneDummy, rMove);
//   for( int i=0; i<ndene; i++ )
//     std::cout << dene[i] << ' ';
//   std::cout << "\n\n\n";
// }
//
// sc.sortX(dene, shift, nshift, ndene);
// for( int i=0; i<ndene; i++ )
//   std::cout << dene[i] << ' ';
// std::cout << "\n\n\n";