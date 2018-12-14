// getCost: returns the nonlinear cost given state matrix X
double SCvx::getCost(const stateVecThread X, const ctrlVecThread U, Eigen::VectorXd po_d, double *w)
{
    Eigen::VectorXd po_f(6);
    // terminal cost
    for( int i=0; i<6; i++ )
    {
        po_f[i] = X[NTS](i);
    }
    double Jt = 0.5*(w[0]*(po_d.block<2,1>(0,0)-po_f.block<2,1>(0,0)).squaredNorm()+
                     w[1]*(po_d.block<4,1>(2,0)-po_f.block<4,1>(2,0)).squaredNorm());

    // integrated cost
    // kcon terms
    Eigen::Matrix<double, NTS, 1> k; k.setZero();
    for( int i=0; i<NTS; i++ )
    {
        k[i] = U[i](ndof);
    }
    double Ji = 0.5*w[2]*k.squaredNorm();

    // total cost
    double J = Jt + Ji;

    return J;
}

