// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONT class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

#ifndef CITO_CONT_H
#define CITO_CONT_H

class CitoCont
{
private:
  // control variables
  double kcon[npair];
  Eigen::Matrix<double, 6*nfree, 1> hcon;
  // contact model variables
  double phi_e, phi_n, zeta, phi_c, fn;
  Eigen::Matrix<double, 3, 1> p_sr, p_se, p_fb, n_cs, v_re, v_ef, lambda;
  Eigen::Matrix<double, 6*nfree, 1> h;

public:
  void                              setControl(mjData* d, const ctrlVec_t u);
  Eigen::Matrix<double, 6*nfree, 1> contactModel(const mjData* d, double[npair] kcon);
};

#endif
