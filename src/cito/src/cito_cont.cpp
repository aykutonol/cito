// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONT class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

// ***** CLASS TYPE ************************************************************
// Robot and environment specific

#include "cito_cont.h"

// intermediate control variables **********************************************


//setControl: sets generalized forces on joints and bodies
void SCvx::setControl(const mjModel* m, mjData* d, const ctrlVec_t u)
{
    // set control
    for( int i=0; i<ndof; i++ )
    {
      d->ctrl[i] = u[i] + d->qfrc_bias[u_id[i]];
    }
    // get kcon
    double kcon = u[ndof];
    // contact model
    ho = this->contactModel(m, d, kcon);
    // set external force on the torso
    for( int i=0; i<6; i++ )
    {
      d->xfrc_applied[mj_name2id(m, mjOBJ_BODY, "torso")*6+i] = ho[i];
    }
}
