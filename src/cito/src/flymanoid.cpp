// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include <iostream>
#include "mujoco.h"
#include "cito_cont.h"

// ********* mujoco model & data **********************************************/
mjModel*  m = NULL;
mjData*   d = NULL;
//==============================================================================
int main(int argc, char const *argv[]) {
  // ********* mujoco initialization ******************************************/
  // activate mujoco
  mj_activate("/home/aykut/Development/cito/src/cito/bin/mjkey.txt");
  // load xml model
  m = mj_loadXML("/home/aykut/Development/cito/src/cito/model/flymanoid.xml", 0, 0, 0);
  if( !m )
      mju_error("Cannot load the model");
  // create data
  d = mj_makeData(m);


  // ********* mujoco shut down ***********************************************/
  mj_deleteData(d);
  mj_deleteModel(m);
  mj_deactivate();

  return 0;
}
