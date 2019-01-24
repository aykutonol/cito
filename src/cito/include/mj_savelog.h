// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// MJ_SAVE class consists of functions that save MuJoCo logs.

// ***** CLASS TYPE ************************************************************
// Physics engine specific

#ifndef MJ_SAVELOG_H
#define MJ_SAVELOG_H

#include <stdio.h>
#include <string.h>
#include "mujoco.h"
#include "cito_params.h"

class MjSaveLog{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
//    MjSaveLog();
    MjSaveLog(const mjModel* model);
    ~MjSaveLog();
    // ***** FUNCTIONS *************************************************************
    void writeData(const mjData *d);
    // ***** PARAMETERS ************************************************************
    const mjModel* m;
    FILE* printFile;
    int header[6];
    int dataSize;
private:
};

#endif //MJ_SAVELOG_H
