/*! Save Log */
/**
 *  \brief MjSaveLog saves MuJoCo log
 *
 *  This class consists of functions for saving data from a MuJoCo simulation.
 *
 *  \author Aykut Onol
 */

#ifndef MJ_SAVELOG_H
#define MJ_SAVELOG_H

#include "cito_params.h"
#include <fstream>

class MjSaveLog{
public:
    /// Constructor
    MjSaveLog(const mjModel* model);
    /// Destructor
    ~MjSaveLog();
    /// writes simulation data to the print file
    void writeData(const mjData *d);
private:
    /// MuJoCo model file
    const mjModel* m;
    /// Log file for playback
    FILE* logFile;
    /// Data header
    int header[6];
    /// Data size
    int dataSize;
    /// Trajectory file for execution
    std::ofstream trajFile;
    /// Objects
    CitoParams cp;
};

#endif //MJ_SAVELOG_H
