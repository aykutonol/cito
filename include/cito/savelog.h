/*! Save Log */
/**
 *  \brief SaveLog saves MuJoCo log
 *
 *  This class consists of methods for saving data from a MuJoCo simulation.
 *
 *  \author Aykut Onol
 */

#ifndef SAVELOG_H
#define SAVELOG_H

#include "cito/params.h"
#include <fstream>

class SaveLog
{
public:
    /// Constructor
    SaveLog(const mjModel *m_, Params *cp_);
    /// Destructor
    ~SaveLog();
    /// writes simulation data to the print file
    void writeData(const mjData *d, int save);

private:
    /// MuJoCo model file
    const mjModel *m;
    /// Log file for playback
    FILE *logFile;
    /// Data header
    int header[6];
    /// Data size
    int dataSize;
    /// Trajectory file for execution
    std::ofstream trajFile;
    /// Objects
    Params *cp;
};

#endif //MJ_SAVELOG_H
