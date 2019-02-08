// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// MJ_SAVE class consists of functions that save MuJoCo logs.

// ***** CLASS TYPE ************************************************************
// Physics engine specific

#include "mj_savelog.h"
#include <iostream>
// ***** CONSTRUCTOR ***********************************************************
MjSaveLog::MjSaveLog(const mjModel* model) : m(model)
{
    // open the log file
    printFile = fopen(paths::logFile, "wb");
    if( printFile == NULL ) { mju_error("Unable to open the log file."); }
    // create and write the header
    header[0] = m->nq;
    header[1] = m->nv;
    header[2] = m->nu;
    header[3] = m->nmocap;
    header[4] = m->nsensordata;
    header[5] = strlen(m->names);
    fwrite(header, sizeof(int), 6, printFile);
    this->dataSize = 1 + m->nq + m->nv + m->nu + 7*m->nmocap + m->nsensordata;
}
MjSaveLog::~MjSaveLog()
{
    // close the log file
    fclose(printFile);
}

// ***** FUNCTIONS *************************************************************
// writeData: writes simulation data to the print file
void MjSaveLog::writeData(const mjData* d)
{
    // create data for recording
    float data[dataSize];
    // copy data to the record
    data[0] = (float) d->time;
    mju_n2f(data+1, d->qpos, m->nq);
    mju_n2f(data+1+m->nq, d->qvel, m->nv);
    mju_n2f(data+1+m->nq+m->nv, d->ctrl, m->nu);
    mju_n2f(data+1+m->nq+m->nv+m->nu, d->mocap_pos, 3*m->nmocap);
    mju_n2f(data+1+m->nq+m->nv+m->nu+3*m->nmocap, d->mocap_quat, 4*m->nmocap);
    mju_n2f(data+1+m->nq+m->nv+m->nu+7*m->nmocap, d->sensordata, m->nsensordata);
	// read data, print size
	size_t nn = fwrite(data, this->dataSize*sizeof(float), 1, printFile);
}