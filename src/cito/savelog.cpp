// ***** DESCRIPTION ***********************************************************
// SaveLog consists of methods for saving data from a MuJoCo simulation.

#include "cito/savelog.h"

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
SaveLog::SaveLog(const mjModel *m_, Params *cp_) : m(m_), cp(cp_)
{
    // open the log file
    YAML::Node params = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/params.yaml");
    std::string modelName = params["model"].as<std::string>();
    modelName.erase(modelName.end() - 4, modelName.end());
    std::string logPathStr = paths::workspaceDir + "/logs/mjLog_" + modelName + ".log";
    const char *logPath = logPathStr.c_str();
    logFile = fopen(logPath, "wb");
    if (logFile == NULL)
    {
        mju_error("Unable to open the log file.");
    }
    // create and write the header
    header[0] = m->nq;
    header[1] = m->nv;
    header[2] = m->nu;
    header[3] = m->nmocap;
    header[4] = m->nsensordata;
    header[5] = strlen(m->names);
    fwrite(header, sizeof(int), 6, logFile);
    this->dataSize = 1 + m->nq + m->nv + m->nu + 7 * m->nmocap + m->nsensordata;
    // open the trajectory file
    std::string trajPathStr = paths::workspaceDir + "/logs/traj_" + modelName + ".txt";
    const char *trajPath = trajPathStr.c_str();
    trajFile.open(trajPath);
    /// write header to outFile
    if (trajFile.is_open())
    {
        trajFile << "time,positions,velocities,accelerations,controls\n";
    }
    else
    {
        std::cout << "WARNING: Unable to open the trajectory file.\n";
    }
}
SaveLog::~SaveLog()
{
    // close the log file
    fclose(logFile);
    // close the trajectory file
    trajFile.close();
}

// ***** FUNCTIONS *************************************************************
// writeData: writes simulation data to logFile
void SaveLog::writeData(const mjData *d, int save)
{
    /// logFile
    // create data for the record
    float *data = new float[dataSize];
    // write to log file
    data[0] = (float)d->time;
    mju_n2f(data + 1, d->qpos, m->nq);
    mju_n2f(data + 1 + m->nq, d->qvel, m->nv);
    mju_n2f(data + 1 + m->nq + m->nv, d->ctrl, m->nu);
    mju_n2f(data + 1 + m->nq + m->nv + m->nu, d->mocap_pos, 3 * m->nmocap);
    mju_n2f(data + 1 + m->nq + m->nv + m->nu + 3 * m->nmocap, d->mocap_quat, 4 * m->nmocap);
    mju_n2f(data + 1 + m->nq + m->nv + m->nu + 7 * m->nmocap, d->sensordata, m->nsensordata);
    // read data, print size
    size_t nn = fwrite(data, this->dataSize * sizeof(float), 1, logFile);
    /// trajFile
    // write to trajectory file
    if (trajFile.is_open() && save > 1)
    {
        trajFile << d->time << ",";
        for (int i = 0; i < m->nu; i++)
        {
            trajFile << d->qpos[cp->pAct[i]] << ",";
        }
        for (int i = 0; i < m->nu; i++)
        {
            trajFile << d->qvel[cp->dAct[i]] << ",";
        }
        for (int i = 0; i < m->nu; i++)
        {
            trajFile << d->qacc[cp->dAct[i]] << ",";
        }
        for (int i = 0; i < m->nu; i++)
        {
            trajFile << d->ctrl[i] << ",";
        }
        trajFile << "\n";
    }
    delete[] data;
}
