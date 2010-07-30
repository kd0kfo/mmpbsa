/* 
 * File:   SanderInterface.h
 * Author: dcoss
 *
 * Created on July 29, 2010, 9:30 AM
 */

#ifndef SANDERINTERFACE_H
#define	SANDERINTERFACE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mmpbsa_utils.h"

#ifdef __USE_BOINC__
#include "boinc_api.h"
#include "filesys.h"
#include "error_numbers.h"
#include "util.h"
#define HAVE_STRCASESTR 1
#include "str_replace.h"
#include "str_util.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <signal.h>
#else
#define ERR_FORK -147
#define PROCESS_IDLE_PRIORITY 19
#define ERR_EXEC -148
#endif

namespace mmpbsa{

class SanderInterface {
public:
    SanderInterface();
    SanderInterface(const SanderInterface& orig);
    virtual ~SanderInterface();

    int start(const double& start_time = 0.0);
    bool poll(int& status);
    void kill();
    void stop();
    void resume();
    void poll_boinc_messages();
    

    const double& start_time()const{return starting_cpu;}
    const bool& isSuspended()const{return suspended;}
    const int& getPID()const{return pid;}
    const double netRuntime()const{return starting_cpu + final_cpu_time;}
    double cpu_time()const;

    void setPID(const int& newPID){pid = newPID;}
    void set_start_time(const double& new_start_time){starting_cpu = new_start_time;}

    std::string mdinFilename;
    std::string mdoutFilename;
    std::string prmtopFilename;
    std::string inpcrdFilename;
    std::string restartFilename;
    static const double pollPeriod = 1.0;
    bool completed;
    
private:
    int pid;
    double wall_cpu_time;
    double starting_cpu;
    double final_cpu_time;
    bool suspended;
};


}//end namespace mmpbsa

#endif	/* SANDERINTERFACE_H */

