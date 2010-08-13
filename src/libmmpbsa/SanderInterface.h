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

#if defined(_WIN32) || defined(__MINGW_WIN32__)
#include <windows.h>
#include <winbase.h>
#include <tlhelp32.h>
#endif

#ifdef __USE_BOINC__
#if defined(_WIN32) || defined(__MINGW_WIN32__)
#include "boinc_win.h"
#include "win_util.h"
#endif

#include "boinc_api.h"
#include "filesys.h"
#include "error_numbers.h"
#include "util.h"

#define HAVE_STRCASESTR 1
#include "str_replace.h"
#include "str_util.h"
#else
extern size_t strlcpy(char *dst, const char *src, size_t size);
#define ERR_FORK -147
#define PROCESS_IDLE_PRIORITY 19
#define ERR_EXEC -148
#define ERR_FOPEN -108
#endif

#ifdef __MINGW_WIN32__
#include <windows.h>
#include <winbase.h>
#include <tlhelp32.h>
#endif

#if !defined(_WIN32)
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <signal.h>
#endif



namespace mmpbsa{

class SanderInterface {
public:
    /**
     * SanderInterface is used to create and monitor a Sander process.
     */
    SanderInterface();
    SanderInterface(const SanderInterface& orig);
    virtual ~SanderInterface();

    /**
     * Starts a Sander process. If an optional start_time is provided, it is assumed
     * by the sander process that the process began at start_time seconds when
     * reporting the amount of time the process has be running, when netRuntime()
     * is called.
     * 
     * @param start_time
     * @return
     */
    int start(const double& start_time = 0.0);

    /**
     * Checks on the status of the sander process. Returns true when the process
     * ends, either successfully or with an error. The return value of the process
     * is reported by the status variable provided to poll.
     * 
     * @param status
     * @return
     */
    bool poll(int& status);

    /**
     * Terminates the Sander process. In linux this is done using SIGKILL.
     * 
     */
    void kill();

    /**
     * Suspends the Sander process. In linux this is done using SIGSTOP.
     *
     */
    void stop();

    /**
     * Resumes a suspended Sander process. In this this is done using SIGCONT.
     */
    void resume();

    /**
     * Returns the starting time (in seconds) of the sander process. Generally,
     * this is zero.
     * 
     * @return
     */
    const double& start_time()const{return starting_cpu;}

    /**
     * Indicates whether or not the Sander process has been suspended.
     * 
     * @return
     */
    const bool& isSuspended()const{return suspended;}

    /**
     * Returns the process id (pid) of the Sander process.
     * 
     * @return
     */
#if defined(_WIN32) || defined(__MINGW_WIN32__)
    const DWORD& getPID()const{return pid;}
#else
    const int& getPID()const{return pid;}
#endif
    
    /**
     * Returns the total amount of time (in seconds) it took for the Sander
     * process to finish. This should only be called after the Sander process
     * finishes. To get the amount of cpu time of a current running process,
     * use cpu_time()
     * 
     * @return
     */
    const double netRuntime()const{return starting_cpu + final_cpu_time;}

    /**
     * Returns the amount of time the process has been running (in seconds).
     * This amount does not include the start time offset, if one was provided.
     * 
     * @return
     */
    double cpu_time()const;

    /**
     * Sets the process id of the Sander process to be monitored. This would be
     * used to re-attach to a Sander process.
     * 
     * @param newPID
     */
    void setPID(const int& newPID){pid = newPID;}

    /**
     * Sets the starting time (in seconds) of the Sander process.
     * 
     * @param new_start_time
     */
    void set_start_time(const double& new_start_time){starting_cpu = new_start_time;}

    std::string mdinFilename;
    std::string mdoutFilename;
    std::string mdcrdFilename;
    std::string prmtopFilename;
    std::string inpcrdFilename;
    std::string restartFilename;
    static const double pollPeriod = 1.0;///<Amount of time to wait between checking on the Sander process using poll()
    bool completed;
    
private:
    double wall_cpu_time;
    double starting_cpu;
    double final_cpu_time;
    bool suspended;

#if defined(_WIN32) || defined(__MINGW_WIN32__)
    HANDLE pid_handle;
    DWORD pid;
    HANDLE thread_handle;
    struct _stat last_stat;
#else
    int pid;
#endif
};


}//end namespace mmpbsa

#endif	/* SANDERINTERFACE_H */

