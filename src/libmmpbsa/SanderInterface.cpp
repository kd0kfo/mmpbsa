#include "SanderInterface.h"

mmpbsa::SanderInterface::SanderInterface() {
    pid = 0;
    double wall_cpu_time = 0;
    double final_cpu_time = 0;
    bool suspended = false;
}

mmpbsa::SanderInterface::SanderInterface(const mmpbsa::SanderInterface& orig)
{
    mdinFilename = orig.mdinFilename;
    mdoutFilename = orig.mdoutFilename;
    prmtopFilename = orig.prmtopFilename;
    inpcrdFilename = orig.inpcrdFilename;
    restartFilename = orig.restartFilename;
    pid = orig.pid;
    wall_cpu_time = orig.wall_cpu_time;
    final_cpu_time = orig.final_cpu_time;
    suspended = orig.suspended;
}

mmpbsa::SanderInterface::~SanderInterface() {
}

int mmpbsa::SanderInterface::start() {
    using std::string;
    string stdout_path = "sander-stdout.txt";
    string stdin_path = "sander-stdin.txt";
    string stderr_path = "sander-stdin.txt";
    string fraction_done_filename = "sander-pace.txt";

#ifdef __USE_BOINC__/* Clear files to be used by this calling of start()*/
    if (fraction_done_filename.size()) {
        boinc_delete_file(fraction_done_filename.c_str());
    }
#endif
    
    string application = "moldyn";
#ifdef __USE_BOINC__
    char buff[256];
    boinc_resolve_filename_s("moldyn", application);
#endif
    
    std::string command_line = "-O  -i " + mdinFilename + " -o  " + mdoutFilename
            + "-c " + inpcrdFilename + " -p " + prmtopFilename + " -r " + restartFilename;

    fprintf(stderr, "%s running with arguments: %s\n",
        application.c_str(), command_line.c_str()
    );

    

#ifdef _WIN32
    PROCESS_INFORMATION process_info;
    STARTUPINFO startup_info;
    string command;

    memset(&process_info, 0, sizeof(process_info));
    memset(&startup_info, 0, sizeof(startup_info));
    command = string("\"") + application + string("\" ") + command_line;

    // pass std handles to app
    //
    startup_info.dwFlags = STARTF_USESTDHANDLES;
	if (stdout_filename != "") {
		boinc_resolve_filename_s(stdout_filename.c_str(), stdout_path);
		startup_info.hStdOutput = win_fopen(stdout_path.c_str(), "a");
	}
	if (stdin_filename != "") {
		boinc_resolve_filename_s(stdin_filename.c_str(), stdin_path);
		startup_info.hStdInput = win_fopen(stdin_path.c_str(), "r");
	}
    if (stderr_filename != "") {
        boinc_resolve_filename_s(stderr_filename.c_str(), stderr_path);
        startup_info.hStdError = win_fopen(stderr_path.c_str(), "a");
    } else {
        startup_info.hStdError = win_fopen(STDERR_FILE, "a");
    }

    if (!CreateProcess(
        app_path,
        (LPSTR)command.c_str(),
        NULL,
        NULL,
        TRUE,		// bInheritHandles
        CREATE_NO_WINDOW|IDLE_PRIORITY_CLASS,
        NULL,
        NULL,
        &startup_info,
        &process_info
    )) {
        char error_msg[1024];
        windows_error_string(error_msg, sizeof(error_msg));
        fprintf(stderr, "can't run app: %s\n", error_msg);
        return ERR_EXEC;
    }
    pid_handle = process_info.hProcess;
    pid = process_info.dwProcessId;
    thread_handle = process_info.hThread;
    SetThreadPriority(thread_handle, THREAD_PRIORITY_IDLE);
#else
    int retval;
    

    pid = fork();
    if (pid == -1)
    {
        perror("fork(): ");
        return ERR_FORK;
    }
    if (pid == 0) 
    {

        // we're in the child process here
        //
        // open stdout, stdin if file names are given
        // NOTE: if the application is restartable,
        // we should deal with atomicity somehow
        //
        FILE* stdout_file;
        FILE* stdin_file;
        FILE* stderr_file;
        stdout_file = freopen(stdout_path.c_str(), "a", stdout);
        if (!stdout_file) return ERR_FOPEN;

        stdin_file = freopen(stdin_path.c_str(), "r", stdin);
        if (!stdin_file) return ERR_FOPEN;

        stderr_file = freopen(stderr_path.c_str(), "a", stderr);
        if (!stderr_file) return ERR_FOPEN;
		// construct argv
        // TODO: use malloc instead of stack var
        //
        string argv = application + " " + command_line;
        setpriority(PRIO_PROCESS, 0, PROCESS_IDLE_PRIORITY);
        char arglist[4096];
        ::strlcpy(arglist,argv.c_str(),sizeof(arglist));
        char* argvs[256];
        argvs[0] = "";
        strcat(argvs[0],application.c_str());
        int argc = parse_command_line(arglist, argvs+1);
        retval = execv(application.c_str(), argvs);
        perror("execv() failed: ");
        exit(ERR_EXEC);
    }
#endif
    wall_cpu_time = 0;
    suspended = false;
    return 0;
}


bool mmpbsa::SanderInterface::poll(int& status) {
    if (!suspended) wall_cpu_time += pollPeriod;
#ifdef _WIN32
    unsigned long exit_code;
    if (GetExitCodeProcess(pid_handle, &exit_code)) {
        if (exit_code != STILL_ACTIVE) {
            status = exit_code;
            final_cpu_time = cpu_time();
            return true;
        }
    }
#else
    int wpid, stat;
    struct rusage ru;

    wpid = wait4(pid, &status, WNOHANG, &ru);
    if (wpid) {
        final_cpu_time = (float)ru.ru_utime.tv_sec + ((float)ru.ru_utime.tv_usec)/1e+6;
        return true;
    }
#endif
    return false;
}

void mmpbsa::SanderInterface::kill() {
#ifdef _WIN32
    TerminateProcess(pid_handle, -1);
#else
    ::kill(pid, SIGKILL);
#endif
}

void mmpbsa::SanderInterface::stop() {
#ifdef _WIN32
    suspend_or_resume_threads(pid, 0, false);
#else
    ::kill(pid, SIGSTOP);
#endif
    suspended = true;
}

void mmpbsa::SanderInterface::resume() {
#ifdef _WIN32
    suspend_or_resume_threads(pid, 0, true);
#else
    ::kill(pid, SIGCONT);
#endif
    suspended = false;
}


void mmpbsa::SanderInterface::poll_boinc_messages()
{
    BOINC_STATUS status;
    boinc_get_status(&status);
    if (status.no_heartbeat) {
        this->kill();
    }
    if (status.quit_request) {
        this->kill();
    }
    if (status.abort_request) {
        this->kill();
    }
    if (status.suspended) {
        if (!suspended) {
            this->stop();
        }
    } else {
        if (this->suspended) {
            this->resume();
        }
    }
}

double mmpbsa::SanderInterface::cpu_time()const
{
#ifdef _WIN32
    double x;
    int retval = boinc_process_cpu_time(pid_handle, x);
    if (retval) return wall_cpu_time;
    return x;
#elif defined(__APPLE__)
    // There's no easy way to get another process's CPU time in Mac OS X
    //
    return wall_cpu_time;
#else
    return linux_cpu_time(pid);
#endif
}

