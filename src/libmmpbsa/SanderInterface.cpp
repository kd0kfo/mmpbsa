#include "SanderInterface.h"
#include "mmpbsa_io.h"

mmpbsa::SanderInterface::SanderInterface() {
    pid = 0;
    wall_cpu_time = 0;
    final_cpu_time = 0;
    mdcrdFilename = "";
    suspended = false;
    completed = false;
}

mmpbsa::SanderInterface::SanderInterface(const mmpbsa::SanderInterface& orig)
{
    mdinFilename = orig.mdinFilename;
    mdoutFilename = orig.mdoutFilename;
    prmtopFilename = orig.prmtopFilename;
    inpcrdFilename = orig.inpcrdFilename;
    restartFilename = orig.restartFilename;
    mdcrdFilename = orig.mdcrdFilename;
    pid = orig.pid;
    wall_cpu_time = orig.wall_cpu_time;
    final_cpu_time = orig.final_cpu_time;
    suspended = orig.suspended;
    completed = orig.completed;
}

mmpbsa::SanderInterface::~SanderInterface() {
}


#ifndef __USE_BOINC__
//These are needed for start() below. However, if BOINC is not used, it must
//be provided here.
//Copied from str_util.cpp under the terms of the GNU Lesser General Public License.
//See http://boinc.berkeley.edu or http://www.gnu.org for details.

#define NOT_IN_TOKEN                0
#define IN_SINGLE_QUOTED_TOKEN      1
#define IN_DOUBLE_QUOTED_TOKEN      2
#define IN_UNQUOTED_TOKEN           3

int parse_command_line(char* p, char** argv) {
    int state = NOT_IN_TOKEN;
    int argc=0;

    while (*p) {
        switch(state) {
        case NOT_IN_TOKEN:
            if (isspace(*p)) {
            } else if (*p == '\'') {
                p++;
                argv[argc++] = p;
                state = IN_SINGLE_QUOTED_TOKEN;
                break;
            } else if (*p == '\"') {
                p++;
                argv[argc++] = p;
                state = IN_DOUBLE_QUOTED_TOKEN;
                break;
            } else {
                argv[argc++] = p;
                state = IN_UNQUOTED_TOKEN;
            }
            break;
        case IN_SINGLE_QUOTED_TOKEN:
            if (*p == '\'') {
                *p = 0;
                state = NOT_IN_TOKEN;
            }
            break;
        case IN_DOUBLE_QUOTED_TOKEN:
            if (*p == '\"') {
                *p = 0;
                state = NOT_IN_TOKEN;
            }
            break;
        case IN_UNQUOTED_TOKEN:
            if (isspace(*p)) {
                *p = 0;
                state = NOT_IN_TOKEN;
            }
            break;
        }
        p++;
    }
    argv[argc] = 0;
    return argc;
}

#endif


#if defined(_WIN32) || defined(__MINGW_WIN32__)
//win_fopen(), suspend_or_resume_threads(...) and windows_error_string() have been created based on the code
//from the BOINC project under the terms of the Lesser GNU Public License.
//See http://boinc.berkeley.edu or http://gnu.org for details.
HANDLE win_fopen(const char* path, const char* mode)
{
    SECURITY_ATTRIBUTES sa;
	memset(&sa, 0, sizeof(sa));
	sa.nLength = sizeof(sa);
	sa.bInheritHandle = TRUE;

	if (!strcmp(mode, "r")) {
		return CreateFile(
			path,
			GENERIC_READ,
			FILE_SHARE_READ,
			&sa,
			OPEN_EXISTING,
			0, 0
		);
	} else if (!strcmp(mode, "w")) {
		return CreateFile(
			path,
			GENERIC_WRITE,
			FILE_SHARE_WRITE,
			&sa,
			OPEN_ALWAYS,
			0, 0
		);
	} else if (!strcmp(mode, "a")) {
		HANDLE hAppend = CreateFile(
			path,
			GENERIC_WRITE,
			FILE_SHARE_WRITE,
			&sa,
			OPEN_ALWAYS,
			0, 0
		);
        SetFilePointer(hAppend, 0, NULL, FILE_END);
        return hAppend;
	} else {
		return 0;
	}
}

#ifndef __USE_BOINC__
//provided by boinc/lib/str_util.cpp
// get message for last error
//
char* windows_error_string(char* pszBuf, int iSize) {
    DWORD dwRet;
    LPSTR lpszTemp = NULL;

    dwRet = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_ARGUMENT_ARRAY,
        NULL,
        GetLastError(),
        LANG_NEUTRAL,
        (LPSTR)&lpszTemp,
        0,
        NULL
    );

    // supplied buffer is not long enough
    if ( !dwRet || ( (long)iSize < (long)dwRet+14 ) ) {
        pszBuf[0] = '\0';
    } else {
        lpszTemp[lstrlenA(lpszTemp)-2] = '\0';  //remove cr and newline character
        sprintf ( pszBuf, "%s (0x%x)", lpszTemp, GetLastError() );
    }

    if ( lpszTemp ) {
        LocalFree((HLOCAL) lpszTemp );
    }

    return pszBuf;
}

//provided by boinc/lib/win_util.cpp
// Suspend or resume the threads in a given process.
// The only way to do this on Windows is to enumerate
// all the threads in the entire system,
// and find those belonging to the process (ugh!!)
//

// OpenThread
typedef HANDLE (WINAPI *tOT)(DWORD dwDesiredAccess, BOOL bInheritHandle, DWORD dwThreadId);

int suspend_or_resume_threads(DWORD pid, bool resume) {
    HANDLE threads, thread;
    HMODULE hKernel32Lib = NULL;
    THREADENTRY32 te = {0};
    tOT pOT = NULL;

    // Dynamically link to the proper function pointers.
    hKernel32Lib = GetModuleHandleA("kernel32.dll");
    pOT = (tOT) GetProcAddress( hKernel32Lib, "OpenThread" );

    if (!pOT) {
        return -1;
    }

    threads = CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0);
    if (threads == INVALID_HANDLE_VALUE) return -1;

    te.dwSize = sizeof(THREADENTRY32);
    if (!Thread32First(threads, &te)) {
        CloseHandle(threads);
        return -1;
    }

    do {
        if (te.th32OwnerProcessID == pid) {
            thread = pOT(THREAD_SUSPEND_RESUME, FALSE, te.th32ThreadID);
            resume ?  ResumeThread(thread) : SuspendThread(thread);
            CloseHandle(thread);
        }
    } while (Thread32Next(threads, &te));

    CloseHandle (threads);

    return 0;
}
#endif//use boinc

#endif

int mmpbsa::SanderInterface::start(const double& start_time) {
    using std::string;
    this->starting_cpu = start_time;
    string stdout_path = "sander-stdout.txt";
    string stdin_path = "sander-stdin.txt";
    string stderr_path = "sander-stdin.txt";
    
    
    char buff[] = "./moldyn";
    char application[1024];
    mmpbsa_io::resolve_filename(buff, application,sizeof(application));
    
    std::string command_line = "-O -i " + mdinFilename + " -o  " + mdoutFilename
            + " -c " + inpcrdFilename + " -p " + prmtopFilename + " -r " + restartFilename;

    if(mdcrdFilename.size())
        command_line += " -x " + mdcrdFilename;


    std::cout << application << " running with arguments: " << command_line << std::endl;
    

#if defined(_WIN32)
    PROCESS_INFORMATION process_info;
    STARTUPINFO startup_info;
    string command;

    memset(&process_info, 0, sizeof(process_info));
    memset(&startup_info, 0, sizeof(startup_info));
    command = string("\"") + application + string("\" ") + command_line;

    // pass std handles to app
    //
    startup_info.dwFlags = STARTF_USESTDHANDLES;
	if (stdout_path != "") {
		startup_info.hStdOutput = win_fopen(stdout_path.c_str(), "a");
	}
    else {
        startup_info.hStdError = win_fopen("stderr.txt", "a");
    }

    if (!CreateProcess(
        application,
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
    std::cerr << "PID: " << pid << std::endl;
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
        stdout_file = freopen(stdout_path.c_str(), "a", stdout);
        if (!stdout_file) return ERR_FOPEN;
        //stderr_file = freopen(stderr_path.c_str(), "a", stderr);
        //if (!stderr_file) return ERR_FOPEN;
        	// construct argv
        // TODO: use malloc instead of stack var
        //
        std::cerr << "Generating arguments" << std::endl;
        char* argvs[256];
        argvs[0] = application;
        char arglist[4096];
        strlcpy(arglist,command_line.c_str(),sizeof(arglist));
        int argc = parse_command_line(arglist,argvs+1);
        std::cerr << "Creating fork" << std::endl;
        setpriority(PRIO_PROCESS, 0, PROCESS_IDLE_PRIORITY);
        retval = execv(application,argvs);
        perror("execl() failed: ");
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
    suspend_or_resume_threads(pid, false);
#else
    ::kill(pid, SIGSTOP);
#endif
    suspended = true;
}

void mmpbsa::SanderInterface::resume() {
#ifdef _WIN32
    suspend_or_resume_threads(pid, true);
#else
    ::kill(pid, SIGCONT);
#endif
    suspended = false;
}


#ifndef __USE_BOINC__
//These are needed for start() below. However, if BOINC is not used, it must
//be provided here.
//Copied from str_util.cpp from the BOINC project under the terms of the
//GNU Lesser General Public License.
//See http://boinc.berkeley.edu or http://www.gnu.org for details.
#ifndef _USING_FCGI_
#ifndef _WIN32
// (linux) return current CPU time of the given process
//
double linux_cpu_time(int pid) {
    FILE *file;
    char file_name[24];
    unsigned long utime = 0, stime = 0;
    int n;

    sprintf(file_name,"/proc/%d/stat",pid);
    if ((file = fopen(file_name,"r")) != NULL) {
        n = fscanf(file,"%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lu%lu",&utime,&stime);
        fclose(file);
        if (n != 2) return 0;
    }
    return (double)(utime + stime)/100;
}
#endif
#endif
#endif

double mmpbsa::SanderInterface::cpu_time()const
{
#ifdef _WIN32
    double x;
    int retval = 1;
#ifdef __USE_BOINC__
    retval = boinc_process_cpu_time(pid_handle, x);
#endif
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

