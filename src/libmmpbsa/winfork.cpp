#include "mmpbsa_utils.h"
#include "mmpbsa_io.h"

#include <windows.h> 
#include <winbase.h>
#include <tchar.h>
#include <stdio.h>
#include <sstream>

#define winfork_BUFSIZE 4096 
 
HANDLE g_hChildStd_IN_Rd = NULL;
HANDLE g_hChildStd_IN_Wr = NULL;
HANDLE g_hChildStd_OUT_Rd = NULL;
HANDLE g_hChildStd_OUT_Wr = NULL;

HANDLE g_hInputFile = NULL;
 
void CreateChildProcess(const std::vector<mmpbsa::atom_t>& atoms,
			const std::valarray<mmpbsa::Vector>& crds,
			const std::map<std::string,mead_data_t>& radii,
			int *error_flag);
void WriteToPipe(const std::vector<mmpbsa::atom_t>& atoms,
		 const std::valarray<mmpbsa::Vector>& crds,
		 const std::map<std::string,mead_data_t>& radii,
		 int *error_flag);
mmpbsa_t ReadFromPipe(int *error_flag);
void ErrorExit(PTSTR); 
 
mmpbsa_t molsurf_win32(const std::vector<mmpbsa::atom_t>& atoms,
		       const std::valarray<mmpbsa::Vector>& crds,
		       const std::map<std::string,mead_data_t>& radii,
		       int *error_flag)
{ 
   SECURITY_ATTRIBUTES saAttr; 
   mmpbsa_t area = 0;
   //printf("\n->Start of parent execution.\n");

// Set the bInheritHandle flag so pipe handles are inherited. 
 
   saAttr.nLength = sizeof(SECURITY_ATTRIBUTES); 
   saAttr.bInheritHandle = TRUE; 
   saAttr.lpSecurityDescriptor = NULL; 

// Create a pipe for the child process's STDOUT. 
 
   if ( ! CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0) ) 
      ErrorExit(TEXT("StdoutRd CreatePipe")); 

// Ensure the read handle to the pipe for STDOUT is not inherited.

   if ( ! SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0) )
      ErrorExit(TEXT("Stdout SetHandleInformation")); 

// Create a pipe for the child process's STDIN. 
 
   if (! CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0)) 
      ErrorExit(TEXT("Stdin CreatePipe")); 

// Ensure the write handle to the pipe for STDIN is not inherited. 
 
   if ( ! SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0) )
      ErrorExit(TEXT("Stdin SetHandleInformation")); 
 
// Create the child process. 
   
   CreateChildProcess(atoms,crds,radii,error_flag);

// Get a handle to an input file for the parent. 
// This example assumes a plain text file and uses string output to verify data flow. 
 
   /*if (argc == 1) 
     ErrorExit(TEXT("Please specify an input file.\n")); */

   /*   g_hInputFile = CreateFile(
       argv[1], 
       GENERIC_READ, 
       0, 
       NULL, 
       OPEN_EXISTING, 
       FILE_ATTRIBUTE_READONLY, 
       NULL); 

   if ( g_hInputFile == INVALID_HANDLE_VALUE ) 
   ErrorExit(TEXT("CreateFile")); */
 
// Write to the pipe that is the standard input for a child process. 
// Data is written to the pipe's buffers, so it is not necessary to wait
// until the child process is running before writing data.
 
   WriteToPipe(atoms,crds,radii,error_flag); 
   //printf( "\n->Contents of %s written to child STDIN pipe.\n", argv[1]);
 
// Read from pipe that is the standard output for child process. 
 
   //printf( "\n->Contents of child process STDOUT??:\n\n");
   area = ReadFromPipe(error_flag);

   //printf("\n->End of parent execution.\n");

// The remaining open handles are cleaned up when this process terminates. 
// To avoid resource leaks in a larger application, close handles explicitly. 

   return area; 
} 
 
void CreateChildProcess(const std::vector<mmpbsa::atom_t>& atoms,
						  const std::valarray<mmpbsa::Vector>& crds,
						  const std::map<std::string,mead_data_t>& radii,
						  int *error_flag)
// Create a child process that uses the previously created pipes for STDIN and STDOUT.
{ 
  const size_t max_chars = 4096;
  CHAR cmdstr[max_chars], appname[max_chars];
  PROCESS_INFORMATION piProcInfo; 
  STARTUPINFO siStartInfo;
  BOOL bSuccess = FALSE;

  if(mmpbsa_io::resolve_filename("./xyzr2sas.exe",appname,max_chars) != 0)
    ErrorExit(TEXT("Could not resolve the area application filename"));
  
  snprintf(cmdstr,max_chars,"\"%s\" %lu",appname,crds.size());			      
 
// Set up members of the PROCESS_INFORMATION structure. 
 
   ZeroMemory( &piProcInfo, sizeof(PROCESS_INFORMATION) );
 
// Set up members of the STARTUPINFO structure. 
// This structure specifies the STDIN and STDOUT handles for redirection.
 
   ZeroMemory( &siStartInfo, sizeof(STARTUPINFO) );
   siStartInfo.cb = sizeof(STARTUPINFO); 
   //siStartInfo.hStdError = g_hChildStd_OUT_Wr;
   siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;
   siStartInfo.hStdInput = g_hChildStd_IN_Rd;
   siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
 
// Create the child process. 
   UINT WINAPI orig_err_mask = SetErrorMode(SEM_NOGPFAULTERRORBOX);
   SetErrorMode(orig_err_mask | SEM_NOGPFAULTERRORBOX);

   bSuccess = CreateProcess(appname, 
			    cmdstr,     // command line
      NULL,          // process security attributes 
      NULL,          // primary thread security attributes 
      TRUE,          // handles are inherited 
			    CREATE_NO_WINDOW,			    //original: 0,             // creation flags 
      NULL,          // use parent's environment 
      NULL,          // use parent's current directory 
      &siStartInfo,  // STARTUPINFO pointer 
      &piProcInfo);  // receives PROCESS_INFORMATION 
   
   SetErrorMode(orig_err_mask);

   // If an error occurs, exit the application. 
   if ( ! bSuccess ) 
     {
       if(error_flag != NULL)
	 *error_flag  = 1;
       ErrorExit(TEXT("Could not start child process"));
     }
   else 
   {
      // Close handles to the child process and its primary thread.
      // Some applications might keep these handles to monitor the status
      // of the child process, for example. 
      CloseHandle(piProcInfo.hProcess);
      CloseHandle(piProcInfo.hThread);
   }
}
 
void WriteToPipe(const std::vector<mmpbsa::atom_t>& atoms,
		 const std::valarray<mmpbsa::Vector>& crds,
		 const std::map<std::string,mead_data_t>& radii,
		 int *error_flag) 

// Read from a file and write its contents to the pipe for the child's STDIN.
// Stop when there is no more data. 
{ 
   DWORD dwRead, dwWritten; 
   BOOL bSuccess = FALSE;
   std::vector<mmpbsa::atom_t>::const_iterator curr_atom = atoms.begin();
   size_t atom_counter = 0;
   std::ostringstream buffer;
   static const char delim[] = " ";
   
   for(;curr_atom != atoms.end();curr_atom++,atom_counter++)
     {
       const mmpbsa::Vector& curr_pos = crds[atom_counter];
       buffer.str("");buffer.clear();
       buffer << curr_pos.x() << delim << curr_pos.y() << delim << curr_pos.z() << delim << mmpbsa_utils::lookup_radius(curr_atom->name,radii) << std::endl;
       bSuccess = WriteFile(g_hChildStd_IN_Wr, buffer.str().c_str(), buffer.str().size(), &dwWritten, NULL);
       if ( ! bSuccess )
	 {
	   if(error_flag != NULL)
	     *error_flag = 1;
	   break; 
	 }
     } 
// Close the pipe handle so the child process stops reading. 
   if ( ! CloseHandle(g_hChildStd_IN_Wr) ) 
     fprintf(stderr,"Could not close child input handle.");
} 
 
mmpbsa_t ReadFromPipe(int *error_flag)

// Read output from the child process's pipe for STDOUT
// and write to the parent process's pipe for STDOUT. 
// Stop when there is no more data. 
{ 
   DWORD dwRead, dwWritten; 
   CHAR chBuf[winfork_BUFSIZE],*strtod_result; 
   BOOL bSuccess = FALSE;
   HANDLE hParentStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
   mmpbsa_t childresult = 0;
   int error_val = 0;
// Close the write end of the pipe before reading from the 
// read end of the pipe, to control child process execution.
// The pipe is assumed to have enough buffer space to hold the
// data the child process has already written to it.
 
   if (!CloseHandle(g_hChildStd_OUT_Wr)) 
      ErrorExit(TEXT("StdOutWr CloseHandle")); 
 
   bSuccess = ReadFile( g_hChildStd_OUT_Rd, chBuf, winfork_BUFSIZE, &dwRead, NULL);
   if( ! bSuccess || dwRead == 0 ) 
     {
       if(error_flag != NULL)
	 *error_flag = 1;
       return 0.0;
     }
      
   //if(strchr(chBuf,'\n') != NULL)
   //*strchr(chBuf,'\n') = (char)0;

   if((childresult = strtod(chBuf,&strtod_result)) == 0)
     {
       if(strtod_result == chBuf)
	 {
	   fprintf(stderr,"Parent: could not interpret child data: %s\n",chBuf);
	   if(error_flag != NULL)
	     *error_flag = 1;
	   return 0.0;
	 }
     }
   //printf("Read %f as a double from child (%s).\n",childresult,chBuf);
      
   if(error_flag != NULL)
	   *error_flag = error_val;
   return childresult;
} 
 
void ErrorExit(PTSTR lpszFunction) 

// Format a readable error message, display a message box, 
// and exit from the application.
{ 
    LPVOID lpMsgBuf;
    LPVOID lpDisplayBuf;
    DWORD dw = GetLastError(); 

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL );

    fprintf(stderr,"%s failed with error %d: %s", lpszFunction, dw, lpMsgBuf); 

    LocalFree(lpMsgBuf);
    LocalFree(lpDisplayBuf);
    ExitProcess(1);
}
