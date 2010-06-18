/* 
 * Various exceptions (hopefully not) thrown by mmpbsa classes.
 *
 * Created on June 18, 2010, 8:52 AM
 */
#include <stdexcept>

#ifndef MMPBSA_EXCEPTIONS_H
#define	MMPBSA_EXCEPTIONS_H

//Error types (don't use 0 )
static const int FILE_READ_ERROR = 1;
static const int BROKEN_PRMTOP_FILE = 2;

class SanderIOException : public std::runtime_error {
public:
    SanderIOException( const std::string& error) : runtime_error(error) {errorType = FILE_READ_ERROR;}

    /**
     * Creates an exception, with a specified error type.
     *
     * @param error
     * @param errorType
     */
    SanderIOException(const std::string& error, const int& errorType) : runtime_error(error){this->errorType = errorType;}

    /**
     * Returns the error type, corresponding to the error types listed below.
     * These should be returned if the exception is caught and the program
     * dies gracefully.
     * 
     * @return int error type.
     */
    int getErrType(){return errorType;}

    
private:
    int errorType;

};

#endif	/* MMPBSA_EXCEPTIONS_H */

