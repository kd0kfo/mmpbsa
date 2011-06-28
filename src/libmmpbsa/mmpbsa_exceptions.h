/**
 * @class mmpbsa::MMPBSAException
 * @brief Exception basis class and Error Codes
 *
 * Basis of exceptions (hopefully not) thrown by mmpbsa classes.
 * MMPBSAErrorTypes enumerates various error codes used by 
 * mmpbsa.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */
#include <stdexcept>
#include <sstream>

#ifndef MMPBSA_EXCEPTIONS_H
#define	MMPBSA_EXCEPTIONS_H

namespace mmpbsa{

//MMPBSAErrorTypes could be used as return values. Therefore, zero is not used.
enum MMPBSAErrorTypes {UNKNOWN_ERROR = 1, /*!General (default) error number. Should be used as a last resort as it is non-specific*/
    FILE_IO_ERROR,/*!IO problem reading file.*/
    BROKEN_PRMTOP_FILE,/*!prmtop file is improperly formatted or missing data.*/
    BROKEN_TRAJECTORY_FILE,/*!trajectory file is improperly formatted or missing data.*/
    INVALID_PRMTOP_DATA,/*!Data which was loaded into an array is incorrect based on what is expected.*/
    DATA_FORMAT_ERROR,/*!Use this error, when data *within* the program no longer matches what it should due to formatting problems.*/
    INVALID_ARRAY_SIZE,
    UNEXPECTED_EOF,
    COMMAND_LINE_ERROR,/*!Program supplied an invalid argument in the command line.*/
    BAD_XML_TAG,
    INVALID_XML_REQUEST, /*!Trying to access or modify nodes and node information that does not exist or is broken.*/
    NULL_POINTER,/* :-( */
    MPI_ERROR,
    SYSTEM_ERROR/* OS and/or system call problems */
};

class MMPBSAException : public std::runtime_error
{
    public:
        /**
         * Very General MMPBSAException. Use more specific exception if one is
         * available.
         *
         * @param error
         */
    MMPBSAException( const std::string& error) : runtime_error(error) {errorType = UNKNOWN_ERROR;}

    /**
     * Creates an exception, with a specified error type.
     *
     * @param error
     * @param errorType
     */
    MMPBSAException(const std::string& error, const MMPBSAErrorTypes& errorType) : runtime_error(error){this->errorType = errorType;}
    MMPBSAException(const std::ostringstream& error, const MMPBSAErrorTypes& errorType = UNKNOWN_ERROR) : runtime_error(error.str()){this->errorType = errorType;}
    /**
     * Returns the error type, corresponding to the error types listed below.
     * These should be returned if the exception is caught and the program
     * dies gracefully.
     *
     * @return int error type.
     */
    const MMPBSAErrorTypes& getErrType()const{return errorType;}

    virtual const char* identifier(){return "General MMPBSA Error";}

private:
    MMPBSAErrorTypes errorType;

};

};//end namespace mmpbsa

std::ostream& operator<<(std::ostream& theStream,mmpbsa::MMPBSAException& me);

#endif	/* MMPBSA_EXCEPTIONS_H */

