#include "mmpbsa_exceptions.h"

std::ostream& operator<<(std::ostream& theStream,mmpbsa::MMPBSAException& me)
{
    theStream << me.identifier() << "(" << me.getErrType() << "): " << me.what() << std::endl;
    return theStream;
}

