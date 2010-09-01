#include "mmpbsa_exceptions.h"

std::ostream& mmpbsa::operator<<(std::ostream& theStream,MMPBSAException& me)
{
    theStream << me.identifier() << "(" << me.getErrType() << "): " << me.what() << std::endl;
    return theStream;
}

