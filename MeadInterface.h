/* 
 * File:   MeadInterface.h
 * Author: dcoss
 *
 * Created on July 9, 2010, 11:25 AM
 */

#ifndef MeadInterface_H
#define	MeadInterface_H

//std lib
#include <cmath>
#include <valarray>
#include <algorithm>

//mmpbsa stuff
#include "mmpbsa_utils.h"
#include "mmpbsa_exceptions.h"
#include "Energy.h"
#include "mmpbsa_io.h"

//MEAD
#include "MEAD/Coord.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/CenteringStyle.h"
#include "MEAD/AtomSet.h"
#include "MEAD/Atom.h"
#include "MEAD/AtomChargeSet.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/DielByAtoms.h"
#include "MEAD/ElectrolyteByAtoms.h"
#include "MEAD/FinDiffElstatPot.h"
#include "MEAD/UniformDielectric.h"
#include "MEAD/UniformElectrolyte.h"

class MeadInterface {
public:
    MeadInterface();
    MeadInterface(const MeadInterface& orig);
    virtual ~MeadInterface();

    static FinDiffMethod createFDM(const std::valarray<mmpbsa_t>& complexCrds,
        const std::valarray<mmpbsa_t>& receptorCrds, const std::valarray<mmpbsa_t>& ligandCrds,
        const int& outbox_grid_dim = 41, const mmpbsa_t& fine_grid_spacing = 0.25) throw (MMPBSAException);

    //radii is a pointer so that a null pointer can be sent to indicate the use of default values from a lookup table
    static EMap full_EMap(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mmpbsa_t>* radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (MMPBSAException);

    static mmpbsa_t* pbsa_solvation(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mmpbsa_t>* radii,
        const std::map<std::string,std::string>& residueMap,
        const mmpbsa_t& interactionStrength = 0.0, const mmpbsa_t& exclusionRadius = 2.0) throw (MMPBSAException);

    static mmpbsa_t bondi_lookup(const std::string& atomName);
private:

};

class MeadException : public MMPBSAException {
public:
    /**
     * Exception for when something goes wrong with MMPBSA using MEAD
     * 
     * @param error
     */
    MeadException(const std::string& error) : MMPBSAException( error){}

    MeadException(const std::string& error, const MMPBSAErrorTypes& errorType)
        : MMPBSAException(error,errorType){}

    const char* identifier(){return "MEAD/MMPBSA Error";}
};


#endif	/* MeadInterface_H */

