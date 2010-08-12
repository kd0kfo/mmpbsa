/* 
 * Interface for the MEAD library
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
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
#include "EmpEnerFun.h"
#include "mmpbsa_io.h"
#include "molsurf/molsurf.h"

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

namespace mmpbsa{
class MeadException : public MMPBSAException {
    public:

        /**
         * Exception for when something goes wrong with MMPBSA using MEAD
         *
         * @param error
         */
        MeadException(const std::string& error) 
            : mmpbsa::MMPBSAException(error) {}

        MeadException(const std::string& error, const mmpbsa::MMPBSAErrorTypes& errorType)
            : MMPBSAException(error, errorType) {}

        const char* identifier()
        {
            return "MEAD/MMPBSA Error";
        }
    };

class MeadInterface {
public:
    /**
     * MeadInteraface stores variable values that are used by Mead
     *
     * Default:
     * istrength = 0;
     * surf_tension =  0.00542;// kcal/mol/Ang^2
     * surf_offset = 0.92;// kcal/mol
     */
    MeadInterface();

    MeadInterface(const MeadInterface& orig);
    virtual ~MeadInterface();

    /**
     * Sets up a Finite Difference Method object based on the provided coordinates.
     * 
     * @param complexCrds
     * @param receptorCrds
     * @param ligandCrds
     * @param outbox_grid_dim
     * @param fine_grid_spacing
     * @return
     */
    static FinDiffMethod createFDM(const std::valarray<mmpbsa_t>& complexCrds,
        const std::valarray<mmpbsa_t>& receptorCrds, const std::valarray<mmpbsa_t>& ligandCrds,
        const int& outbox_grid_dim = 41, const mmpbsa_t& fine_grid_spacing = 0.25) throw (mmpbsa::MeadException);

    /**
     * Creates an EMap object which will included Mead and molsurf calculated
     * energies and surface area.
     * 
     * radii is a pointer so that a null pointer can be sent to indicate the 
     * use of default values from a lookup table stored in MeadInterface.
     * 
     * @param efun
     * @param crds
     * @param fdm
     * @param radii
     * @param residueMap
     * @param interactionStrength
     * @param surfTension
     * @param surfOffset
     * @return 
     */
    static EMap full_EMap(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, std::map<std::string,mmpbsa_t>& radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (mmpbsa::MeadException);

    /**
     * Returns an array containing PB Solvation Energy and surface area, respectively.
     * 
     * @param efun
     * @param crds
     * @param fdm
     * @param radii
     * @param residueMap
     * @param interactionStrength
     * @param exclusionRadius
     * @return 
     */
    static mmpbsa_t* pbsa_solvation(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, std::map<std::string,mmpbsa_t>& radii,
        const std::map<std::string,std::string>& residueMap,
        const mmpbsa_t& interactionStrength = 0.0, const mmpbsa_t& exclusionRadius = 2.0) throw (mmpbsa::MeadException);

    //mmpbsa_t bondi_lookup(const std::string& atomName)const;

    mmpbsa_t istrength;
    mmpbsa_t surf_tension;// kcal/mol/Ang^2
    mmpbsa_t surf_offset;// kcal/mol

    std::map<std::string,mmpbsa_t> brad;

    
};



};//end namespace mmpbsa

#endif	/* MeadInterface_H */

