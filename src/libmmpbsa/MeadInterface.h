/**
 * @class mmpbsa::MeadInterface
 * @brief Interface for the MEAD library
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */

#ifndef MeadInterface_H
#define	MeadInterface_H

#include <valarray>
#include <map>
#include "mmpbsa_exceptions.h"
#include "globals.h"
#include "structs.h"

//forward declare mead class
class FinDiffMethod;
class Coord;

namespace mmpbsa{
//forward declarations
class Vector;
class EMap;
class EmpEnerFun;

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
    static FinDiffMethod createFDM(const std::valarray<mmpbsa::Vector>& complexCrds,
        const std::valarray<mmpbsa::Vector>& receptorCrds, const std::valarray<mmpbsa::Vector>& ligandCrds,
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
    static EMap full_EMap(const EmpEnerFun& efun, const std::valarray<mmpbsa::Vector>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,float>& radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (mmpbsa::MeadException);

    static EMap full_EMap(const std::vector<mmpbsa::atom_t>& atoms, const mmpbsa::forcefield_t& ff, const std::valarray<mmpbsa::Vector>& crds,
            const FinDiffMethod& fdm, const std::map<std::string,float>& radii,
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
    static mmpbsa_t* pbsa_solvation(const EmpEnerFun& efun, const std::valarray<mmpbsa::Vector>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,float>& radii,
        const std::map<std::string,std::string>& residueMap,
        const mmpbsa_t& interactionStrength = 0.0, const mmpbsa_t& exclusionRadius = 2.0) throw (mmpbsa::MeadException);

    static mmpbsa_t* pbsa_solvation(const std::vector<mmpbsa::atom_t>& atoms, const mmpbsa::forcefield_t& ff,
    		const std::valarray<mmpbsa::Vector>& crds,
    		const FinDiffMethod& fdm, const std::map<std::string,float>& radii,
    		const std::map<std::string,std::string>& residueMap,
    		const mmpbsa_t& interactionStrength = 0.0, const mmpbsa_t& exclusionRadius = 2.0) throw (mmpbsa::MeadException);


    mmpbsa_t istrength;
    mmpbsa_t surf_tension;///<kcal/mol/Ang^2
    mmpbsa_t surf_offset;///<kcal/mol

    /**
     *
     * offset used to recombine snaplist segments.
     * Snap list are counted in order of appearance in the trajectory file.
     * To ensure that multiple MMPBSA calculation can be combined in the correct order,
     * an offset may be used to be added to the snap shot's index when recombining results.
     */
    int snap_list_offset;

    std::map<std::string,mead_data_t> brad;///<Lookup table for default radius values

    int multithread;///<used to indicate number of threads to be used in MMPBSA calculations. Default = 1
    
};

};//end namespace mmpbsa

Coord ToCoord(const mmpbsa::Vector& v);

#endif	/* MeadInterface_H */

