/**
 * @class mmpbsa::EmpEnerFun
 * @brief Energy Abstraction. 
 * Encapsulates data from Sander Parmtop files which is then used to
 * perform MMPBSA calculations.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */
#ifndef ENERGY_H
#define	ENERGY_H

#include <cmath>
#include <valarray>
#include <vector>
#include <fstream>

#include "mmpbsa_utils.h"
#include "mmpbsa_io.h"
#include "SanderParm.h"
#include "structs.h"
#include "Energy.h"
#include "XMLNode.h"
#include "XMLParser.h"

namespace mmpbsa{

class EmpEnerFun{
public:
    /**
     * Default constructor. Uses default values of scnb, scee and dielc
     * Valarray pointers are set to null.
     * parminfo and resnames are set to null.
     * 
     * @param newparminfo
     */
    EmpEnerFun();

    /**
     * Main Energy Constructor
     * After a Parmtop file has been loaded into a SanderParm object, that object
     * may be used to do energy calculations with the EmpEnerFun class.
     * A pointer to the SanderParm object is stored in EmpEnerFun. Therefore,
     * data is store in one location and shared between SanderParm and EmpEnerFun
     * classes.
     * 
     * @param newparminfo
     * @param scnb
     * @param scee
     * @param dielc
     */
    EmpEnerFun(mmpbsa::SanderParm * newparminfo, const mmpbsa_t& scnb = DEFAULT_SCNB,
            const mmpbsa_t& scee = DEFAULT_SCEE, const mmpbsa_t& dielc = DEFAULT_DIELC);

    /**
     * Copy Constructor
     * 
     * @param orig
     */
    EmpEnerFun(const EmpEnerFun& orig);

    
    /**
     * Destructor. Cleans up pointers used by EmpEnerFun. SanderParm must be 
     * deleted separately, as it might be used elsewhere.
     */
    ~EmpEnerFun();//do not delete parminfo. It is externally made and should be deleted outside of EmpEnerInfo
    

    /**
     * Returns an EmpEnerFun object containg only data need to do energy calculations
     * with the atoms correspoding to "true" values in the boolean valarray.
     * Additionally, a warning can be given if a bond has an atom in it that is not kept,
     * while other atoms are.
     *
     * @param keepers
     * @param dangleWarn
     * @return
     */
    EmpEnerFun stripEnerFun(const std::valarray<bool>& keepers,
        const bool& dangleWarn) const;

    /**
     * Resizes and shuffles data arrays to only use atoms which are flagged to be
     * kept. All data values in the new data array will contain the new index
     * in the stripped Energy Object. The last data column in the array contains
     * pointers to unchanged energy data, eg equilibrium values, and therefore
     * remains unchanged.
     * 
     * @param newIndices
     * @param oldIndices
     * @param slices
     * @param newidx
     * @param keepers
     * @param dangleWarn
     */
    template <class M>
    static void internalConvert(std::valarray<M>& newIndices,
        const std::valarray<M>& oldIndices,const std::valarray<std::slice>& slices,
        const std::valarray<size_t>& newidx,
        const std::valarray<bool>& keepers, const bool& dangleWarn = false);

    /**
     * Shuffles the dihedral masks to correspond to the new locations of the
     * dihedral components in the new dihedral pointer arrays.
     * 
     * @param newMask
     * @param oldMask
     * @param oldIndices
     * @param slices
     * @param keepers
     * @param dangleWarn
     */
    template <class M> static void updatePhiMasks(std::valarray<bool>& newMask,
    const std::valarray<bool>& oldMask,
    const std::valarray<M>& oldIndices, const std::valarray<std::slice>& slices,
        const std::valarray<bool>& keepers, const bool& dangleWarn);


    //Energy Calculations
    /**
     * Calculates the total bond energy for the given snapshot coordinates
     * 
     * @param crds
     * @return 
     */
    mmpbsa_t total_bond_energy(const std::valarray<mmpbsa::Vector>& crds)const;
    mmpbsa::bond_energy_t* extract_bond_structs(std::vector<bond_t>& bonds_with_H,std::vector<bond_t>& bonds_without_H)const;

    /**
     * Calculates the total angle energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_angle_energy(const std::valarray<mmpbsa::Vector>& crds)const;
    bond_energy_t* extract_angle_structs(std::vector<mmpbsa::angle_t>& angles_with_H, std::vector<mmpbsa::angle_t>& angles_without_H)const;
    /**
     * Calculates the total dihedral energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_dihedral_energy(const std::valarray<mmpbsa::Vector>& crds)const;

    mmpbsa::dihedral_energy_t* extract_dihedral_structs(std::vector<mmpbsa::dihedral_t>& dihedrals_with_H,std::vector<mmpbsa::dihedral_t>& dihedrals_without_H)const;
    void extract_atom_structs(std::vector<mmpbsa::atom_t>& atoms)const;

    /**
     * Calculates the total Van der Waals energy between 1-4 pairs for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_vdw14_energy(const std::valarray<mmpbsa::Vector>& crds)const;

    /**
     * Calculates the total Electrostatic energy between 1-4 pairs for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_elstat14_energy(const std::valarray<mmpbsa::Vector>& crds)const;

    /**
     * Calculates the total Van der Waals energy of the system for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_vdwaals_energy(const std::valarray<mmpbsa::Vector>& crds)const;

    /**
     * Calculates the total Electrostatic energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_elstat_energy(const std::valarray<mmpbsa::Vector>& crds)const;

    /**
     * Generates a report of the energy values for the given snapshot using
     * the Parmtop parameters stored in this object.
     *
     * Format:
     *      BOND    =  %12.4f  ANGLE   =  %12.4f  DIHED      =  %12.4f
     *      VDWAALS =  %12.4f  EEL     =  %12.4f  EGB        =  %12.4f
     *      1-4 VDW =  %12.4f  1-4 EEL =  %12.4f  RESTRAINT  =  %12.4f
     *
     * @param crds
     * @return
     */
    std::string ereport(const std::valarray<mmpbsa::Vector>& crds);

    void extract_lj_params(std::vector<mmpbsa::lj_params_t>& lj_params)const;

    void extract_force_field(mmpbsa::forcefield_t& ff)const;

    /**
     * Returns a valarray containing the residue ranges. The array is of the form:
     * [(min,max),(min,max),(min,max),...]
     *
     * @param resptr List of pointers to residues.
     * @return list of ranges for residue pointers.
     */
    static std::valarray<size_t> get_res_ranges(const std::valarray<size_t>& resptr,
            const size_t& natoms);

    /**
     * Returns a valarray containing the molecule ranges. The array is of the form:
     * [(min,max),(min,max),(min,max),...]
     *
     * @param molptrs
     * @return
     */
    static std::valarray<size_t> get_mol_ranges(const std::valarray<size_t>& molptrs);


    //operators
    EmpEnerFun& operator=(const EmpEnerFun& orig);
    

    //data
    mmpbsa::SanderParm * parminfo;//do not delete via EmpEnerFun
    
    int end_solute_atoms;
    int begin_solvent_atoms;
    mmpbsa_t inv_scnb;
    mmpbsa_t inv_scee;
    mmpbsa_t dielc;
    std::slice bond_i;
    std::slice bond_j;
    std::slice bond_k;
    std::slice bond_h_i;
    std::slice bond_h_j;
    std::slice bond_h_k;
    std::slice angle_i;
    std::slice angle_j;
    std::slice angle_k;
    std::slice angle_l;
    std::slice angle_h_i;
    std::slice angle_h_j;
    std::slice angle_h_k;
    std::slice angle_h_l;
    std::slice phi_i;
    std::slice phi_j;
    std::slice phi_k;
    std::slice phi_l;
    std::slice phi_m;
    std::slice phi_h_i;
    std::slice phi_h_j;
    std::slice phi_h_k;
    std::slice phi_h_l;
    std::slice phi_h_m;

    std::vector<std::vector<size_t> > exclst;///<List of Excluded atoms

    std::valarray<mmpbsa_t> LJA;///<Lennard Jones Coefficients
    std::valarray<mmpbsa_t> LJB;///<Lennard Jones Coefficients

    std::valarray<size_t> res_ranges;///<Beginning and Ending indices of residues
    std::valarray<size_t> mol_ranges;///<Beginning and Ending indices of molecules


    
    
private:
    /**
     * Function used in bond walking that increments markers.
     *
     * @param markers
     * @param i
     */
    static void visitor(std::valarray<int>& markers,const size_t& i){markers[i]++;}

};

class BondWalker
{
public:
    /**
     * Bond list creator.
     * When Bond lists are not provided by the parmtop file, BondWalker produces
     * the bond list.
     *
     * @param efun
     */
    BondWalker(EmpEnerFun const * efun);

    void init();

    /**
     * Generate bond list for the specified atom.
     * 
     * @param startAtom
     * @param stoppers
     * @param enerfunMarkers
     * @param func
     * @return
     */
    std::vector<size_t> walk(const size_t& startAtom,
        std::valarray<int> stoppers,std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom));

    const EmpEnerFun * enerfun;//energy function containing energy values
    
private:
    void recwalk(std::vector<size_t>& list, const size_t& atom,
        std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom));

    bool initialized;
    std::valarray<bool> visited;
    
    /**
     * This is a list of lists. Each sublist corresponds to an atom. The contents
     * of that atom's list provide the atoms to which it is bonded.
     */
    std::vector<std::vector<int> > atom_bond_list;
};

};//end namespace mmpbsa

/**
 * Writes Energy data to PDB format (cf http://www.wwpdb.org/documentation/format32/sect9.html)
 */
std::ostream& streamPDB(std::ostream& theStream, const mmpbsa::EmpEnerFun& energy, const std::valarray<mmpbsa::Vector>& crds) throw (mmpbsa::MMPBSAException);

#endif	/* ENERGY_H */


