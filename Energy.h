/*
 * Energy
 *
 * Calculates Energies from Sander parmtop files.
 *
 * Created by David Coss
 */
#ifndef ENERGY_H
#define	ENERGY_H

#include <cmath>
#include <valarray>
#include <vector>
#include <fstream>
extern std::fstream myOutput;

#include "mmpbsa_utils.h"
#include "mmpbsa_io.h"



#define DEFAULT_SCNB 2.0
#define DEFAULT_SCEE 1.2
#define DEFAULT_DIELC 1.0

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
    EmpEnerFun(mmpbsa_io::SanderParm * newparminfo, const mmpbsa_t& scnb = DEFAULT_SCNB,
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
    mmpbsa_t total_bond_energy(const std::valarray<mmpbsa_t>& crds)const{return bond_inc_H(crds)+bond_without_H(crds);}
    mmpbsa_t bond_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t bond_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t bond_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& bondIndices)const;

    /**
     * Calculates the total angle energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_angle_energy(const std::valarray<mmpbsa_t>& crds)const{return angle_inc_H(crds)+angle_without_H(crds);}
    mmpbsa_t angle_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t angle_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t angle_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& angleIndices)const;

    /**
     * Calculates the total dihedral energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_dihedral_energy(const std::valarray<mmpbsa_t>& crds)const{return dihedral_inc_H(crds)+dihedral_without_H(crds);}
    mmpbsa_t dihedral_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t dihedral_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t dihedral_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices)const;

    /**
     * Calculates the total Van der Waals energy between 1-4 pairs for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_vdw14_energy(const std::valarray<mmpbsa_t>& crds)const{return vdw14_inc_H(crds)+vdw14_without_H(crds);}
    mmpbsa_t vdw14_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t vdw14_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t vdw14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices,const std::valarray<bool>& phi_mask)const;

    /**
     * Calculates the total Electrostatic energy between 1-4 pairs for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_elstat14_energy(const std::valarray<mmpbsa_t>& crds)const{return elstat14_inc_H(crds)+elstat14_without_H(crds);}
    mmpbsa_t elstat14_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t elstat14_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t elstat14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices, const std::valarray<bool>& phi_mask)const;

    /**
     * Calculates the total Van der Waals energy of the system for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_vdwaals_energy(const std::valarray<mmpbsa_t>& crds)const;

    /**
     * Calculates the total Electrostatic energy for the given snapshot coordinates
     *
     * @param crds
     * @return
     */
    mmpbsa_t total_elstat_energy(const std::valarray<mmpbsa_t>& crds)const;

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
    std::string ereport(const std::valarray<mmpbsa_t>& crds);


    //operators
    EmpEnerFun& operator=(const EmpEnerFun& orig);
    

    //data
    mmpbsa_io::SanderParm * parminfo;//do not delete via EmpEnerFun
    
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

    std::vector<std::vector<size_t> > * exclst;

    std::valarray<mmpbsa_t> * LJA;
    std::valarray<mmpbsa_t> * LJB;

    std::valarray<size_t> * res_ranges;
    std::valarray<size_t> * mol_ranges;


    
    
private:
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

    /**
     * Function used in bond walking that increments markers.
     *
     * @param markers
     * @param i
     */
    static void visitor(std::valarray<int>& markers,const size_t& i){markers[i]++;}

    void freeMemory();
};

class EMap{
public:
    /**
     * Default Constructor for Energy Map Object
     * @see EMap(const EmpEnerFun* efun,const std::valarray<mmpbsa_t>& crds)
     *
     */
    EMap();

    /**
     * Copy constructor
     *
     * @param orig
     */
    EMap(const EMap& orig);

    /**
     * Creates an Energy Map object for the given Energy Function.
     * Energy maps are used to provide arithmetic operations to Energy data.
     * Operator overloads provided in the EMap class allow one to do arithmetic
     * on correspoding energy values in the Energy object.
     * 
     * @param efun
     * @param crds
     */
    EMap(const EmpEnerFun* efun,const std::valarray<mmpbsa_t>& crds);

    ~EMap(){}

    const mmpbsa_t& set_elstat_solv(const mmpbsa_t& val){return elstat_solv = val;}
    const mmpbsa_t& set_area(const mmpbsa_t& val){return area = val;}
    const mmpbsa_t& set_sasol(const mmpbsa_t& val){return sasol = val;}
    /**
     * ostream output for convienent writing of energy data.
     *
     * Format:
     * BOND
     * ANGLE
     * DIHED
     * VDW14
     * ELE14
     * VDWAALS
     * VACELE
     * 
     * @param rhs
     * @return
     */
    friend std::ostream& operator<<(std::ostream& theStream, const EMap& toWrite);

    //operators
    EMap& operator=(const EMap& rhs);
    EMap operator+(const EMap& rhs);
    EMap& operator+=(const EMap& rhs);
    EMap operator-(const EMap& rhs);
    EMap& operator-=(const EMap& rhs);


private:
    mmpbsa_t bond;
    mmpbsa_t angle;
    mmpbsa_t dihed;
    mmpbsa_t vdw14;
    mmpbsa_t ele14;
    mmpbsa_t vdwaals;
    mmpbsa_t vacele;
    mmpbsa_t elstat_solv;
    mmpbsa_t area;
    mmpbsa_t sasol;
    
};

class BondWalker
{
public:
    BondWalker(EmpEnerFun const * efun);

    void init();

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


#endif	/* ENERGY_H */

