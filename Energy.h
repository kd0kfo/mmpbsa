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

#include "SanderIO.h"
#include "mmpbsa_utils.cpp"



#define DEFAULT_SCNB 2.0
#define DEFAULT_SCEE 1.2
#define DEFAULT_DIELC 1.0

class EmpEnerFun{
public:
    /**
     * Default constructor. Uses default values of scnb, scee and dielc.
     * Initializes the valarrays (of size 0).
     * parminfo and resnames are set to null.
     * 
     * @param newparminfo
     */
    EmpEnerFun();


    /**
     * Creates an Energy Function class with the following defaults:
     * scnb=2.0
     * scee = 1.2
     * dielc = 1.0
     * 
     * @param newparminfo
     */
    //EmpEnerFun(sanderio::SanderParm * newparminfo);
    
    /**
     * Main Energy Constructor
     */
    EmpEnerFun(sanderio::SanderParm * newparminfo, const mmpbsa_t& scnb = DEFAULT_SCNB,
            const mmpbsa_t& scee = DEFAULT_SCEE, const mmpbsa_t& dielc = DEFAULT_DIELC);

    /**
     * Copy Constructor
     * 
     * @param orig
     */
    EmpEnerFun(const EmpEnerFun& orig);

    ~EmpEnerFun();//do not delete parminfo. It is externally made and should be deleted outside of EmpEnerInfo
    

    EmpEnerFun stripEnerFun(const std::valarray<bool>& keepers,
        const bool& dangleWarn) const;

    template <class M>
    static void internalConvert(std::valarray<M>& newIndices,
        const std::valarray<M>& oldIndices,const std::valarray<std::slice>& slices,
        const std::valarray<size_t>& newidx,
        const std::valarray<bool>& keepers, const bool& dangleWarn = false);

    template <class M> static void updatePhiMasks(std::valarray<bool>& newMask,
    const std::valarray<bool>& oldMask,
    const std::valarray<M>& oldIndices, const std::valarray<std::slice>& slices,
        const std::valarray<bool>& keepers, const bool& dangleWarn);


    //Energy Calculations

    mmpbsa_t total_bond_energy(const std::valarray<mmpbsa_t>& crds)const{return bond_inc_H(crds)+bond_without_H(crds);}
    mmpbsa_t bond_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t bond_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t bond_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& bondIndices)const;

    mmpbsa_t total_angle_energy(const std::valarray<mmpbsa_t>& crds)const{return angle_inc_H(crds)+angle_without_H(crds);}
    mmpbsa_t angle_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t angle_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t angle_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& angleIndices)const;

    mmpbsa_t total_dihedral_energy(const std::valarray<mmpbsa_t>& crds)const{return dihedral_inc_H(crds)+dihedral_without_H(crds);}
    mmpbsa_t dihedral_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t dihedral_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t dihedral_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices)const;

    mmpbsa_t total_vdw14_energy(const std::valarray<mmpbsa_t>& crds)const{return vdw14_inc_H(crds)+vdw14_without_H(crds);}
    mmpbsa_t vdw14_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t vdw14_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t vdw14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices,const std::valarray<bool>& phi_mask)const;

    mmpbsa_t total_elstat14_energy(const std::valarray<mmpbsa_t>& crds)const{return elstat14_inc_H(crds)+elstat14_without_H(crds);}
    mmpbsa_t elstat14_inc_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t elstat14_without_H(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t elstat14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices, const std::valarray<bool>& phi_mask)const;

    mmpbsa_t total_vdwaals_energy(const std::valarray<mmpbsa_t>& crds)const;
    mmpbsa_t total_elstat_energy(const std::valarray<mmpbsa_t>& crds)const;
    
    std::string ereport(const std::valarray<mmpbsa_t>& crds);


    //operators
    EmpEnerFun& operator=(const EmpEnerFun& orig);
    

    //data
    sanderio::SanderParm * parminfo;//do not delete via EmpEnerFun
    
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

//    std::valarray<mmpbsa_t> * bond_h_const;
//    std::valarray<mmpbsa_t> * bond_h_eq;
//    std::valarray<mmpbsa_t> * bond_const;
//    std::valarray<mmpbsa_t> * bond_eq;
//    std::valarray<mmpbsa_t> * angle_const;
//    std::valarray<mmpbsa_t> * angle_eq;
//    std::valarray<mmpbsa_t> * angle_h_const;
//    std::valarray<mmpbsa_t> * angle_h_eq;
//    std::valarray<mmpbsa_t> * phi_const;
//    std::valarray<mmpbsa_t> * phi_periodicity;
//    std::valarray<mmpbsa_t> * phi_phase;
//    std::valarray<mmpbsa_t> * phi_cos_phase;
//    std::valarray<mmpbsa_t> * phi_sin_phase;
//    std::valarray<mmpbsa_t> * phi_h_const;
//    std::valarray<mmpbsa_t> * phi_h_periodicity;
//    std::valarray<mmpbsa_t> * phi_h_phase;
//    //really want these?
//    std::valarray<mmpbsa_t> * phi_h_cos_phase;
//    std::valarray<mmpbsa_t> * phi_h_sin_phase;

//    std::valarray<bool> * phi_prd_mask;
//    std::valarray<bool> * phi_ignend_mask;
//    std::valarray<bool> * phi_imp_mask;
//    std::valarray<bool> * phi_mask;
//    std::valarray<bool> * phi_h_prd_mask;
//    std::valarray<bool> * phi_h_ignend_mask;
//    std::valarray<bool> * phi_h_imp_mask;
//    std::valarray<bool> * phi_h_mask;

    
    
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

    static std::valarray<size_t> get_mol_ranges(const std::valarray<size_t>& molptrs);

    static void visitor(std::valarray<int>& markers,const size_t& i){markers[i]++;}

    void freeMemory();
};

class EMap{
public:
    EMap();
    EMap(const EMap& orig);
    EMap(const EmpEnerFun* efun,const std::valarray<mmpbsa_t>& crds);
    ~EMap(){}

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

