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
    EmpEnerFun(sanderio::SanderParm * newparminfo);
    
    /**
     * Main Energy Constructor
     */
    EmpEnerFun(sanderio::SanderParm * newparminfo, const mmpbsa_t& scnb,
            const mmpbsa_t& scee, const mmpbsa_t& dielc);

    /**
     * Copy Constructor
     * 
     * @param orig
     */
    EmpEnerFun(const EmpEnerFun& orig);

    ~EmpEnerFun(){//do not delete parminfo. It is externally made and should be deleted outside of EmpEnerInfo
    }

    EmpEnerFun stripEnerFun(const std::valarray<bool>& keepers,
        const bool& dangleWarn) const;

    template <class M, class N>
    static void internalConvert(std::valarray<std::valarray<M> >& newIndices,std::valarray<std::valarray<N> >& newTerms,
        const std::valarray<M>& oldIndices,const std::valarray<std::slice>& slices,
        const std::valarray<const std::valarray<N>* >& oldTerms,
        const std::valarray<size_t>& newidx,
        const std::valarray<bool> keepers, const bool& dangleWarn = false);

    sanderio::SanderParm * parminfo;//do not delete via EmpEnerFun
    const std::valarray<std::string> * resnames;//do not delete. Can you?

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

    std::vector<std::vector<size_t> > exclst;

    std::valarray<size_t> res_ranges;
    std::valarray<size_t> mol_ranges;

    std::valarray<mmpbsa_t> bond_h_const;
    std::valarray<mmpbsa_t> bond_h_eq;
    std::valarray<mmpbsa_t> bond_const;
    std::valarray<mmpbsa_t> bond_eq;
    std::valarray<mmpbsa_t> angle_const;
    std::valarray<mmpbsa_t> angle_eq;
    std::valarray<mmpbsa_t> angle_h_const;
    std::valarray<mmpbsa_t> angle_h_eq;
    std::valarray<mmpbsa_t> phi_const;
    std::valarray<mmpbsa_t> phi_periodicity;
    std::valarray<mmpbsa_t> phi_phase;
    std::valarray<mmpbsa_t> phi_cos_phase;
    std::valarray<mmpbsa_t> phi_sin_phase;
    std::valarray<mmpbsa_t> phi_h_const;
    std::valarray<mmpbsa_t> phi_h_periodicity;
    std::valarray<mmpbsa_t> phi_h_phase;
    //really want these?
    std::valarray<mmpbsa_t> phi_h_cos_phase;
    std::valarray<mmpbsa_t> phi_h_sin_phase;

    std::valarray<bool> phi_prd_mask;
    std::valarray<bool> phi_ignend_mask;
    std::valarray<bool> phi_imp_mask;
    std::valarray<bool> phi_mask;
    std::valarray<bool> phi_h_prd_mask;
    std::valarray<bool> phi_h_ignend_mask;
    std::valarray<bool> phi_h_imp_mask;
    std::valarray<bool> phi_h_mask;

    
    
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

