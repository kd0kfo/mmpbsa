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




class EmpEnerFun{
public:
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
    EmpEnerFun(sanderio::SanderParm * newparminfo, const mmpbsa_t& scn,
            const mmpbsa_t& scee, const mmpbsa_t& dielc);

    /**
     * Copy Constructor
     * 
     * @param orig
     */
    EmpEnerFun(const EmpEnerFun& orig);

    ~EmpEnerFun(){//do not delete parminfo. It is externally made and should be deleted outside of EmpEnerInfo
    }

    sanderio::SanderParm * parminfo;//do not delete via EmpEnerFun
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



    static void visitor(std::valarray<int>& markers,const size_t& i){markers[i]++;}
    
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

