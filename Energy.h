/*
 * Energy
 *
 * Calculates Energies from Sander parmtop files.
 *
 * Created by David Coss
 */

#ifndef ENERGY_H
#define	ENERGY_H

#include <valarray>
#include <vector>

#include "SanderIO.h"
#include "mmpbsa_utils.h"


class EmpEnerFun{
public:
    /**
     * Main Energy Constructor
     */
    EmpEnerFun(const sanderio::SanderParm& newparminfo, const mmpbsa_t& scnb=2.0,
            const mmpbsa_t& scee = 1.2, const mmpbsa_t& dielc = 1.0);


private:
    /**
     * Returns a valarray containing the residue ranges. The array is of the form:
     * [(min,max),(min,max),(min,max),...]
     * 
     * @param resptr List of pointers to residues.
     * @return list of ranges for residue pointers.
     */
    static std::valarray<int> get_res_ranges(const std::valarray<int>& resptr);

    static std::valarray<int> get_mol_ranges(const std::valarray<int>& molptrs);

    sanderio::SanderParm parminfo;///<Contains Paramters from the Sander prmtop file.

    static void visitor(int& i){return i++;}
    
};

class BondWalker
{
public:
    BondWalker();

private:
    bool initiated;
    std::valarray<bool> visited;
    
};


#endif	/* ENERGY_H */

