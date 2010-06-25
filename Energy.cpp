#include "Energy.h"

EmpEnerFun::EmpEnerFun(sanderio::SanderParm * newparminfo, const mmpbsa_t& scnb,
        const mmpbsa_t& scee, const mmpbsa_t& dielc)
{
    using namespace mmpbsa_utils;
    using namespace std;

    parminfo = newparminfo;
    mmpbsa_t inv_scnb = 1.0 / scnb;
    mmpbsa_t inv_scee = 1.0 / scee;
    this->dielc = dielc;
    const int& ntypes = parminfo->ntypes;
    const int& natom = parminfo->natom;
    valarray<size_t> resptr = parminfo->residue_pointers - size_t(1); //sander file pointers are 1-indexed
    valarray<size_t> res_ranges = get_res_ranges(resptr,natom); //getranges should be of the type (min,max),(min,max),...


    //Non-bonded stuff
    const valarray<size_t>& nb_parm_idx = parminfo->nonbonded_parm_indices;
    const valarray<mmpbsa_t>& CA = parminfo->lennard_jones_acoefs;
    const valarray<mmpbsa_t>& CB = parminfo->lennard_jones_bcoefs;
    const valarray<mmpbsa_t>& CHA = parminfo->hbond_acoefs;
    const valarray<mmpbsa_t>& CHB = parminfo->hbond_bcoefs;
    valarray<mmpbsa_t> LJA(ntypes * ntypes);
    valarray<mmpbsa_t> LJB(ntypes * ntypes);
    valarray<mmpbsa_t> LJHA(ntypes * ntypes);
    valarray<mmpbsa_t> LJHB(ntypes, ntypes);

    for (size_t i = 0; i < ntypes; i++) {
        for (size_t j = i; j < ntypes; j++) {
            size_t ico = nb_parm_idx[i+j*ntypes];
            bool isTenTwelvePair = parminfo->nonbonded_parm_mask[i+j*ntypes];
            if (!isTenTwelvePair) {
                LJA[i+j*ntypes] = CA[ico - 1];
                LJB[i+j*ntypes] = CB[ico - 1];
            } else {
                LJHA[i+j*ntypes] = CHA[ico - 1];
                LJHB[i+j*ntypes] = CHB[ico - 1];
            }
        }
    }

    //# Fill in the other triangle
    for (size_t i = 0; i < ntypes; i++) {
        for (size_t j = i + 1; j < ntypes; j++) {
            LJA[j+i*ntypes] = LJA[i+j*ntypes];
            LJB[j+i*ntypes] = LJB[i+j*ntypes];
            LJHA[j+i*ntypes] = LJHA[i+j*ntypes];
            LJHB[j+i*ntypes] = LJHB[i+j*ntypes];
        }
    }

    //# excluded atom stuff
    const valarray<size_t>& numex = parminfo->number_excluded_atoms;
    valarray<size_t> natex = parminfo->excluded_atoms_list - 1;//1-indexed
    vector<size_t> exclst;
    size_t num = 0;
    for (size_t i = 0; i < natom; i++) {
        //# This compress is because a single zero element seems to
        //# be the parmtop way of indicating an empty list
        vector<size_t> irow = compress_ge(natex,slice(num, num + numex[i], 1), size_t(0));

        //# shift indices to refer to list of atoms AFTER i:
        for (size_t j = 0; j < irow.size(); j++) {
            irow[j] -= i + 1;
            exclst.push_back(irow[j]);
        }
        num += numex[i];
    }

    //temporary loading variables
    size_t columnWidth;

    //# bond stuff
    valarray<size_t>& bondcodes = parminfo->bonds_inc_hydrogen;
    columnWidth = 3;
    //divide for (i,j,k) bondcodes, divide i and j by 3 and decrement k
    bond_h_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    bond_h_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    bond_h_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    
    valarray<mmpbsa_t> bond_h_const = parminfo->bond_force_constants[bondcodes[bond_h_k]];
    valarray<mmpbsa_t> bond_h_eq = parminfo->bond_equil_values[bondcodes[bond_h_k]];

    bondcodes = parminfo->bonds_without_hydrogen;
    bond_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    bond_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    bond_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    valarray<mmpbsa_t> bond_const = parminfo->bond_force_constants[bondcodes[bond_k]];
    valarray<mmpbsa_t> bond_eq = parminfo->bond_equil_values[bondcodes[bond_k]];

    //# angle stuff
    bondcodes = parminfo->angles_without_hydrogen;
    columnWidth = 4;
    angle_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    angle_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    angle_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    angle_l = slice(3, bondcodes.size() / columnWidth, columnWidth);
    valarray<mmpbsa_t> angle_const = parminfo->angle_force_constants[bondcodes[angle_l]];
    valarray<mmpbsa_t> angle_eq = parminfo->angle_equil_values[bondcodes[angle_l]];

    bondcodes = parminfo->angles_inc_hydrogen;
    angle_h_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    angle_h_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    angle_h_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    angle_h_l = slice(3, bondcodes.size() / columnWidth, columnWidth);
    valarray<mmpbsa_t> angle_h_const = parminfo->angle_force_constants[bondcodes[angle_h_l]];
    valarray<mmpbsa_t> angle_h_eq = parminfo->angle_equil_values[bondcodes[angle_h_l]];

    //# dihedral stuff
    bondcodes = parminfo->dihedrals_without_hydrogen;
    columnWidth = 5;
    phi_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    phi_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    phi_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    phi_l = slice(3, bondcodes.size() / columnWidth, columnWidth);
    phi_m = slice(4, bondcodes.size() / columnWidth, columnWidth);
    valarray<bool> phi_ignend_mask = parminfo->dihedral_mask[phi_k];
    valarray<bool> phi_imp_mask = parminfo->dihedral_mask[phi_l];
    valarray<mmpbsa_t> phi_const = parminfo->dihedral_force_constants[bondcodes[phi_m]];
    valarray<mmpbsa_t> phi_periodicity = parminfo->dihedral_periodicities[bondcodes[phi_m]];
    valarray<bool> phi_prd_mask = phi_periodicity < 0.0;
    phi_periodicity = abs(phi_periodicity);

    valarray<bool> phi_mask = mmpbsa_utils::logical_not(phi_imp_mask) && mmpbsa_utils::logical_not(phi_ignend_mask)
            && mmpbsa_utils::logical_not(phi_prd_mask);

    valarray<mmpbsa_t> phi_phase = parminfo->dihedral_phases[parminfo->dihedrals_without_hydrogen[phi_m]];
    //really want these???
    valarray<mmpbsa_t> phi_cos_phase = cos(phi_phase);
    valarray<mmpbsa_t> phi_sin_phase = sin(phi_phase);

    bondcodes = parminfo->dihedrals_inc_hydrogen;
    phi_h_i = slice(0, bondcodes.size() / columnWidth, columnWidth);
    phi_h_j = slice(1, bondcodes.size() / columnWidth, columnWidth);
    phi_h_k = slice(2, bondcodes.size() / columnWidth, columnWidth);
    phi_h_l = slice(3, bondcodes.size() / columnWidth, columnWidth);
    phi_h_m = slice(4, bondcodes.size() / columnWidth, columnWidth);
    valarray<bool> phi_h_ignend_mask = parminfo->dihedral_h_mask[phi_h_k];
    valarray<bool> phi_h_imp_mask = parminfo->dihedral_h_mask[phi_h_l];
    valarray<mmpbsa_t> phi_h_const = parminfo->dihedral_force_constants[bondcodes[phi_h_m]];
    valarray<mmpbsa_t> phi_h_periodicity = parminfo->dihedral_periodicities[bondcodes[phi_h_m]];
    valarray<mmpbsa_t> phi_h_phase = parminfo->dihedral_phases[bondcodes[phi_h_m]];
    //really want these?
    valarray<mmpbsa_t> phi_h_cos_phase = cos(phi_h_phase);
    valarray<mmpbsa_t> phi_h_sin_phase = sin(phi_h_phase);

    valarray<bool> phi_h_prd_mask = phi_h_periodicity < 0.0;
    phi_h_periodicity = abs(phi_h_periodicity);

    valarray<bool> phi_h_mask = mmpbsa_utils::logical_not(phi_h_imp_mask)
            && mmpbsa_utils::logical_not(phi_h_ignend_mask) && mmpbsa_utils::logical_not(phi_h_prd_mask);

    valarray<size_t> mol_ranges;
    if (parminfo->ifbox) //# We have solvent pointer info
    {
        valarray<size_t> molptrs = cumsum(parminfo->atoms_per_molecule);
        if (molptrs[molptrs.size() - 1] != natom)
            throw "The last column sum of atoms per molecule should equal "
            "the number of atoms";

        valarray<size_t> mol_ranges = get_mol_ranges(molptrs); //getranges should be of the type (min,max),(min,max),...
        if (mol_ranges.size() != parminfo->nspm)
            throw "The number of ranges must match the total number of molecules";

        int end_solute_atoms = res_ranges[1+(parminfo->iptres-1)*res_ranges.size()];
        int begin_solvent_atoms = mol_ranges[mol_ranges.size()*(parminfo->nspsol-1)];
        //# TODO Make sure stripEnerFun fixes the above begin and end ptrs!
    }
    else//# Darn, no ATOMS_PER_MOLECULE INFO!  Walk the bond list
    {
        BondWalker bw(this);
        vector<size_t> beg_ptrs;
        valarray<int> markers(natom);
        while(markers.min() == 0)
        {
            size_t first_unmarked = find_first(markers,0);
            throw "add assert";
            beg_ptrs.push_back(first_unmarked);
            bw.walk(first_unmarked,valarray<int>(),markers,visitor);
        }
        valarray<size_t> mol_max = mmpbsa_utils::cshift(beg_ptrs,1);
        mol_max[mol_max.size()-1] = natom;
        mol_ranges = zip(beg_ptrs,mol_max);
    }

}/*end of constructor*/

EmpEnerFun::EmpEnerFun(sanderio::SanderParm * newparminfo)
{
    EmpEnerFun tmp(newparminfo, 2.0, 1.2,1.0);
    *this = tmp;
}/**/

EmpEnerFun::EmpEnerFun(const EmpEnerFun& orig)
{

}

std::valarray<size_t> EmpEnerFun::get_res_ranges(const std::valarray<size_t>& resptr,
        const size_t& natoms) {
    std::valarray<size_t> res_ranges(2 * resptr.size());

    for (int i = 0; i < resptr.size() - 2; i += 2) {
        res_ranges[i] = resptr[i];
        res_ranges[i + 1] = resptr[i + 1];
    }
    res_ranges[resptr.size() - 2] = resptr[resptr.size() - 1];
    res_ranges[resptr.size() - 1] = natoms;

    return res_ranges;
}

std::valarray<size_t> EmpEnerFun::get_mol_ranges(const std::valarray<size_t>& molptrs) {
    //[0] + molptrs[:-1], molptrs
    std::valarray<size_t> mol_ranges(2 * molptrs.size());

    mol_ranges[0] = molptrs[molptrs.size() - 1];
    mol_ranges[1] = molptrs[0];

    for (int i = 1; i < molptrs.size(); i += 2) {
        mol_ranges[i + 1] = molptrs[i - 1];
        mol_ranges[i + 2] = molptrs[i];
    }

    return mol_ranges;
}

BondWalker::BondWalker(EmpEnerFun const * efun)
{
    enerfun = efun;
    initialized = false;
    visited;
    atom_bond_list;
}

void BondWalker::init()
{
    using std::valarray;
    using std::slice;
    using std::slice_array;

    sanderio::SanderParm const * const parminfo = enerfun->parminfo;
    atom_bond_list.resize(parminfo->natom);
    size_t numBondPairs = size_t( parminfo->bonds_without_hydrogen.size()/3
        + parminfo->bonds_inc_hydrogen.size()/3 );//divide by 3 because each bond has an i,j pair plus force info

    //multiple by 2 because this contains pairs of atoms in each bond.
    //This array will have the form (bond, bond, ...) or ((i,j), (i,j), ...)
    valarray<size_t> bondPairs(2*numBondPairs);
    bondPairs[slice(0,numBondPairs,1)] = mmpbsa_utils::zip<size_t>(parminfo->bonds_inc_hydrogen,enerfun->bond_i,
            parminfo->bonds_inc_hydrogen,enerfun->bond_j);
    bondPairs[slice(numBondPairs,2*numBondPairs,1)] = mmpbsa_utils::zip(parminfo->bonds_without_hydrogen,enerfun->bond_h_i,
            parminfo->bonds_without_hydrogen,enerfun->bond_h_j);

    if(bondPairs.max() > parminfo->natom)
        throw MMPBSAException("An atom index in the bonds pairs exceed the "
                "maximum value.",DATA_FORMAT_ERROR);

    size_t atom_i,atom_j;
    for(size_t i = 0;i<numBondPairs;i++)
    {
        atom_i = bondPairs[2*i];
        atom_j = bondPairs[2*i+1];
        atom_bond_list[atom_i].push_back(atom_j);
        atom_bond_list[atom_j].push_back(atom_i);
    }

    visited.resize(parminfo->natom,false);
    initialized = true;

}

std::vector<size_t> BondWalker::walk(const size_t& startAtom,
        std::valarray<int> stoppers,std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom))
{
    if(initialized)
        visited != visited;
    else
        init();

    if(stoppers.size() > 0)
        stoppers = 1;

    std::vector<size_t> returnMe;
    recwalk(returnMe,startAtom,enerfunMarkers,func);

    return returnMe;
}

void BondWalker::recwalk(std::vector<size_t>& list, const size_t& atom,
        std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom))
{
    visited[atom] = true;
    for(size_t i = 0;i<atom_bond_list[atom].size();i++)
    {
        size_t bondedAtom = atom_bond_list[atom][i];
        if(visited[bondedAtom])
            continue;
        recwalk(list,bondedAtom,enerfunMarkers,func);
    }

}
