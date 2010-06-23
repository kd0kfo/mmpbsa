#include "Energy.h"

EmpEnerFun::EmpEnerFun(const sanderio::SanderParm& newparminfo, const mmpbsa_t& scnb,
            const mmpbsa_t& scee, const mmpbsa_t& dielc)
{
    using namespace mmpbsa_utils;
    using namespace std;

    parminfo = newparminfo;
        mmpbsa_t inv_scnb = 1.0/scnb;
        mmpbsa_t inv_scee = 1.0/scee;
        mmpbsa_t dielc = dielc;
        const int& ntypes = parminfo.ntypes;
        const int& natom = parminfo.natom;
        valarray<int> resptr = parminfo.residue_pointers - 1;//sander file pointers are 1-indexed
        valarray<int> res_ranges = get_res_ranges(resptr);//getranges should be of the type (min,max),(min,max),...


        //Non-bonded stuff
        const valarray<int>& nb_parm_idx = parminfo.nonbonded_parm_indices;
        const valarray<mmpbsa_t>& CA = parminfo.lennard_jones_acoefs;
        const valarray<mmpbsa_t>& CB = parminfo.lennard_jones_bcoefs;
        const valarray<mmpbsa_t>& CHA = parminfo.hbond_acoefs;
        const valarray<mmpbsa_t>& CHB = parminfo.hbond_bcoefs;
        valarray<mmpbsa_t> LJA(ntypes*ntypes);
        valarray<mmpbsa_t> LJB(ntypes*ntypes);
        valarray<mmpbsa_t> LJHA(ntypes*ntypes);
        valarray<mmpbsa_t> LJHB(ntypes, ntypes);

        for(int i = 0;i<ntypes;i++)
        {
            for(int j = i;j<ntypes;j++)
            {
                int ico = nb_parm_idx[std::slice(j,ntypes,ntypes)][i];
                if(ico > 0)
                {
                    LJA[std::slice(j,ntypes,ntypes)][i] = CA[ico-1];
                    LJB[std::slice(j,ntypes,ntypes)][i] = CB[ico-1];
                }
                else
                {
                    LJHA[std::slice(j,ntypes,ntypes)][i] = CHA[-ico - 1];
                    LJHB[std::slice(j,ntypes,ntypes)][i] = CHB[-ico - 1];
                }
            }
        }

        //# Fill in the other triangle
        for(int i= 0;i<ntypes;i++)
        {
            for(int j = i+1;j<ntypes;j++)
            {
                LJA[std::slice(i,ntypes,ntypes)][j] = LJA[std::slice(j,ntypes,ntypes)][i];
                LJB[std::slice(i,ntypes,ntypes)][j] = LJB[std::slice(j,ntypes,ntypes)][i];
                LJHA[std::slice(i,ntypes,ntypes)][j] = LJHA[std::slice(j,ntypes,ntypes)][i];
                LJHB[std::slice(i,ntypes,ntypes)][j] = LJHB[std::slice(j,ntypes,ntypes)][i];
            }
        }

        //# excluded atom stuff
        const valarray<int>& numex = parminfo.number_excluded_atoms;
        valarray<int>& natex = parminfo.excluded_atoms_list-1;//1-indexed
        vector<int> exclst;
        int num = 0;
        for(int i = 0;i<natom;i++)
        {
            std::slice_array<int> pre_irow = natex[slice(num,num+numex[i],1)];
            //# This compress is because a single zero element seems to
            //# be the parmtop way of indicating an empty list
            vector<int> irow = compress_ge(pre_irow, 0);

            //# shift indices to refer to list of atoms AFTER i:
            for(int j=0;j<irow.size();j++)
            {
                irow[j] -= i+1;
                exclst.push_back(irow[j]);
            }
            num += numex[i];
        }

        //temporary loading variables
        const valarray<int>& bondcodes;
        int columnWidth;
        
        //# bond stuff
        bondcodes = parminfo.bonds_inc_hydrogen;
        columnWidth = 3;
        //pyamber: bondcodes.shape = (-1,3)

        //divide for (i,j,k) bondcodes, divide i and j by 3 and decrement k
        slice_array<int> bond_h_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> bond_h_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> bond_h_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        bond_h_i /= 3;
        bond_h_j /= 3;
        bond_h_k -= 1;

        valarray<mmpbsa_t> bond_h_const = take(parminfo.bond_force_constants, bond_h_k);
        valarray<mmpbsa_t> bond_h_eq = take(parminfo.bond_equil_values, bond_h_k);

        bondcodes = parminfo.bonds_without_hydrogen;
        slice_array<int> bond_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> bond_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> bond_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        bond_i /= 3;
        bond_j /= 3;
        bond_k -= 1;
        valarray<mmpbsa_t> bond_const = take(parminfo.bond_force_constants, bond_h_k);
        valarray<mmpbsa_t> bond_eq = take(parminfo.bond_equil_values, bond_h_k);

        //# angle stuff
        bondcodes = parminfo.angles_without_hydrogen;
        columnWidth = 4;
        slice_array<int> angle_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_l = bondcodes[slice(3,bondcodes.size()/columnWidth,columnWidth)];
        angle_i /= 3;
        angle_j /= 3;
        angle_k /= 3;
        angle_l -= 1;
        valarray<mmpbsa_t> angle_const = take(parminfo.angle_force_constants, angle_l);
        valarray<mmpbsa_t> angle_eq = take(parminfo.angle_equil_values, angle_l);

        bondcodes = parminfo.angles_inc_hydrogen;
        slice_array<int> angle_h_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_h_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_h_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> angle_h_l = bondcodes[slice(3,bondcodes.size()/columnWidth,columnWidth)];
        angle_h_i /= 3;
        angle_h_j /= 3;
        angle_h_k /= 3;
        angle_h_l -= 1;
        valarray<mmpbsa_t> angle_h_const = take(parminfo.angle_force_constants, angle_h_l);
        valarray<mmpbsa_t> angle_h_eq = take(parminfo.angle_equil_values, angle_h_l)

        //# dihedral stuff
        bondcodes = parminfo.dihedrals_without_hydrogen;
        columnWidth = 5;
        slice_array<int> phi_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_l = bondcodes[slice(3,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_m = bondcodes[slice(4,bondcodes.size()/columnWidth,columnWidth)];
        phi_i /= 3;
        phi_j /= 3;
        phi_k = abs(phi_k)/3;
        phi_l = abs(phi_l)/3;
        phi_m -= 1;
        valarray<bool> phi_imp_mask = phi_l < 0;
        valarray<bool> phi_ignend_mask = phi_k < 0;
        valarray<mmpbsa_t> phi_const = take(parminfo.dihedral_force_constants,
                              phi_m);
        valarray<mmpbsa_t> phi_periodicity = take(parminfo.dihedral_periodicities,
                                    phi_m);
        
        valarray<bool> phi_prd_mask = phi_periodicity < 0;
        phi_periodicity = abs(phi_periodicity);

        throw "test next line: phi_mask";
        valarray<bool> phi_mask = logical_not(phi_imp_mask) && logical_not(phi_ignend_mask)
                         && logical_not(phi_prd_mask);

        valarray<mmpbsa_t> phi_phase = take(parminfo.dihedral_phases, phi_m);
        //really want these???
        valarray<mmpbsa_t> phi_cos_phase = cos(phi_phase);
        valarray<mmpbsa_t> phi_sin_phase = sin(phi_phase);

        bondcodes = parminfo.dihedrals_inc_hydrogen;
        slice_array<int> phi_h_i = bondcodes[slice(0,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_h_j = bondcodes[slice(1,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_h_k = bondcodes[slice(2,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_h_l = bondcodes[slice(3,bondcodes.size()/columnWidth,columnWidth)];
        slice_array<int> phi_h_m = bondcodes[slice(4,bondcodes.size()/columnWidth,columnWidth)];
        phi_h_i /= 3;
        phi_h_j /= 3;
        phi_h_k = abs(phi_h_k)/3;
        phi_h_l = abs(phi_h_l)/3;
        valarray<bool> phi_h_imp_mask = phi_h_l < 0;
        valarray<bool> phi_h_ignend_mask = phi_h_k < 0;
        valarray<mmpbsa_t> phi_h_const = take(parminfo.dihedral_force_constants, phi_h_m);
        valarray<mmpbsa_t> phi_h_periodicity = take(parminfo.dihedral_periodicities,
                                      phi_h_m);
        valarray<mmpbsa_t> phi_h_phase = take(parminfo.dihedral_phases, phi_h_m);
        //really want these?
        valarray<mmpbsa_t> phi_h_cos_phase = cos(phi_h_phase);
        valarray<mmpbsa_t> phi_h_sin_phase = sin(phi_h_phase);

        valarray<mmpbsa_t> phi_h_periodicity = take(parminfo.dihedral_periodicities,
                                      phi_h_m);
        valarray<bool> phi_h_prd_mask = phi_h_periodicity < 0 ;
        phi_h_periodicity = abs(phi_h_periodicity);

        valarray<bool> phi_h_mask = logical_not(phi_h_imp_mask)
            && logical_not(phi_h_ignend_mask) && logical_not(phi_h_prd_mask);

        valarray<int> mol_ranges;
        if(parminfo.ifbox) //# We have solvent pointer info
        {
            valarray<int> molptrs = cumsum(parminfo.atoms_per_molecule);
            if(molptrs[molptrs.size()-1] != natom)
                throw "The last column sum of atoms per molecule should equal "
                        "the number of atoms";
            
            valarray<int> mol_ranges = get_mol_ranges(molptrs);//getranges should be of the type (min,max),(min,max),...
            if(mol_ranges.size() != parminfo.nspm)
                throw "The number of ranges must match the total number of molecules";
            
            int end_solute_atoms = res_ranges[slice(parminfo.iptres-1,res_ranges.size()/2,2)][1];
            int begin_solvent_atoms = mol_ranges[slice(parminfo.nspsol - 1,mol_ranges.size()/2,2)][0];
            //# TODO Make sure stripEnerFun fixes the above begin and end ptrs!
        }
        else//# Darn, no ATOMS_PER_MOLECULE INFO!  Walk the bond list
        {
            BondWalker bw();
            vector<int> beg_ptrs;
            valarray<int> markers(natom);
            throw "Check Next line: max() and just below in if...then throw";
            while((markers > 0).max())
            {
                int first_unmarked = find_first(markers,0);//find the index first then get the data
                if(find_first == markers.size())
                    throw MMPBSAException("find first did not find the trivial marker",DATA_FORMAT_ERROR);
                
                first_unmarked = markers[find_first(markers,0)];
                //# we insist that atoms of a molecule are all together
                //# in the list.
                valarray<bool> unmarked = (markers[first_unmarked:] > 0);
                if(unmarked.max())
                    throw MMPBSAException("All other markers should have been "
                            "zero", DATA_FORMAT_ERROR);
                
                beg_ptrs.push_back(first_unmarked);
                bw.walk(first_unmarked, [], visitor);
                //# bondwalk(self, first_unmarked, [], visitor)
            }
            valarray<int> shifted_beg_ptrs = cshift(beg_ptrs,1);
            shifted_beg_ptrs[shifted_beg_ptrs.size()-1] = natom;
            mol_ranges = zip(beg_ptrs, shifted_beg_ptrs);
        }
        
}/*end of constructor*/


std::valarray<int> EmpEnerFun::get_res_ranges(const std::valarray<int>& resptr, const int& natoms)
{
    std::valarray<int> res_ranges(2*resptr.size());

    for(int i = 0;i<resptr.size()-2;i+=2)
    {
        res_ranges[i] = resptr[i];
        res_ranges[i+1] = resptr[i+1];
    }
    res_ranges[resptr.size()-2] = resptr[resptr.size()-1];
    res_ranges[resptr.size()-1] = natoms;

    return res_ranges;
}

std::valarray<int> EmpEnerFun::get_mol_ranges(const std::valarray<int>& molptrs)
{
    //[0] + molptrs[:-1], molptrs
    std::valarray<int> mol_ranges(2*molptrs.size());

    mol_ranges[0] = molptrs[molptrs.size()-1];
    mol_ranges[1] = molptrs[0];

    for(int i = 1;i<molptrs.size();i+=2)
    {
        mol_ranges[i+1] = molptrs[i-1];
        mol_ranges[i+2] = molptrs[i];
    }

    return mol_ranges;
}


