#include "Energy.h"

EmpEnerFun::EmpEnerFun()
{
    parminfo = 0;
    resnames = 0;

    end_solute_atoms = -1;
    int begin_solvent_atoms = -1;
    inv_scnb = DEFAULT_SCNB;
    inv_scee = DEFAULT_SCEE;
    dielc = DEFAULT_DIELC;
    bond_i;
    bond_j;
    bond_k;
    bond_h_i;
    bond_h_j;
    bond_h_k;
    angle_i;
    angle_j;
    angle_k;
    angle_l;
    angle_h_i;
    angle_h_j;
    angle_h_k;
    angle_h_l;
    phi_i;
    phi_j;
    phi_k;
    phi_l;
    phi_m;
    phi_h_i;
    phi_h_j;
    phi_h_k;
    phi_h_l;
    phi_h_m;

    exclst;

    res_ranges;
    mol_ranges;

    bond_h_const;
    bond_h_eq;
    bond_const;
    bond_eq;
    angle_const;
    angle_eq;
    angle_h_const;
    angle_h_eq;
    phi_const;
    phi_periodicity;
    phi_phase;
    phi_cos_phase;
    phi_sin_phase;
    phi_h_const;
    phi_h_periodicity;
    phi_h_phase;
    //really want these?
    phi_h_cos_phase;
    phi_h_sin_phase;

    phi_prd_mask;
    phi_ignend_mask;
    phi_imp_mask;
    phi_mask;
    phi_h_prd_mask;
    phi_h_ignend_mask;
    phi_h_imp_mask;
    phi_h_mask;
}

EmpEnerFun::EmpEnerFun(sanderio::SanderParm * newparminfo, const mmpbsa_t& scnb,
        const mmpbsa_t& scee, const mmpbsa_t& dielc)
{
    using namespace mmpbsa_utils;
    using namespace std;

    parminfo = newparminfo;
    mmpbsa_t inv_scnb = 1.0 / scnb;
    mmpbsa_t inv_scee = 1.0 / scee;
    this->dielc = dielc;
    resnames = &(parminfo->residue_labels);
    end_solute_atoms = -1;//only needed when there is a solvent. Otherwise these values are -1 to indicate this is not the case.
    begin_solvent_atoms = -1;

    const int& ntypes = parminfo->ntypes;
    const int& natom = parminfo->natom;
    valarray<size_t> resptr = parminfo->residue_pointers - size_t(1); //sander file pointers are 1-indexed
    res_ranges = get_res_ranges(resptr,natom); //ranges is of the type (min,max),(min,max),...


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
    exclst;
    size_t num = 0;
    for (size_t i = 0; i < natom; i++) {
        size_t currentExclCount = numex[i];
        //# This compress is because a single zero element seems to
        //# be the parmtop way of indicating an empty list
        vector<size_t> irow;
        irow.reserve(currentExclCount);
        for(size_t j = 0;j<currentExclCount;j++)
            if(size_t curr = parminfo->excluded_atoms_list[num+j])
                irow.push_back(curr);
        

        //# shift indices to refer to list of atoms AFTER i:
        for (size_t j = 0; j < irow.size(); j++) {
            irow.at(j) -= i+1;//excluded_atoms_list is 0-indexed once read into SanderParm
        }
        exclst.push_back(irow);
        num += numex[i];
    }

    //temporary loading variables
    size_t columnWidth;

    //# bond stuff
    valarray<size_t>* bondcodes = &(parminfo->bonds_inc_hydrogen);
    columnWidth = 3;
    //divide for (i,j,k) bondcodes, divide i and j by 3 and decrement k
    bond_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    bond_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    bond_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    
    bond_h_const = parminfo->bond_force_constants[(*bondcodes)[bond_h_k]];
    bond_h_eq = parminfo->bond_equil_values[(*bondcodes)[bond_h_k]];

    bondcodes = &(parminfo->bonds_without_hydrogen);
    bond_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    bond_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    bond_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    bond_const = parminfo->bond_force_constants[(*bondcodes)[bond_k]];
    bond_eq = parminfo->bond_equil_values[(*bondcodes)[bond_k]];

    //# angle stuff
    bondcodes = &(parminfo->angles_without_hydrogen);
    columnWidth = 4;
    angle_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    angle_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    angle_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    angle_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    angle_const = parminfo->angle_force_constants[(*bondcodes)[angle_l]];
    angle_eq = parminfo->angle_equil_values[(*bondcodes)[angle_l]];

    bondcodes = &(parminfo->angles_inc_hydrogen);
    angle_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    angle_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    angle_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    angle_h_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    angle_h_const = parminfo->angle_force_constants[(*bondcodes)[angle_h_l]];
    angle_h_eq = parminfo->angle_equil_values[(*bondcodes)[angle_h_l]];

    //# dihedral stuff
    bondcodes = &(parminfo->dihedrals_without_hydrogen);
    columnWidth = 5;
    phi_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    phi_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    phi_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    phi_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    phi_m = slice(4, bondcodes->size() / columnWidth, columnWidth);
    phi_const = parminfo->dihedral_force_constants[(*bondcodes)[phi_m]];
    phi_periodicity = parminfo->dihedral_periodicities[(*bondcodes)[phi_m]];
    phi_prd_mask = phi_periodicity < 0.0;
    phi_ignend_mask = parminfo->dihedral_mask[phi_k];
    phi_imp_mask = parminfo->dihedral_mask[phi_l];
    phi_periodicity = abs(phi_periodicity);

    phi_mask = (!phi_imp_mask) && (!phi_ignend_mask) && (!phi_prd_mask);

    phi_phase = parminfo->dihedral_phases[parminfo->dihedrals_without_hydrogen[phi_m]];
    //really want these???
    phi_cos_phase = cos(phi_phase);
    phi_sin_phase = sin(phi_phase);

    bondcodes = &(parminfo->dihedrals_inc_hydrogen);
    phi_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    phi_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    phi_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    phi_h_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    phi_h_m = slice(4, bondcodes->size() / columnWidth, columnWidth);
    phi_h_ignend_mask = parminfo->dihedral_h_mask[phi_h_k];
    phi_h_imp_mask = parminfo->dihedral_h_mask[phi_h_l];
    phi_h_const = parminfo->dihedral_force_constants[(*bondcodes)[phi_h_m]];
    phi_h_periodicity = parminfo->dihedral_periodicities[(*bondcodes)[phi_h_m]];
    phi_h_phase = parminfo->dihedral_phases[(*bondcodes)[phi_h_m]];
    //really want these?
    phi_h_cos_phase = cos(phi_h_phase);
    phi_h_sin_phase = sin(phi_h_phase);

    phi_h_prd_mask = phi_h_periodicity < 0.0;
    phi_h_periodicity = abs(phi_h_periodicity);

    phi_h_mask = (!phi_h_imp_mask) && (!phi_h_ignend_mask) && (!phi_h_prd_mask);

    mol_ranges;
    if (parminfo->ifbox) //# We have solvent pointer info
    {
        valarray<size_t> molptrs = cumsum(parminfo->atoms_per_molecule);
        if (molptrs[molptrs.size() - 1] != natom)
        {
            char * error;
            int lastSum = int(molptrs[molptrs.size()-1]);
            sprintf(error,"The last column sum of atoms per molecule should equal "
            "the number of atoms (%d) but instead it equaled %d",natom,lastSum);
            throw MMPBSAException(error,DATA_FORMAT_ERROR);
        }

        valarray<size_t> mol_ranges = get_mol_ranges(molptrs); //getranges should be of the type (min,max),(min,max),...
        if (mol_ranges.size()/2 != parminfo->nspm)
            throw "The number of ranges must match the total number of molecules";

        end_solute_atoms = res_ranges[1+2*(parminfo->iptres - 1)];
        begin_solvent_atoms = mol_ranges[2*(parminfo->nspsol - 1)];
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
    EmpEnerFun tmp(newparminfo, DEFAULT_SCNB, DEFAULT_SCEE, DEFAULT_DIELC);
    *this = tmp;
}/**/

EmpEnerFun::EmpEnerFun(const EmpEnerFun& orig)
{

}

std::valarray<size_t> EmpEnerFun::get_res_ranges(const std::valarray<size_t>& resptr,
        const size_t& natoms) {
    std::valarray<size_t> res_ranges(2 * resptr.size());

    size_t res_range_index = 0;
    for(size_t i = 0;i<resptr.size()-2;i++)
    {
        res_ranges[res_range_index++] = resptr[i];
        res_ranges[res_range_index++] = resptr[i+1];
    }
    res_ranges[resptr.size() - 2] = resptr[resptr.size() - 1];
    res_ranges[resptr.size() - 1] = natoms;

    return res_ranges;
}

std::valarray<size_t> EmpEnerFun::get_mol_ranges(const std::valarray<size_t>& molptrs) {
    //[0] + molptrs[:-1], molptrs
    std::valarray<size_t> mol_ranges(2 * molptrs.size());

    size_t mol_range_index = 0;
    mol_ranges[mol_range_index++] = 0;
    mol_ranges[mol_range_index++] = molptrs[0];

    for (size_t i = 1; i < molptrs.size(); i++)
    {
        mol_ranges[mol_range_index++] = molptrs[i];
        mol_ranges[mol_range_index++] = molptrs[i+1];
    }

    return mol_ranges;
}

EmpEnerFun EmpEnerFun::stripEnerFun(const std::valarray<bool>& keepers,
        const bool& dangleWarn) const
{
    using std::valarray;
    using std::vector;
    using std::slice;
    using namespace sanderio;
    
    if(keepers.size() != parminfo->natom)
    {
        char* error;
        sprintf(error,"Energy function has a size %d. The size of the array used"
            "to strip the energy function is %d .",parminfo->natom,keepers.size());
        throw MMPBSAException(error,DATA_FORMAT_ERROR);
    }

    //Will belong to the new EmpEnerFun (therefore, do not delete at the end)
    SanderParm * sp = new SanderParm(*parminfo);
    //new EmpEnerFun 
    EmpEnerFun returnMe;

    //create new atom index list
    size_t keepersSize = keepers.size();
    sp->natom = 0;
    valarray<size_t> newidx = mmpbsa_utils::cumBoolSum<bool>(keepers);
    for(size_t i = 0;i<keepersSize;i++)
    {
        if(keepers[i])
        {
            newidx[i]--;
        }
        else
        {
            newidx[i] = keepersSize;
        }
    }
    sp->natom = keepersSize;
    sp->atom_names = parminfo->atom_names[keepers];
    sp->charges = parminfo->charges[keepers];
    sp->masses = parminfo->masses[keepers];
    sp->atom_type_indices = parminfo->atom_type_indices[keepers];

    returnMe.dielc = dielc;
    returnMe.inv_scnb = inv_scnb;
    returnMe.inv_scee = inv_scee;

    //update residue pointers and names
    size_t apos = 0;
    vector<size_t> new_res_ranges;
    for(size_t i = 0;i<resnames->size()-1;i++)
    {
        size_t rbegin = res_ranges[2*i];
        size_t rend = res_ranges[2*i+1];
        size_t rnats;
        if(rbegin > 0)
            rnats = newidx[rend]-newidx[rbegin-1]+2;//plus two because newidx = cumsum(keepers[0:position])-1
        else
            rnats = newidx[rend];

        if(rnats)
        {
            new_res_ranges.push_back(apos);
            apos += rnats;
            new_res_ranges.push_back(apos);

        }
    }
    
    //update mol ranges
    size_t molpos = 0;
    vector<size_t> new_mol_ranges;
    for(size_t i = 0;i<mol_ranges.size()/2;i++)
    {
        size_t molbegin = mol_ranges[2*i];
        size_t molend = mol_ranges[2*i+1];
        size_t molnats;
        if(molbegin > 0)
            molnats = newidx[molend]-newidx[molbegin-1]+2;//plus two because newidx = cumsum(keepers[0:position])-1
        else
            molnats = newidx[molend];

        if(molnats)
        {
            new_mol_ranges.push_back(molpos);
            molpos += molnats;
            new_mol_ranges.push_back(molpos);

        }
    }

    //check to see if there are solvent pointers to be kept. If so, clean them.
    if(begin_solvent_atoms != -1)//-1 means there is no solvent
    {
        if(newidx[newidx.size()-1] - newidx[size_t(begin_solvent_atoms)] > 0)
        {
            returnMe.begin_solvent_atoms = newidx[size_t(begin_solvent_atoms)]+1;
            returnMe.end_solute_atoms = newidx[size_t(end_solute_atoms)]+1;
        }
        else
        {
            returnMe.begin_solvent_atoms = -1;
            returnMe.end_solute_atoms = -1;
        }
    }

    //rebuild excluded atoms list
    //sp->excluded_atoms_list.resize();
    returnMe.exclst.resize(0);
    for(size_t i = 0;i<parminfo->natom;i++)
        if(keepers[i])
        {
            valarray<size_t> oldirow(exclst[i].size());
            for(size_t j = 0;j<oldirow.size();j++)
                oldirow[j] = exclst[i][j] + i + 1;
            
            valarray<bool> keepo = keepers[oldirow];
            if(dangleWarn && keepo.min() == 0)
                std::cerr << "some danglers in excl list of old atom" << i << std::endl;

            vector<size_t> irow = mmpbsa_utils::take(newidx,oldirow[keepo]);
            size_t shiftNewValues = returnMe.exclst.size()+1;
            for(size_t j = 0;j<irow.size();j++)
                irow[j] -= shiftNewValues;
            returnMe.exclst.push_back(irow);
        }
        
    //convert energy value tables and their indices
    valarray<slice> bond_slices(3);
    valarray<const valarray<mmpbsa_t>* > bond_values(2);
    valarray<valarray<size_t> > newIndices(2);
    valarray<valarray<mmpbsa_t> > newTerms(2);

    //convert H-bonds
    bond_slices[0] = bond_h_i;bond_slices[1] = bond_h_j;bond_slices[2] = bond_h_k;
    bond_values[0] = &bond_h_const;bond_values[1] = &bond_h_eq;
    EmpEnerFun::internalConvert<size_t,mmpbsa_t>(
        newIndices,newTerms,parminfo->bonds_inc_hydrogen,
            bond_slices,bond_values,newidx,keepers,false);
    sp->bonds_inc_hydrogen = mmpbsa_utils::zip(newIndices);
    returnMe.bond_h_i = slice(0,sp->bonds_inc_hydrogen.size()/3,3);
    returnMe.bond_h_j= slice(1,sp->bonds_inc_hydrogen.size()/3,3);
    returnMe.bond_h_k = slice(2,sp->bonds_inc_hydrogen.size()/3,3);
    returnMe.bond_h_const = newTerms[0];
    returnMe.bond_h_eq = newTerms[1];

    //convert non-H-bonds
    bond_slices[0] = bond_i;bond_slices[1] = bond_j;bond_slices[2] = bond_k;
    bond_values[0] = &bond_const;bond_values[1] = &bond_eq;
    internalConvert(newIndices,newTerms,parminfo->bonds_without_hydrogen,bond_slices,bond_values,newidx,keepers,false);
    sp->bonds_without_hydrogen = mmpbsa_utils::zip(newIndices);//same as with hydrogen.
    returnMe.bond_i = slice(0,sp->bonds_without_hydrogen.size()/3,3);
    returnMe.bond_j= slice(1,sp->bonds_without_hydrogen.size()/3,3);
    returnMe.bond_k = slice(2,sp->bonds_without_hydrogen.size()/3,3);
    returnMe.bond_const = newTerms[0];
    returnMe.bond_eq = newTerms[1];

    //convert non-H-angles
    bond_slices.resize(4);
    bond_slices[0] = angle_i;bond_slices[1] = angle_j;bond_slices[2] = angle_k;bond_slices[3] = angle_l;
    bond_values[0] = &angle_const;bond_values[1] = &angle_eq;
    internalConvert(newIndices,newTerms,parminfo->angles_without_hydrogen,bond_slices,bond_values,newidx,keepers,false);
    sp->angles_without_hydrogen = mmpbsa_utils::zip(newIndices);//same as with hydrogen.
    returnMe.angle_i = slice(0,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_j = slice(1,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_k = slice(2,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_l = slice(3,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_const = newTerms[0];
    returnMe.angle_eq = newTerms[1];

    //convert H-angles
    bond_slices[0] = angle_h_i;bond_slices[1] = angle_h_j;bond_slices[2] = angle_h_k;bond_slices[3] = angle_h_l;
    bond_values[0] = &angle_h_const;bond_values[1] = &angle_h_eq;
    internalConvert(newIndices,newTerms,parminfo->angles_inc_hydrogen,bond_slices,bond_values,newidx,keepers,false);
    sp->angles_inc_hydrogen = mmpbsa_utils::zip(newIndices);//same as with hydrogen.
    returnMe.angle_h_i = slice(0,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_j = slice(1,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_k = slice(2,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_l = slice(3,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_const = newTerms[0];
    returnMe.angle_eq = newTerms[1];

    //convert non-H-dihedrals
    bond_slices.resize(5);
    bond_slices[0] = phi_i;bond_slices[1] = phi_j;bond_slices[2] = phi_k;bond_slices[3] = phi_l;bond_slices[4] = phi_m;
    bond_values.resize(5);
    bond_values[0] = &phi_const;bond_values[1] = &phi_periodicity;bond_values[2] = &phi_phase;
    bond_values[3] = &phi_cos_phase;bond_values[4] = &phi_sin_phase;
    internalConvert(newIndices,newTerms,parminfo->dihedrals_without_hydrogen,bond_slices,bond_values,newidx,keepers,false);
    sp->dihedrals_without_hydrogen = mmpbsa_utils::zip(newIndices);
    returnMe.phi_i = slice(0,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_j = slice(1,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_k = slice(2,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_l = slice(3,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_m = slice(4,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_const = newTerms[0];returnMe.phi_periodicity = newTerms[1];
    returnMe.phi_phase = newTerms[2];returnMe.phi_cos_phase = newTerms[3];
    returnMe.phi_sin_phase = newTerms[4];
    //create masks
    returnMe.phi_prd_mask = returnMe.phi_periodicity < 0.0;
    returnMe.phi_ignend_mask = sp->dihedral_mask[phi_k];
    returnMe.phi_imp_mask = sp->dihedral_mask[phi_l];
    returnMe.phi_periodicity = abs(phi_periodicity);
    returnMe.phi_mask = (!returnMe.phi_imp_mask) && (!returnMe.phi_ignend_mask)
            && (!returnMe.phi_prd_mask);
    
    
    //convert H-dihedrals
    bond_slices[0] = phi_h_i;bond_slices[1] = phi_h_j;bond_slices[2] = phi_h_k;bond_slices[3] = phi_h_l;bond_slices[4] = phi_h_m;
    bond_values[0] = &phi_h_const;bond_values[1] = &phi_h_periodicity;bond_values[2] = &phi_h_phase;
    bond_values[3] = &phi_h_cos_phase;bond_values[4] = &phi_h_sin_phase;
    internalConvert(newIndices,newTerms,parminfo->dihedrals_inc_hydrogen,bond_slices,bond_values,newidx,keepers,false);
    sp->dihedrals_inc_hydrogen = mmpbsa_utils::zip(newIndices);
    returnMe.phi_h_i = slice(0,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_j = slice(1,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_k = slice(2,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_l = slice(3,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_m = slice(4,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_const = newTerms[0];returnMe.phi_h_periodicity = newTerms[1];
    returnMe.phi_h_phase = newTerms[2];returnMe.phi_h_cos_phase = newTerms[3];
    returnMe.phi_h_sin_phase = newTerms[4];
    //create masks
    returnMe.phi_h_prd_mask = returnMe.phi_h_periodicity < 0.0;
    returnMe.phi_h_ignend_mask = sp->dihedral_h_mask[phi_k];
    returnMe.phi_h_imp_mask = sp->dihedral_h_mask[phi_l];
    returnMe.phi_h_periodicity = abs(phi_h_periodicity);
    returnMe.phi_h_mask = (!returnMe.phi_h_imp_mask) && (!returnMe.phi_h_ignend_mask)
            && (!returnMe.phi_h_prd_mask);

    returnMe.parminfo = sp;//thus, don't delete sp.

}//end stripEnerFun(...)

template <class M, class N> void EmpEnerFun::internalConvert(
    std::valarray<std::valarray<M> >& newIndices,std::valarray<std::valarray<N> >& newTerms,
        const std::valarray<M>& oldIndices,const std::valarray<std::slice>& slices,
        const std::valarray<const std::valarray<N>* >& oldTerms,
        const std::valarray<size_t>& newidx,
        const std::valarray<bool> keepers, const bool& dangleWarn)
{
    using std::vector;
    using std::valarray;
    throw "Verify array size match.";
    vector<valarray<size_t> > proto_newidx;
    valarray<bool> dangleAnd;
    valarray<bool> dangleOr;

    //dangle test dangleOr is false if an element is ALWAYS zero,
    //dangleAnd is false if it's been zero at least one time.
    //TEST:
    //!(dangleAnd || !dangleOr) is a dangle
    if(dangleWarn)
    {
        dangleAnd.resize(keepers.size(),false);
        dangleOr.resize(keepers.size(),false);
    }

    for(size_t i = 0;i<slices.size();i++)
    {
        std::slice currSlice = slices[i];
        valarray<M> keep_array = (oldIndices[currSlice])[keepers];
        proto_newidx.push_back(newidx[keep_array]);
        if(i != 0)
        {
            dangleAnd = dangleAnd && (keep_array == M(0));
            dangleOr = dangleOr || (keep_array == M(0));
        }
        else
        {
            dangleAnd = dangleOr = (keep_array == M(0));
        }
        
    }

    valarray<bool> bkeep = dangleAnd;//keep only atoms for which the lists have always had values.

    //Pack new arrays
    newIndices.resize(proto_newidx.size());
    for(size_t i = 0;i<proto_newidx.size();i++)
        newIndices[i] = (proto_newidx[i][bkeep]);

    newTerms.resize(oldTerms.size());
    for(size_t i = 0;i<oldTerms.size();i++)
        newTerms[i] = (*oldTerms[i])[bkeep];

    if(dangleWarn)
    {
        valarray<bool> danglers = !(bkeep);
        if(danglers.max())
            std::cerr << "Danglers were found.\n" << std::endl;
    }

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
