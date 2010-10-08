#include "EmpEnerFun.h"

mmpbsa::EmpEnerFun::EmpEnerFun()
{
    parminfo = 0;

    end_solute_atoms = -1;
    begin_solvent_atoms = -1;
    inv_scnb = 1.0/DEFAULT_SCNB;
    inv_scee =  1.0/DEFAULT_SCEE;
    dielc =  DEFAULT_DIELC;
    
}

mmpbsa::EmpEnerFun::EmpEnerFun(mmpbsa::SanderParm * newparminfo, const mmpbsa_t& scnb,
        const mmpbsa_t& scee, const mmpbsa_t& dielc)
{
    using namespace mmpbsa_utils;
    using namespace std;

    parminfo = newparminfo;
    this->inv_scnb = 1.0 / scnb;
    this->inv_scee = 1.0 / scee;
    this->dielc = dielc;
    end_solute_atoms = -1;//only needed when there is a solvent. Otherwise these values are -1 to indicate this is not the case.
    begin_solvent_atoms = -1;

    const size_t& ntypes = parminfo->ntypes;
    const size_t& natom = parminfo->natom;
    valarray<size_t> resptr = parminfo->residue_pointers - size_t(1); //sander file pointers are 1-indexed
    res_ranges.resize(2*resptr.size());
    res_ranges = get_res_ranges(resptr,natom); //ranges is of the type (min,max),(min,max),...


    //Non-bonded stuff
    const valarray<size_t>& nb_parm_idx = parminfo->nonbonded_parm_indices;
    const valarray<mmpbsa_t>& CA = parminfo->lennard_jones_acoefs;
    const valarray<mmpbsa_t>& CB = parminfo->lennard_jones_bcoefs;
    //const valarray<mmpbsa_t>& CHA = parminfo->hbond_acoefs;
    //const valarray<mmpbsa_t>& CHB = parminfo->hbond_bcoefs;
    LJA.resize(ntypes * ntypes);
    LJB.resize(ntypes * ntypes);
    //valarray<mmpbsa_t> LJHA(ntypes * ntypes);//according to grep in pyamber, unused.
    //valarray<mmpbsa_t> LJHB(ntypes, ntypes);

    for (size_t i = 0; i < ntypes; i++) {
        for (size_t j = i; j < ntypes; j++) {
            size_t ico = nb_parm_idx[j+i*ntypes];
            bool isTenTwelvePair = parminfo->nonbonded_parm_mask[i+j*ntypes];
            if (!isTenTwelvePair) {
                LJA[j+i*ntypes] = CA[ico - 1];
                LJB[j+i*ntypes] = CB[ico - 1];
            } /*/else {
                LJHA[j+i*ntypes] = CHA[ico - 1];
                LJHB[j+i*ntypes] = CHB[ico - 1];
            } see above note about LJHA */
        }
    }

    //# Fill in the other triangle
    for (size_t i = 0; i < ntypes; i++) {
        for (size_t j = i + 1; j < ntypes; j++) {
            LJA[i+j*ntypes] = LJA[j+i*ntypes];
            LJB[i+j*ntypes] = LJB[j+i*ntypes];
            //LJHA[i+j*ntypes] = LJHA[j+i*ntypes];/* see above note about LJHA */
            //LJHB[i+j*ntypes] = LJHB[j+i*ntypes];
        }
    }

    //# excluded atom stuff
    //const valarray<size_t>& numex = parminfo->number_excluded_atoms;
    exclst.resize(natom);
    size_t num = 0;
    for (size_t i = 0; i < natom; i++) {
        exclst.at(i).clear();
        size_t currentExclCount = parminfo->number_excluded_atoms[i];
        //# This compress is because a single zero element seems to
        //# be the parmtop way of indicating an empty list
        for(size_t j = 0;j<currentExclCount;j++)
        {
            //A signal zero-valued element in the parmtop file means an empty list.
            if(parminfo->excluded_atoms_list[num+j] > 0)
            {
                size_t curr = parminfo->excluded_atoms_list[num+j] - size_t(1);//Re: offset, see Amber8 Manual Appendix C
                exclst.at(i).push_back(curr);
            }
        }
        

        //# shift indices to refer to list of atoms AFTER i:
        for (size_t j = 0; j < exclst.at(i).size(); j++) {
            exclst.at(i).at(j) -= i+1;//excluded_atoms_list is 0-indexed once read into SanderParm
        }
        num += currentExclCount;
    }

    //temporary loading variables
    size_t columnWidth;

    //# bond stuff
    //Bond Codes are stored in SanderParm. Slices are provided here for easy
    //access to i,j,k components. However, there is no need to duplicate the data;
    //therefore it remains in SanderParm.
    valarray<size_t>* bondcodes = &(parminfo->bonds_inc_hydrogen);
    columnWidth = 3;
    //divide for (i,j,k) bondcodes, divide i and j by 3 and decrement k
    bond_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    bond_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    bond_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    
    bondcodes = &(parminfo->bonds_without_hydrogen);
    bond_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    bond_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    bond_k = slice(2, bondcodes->size() / columnWidth, columnWidth);

    //# angle stuff
    bondcodes = &(parminfo->angles_without_hydrogen);
    columnWidth = 4;
    angle_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    angle_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    angle_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    angle_l = slice(3, bondcodes->size() / columnWidth, columnWidth);

    bondcodes = &(parminfo->angles_inc_hydrogen);
    angle_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    angle_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    angle_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    angle_h_l = slice(3, bondcodes->size() / columnWidth, columnWidth);

    //# dihedral stuff
    bondcodes = &(parminfo->dihedrals_without_hydrogen);
    columnWidth = 5;
    phi_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    phi_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    phi_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    phi_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    phi_m = slice(4, bondcodes->size() / columnWidth, columnWidth);

    bondcodes = &(parminfo->dihedrals_inc_hydrogen);
    phi_h_i = slice(0, bondcodes->size() / columnWidth, columnWidth);
    phi_h_j = slice(1, bondcodes->size() / columnWidth, columnWidth);
    phi_h_k = slice(2, bondcodes->size() / columnWidth, columnWidth);
    phi_h_l = slice(3, bondcodes->size() / columnWidth, columnWidth);
    phi_h_m = slice(4, bondcodes->size() / columnWidth, columnWidth);

    if (parminfo->ifbox) //# We have solvent pointer info
    {
        valarray<size_t> molptrs = cumsum(parminfo->atoms_per_molecule);
        if (molptrs[molptrs.size() - 1] != natom)
        {
            std::ostringstream error;
            int lastSum = int(molptrs[molptrs.size()-1]);
            error << "The last column sum of atoms per molecule should equal "
            "the number of atoms (" << natom << ") but instead it equaled " << lastSum;
            throw MMPBSAException(error,DATA_FORMAT_ERROR);
        }

        mol_ranges.resize(2*molptrs.size());
        mol_ranges = get_mol_ranges(molptrs); //getranges should be of the type (min,max),(min,max),...
        if (mol_ranges.size()/2 != parminfo->nspm)
            throw "The number of ranges must match the total number of molecules";

        end_solute_atoms = res_ranges[1+2*(parminfo->iptres - 1)];
        begin_solvent_atoms = mol_ranges[2*(parminfo->nspsol - 1)];
        //# TODO Make sure stripEnerFun fixes the above begin and end ptrs!
    }
    else//# Darn, no ATOMS_PER_MOLECULE INFO!  Walk the bond list
    {
        BondWalker bw(this);//NEED to test BondWalker
        vector<size_t> beg_ptrs;
        valarray<int> markers(natom);
        while(markers.min() == 0)
        {
            size_t first_unmarked = find_first(markers,0);
            if(first_unmarked == markers.size())//this means no value of zero was found. But that shouldn't happen in this scope.
            	throw MMPBSAException("find_first gave unexpected results in EmpEnerFun::stripEnerFun.",DATA_FORMAT_ERROR);
            for(size_t i = first_unmarked;i<markers.size();i++)
                if(markers[i])
                    throw MMPBSAException("Atoms of a particular molecule must be stored "
                            "together",DATA_FORMAT_ERROR);
            beg_ptrs.push_back(first_unmarked);
            bw.walk(first_unmarked,valarray<int>(0),markers,visitor);
        }
        beg_ptrs.push_back(parminfo->natom);
        mol_ranges.resize(size_t(2*beg_ptrs.size()-2));
        size_t mol_ranges_index = 0;
        for(vector<size_t>::iterator it = beg_ptrs.begin();it != beg_ptrs.end()-1;it++)
        {
            mol_ranges[mol_ranges_index++] = *it;
            mol_ranges[mol_ranges_index++] = *(it+1);
        }
    }

}/*end of constructor*/

mmpbsa::EmpEnerFun::EmpEnerFun(const mmpbsa::EmpEnerFun& orig)
{
    using std::vector;
    using std::valarray;
    
    parminfo = orig.parminfo;
    
    end_solute_atoms = orig.end_solute_atoms;
    begin_solvent_atoms = orig.begin_solvent_atoms;
    inv_scnb = orig.inv_scnb;
    inv_scee = orig.inv_scee;
    dielc = orig.dielc;
    bond_i = orig.bond_i;
    bond_j = orig.bond_j;
    bond_k = orig.bond_k;
    bond_h_i = orig.bond_h_i;
    bond_h_j = orig.bond_h_j;
    bond_h_k = orig.bond_h_k;
    angle_i = orig.angle_i;
    angle_j = orig.angle_j;
    angle_k = orig.angle_k;
    angle_l = orig.angle_l;
    angle_h_i = orig.angle_h_i;
    angle_h_j = orig.angle_h_j;
    angle_h_k = orig.angle_h_k;
    angle_h_l = orig.angle_h_l;
    phi_i = orig.phi_i;
    phi_j = orig.phi_j;
    phi_k = orig.phi_k;
    phi_l = orig.phi_l;
    phi_m = orig.phi_m;
    phi_h_i = orig.phi_h_i;
    phi_h_j = orig.phi_h_j;
    phi_h_k = orig.phi_h_k;
    phi_h_l = orig.phi_h_l;
    phi_h_m = orig.phi_h_m;

    exclst.resize(orig.exclst.size());
    exclst = orig.exclst;

    res_ranges.resize(orig.res_ranges.size());
    res_ranges = orig.res_ranges;

    mol_ranges.resize(orig.mol_ranges.size());
    mol_ranges = orig.mol_ranges;

    LJA.resize(orig.LJA.size());
    LJA = orig.LJA;
    LJB.resize(orig.LJB.size());
    LJB = orig.LJB;
   
}

//do not delete parminfo. It is externally made and should be deleted outside of EmpEnerInfo
mmpbsa::EmpEnerFun::~EmpEnerFun()
{
    
}

mmpbsa::EmpEnerFun& mmpbsa::EmpEnerFun::operator=(const mmpbsa::EmpEnerFun& orig)
{
    using std::vector;
    using std::valarray;

    if(this == &orig)
        return *this;

    parminfo = orig.parminfo;
    
    end_solute_atoms = orig.end_solute_atoms;
    begin_solvent_atoms = orig.begin_solvent_atoms;
    inv_scnb = orig.inv_scnb;
    inv_scee = orig.inv_scee;
    dielc = orig.dielc;
    bond_i = orig.bond_i;
    bond_j = orig.bond_j;
    bond_k = orig.bond_k;
    bond_h_i = orig.bond_h_i;
    bond_h_j = orig.bond_h_j;
    bond_h_k = orig.bond_h_k;
    angle_i = orig.angle_i;
    angle_j = orig.angle_j;
    angle_k = orig.angle_k;
    angle_l = orig.angle_l;
    angle_h_i = orig.angle_h_i;
    angle_h_j = orig.angle_h_j;
    angle_h_k = orig.angle_h_k;
    angle_h_l = orig.angle_h_l;
    phi_i = orig.phi_i;
    phi_j = orig.phi_j;
    phi_k = orig.phi_k;
    phi_l = orig.phi_l;
    phi_m = orig.phi_m;
    phi_h_i = orig.phi_h_i;
    phi_h_j = orig.phi_h_j;
    phi_h_k = orig.phi_h_k;
    phi_h_l = orig.phi_h_l;
    phi_h_m = orig.phi_h_m;

    exclst.resize(orig.exclst.size());
    exclst = orig.exclst;

    res_ranges.resize(orig.res_ranges.size());
    res_ranges = orig.res_ranges;

    mol_ranges.resize(orig.mol_ranges.size());
    mol_ranges = orig.mol_ranges;

    LJA.resize(orig.LJA.size());
    LJA = orig.LJA;
    LJB.resize(orig.LJB.size());
    LJB = orig.LJB;

    return *this;
}

std::string mmpbsa::EmpEnerFun::ereport(const std::valarray<mmpbsa_t>& crds)
{
    char ereport[512];
    mmpbsa_t bon,ang,dihe,vdw14,ele14,vdw,ele;
    bon = total_bond_energy(crds);
    ang = total_angle_energy(crds);
    dihe = total_dihedral_energy(crds);
    vdw14 = total_vdw14_energy(crds);
    ele14 = total_elstat14_energy(crds);
    vdw = total_vdwaals_energy(crds);
    ele = total_elstat_energy(crds);

    std::sprintf(ereport,"\n"
            "BOND    =  %12.4f  ANGLE   =  %12.4f  DIHED      =  %12.4f\n"
            "VDWAALS =  %12.4f  EEL     =  %12.4f  EGB        =  %12.4f\n"
            "1-4 VDW =  %12.4f  1-4 EEL =  %12.4f  RESTRAINT  =  %12.4f\n",
            bon, ang, dihe, 
            vdw, ele, 0.0,
            vdw14, ele14, 0.0 );

    std::string returnMe(ereport);
    return returnMe;
}

std::valarray<size_t> mmpbsa::EmpEnerFun::get_res_ranges(const std::valarray<size_t>& resptr,
        const size_t& natoms) {
    std::valarray<size_t> res_ranges(2 * resptr.size());

    size_t res_range_index = 0;
    for(size_t i = 0;i<resptr.size()-1;i++)
    {
        res_ranges[res_range_index++] = resptr[i];
        res_ranges[res_range_index++] = resptr[i+1];
    }
    res_ranges[resptr.size() - 2] = resptr[resptr.size() - 1];
    res_ranges[resptr.size() - 1] = natoms;

    return res_ranges;
}

std::valarray<size_t> mmpbsa::EmpEnerFun::get_mol_ranges(const std::valarray<size_t>& molptrs) {
    std::valarray<size_t> mol_ranges(2 * molptrs.size());

    size_t mol_range_index = 0;
    mol_ranges[mol_range_index++] = 0;
    mol_ranges[mol_range_index++] = molptrs[0];

    for (size_t i = 0; i < molptrs.size()-1; i++)
    {
        mol_ranges[mol_range_index++] = molptrs[i];
        mol_ranges[mol_range_index++] = molptrs[i+1];
    }

    return mol_ranges;
}

mmpbsa::EmpEnerFun mmpbsa::EmpEnerFun::stripEnerFun(const std::valarray<bool>& keepers,
        const bool& dangleWarn) const
{
    using std::valarray;
    using std::vector;
    using std::slice;
    using namespace mmpbsa;
    
    if(keepers.size() != parminfo->natom)
    {
        std::ostringstream error;
        error << "Energy function has a size " << parminfo->natom << ". The size of the array used"
            "to strip the energy function is " << keepers.size() << " .";
        throw MMPBSAException(error,DATA_FORMAT_ERROR);
    }

    //Will belong to the new EmpEnerFun (therefore, do not delete at the end)
    SanderParm * sp = new SanderParm(*parminfo);
    //new EmpEnerFun 
    EmpEnerFun returnMe;
    returnMe.dielc = dielc;
    returnMe.inv_scnb = inv_scnb;
    returnMe.inv_scee = inv_scee;

    //create new atom index list
    size_t keepersSize = keepers.size();
    sp->natom = 0;
    valarray<size_t> newidx = mmpbsa_utils::cumBoolSum<bool>(keepers);
    for(size_t i = 0;i<keepersSize;i++)
    {
        if(keepers[i])
        {
            newidx[i]--;//newidx refers to the index the atom will have in *new* data structures.
            sp->natom++;
        }
        else
        {
            newidx[i] = keepersSize;
        }
    }

    //copy to the new Energy object only values that correspond to kept atoms.
    sp->atom_names.resize(sp->natom);
    sp->charges.resize(sp->natom);
    sp->masses.resize(sp->natom);
    sp->atom_type_indices.resize(sp->natom);
    size_t newAtomNameIndex = 0;
    for(size_t i = 0;i<parminfo->natom;i++)
        if(keepers[i])
        {
            sp->atom_names[newAtomNameIndex] = parminfo->atom_names[i];
            sp->charges[newAtomNameIndex] = parminfo->charges[i];
            sp->masses[newAtomNameIndex] = parminfo->masses[i];
            sp->atom_type_indices[newAtomNameIndex] = parminfo->atom_type_indices[i];
            newAtomNameIndex++;
        }

    
    //update residue pointers and names
    size_t apos = 0;
    vector<size_t> new_res_ranges;
    size_t residue_lable_index = 0;
    for(size_t i = 0;i<res_ranges.size();i+=2)
    {
        size_t rbegin = res_ranges[i];
        size_t rend = res_ranges[i+1];
        size_t rnats = 0;
        for(size_t j = rbegin;j<rend;j++)
            if(keepers[j])
                rnats++;
        
        if(rnats)
        {
            new_res_ranges.push_back(apos);
            apos += rnats;
            new_res_ranges.push_back(apos);
            new_res_ranges.push_back(residue_lable_index);//pointer to residue_label used below.
        }
        residue_lable_index++;
    }

    size_t rangeLength = size_t(new_res_ranges.size()*2/3);
    returnMe.res_ranges.resize(rangeLength);
    sp->residue_labels.resize(size_t(new_res_ranges.size()/3));
    residue_lable_index = 0;
    size_t residue_range_index = 0;
    for(size_t i = 0;i<new_res_ranges.size();i+=3)
    {
        returnMe.res_ranges[residue_range_index++] = new_res_ranges[i];
        returnMe.res_ranges[residue_range_index++] = new_res_ranges[i+1];
        sp->residue_labels[residue_lable_index++] = parminfo->residue_labels[new_res_ranges[i+2]];
    }
    
    //update mol ranges
    size_t molpos = 0;
    vector<size_t> new_mol_ranges;
    for(size_t i = 0;i<mol_ranges.size();i+=2)
    {
        size_t molbegin = mol_ranges[i];
        size_t molend = mol_ranges[i+1];
        size_t molnats = 0;
        for(size_t j = molbegin;j<molend;j++)
            if(keepers[j])
                molnats++;

        if(molnats)
        {
            new_mol_ranges.push_back(molpos);
            molpos += molnats;
            new_mol_ranges.push_back(molpos);
        }
    }
    returnMe.mol_ranges.resize(new_mol_ranges.size());
    for(size_t i = 0;i<new_mol_ranges.size();i++)
        returnMe.mol_ranges[i] = new_mol_ranges[i];


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
    for(size_t i = 0;i<parminfo->natom;i++)
        if(keepers[i])
        {
            valarray<size_t> oldirow(exclst.at(i).size());
            valarray<bool> keepo(oldirow.size());
            bool haveDangler = false;
            for(size_t j = 0;j<oldirow.size();j++)
            {
                oldirow[j] = exclst.at(i).at(j) + i + 1;
                keepo[j] = keepers[oldirow[j]];
                haveDangler |= !(keepo[j]);
            }
            
            if(dangleWarn && haveDangler)
                std::cerr << "some danglers in excl list of old atom whose id = " << i << std::endl;
            
            vector<size_t> irow = mmpbsa_utils::take(newidx,oldirow,keepo);
            size_t shiftNewValues = returnMe.exclst.size()+1;
            for(size_t j = 0;j<irow.size();j++)
                irow[j] -= shiftNewValues;
            returnMe.exclst.push_back(irow);
        }

        
    //convert energy value tables and their indices to correspond to new atom indices
    valarray<slice> bond_slices(3);
    valarray<const valarray<mmpbsa_t>* > bond_values(2);
    
    //convert H-bonds
    bond_slices[0] = bond_h_i;bond_slices[1] = bond_h_j;bond_slices[2] = bond_h_k;
    EmpEnerFun::internalConvert(sp->bonds_inc_hydrogen,parminfo->bonds_inc_hydrogen,
            bond_slices,newidx,keepers,false);
    returnMe.bond_h_i = slice(0,sp->bonds_inc_hydrogen.size()/3,3);
    returnMe.bond_h_j= slice(1,sp->bonds_inc_hydrogen.size()/3,3);
    returnMe.bond_h_k = slice(2,sp->bonds_inc_hydrogen.size()/3,3);

    //convert non-H-bonds
    bond_slices[0] = bond_i;bond_slices[1] = bond_j;bond_slices[2] = bond_k;
    internalConvert(sp->bonds_without_hydrogen,parminfo->bonds_without_hydrogen,bond_slices,newidx,keepers,false);
    returnMe.bond_i = slice(0,sp->bonds_without_hydrogen.size()/3,3);
    returnMe.bond_j= slice(1,sp->bonds_without_hydrogen.size()/3,3);
    returnMe.bond_k = slice(2,sp->bonds_without_hydrogen.size()/3,3);

    //convert non-H-angles
    bond_slices.resize(4);
    bond_slices[0] = angle_i;bond_slices[1] = angle_j;bond_slices[2] = angle_k;bond_slices[3] = angle_l;
    internalConvert(sp->angles_without_hydrogen,parminfo->angles_without_hydrogen,bond_slices,newidx,keepers,false);
    returnMe.angle_i = slice(0,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_j = slice(1,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_k = slice(2,sp->angles_without_hydrogen.size()/4,4);
    returnMe.angle_l = slice(3,sp->angles_without_hydrogen.size()/4,4);

    //convert H-angles
    bond_slices[0] = angle_h_i;bond_slices[1] = angle_h_j;bond_slices[2] = angle_h_k;bond_slices[3] = angle_h_l;
    internalConvert(sp->angles_inc_hydrogen,parminfo->angles_inc_hydrogen,bond_slices,newidx,keepers,false);
    returnMe.angle_h_i = slice(0,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_j = slice(1,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_k = slice(2,sp->angles_inc_hydrogen.size()/4,4);
    returnMe.angle_h_l = slice(3,sp->angles_inc_hydrogen.size()/4,4);


    //convert non-H-dihedrals
    bond_slices.resize(5);
    bond_slices[0] = phi_i;bond_slices[1] = phi_j;bond_slices[2] = phi_k;bond_slices[3] = phi_l;bond_slices[4] = phi_m;
    internalConvert(sp->dihedrals_without_hydrogen,parminfo->dihedrals_without_hydrogen,bond_slices,newidx,keepers,false);
    returnMe.phi_i = slice(0,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_j = slice(1,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_k = slice(2,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_l = slice(3,sp->dihedrals_without_hydrogen.size()/5,5);
    returnMe.phi_m = slice(4,sp->dihedrals_without_hydrogen.size()/5,5);
    //update masks
    //While there is no need to update energy constants and equilibrium values,
    //the dihedral masks correspond to the locations of the dihedral elements in
    //dihedral_without_hydrogen. Therefore, the location of mask values must be
    //shuffled to new positions to correspond to the new positions in 
    //dihedrals_without_hydrogen
    updatePhiMasks(sp->dihedral_mask,parminfo->dihedral_mask,
        parminfo->dihedrals_without_hydrogen, bond_slices,
            keepers, false);
    
    //convert H-dihedrals
    bond_slices.resize(5);
    bond_slices[0] = phi_h_i;bond_slices[1] = phi_h_j;bond_slices[2] = phi_h_k;
    bond_slices[3] = phi_h_l;bond_slices[4] = phi_h_m;
    internalConvert(sp->dihedrals_inc_hydrogen,parminfo->dihedrals_inc_hydrogen,bond_slices,newidx,keepers,false);
    returnMe.phi_h_i = slice(0,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_j = slice(1,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_k = slice(2,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_l = slice(3,sp->dihedrals_inc_hydrogen.size()/5,5);
    returnMe.phi_h_m = slice(4,sp->dihedrals_inc_hydrogen.size()/5,5);
    //update masks
    updatePhiMasks(sp->dihedral_h_mask,parminfo->dihedral_h_mask,
        parminfo->dihedrals_inc_hydrogen, bond_slices,
            keepers, false);

    returnMe.LJA.resize(LJA.size());
    returnMe.LJA = LJA;
    returnMe.LJB.resize(LJB.size());
    returnMe.LJB = LJB;

    returnMe.parminfo = sp;//thus, don't delete sp.

    return returnMe;
}//end stripEnerFun(...)

//Energy Calculations
mmpbsa_t mmpbsa::EmpEnerFun::bond_inc_H(const std::valarray<mmpbsa_t>& crds)const
{
    return bond_energy_calc(crds,parminfo->bonds_inc_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::bond_without_H(const std::valarray<mmpbsa_t>& crds)const
{
    return bond_energy_calc(crds,parminfo->bonds_without_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::bond_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& bondIndices)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",INVALID_ARRAY_SIZE);
    
    size_t numBonds = size_t(bondIndices.size()/3);
    
    mmpbsa_t totalEnergy = 0;
    size_t bndi,bndj,bondid;
    mmpbsa_t ix,iy,iz,jx,jy,jz,bconst,beq,disp;//disp = displacement
    for(size_t i = 0;i<numBonds;i++)
    {
        bndi = bondIndices[3*i];bndj = bondIndices[3*i+1];
        bondid = bondIndices[3*i+2];
        ix = crds[3*bndi];iy = crds[3*bndi+1];iz = crds[3*bndi+2];
        jx = crds[3*bndj];jy = crds[3*bndj+1];jz = crds[3*bndj+2];
        bconst = parminfo->bond_force_constants[bondid];
        beq = parminfo->bond_equil_values[bondid];
        disp = sqrt(pow(ix-jx,2)+pow(iy-jy,2)+pow(iz-jz,2))-beq;
        totalEnergy += bconst*disp*disp;
    }

    return totalEnergy;

}

mmpbsa_t mmpbsa::EmpEnerFun::angle_inc_H(const std::valarray<mmpbsa_t>& crds)const
{
    return angle_energy_calc(crds,parminfo->angles_inc_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::angle_without_H(const std::valarray<mmpbsa_t>& crds)const
{
    return angle_energy_calc(crds,parminfo->angles_without_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::angle_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& angleIndices)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);
    
    size_t numAngles = size_t(angleIndices.size()/4);
    mmpbsa_t totalEnergy = 0;//total energy = \sum^N_m (const_m*(angle_m-eq_m)^2)
    size_t angle_id,atom_i,atom_j,atom_k;
    mmpbsa_t ix,iy,iz,jx,jy,jz,kx,ky,kz,r_ij,r_jk,dotprod,net_angle,angle_eq,angle_const;
    for(size_t i = 0;i<numAngles;i++)
    {
        atom_i = angleIndices[4*i];atom_j = angleIndices[4*i+1];atom_k = angleIndices[4*i+2];
        angle_id = angleIndices[4*i+3];
        ix = crds[3*atom_i];iy = crds[3*atom_i+1];iz = crds[3*atom_i+2];
        jx = crds[3*atom_j];jy = crds[3*atom_j+1];jz = crds[3*atom_j+2];
        kx = crds[3*atom_k];ky = crds[3*atom_k+1];kz = crds[3*atom_k+2];
        r_ij = sqrt(pow(ix-jx,2)+pow(iy-jy,2)+pow(iz-jz,2));
        r_jk = sqrt(pow(kx-jx,2)+pow(ky-jy,2)+pow(kz-jz,2));
        dotprod = (ix-jx)*(kx-jx)+(iy-jy)*(ky-jy)+(iz-jz)*(kz-jz);
        angle_eq = parminfo->angle_equil_values[angle_id];
        net_angle = acos(dotprod/(r_ij*r_jk)) - angle_eq;
        angle_const = parminfo->angle_force_constants[angle_id];
        totalEnergy += angle_const*net_angle*net_angle;
    }

    return totalEnergy;
}

mmpbsa_t mmpbsa::EmpEnerFun::dihedral_inc_H(const std::valarray<mmpbsa_t>& crds)const
{
    return dihedral_energy_calc(crds,parminfo->dihedrals_inc_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::dihedral_without_H(const std::valarray<mmpbsa_t>& crds)const
{
    return dihedral_energy_calc(crds,parminfo->dihedrals_without_hydrogen);
}

mmpbsa_t mmpbsa::EmpEnerFun::dihedral_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    using namespace mmpbsa_utils;
    mmpbsa_t * r_ij = new mmpbsa_t[3];//simply easier to put coordinates into an array
    mmpbsa_t * r_kj = new mmpbsa_t[3];//and use a separate method to calculate
    mmpbsa_t * r_kl = new mmpbsa_t[3];//cross product

    mmpbsa_t * d;//for use with cross products
    mmpbsa_t * g;//for use with cross products
    mmpbsa_t * s;
    mmpbsa_t totalEnergy,nphi,phi,ap0,ct1,gmag,dmag,dotprod,periodicity,dihedral_constant,phase;
    size_t d_i,d_j,d_l,d_k,d_id;
    size_t numDihedrals = size_t(dihedralIndices.size()/5);
    totalEnergy = 0;
    for(size_t i = 0;i<numDihedrals;i++)
    {
        d_i = dihedralIndices[5*i];d_j = dihedralIndices[5*i+1];d_k = dihedralIndices[5*i+2];d_l = dihedralIndices[5*i+3];
        d_id = dihedralIndices[5*i+4];
        for(size_t j = 0;j<3;j++)//fill vectors
        {
            r_ij[j] = crds[3*d_i+j]-crds[3*d_j+j];
            r_kj[j] = crds[3*d_k+j]-crds[3*d_j+j];
            r_kl[j] = crds[3*d_k+j]-crds[3*d_l+j];
        }
        d = cross_product(r_ij,r_kj,3);
        g = cross_product(r_kl,r_kj,3);
        dotprod = d[0]*g[0]+d[1]*g[1]+d[2]*g[2];
        dmag = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        gmag = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
        ct1 = dotprod/(dmag*gmag);
        
        if(std::abs(ct1) > 1.1)
        {
            std::ostringstream error;
            error << "dihedral routine fails on ";
            for(size_t concatIndices = 5*i;concatIndices<5*i+4;concatIndices++)
                    error << dihedralIndices[concatIndices] << ", ";
            error << dihedralIndices[5*i+4] << ": ct1 = " << ct1;
            throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
        }

        if(ct1 > 1)
        {
            ap0 = 0;
            std::cerr << "Warning: In dihedral";
            for(size_t concatIndices = 5*i;concatIndices < 5*i+4;concatIndices++)
                std::cerr << dihedralIndices[concatIndices] << " ";
            std::cerr << dihedralIndices[5*i+4]
                    <<" cosine comes to " << ct1
                    << ". By default, arccos is set to zero." << std::endl;
        }
        else if(ct1 < -1)
        {
            ap0 = MMPBSA_PI;
            std::cerr << "Warning: In dihedral"<< dihedralIndices[5*i]<< " " <<
                    dihedralIndices[5*i+1]<< " " << dihedralIndices[5*i+2] <<
                    " " << dihedralIndices[5*i+3] << " " << dihedralIndices[5*i+4] <<
                    " cos comes to " << ct1 << ", taking arccos as " << MMPBSA_PI << std::endl;
        }
        else
            ap0 = acos(ct1);

        s = cross_product(g,d,3);
        s[0] *= r_kj[0];s[1] *= r_kj[1];s[2] *= r_kj[2];
        if(s[0]+s[1]+s[2] < 0)
            phi = MMPBSA_PI+ap0;
        else
            phi = MMPBSA_PI-ap0;
        periodicity = std::abs(parminfo->dihedral_periodicities[d_id]);
        nphi = periodicity * phi;
        dihedral_constant = parminfo->dihedral_force_constants[d_id];
        phase = parminfo->dihedral_phases[d_id];
        totalEnergy += dihedral_constant * (1+cos(nphi)*cos(phase)+sin(nphi)*sin(phase));

        delete [] s;
        delete [] d;
        delete [] g;
    }

    delete [] r_ij;
    delete [] r_kj;
    delete [] r_kl;

    return totalEnergy;
}

mmpbsa_t mmpbsa::EmpEnerFun::vdw14_inc_H(const std::valarray<mmpbsa_t>& crds)const
{
    return vdw14_energy_calc(crds,parminfo->dihedrals_inc_hydrogen,parminfo->dihedral_h_mask);
}

mmpbsa_t mmpbsa::EmpEnerFun::vdw14_without_H(const std::valarray<mmpbsa_t>& crds)const
{
    return vdw14_energy_calc(crds,parminfo->dihedrals_without_hydrogen,parminfo->dihedral_mask);
}


mmpbsa_t mmpbsa::EmpEnerFun::vdw14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices,const std::valarray<bool>& phi_mask)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    const size_t& ntypes = parminfo->ntypes;
    size_t numDihedrals = size_t(dihedralIndices.size()/5);
    size_t d_i,d_l,d_m,type_i,type_l,flindex;
    mmpbsa_t rsqrd,a,b,inv_r6,inv_r12;
    mmpbsa_t totalEnergy = 0;
    bool imp_mask,igend_mask,period_mask,mask;
    for(size_t i = 0;i<numDihedrals;i++)
    {
        d_i = dihedralIndices[5*i];
        d_l = dihedralIndices[5*i+3];
        d_m = dihedralIndices[5*i+4];
        imp_mask = phi_mask[5*i+3];
        igend_mask = phi_mask[5*i+2];
        period_mask = parminfo->dihedral_periodicities[d_m] < 0;
        mask = imp_mask || igend_mask || period_mask;
        if (!mask)
        {
            rsqrd = 0;
            for (size_t j = 0; j < 3; j++)
                rsqrd += pow(crds[3*d_i+j] - crds[3*d_l+j], 2);
            type_i = parminfo->atom_type_indices[d_i] - 1;//There is an offset. See Amber 8 manual appendix C
            type_l = parminfo->atom_type_indices[d_l] - 1;
            flindex = type_i*ntypes + type_l; 
            a = LJA[flindex];
            b = LJB[flindex];
            inv_r6 = 1/pow(rsqrd, 3);
            inv_r12 = inv_r6*inv_r6;
            totalEnergy += a*inv_r12 - b*inv_r6;
            
        }
    }
    return totalEnergy*inv_scnb;

}

mmpbsa_t mmpbsa::EmpEnerFun::elstat14_inc_H(const std::valarray<mmpbsa_t>& crds)const
{
    return elstat14_energy_calc(crds,parminfo->dihedrals_inc_hydrogen,parminfo->dihedral_h_mask);
}

mmpbsa_t mmpbsa::EmpEnerFun::elstat14_without_H(const std::valarray<mmpbsa_t>& crds)const
{
    return elstat14_energy_calc(crds,parminfo->dihedrals_without_hydrogen,parminfo->dihedral_mask);
}

mmpbsa_t mmpbsa::EmpEnerFun::elstat14_energy_calc(const std::valarray<mmpbsa_t>& crds,
        const std::valarray<size_t>& dihedralIndices, const std::valarray<bool>& phi_mask)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);


    size_t numDihedrals = size_t(dihedralIndices.size()/5);
    size_t d_i,d_l;
    mmpbsa_t rsqrd,q_i,q_l;
    mmpbsa_t totalEnergy = 0;
    bool imp_mask,igend_mask,period_mask,mask;
    for(size_t i = 0;i<numDihedrals;i++)
    {
        imp_mask = phi_mask[5*i+3];
        igend_mask = phi_mask[5*i+2];
        period_mask = parminfo->dihedral_periodicities[dihedralIndices[5*i+4]] < 0;
        mask = imp_mask || igend_mask || period_mask;
        if(!mask)
        {
            rsqrd = 0;
            d_i = dihedralIndices[5*i];
            d_l = dihedralIndices[5*i+3];
            for(size_t j = 0; j < 3; j++)
                rsqrd += pow(crds[3*d_i+j] - crds[3*d_l+j], 2);
            q_i = parminfo->charges[d_i];
            q_l = parminfo->charges[d_l];
            totalEnergy += (inv_scee/dielc)*q_i*q_l/sqrt(rsqrd);//Ah, Coulomb's law :-)
        }
    }
    return totalEnergy;
}

mmpbsa_t mmpbsa::EmpEnerFun::total_vdwaals_energy(const std::valarray<mmpbsa_t>& crds)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    mmpbsa_t totalEnergy = 0;
    size_t natom = parminfo->natom;
    size_t ntypes = parminfo->ntypes;
    size_t type,type_2;
    mmpbsa_t x,y,z,rsqrd,a,b,atomEnergy;
    
    for(size_t i = 0;i<natom;i++)
    {
        x = crds[3*i];y = crds[3*i+1];z = crds[3*i+2];
        type = parminfo->atom_type_indices[i] - 1;//index is offset in parmtop file
        atomEnergy = 0;
        for(size_t j = i+1;j<natom;j++)//sum over all other atoms after the i-th atom
        {
            if(!mmpbsa_utils::contains(exclst.at(i),j-i-1))
            {
                type_2 = parminfo->atom_type_indices[j]-1;
                a = LJA[type*ntypes+type_2];b = LJB[type*ntypes+type_2];
                rsqrd = pow(x-crds[3*j],2) + pow(y-crds[3*j+1],2) + pow(z-crds[3*j+2],2);
                atomEnergy += a/pow(rsqrd,6) - b/pow(rsqrd,3);
            }
                
        }

        totalEnergy += atomEnergy;
    }
        return totalEnergy;
}


mmpbsa_t mmpbsa::EmpEnerFun::total_elstat_energy(const std::valarray<mmpbsa_t>& crds)const
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    mmpbsa_t totalEnergy = 0;
    size_t natom = parminfo->natom;
    mmpbsa_t q_i,q_j,x,y,z,r;
    for(size_t i = 0;i<natom;i++)
    {
        q_i = parminfo->charges[i];
        x = crds[3*i];y = crds[3*i+1];z = crds[3*i+2];

        for(size_t j = i+1;j<natom;j++)//sum over all other atoms after the i-th atom
        {
            if(!mmpbsa_utils::contains(exclst.at(i),j-i-1))
            {
                q_j = parminfo->charges[j];
                r = sqrt( pow(x-crds[3*j],2) + pow(y-crds[3*j+1],2) + pow(z-crds[3*j+2],2) );
                totalEnergy += q_i*q_j/r;
            }
        }
    }

    return totalEnergy;//Amber's charge units give kcal/mol without extra factor.
}


template <class M> void mmpbsa::EmpEnerFun::internalConvert(
    std::valarray<M>& newIndices,
        const std::valarray<M>& oldIndices,const std::valarray<std::slice>& slices,
        const std::valarray<size_t>& newidx,
        const std::valarray<bool>& keepers, const bool& dangleWarn)
{
    using std::vector;
    using std::valarray;

    if(slices.size() == 0)
        throw mmpbsa::MMPBSAException("At least one slice is needed for "
                "EmpEnerFun::internalConvert",mmpbsa::UNKNOWN_ERROR);

    size_t sliceSize = slices[0].size();
    for(size_t i = 1;i<slices.size();i++)
        if(sliceSize != slices[i].size())
            throw mmpbsa::MMPBSAException("Slices used by EmpEnerFun::internalConvert "
                    "must have the same size",mmpbsa::UNKNOWN_ERROR);
    size_t numKeptItems = 0;

    //An atom is dangling if it is not kept, but other are.
    //These tests are to see if one atom is not kept, are any others?
    for(size_t i = 0;i<oldIndices.size();i+=slices[0].stride())
    {
        size_t currAtom = oldIndices[i];
        bool dangleAlways = !keepers[currAtom];
        bool dangleOnce = !keepers[currAtom];
        for(size_t j = 1;j<slices[0].stride()-1;j++)
        {
            currAtom = oldIndices[i+j];
            dangleOnce |= !keepers[currAtom];
            dangleAlways &= !keepers[currAtom];
                
        }
        if(dangleOnce && !dangleAlways && dangleWarn)
        {
            std::cerr << "Dangler found at " << i << std::endl;
        }

        if(!dangleAlways)
            numKeptItems++;
        
    }

    //resize newIndices if necessary
    if(newIndices.size() != numKeptItems)
        newIndices.resize(numKeptItems*slices.size());

    //fill newIndices with the NEW INDEX of the kept atoms
    size_t newIndicesIndex = 0;
    for(size_t i =0;i<oldIndices.size();i+=slices[0].stride())
    {
        size_t currAtom = oldIndices[i];
        bool dangleAlways = !keepers[currAtom];
        bool dangleOnce = !keepers[currAtom];
        for(size_t j = 1;j<slices[0].stride()-1;j++)
        {
            currAtom = oldIndices[i+j];
            dangleOnce |= !keepers[currAtom];
            dangleAlways &= !keepers[currAtom];

        }

        if(!dangleAlways)
        {
            for(size_t j = 0;j<slices[0].stride() - 1;j++)
                newIndices[newIndicesIndex++] = newidx[oldIndices[i+j]];

            //last element should remain unchanged, since the bond data value
            //tables are unchanged as well
            newIndices[newIndicesIndex++] = oldIndices[i+slices[0].stride() - 1];
        }
    }


}//end internalConvert

template <class M> void mmpbsa::EmpEnerFun::updatePhiMasks(std::valarray<bool>& newMask,
    const std::valarray<bool>& oldMask,
    const std::valarray<M>& oldIndices, const std::valarray<std::slice>& slices,
        const std::valarray<bool>& keepers, const bool& dangleWarn)
{
    using std::vector;
    using std::valarray;

    if(slices.size() == 0)
        throw mmpbsa::MMPBSAException("At least one slice is needed for "
                "EmpEnerFun::internalConvert",mmpbsa::UNKNOWN_ERROR);

    size_t sliceSize = slices[0].size();
    for(size_t i = 1;i<slices.size();i++)
        if(sliceSize != slices[i].size())
            throw mmpbsa::MMPBSAException("Slices used by EmpEnerFun::internalConvert "
                    "must have the same size",mmpbsa::UNKNOWN_ERROR);
    size_t numKeptItems = 0;

    //An atom is dangling if it is not kept, but other are.
    //These tests are to see if one atom is not kept, are any others?
    for(size_t i = 0;i<oldIndices.size();i+=slices[0].stride())
    {
        size_t currAtom = oldIndices[i];
        bool dangleAlways = !keepers[currAtom];
        bool dangleOnce = !keepers[currAtom];
        for(size_t j = 1;j<slices[0].stride()-1;j++)
        {
            currAtom = oldIndices[i+j];
            dangleOnce |= !keepers[currAtom];
            dangleAlways &= !keepers[currAtom];

        }
        if(dangleOnce && !dangleAlways && dangleWarn)
        {
            std::cerr << "Dangler found at " << i << std::endl;
        }

        if(!dangleAlways)
            numKeptItems++;

    }

    //resize newIndices if necessary
    if(newMask.size() != numKeptItems*slices.size())
        newMask.resize(numKeptItems*slices.size());

    //fill newIndices with the NEW INDEX of the kept atoms
    size_t newMaskIndex = 0;
    for(size_t i =0;i<oldIndices.size();i+=slices[0].stride())
    {
        size_t currAtom = oldIndices[i];
        bool dangleAlways = !keepers[currAtom];
        bool dangleOnce = !keepers[currAtom];
        for(size_t j = 1;j<slices[0].stride()-1;j++)
        {
            currAtom = oldIndices[i+j];
            dangleOnce |= !keepers[currAtom];
            dangleAlways &= !keepers[currAtom];

        }

        if(!dangleAlways)
        {
            for(size_t j = 0;j<slices[0].stride();j++)
                newMask[newMaskIndex++] = oldMask[i+j];
        }
    }
}//end updatePhiMasks


std::string pdbPad(const int& neededDigits,const int& currentNumber)
{
	std::string returnMe = "";
	for(size_t serialPad = 0;serialPad < neededDigits - floor(log(currentNumber)/log(10.0));serialPad++)
					returnMe += " ";
	return returnMe;
}

std::ostream& streamPDB(std::ostream& theStream, const mmpbsa::EmpEnerFun& energy, const std::valarray<mmpbsa_t>& crds) throw (mmpbsa::MMPBSAException)
{
	using namespace mmpbsa;
	if(!theStream.good())
		throw MMPBSAException("Could not write to the stream provided to streamPDB",FILE_IO_ERROR);

	const mmpbsa::SanderParm* parm = energy.parminfo;

	//get molecule starting and stopping points.
	const std::valarray<size_t>& mol_ranges = energy.mol_ranges;

	//make sure mol_ranges has both starting and stopping indices.
	if(mol_ranges.size() % 2 != 0 || mol_ranges.size() == 0)
	{
		std::ostringstream error;
		error << "mol_ranges supplied to mmpbsa::streamPDB is corrupt."
				<< "Expected an array with an even number of molecule pointers,"
				<< " but was given an array of size " << mol_ranges.size();
		throw MMPBSAException(error,DATA_FORMAT_ERROR);
	}

	//iterate through list. output using PDB format.
	size_t first,last,currAtom,currRes;
	bool haveFormatWarned = false;
	size_t serialNumber = 1;
	for(size_t i = 0;i<mol_ranges.size();i+=2)
	{
		first = mol_ranges[i];
		last = mol_ranges[i+1];
		for(currAtom = first;currAtom<last;currAtom++)
		{
			//record name (cols 1-6)
			theStream << "ATOM  ";

			//serial number (cols 7-11)
			theStream << pdbPad(4,serialNumber);
			theStream << serialNumber++;

			//empty space in column 12
			theStream << " ";

			//atom name (from parm top file) (cols 13-16)
			theStream << parm->atom_names[currAtom];

			//alt location indicator (col 17)
			theStream << " ";

			//Residue name(cols 18-20)
			std::string resName = mmpbsa::getResidueLabel(currAtom,energy.res_ranges,parm,&currRes);
			theStream << resName.substr(0,3);

			//empty space column 21 chain id column 22
			theStream << "  ";

			//residue sequence (cols 23 - 26)
			theStream << pdbPad(4,currRes);
			theStream << currRes;

			//icode (col 27) with columns 28-30 empty
			theStream << "    ";

			//Coordinates (cols 31 - 54)
			std::ostringstream coordinateBuffer;
			coordinateBuffer.setf(std::ios::fixed,std::ios_base::floatfield);
			coordinateBuffer.precision(3);
			for(size_t coordIndex = 0;coordIndex < 3;coordIndex++)
			{
				coordinateBuffer.width(8);
				coordinateBuffer << crds[3*currAtom + coordIndex];
				theStream << coordinateBuffer.str();
				coordinateBuffer.str("");
			}

			//occupancy & temp factor
			theStream << "  1.00" << "  0.00";

			theStream << std::endl;
		}
	}

	return theStream;
}


mmpbsa::BondWalker::BondWalker(mmpbsa::EmpEnerFun const * efun)
{
    enerfun = efun;
    initialized = false;
}

void mmpbsa::BondWalker::init()
{
    using std::valarray;
    using std::slice;
    using std::slice_array;

    mmpbsa::SanderParm const * const parminfo = enerfun->parminfo;
    atom_bond_list.resize(parminfo->natom);
    size_t numBondPairs = size_t( parminfo->bonds_without_hydrogen.size()/3
        + parminfo->bonds_inc_hydrogen.size()/3 );//divide by 3 because each bond has an i,j pair plus force info

    //multiply by 2 because this contains pairs of atoms in each bond.
    //This array will have the form (bond, bond, ...) or ((i,j), (i,j), ...)
    valarray<size_t> bondPairs(2*numBondPairs);
    size_t withHydrogenOffset = size_t( parminfo->bonds_inc_hydrogen.size()/3);//first indicates the number bonds.
    for(size_t i = 0;i<withHydrogenOffset;i++)
    {
        bondPairs[2*i] = parminfo->bonds_inc_hydrogen[3*i];
        bondPairs[2*i+1] = parminfo->bonds_inc_hydrogen[3*i+1];
    }
    withHydrogenOffset *= 2;//now it indicates how far to skip to the without hydrogen section.
    for(size_t i = 0;i<size_t( parminfo->bonds_without_hydrogen.size()/3);i++)
    {
        bondPairs[2*i+withHydrogenOffset] = parminfo->bonds_without_hydrogen[3*i];
        bondPairs[2*i+1+withHydrogenOffset] = parminfo->bonds_without_hydrogen[3*i+1];
    }
    if(bondPairs.max() > parminfo->natom)
    {
        std::ostringstream error;
        error << "An atom index of " << bondPairs.max() << " in the bonds pairs exceed the "
                "maximum value of " << parminfo->natom << ".";
        std::cerr << "Max inc: " << parminfo->bonds_inc_hydrogen.max()
                << ", max without: " << parminfo->bonds_without_hydrogen.max() << std::endl;
        throw MMPBSAException(error,DATA_FORMAT_ERROR);
    }

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

std::vector<size_t> mmpbsa::BondWalker::walk(const size_t& startAtom,
        std::valarray<int> stoppers,std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom))
{
    if(initialized)
        visited != visited;
    else
        init();

    for(size_t i = 0;i<stoppers.size();i++)
        visited[stoppers[i]] = 1;

    std::vector<size_t> returnMe;
    recwalk(returnMe,startAtom,enerfunMarkers,func);

    return returnMe;
}

void mmpbsa::BondWalker::recwalk(std::vector<size_t>& list, const size_t& atom,
        std::valarray<int>& enerfunMarkers,
        void func(std::valarray<int>& markers, const size_t& funcAtom))
{
    func(enerfunMarkers,atom);
    visited[atom] = 1;
    for(size_t i = 0;i<atom_bond_list[atom].size();i++)
    {
        size_t bondedAtom = atom_bond_list[atom][i];
        if(visited[bondedAtom])
            continue;
        recwalk(list,bondedAtom,enerfunMarkers,func);
    }

}





