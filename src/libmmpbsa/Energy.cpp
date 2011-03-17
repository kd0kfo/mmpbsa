#include "Energy.h"

//Energy Calculations
mmpbsa_t mmpbsa::bond_energy_calc(const std::vector<bond_t>& bonds,const std::valarray<mmpbsa::Vector>& crds)
{
	using namespace mmpbsa;

    mmpbsa_t totalEnergy = 0,distance;
    Vector disp;
    std::vector<bond_t>::const_iterator bond;
   for(bond = bonds.begin();bond != bonds.end();bond++)
    {
	   const Vector& c_i = crds[bond->atom_i];
	   const Vector& c_j = crds[bond->atom_j];
	   const mmpbsa_t& bconst = bond->bond_energy->energy_const;
	   const mmpbsa_t& beq = bond->bond_energy->eq_distance;
	   disp = c_i - c_j;
	   distance = disp.modulus()-beq;
	   totalEnergy += bconst*distance*distance;
    }
    return totalEnergy;

}

mmpbsa_t mmpbsa::angle_energy_calc(const std::vector<angle_t>& angles,
		const std::valarray<mmpbsa::Vector>& crds)
{
	using mmpbsa::Vector;

	mmpbsa_t totalEnergy = 0;//total energy = \sum^N_m (const_m*(angle_m-eq_m)^2)
    mmpbsa_t dotprod,net_angle;
    Vector r_ij,r_jk;
    std::vector<angle_t>::const_iterator angle;
    for(angle = angles.begin();angle != angles.end();angle++)
    {
    	const mmpbsa_t& angle_eq = angle->angle_energy->eq_distance;
    	const mmpbsa_t& angle_const = angle->angle_energy->energy_const;
    	const Vector& c_i = crds[angle->atom_i];
    	const Vector& c_j = crds[angle->atom_j];
    	const Vector& c_k = crds[angle->atom_k];
    	r_ij = c_i - c_j;r_ij /= r_ij.modulus();
    	r_jk = c_k - c_j;r_jk /= r_jk.modulus();
    	net_angle = acos(r_ij*r_jk) - angle_eq;
        totalEnergy += angle_const*net_angle*net_angle;
    }
    return totalEnergy;
}

mmpbsa_t mmpbsa::dihedral_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::valarray<mmpbsa::Vector>& crds)
{
	using namespace mmpbsa_utils;
    using mmpbsa::Vector;
    Vector r_ij, r_kj, r_kl,s;//Interatom vectors

    Vector d, g;//vectors normal to the dihedral planes
    mmpbsa_t totalEnergy,nphi,phi,ap0;

    totalEnergy = 0;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)
    {
    	if(dihedral->dihedral_energy == 0)//this could happen, for example, if there is only LJ14 interaction, but no dihedral.
    		continue;

        r_ij = crds[dihedral->atom_i]-crds[dihedral->atom_j];
        r_kj = crds[dihedral->atom_k]-crds[dihedral->atom_j];
        r_kl = crds[dihedral->atom_k]-crds[dihedral->atom_l];
        d = cross(r_ij,r_kj);
        g = cross(r_kl,r_kj);
        ap0 = dihedral_angle(d,g);

        if(isnan(ap0))
        {
            std::cerr << "Warning: For dihedral"
            		<< dihedral->atom_i << ", " << dihedral->atom_j << ", " << dihedral->atom_k << ", " << dihedral->atom_l
            		<< ", angle is undefinied. " << std::endl;
        }

        s = cross(g,d);
        for(size_t i = 0;i<3;i++)
        	s.at(i) *= r_kj.at(i);
        if(s[0]+s[1]+s[2] < 0)
            phi = MMPBSA_PI+ap0;
        else
            phi = MMPBSA_PI-ap0;
        nphi = dihedral->dihedral_energy->periodicity * phi;
        const mmpbsa_t& dihedral_constant = dihedral->dihedral_energy->energy_const;
        const mmpbsa_t& phase = dihedral->dihedral_energy->phase;
        totalEnergy += dihedral_constant * (1+cos(nphi)*cos(phase)+sin(nphi)*sin(phase));
    }

    return totalEnergy;
}

mmpbsa_t mmpbsa::vdw14_energy_calc(const std::vector<dihedral_t>& dihedrals,const std::valarray<mmpbsa::Vector>& crds,const mmpbsa_t& inv_scnb)
{
    mmpbsa_t rsqrd,inv_r6,inv_r12;
    mmpbsa_t totalEnergy = 0;
    bool period_mask,mask;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)//for(size_t i = 0;i<numDihedrals;i++)
    {
        period_mask = (dihedral->dihedral_energy != 0 && dihedral->dihedral_energy->periodicity < 0);
        mask = dihedral->nonbonded_masks.is_improper || dihedral->nonbonded_masks.should_ignore_end_grp || period_mask;
        if(!mask)
        {
            rsqrd = (crds[dihedral->atom_i]-crds[dihedral->atom_l])*(crds[dihedral->atom_i]-crds[dihedral->atom_l]);
            inv_r6 = 1/pow(rsqrd, 3);
            inv_r12 = inv_r6*inv_r6;
            totalEnergy += dihedral->lj.c12*inv_r12 - dihedral->lj.c6*inv_r6;
        }
    }
    return totalEnergy*inv_scnb;

}

mmpbsa_t mmpbsa::elstat14_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::vector<atom_t>& atoms, const std::valarray<mmpbsa::Vector>& crds, const mmpbsa_t& inv_scee, const mmpbsa_t& dielc)
{
    mmpbsa_t rsqrd,q_i,q_l;
    mmpbsa_t totalEnergy = 0;
    bool mask,period_mask;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)
    {
    	period_mask = (dihedral->dihedral_energy != 0 && dihedral->dihedral_energy->periodicity < 0);
    	mask = dihedral->nonbonded_masks.is_improper || dihedral->nonbonded_masks.should_ignore_end_grp || period_mask;
    	if(!mask)
        {
    		rsqrd = (crds[dihedral->atom_i]-crds[dihedral->atom_l])*(crds[dihedral->atom_i]-crds[dihedral->atom_l]);
    		q_i = atoms.at(dihedral->atom_i).charge;
    		q_l = atoms.at(dihedral->atom_l).charge;
    		totalEnergy += (inv_scee/dielc)*q_i*q_l/sqrt(rsqrd);//Ah, Coulomb's law :-)
        }
    }
    return totalEnergy;
}

mmpbsa_t mmpbsa::vdwaals_energy(const std::vector<atom_t>& atoms, const std::vector<lj_params_t>& lj_params,const std::valarray<mmpbsa::Vector>& crds)
{
	using mmpbsa::Vector;
	mmpbsa_t totalEnergy = 0;
	size_t natom,ntypes,type_row;
	mmpbsa_t rsqrd;

	if(floor(sqrt(lj_params.size())) != ceil(sqrt(lj_params.size())))
		throw mmpbsa::MMPBSAException("mmpbsa::vdwaals_energy: Lennard Jones parameter matrix must be a square matrix.",mmpbsa::DATA_FORMAT_ERROR);

 	natom = atoms.size();
	ntypes = floor(sqrt(lj_params.size()));
	for(size_t i = 0;i<natom;i++)
	{
		const atom_t& atom = atoms.at(i);
		const Vector& c_i = crds[i];
		type_row = atom.atom_type*ntypes;
		for(size_t j = i+1;j<natom;j++)//sum over all other atoms after the i-th atom
		{
			if(atom.exclusion_list.find(j) == atom.exclusion_list.end())
			{
				const lj_params_t& lj = lj_params.at(type_row + atoms.at(j).atom_type);
				rsqrd =(c_i-crds[j])*(c_i-crds[j]);
				totalEnergy += lj.c12/pow(rsqrd,6) - lj.c6/pow(rsqrd,3);
			}
		}
	}
	return totalEnergy;
}


mmpbsa_t mmpbsa::total_elstat_energy(const std::vector<mmpbsa::atom_t>& atoms, const std::valarray<mmpbsa::Vector>& crds, const mmpbsa_t& coulomb_const)
{
	using mmpbsa::Vector;
    mmpbsa_t totalEnergy = 0,atom_potential;
    mmpbsa_t r;

    size_t natom = atoms.size();
    for(size_t i = 0;i<natom;i++)
    {
    	const std::set<size_t>& exclusion_list = atoms.at(i).exclusion_list;
    	const Vector& c_i = crds[i];
    	atom_potential = 0;
    	for(size_t j = i+1;j<natom;j++)//sum over all other atoms after the i-th atom
    	{
    		if(exclusion_list.find(j) == exclusion_list.end())
    		{
    			const atom_t& atom = atoms.at(j);
    			r = (c_i-crds[j]).modulus();
    			atom_potential += atom.charge/r;
    		}
    	}
    	totalEnergy += atom_potential * atoms.at(i).charge;
    }
    return coulomb_const * totalEnergy;
}

