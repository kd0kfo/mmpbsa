#include "Energy.h"

//Energy Calculations

mmpbsa_t mmpbsa::bond_energy_calc(const std::vector<bond_t>& bonds,const std::valarray<mmpbsa_t>& crds)
{
	using namespace mmpbsa;
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("bond_energy_calc: Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",INVALID_ARRAY_SIZE);

    mmpbsa_t totalEnergy = 0;
    size_t bndi,bndj,bondid;
    mmpbsa_t ix,iy,iz,jx,jy,jz,disp;//disp = displacement
    std::vector<bond_t>::const_iterator bond;
    for(bond = bonds.begin();bond != bonds.end();bond++)
    {
        ix = crds[3*bond->atom_i];iy = crds[3*bond->atom_i+1];iz = crds[3*bond->atom_i+2];
        jx = crds[3*bond->atom_j];jy = crds[3*bond->atom_j+1];jz = crds[3*bond->atom_j+2];
        const mmpbsa_t& bconst = bond->bond_energy->energy_const;
        const mmpbsa_t& beq = bond->bond_energy->eq_distance;
        disp = sqrt(pow(ix-jx,2)+pow(iy-jy,2)+pow(iz-jz,2))-beq;
        totalEnergy += bconst*disp*disp;
    }

    return totalEnergy;

}

mmpbsa_t mmpbsa::angle_energy_calc(const std::vector<angle_t>& angles,
		const std::valarray<mmpbsa_t>& crds)
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("mmpbsa::angle_energy_calc: Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    mmpbsa_t totalEnergy = 0;//total energy = \sum^N_m (const_m*(angle_m-eq_m)^2)
    mmpbsa_t ix,iy,iz,jx,jy,jz,kx,ky,kz,r_ij,r_jk,dotprod,net_angle;
    std::vector<angle_t>::const_iterator angle;
    for(angle = angles.begin();angle != angles.end();angle++)
    {
    	const mmpbsa_t& angle_eq = angle->angle_energy->eq_distance;
    	const mmpbsa_t& angle_const = angle->angle_energy->energy_const;
    	ix = crds[3*angle->atom_i];iy = crds[3*angle->atom_i+1];iz = crds[3*angle->atom_i+2];
        jx = crds[3*angle->atom_j];jy = crds[3*angle->atom_j+1];jz = crds[3*angle->atom_j+2];
        kx = crds[3*angle->atom_k];ky = crds[3*angle->atom_k+1];kz = crds[3*angle->atom_k+2];
        r_ij = sqrt(pow(ix-jx,2)+pow(iy-jy,2)+pow(iz-jz,2));
        r_jk = sqrt(pow(kx-jx,2)+pow(ky-jy,2)+pow(kz-jz,2));
        dotprod = (ix-jx)*(kx-jx)+(iy-jy)*(ky-jy)+(iz-jz)*(kz-jz);
        net_angle = acos(dotprod/(r_ij*r_jk)) - angle_eq;
        totalEnergy += angle_const*net_angle*net_angle;
    }

    return totalEnergy;
}

mmpbsa_t mmpbsa::dihedral_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::valarray<mmpbsa_t>& crds)
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("mmpbsa::dihedral_energy_calc: Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    using namespace mmpbsa_utils;
    mmpbsa_t r_ij[3];//simply easier to put coordinates into an array
    mmpbsa_t r_kj[3];//and use a separate method to calculate
    mmpbsa_t r_kl[3];//cross product

    mmpbsa_t *d;//for use with cross products
    mmpbsa_t *g;//for use with cross products
    mmpbsa_t *s;
    mmpbsa_t totalEnergy,nphi,phi,ap0,ct1,gmag,dmag,dotprod,periodicity,dihedral_constant,phase;

    totalEnergy = 0;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)
    {
        for(size_t j = 0;j<3;j++)//fill vectors
        {
            r_ij[j] = crds[3*dihedral->atom_i+j]-crds[3*dihedral->atom_j+j];
            r_kj[j] = crds[3*dihedral->atom_k+j]-crds[3*dihedral->atom_j+j];
            r_kl[j] = crds[3*dihedral->atom_k+j]-crds[3*dihedral->atom_l+j];
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
            error << dihedral->atom_i << ", " << dihedral->atom_j << ", " << dihedral->atom_k << ", " << dihedral->atom_l;
            error << ": ct1 = " << ct1;
            throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
        }

        if(ct1 > 1)
        {
            ap0 = 0;
            std::cerr << "Warning: In dihedral";
            std::cerr << dihedral->atom_i << ", " << dihedral->atom_j << ", " << dihedral->atom_k << ", " << dihedral->atom_l;
            std::cerr <<" cosine comes to " << ct1
                    << ". By default, arccos is set to zero." << std::endl;
        }
        else if(ct1 < -1)
        {
            ap0 = MMPBSA_PI;
            std::cerr << "Warning: In dihedral";
            std::cerr << dihedral->atom_i << ", " << dihedral->atom_j << ", " << dihedral->atom_k << ", " << dihedral->atom_l;
            std::cerr << " cos comes to " << ct1 << ", taking arccos as " << MMPBSA_PI << std::endl;
        }
        else
            ap0 = acos(ct1);

        s = cross_product(g,d,3);
        s[0] *= r_kj[0];s[1] *= r_kj[1];s[2] *= r_kj[2];
        if(s[0]+s[1]+s[2] < 0)
            phi = MMPBSA_PI+ap0;
        else
            phi = MMPBSA_PI-ap0;
        nphi = dihedral->dihedral_energy->periodicity * phi;
        const mmpbsa_t& dihedral_constant = dihedral->dihedral_energy->energy_const;
        const mmpbsa_t& phase = dihedral->dihedral_energy->phase;
        totalEnergy += dihedral_constant * (1+cos(nphi)*cos(phase)+sin(nphi)*sin(phase));

        delete [] d;
        delete [] g;
        delete [] s;

    }

    return totalEnergy;
}

mmpbsa_t mmpbsa::vdw14_energy_calc(const std::vector<dihedral_t>& dihedrals,const std::valarray<mmpbsa_t>& crds,const mmpbsa_t& inv_scnb)
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("mmpbsa::vdw14_energy_calc: Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

    mmpbsa_t rsqrd,a,b,inv_r6,inv_r12;
    mmpbsa_t totalEnergy = 0;
    bool period_mask,mask;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)//for(size_t i = 0;i<numDihedrals;i++)
    {
        period_mask = (dihedral->dihedral_energy->periodicity < 0);
        mask = dihedral->nonbonded_masks.is_improper || dihedral->nonbonded_masks.should_ignore_end_grp || period_mask;
        if(!mask)
        {
            rsqrd = 0;
            for (size_t j = 0; j < 3; j++)
                rsqrd += pow(crds[3*dihedral->atom_i+j] - crds[3*dihedral->atom_l+j], 2);
            inv_r6 = 1/pow(rsqrd, 3);
            inv_r12 = inv_r6*inv_r6;
            totalEnergy += dihedral->lj.c12*inv_r12 - dihedral->lj.c6*inv_r6;
        }
    }
    return totalEnergy*inv_scnb;

}

mmpbsa_t mmpbsa::elstat14_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::vector<atom_t>& atoms, const std::valarray<mmpbsa_t>& crds, const mmpbsa_t& inv_scee, const mmpbsa_t& dielc)
{
    if(crds.size() % 3 != 0)
        throw mmpbsa::MMPBSAException("mmpbsa::elstat14_energy_calc: Coordinate arrays must be a multiple of 3. "
                "bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);


    mmpbsa_t rsqrd,q_i,q_l;
    mmpbsa_t totalEnergy = 0;
    bool mask,period_mask;

    std::vector<dihedral_t>::const_iterator dihedral;
    for(dihedral = dihedrals.begin();dihedral != dihedrals.end();dihedral++)
    {
    	period_mask = (dihedral->dihedral_energy->periodicity < 0);
    	mask = dihedral->nonbonded_masks.is_improper || dihedral->nonbonded_masks.should_ignore_end_grp || period_mask;
    	if(!mask)
        {
            rsqrd = 0;
            for(size_t j = 0; j < 3; j++)
                rsqrd += pow(crds[3*dihedral->atom_i+j] - crds[3*dihedral->atom_l+j], 2);
            q_i = atoms.at(dihedral->atom_i).charge;
            q_l = atoms.at(dihedral->atom_l).charge;
            totalEnergy += (inv_scee/dielc)*q_i*q_l/sqrt(rsqrd);//Ah, Coulomb's law :-)
        }
    }
    return totalEnergy;
}

mmpbsa_t mmpbsa::vdwaals_energy(const std::vector<atom_t>& atoms, const std::vector<lj_params_t>& lj_params,const std::valarray<mmpbsa_t>& crds)
{
	if(crds.size() % 3 != 0)
		throw mmpbsa::MMPBSAException("Coordinate arrays must be a multiple of 3. "
				"bond_energy_calc was given one that was not.",mmpbsa::INVALID_ARRAY_SIZE);

	mmpbsa_t totalEnergy = 0;
	size_t natom,ntypes,type_row,type_2;
	mmpbsa_t x,y,z,rsqrd,a,b,atomEnergy;

	natom = atoms.size();
	ntypes = sqrt(lj_params.size());
	for(size_t i = 0;i<natom;i++)
	{
		const atom_t& atom = atoms.at(i);
		x = crds[3*i];y = crds[3*i+1];z = crds[3*i+2];
		atomEnergy = 0;
		type_row = atom.atom_type*ntypes;
		for(size_t j = i+1;j<natom;j++)//sum over all other atoms after the i-th atom
		{
			if(atom.exclusion_list.find(j) == atom.exclusion_list.end())
			{
				const lj_params_t& lj = lj_params.at(type_row + atoms.at(j).atom_type);
				rsqrd = pow(x-crds[3*j],2) + pow(y-crds[3*j+1],2) + pow(z-crds[3*j+2],2);
				atomEnergy += lj.c12/pow(rsqrd,6) - lj.c6/pow(rsqrd,3);
			}

		}

		totalEnergy += atomEnergy;
	}
	return totalEnergy;
}





