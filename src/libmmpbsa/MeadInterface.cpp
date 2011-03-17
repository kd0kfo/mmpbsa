#include "MeadInterface.h"

mmpbsa::MeadInterface::MeadInterface() {
    brad["N"] = 1.550;
    brad["H"] = 1.200;
    brad["C"] = 1.700;
    brad["O"] = 1.500;
    brad["P"] = 1.800;
    brad["S"] = 1.800;
    brad["FE"] = 1.300;
    brad["Na+"] = 1.200;
    brad["Cl-"] = 1.700;
    brad["MG"] = 1.180;
    istrength = 0;
    surf_tension =  0.00542;// kcal/mol/Ang^2
    surf_offset = 0.92;// kcal/mol
    multithread = 0;
    snap_list_offset = 0;
}

mmpbsa::MeadInterface::MeadInterface(const mmpbsa::MeadInterface& orig) {
    brad = orig.brad;
    istrength = orig.istrength;
    surf_offset = orig.surf_offset;
    surf_tension = orig.surf_tension;
    multithread = orig.multithread;
    snap_list_offset = orig.snap_list_offset;
}

mmpbsa::MeadInterface::~MeadInterface() {
    brad.clear();
}

FinDiffMethod mmpbsa::MeadInterface::createFDM(const std::valarray<mmpbsa::Vector>& complexCrds,
        const std::valarray<mmpbsa::Vector>& receptorCrds, const std::valarray<mmpbsa::Vector>& ligandCrds,
        const int& outbox_grid_dim, const mmpbsa_t& fine_grid_spacing) throw (mmpbsa::MeadException)
{
    using std::max;
    using mmpbsa::Vector;
    
    if(complexCrds.size() == 0 || receptorCrds.size() == 0 || ligandCrds.size() == 0)
        throw MeadException("Trivial coordinates supplied to createFDM must be 3-D.",DATA_FORMAT_ERROR);
    
    //Obtain Complex dimensions and size
    Vector comMax = complexCrds[0];
    Vector comMin = complexCrds[0];
    mmpbsa_t comSize[3];
    mmpbsa_t maxComSize;
    Vector geoCenter;
    for(size_t i = 1;i<complexCrds.size();i++)
    {
    	std::cout << "Coord: " << complexCrds[i] << std::endl;
        for(size_t j = 0;j<3;j++)
        {
            if(complexCrds[i].at(j) > comMax.at(j))
                comMax.at(j) = complexCrds[i].at(j);
            if(complexCrds[i].at(j) < comMin.at(j))
                comMin.at(j) = complexCrds[i].at(j);
        }
    }

    for(size_t i = 0;i<3;i++)
    	comSize[i] = comMax.at(i) - comMin.at(i);
    geoCenter = (comMax+comMin)/2;
    
    maxComSize = comSize[0];
    if(comSize[1] > maxComSize)
        maxComSize = comSize[1];
    if(comSize[2] > maxComSize)
        maxComSize = comSize[2];

    //initialize FinDiffMethod Object
    FinDiffMethod fdm;
    
    //determine outbox grid size and add level(s) accordingly
    mmpbsa_t outboxlen = maxComSize + 30.0;
    mmpbsa_t outbox_grid_spacing = outboxlen/(outbox_grid_dim - 1);
    fdm.add_level(outbox_grid_dim,outbox_grid_spacing,ON_GEOM_CENT);
    if(outbox_grid_dim > 1.0)
    {
        mmpbsa_t med_grid_spacing = outbox_grid_spacing/2;
        mmpbsa_t flen = maxComSize+10.0;
        mmpbsa_t fdim = flen/med_grid_spacing + med_grid_spacing;
        //now how to do rounding up, etc?
        int med_grid_dim = int(floor(fdim) + 1);
        if(med_grid_dim % 2 == 0)
            med_grid_dim++;

        fdm.add_level(int(med_grid_dim), med_grid_spacing, ON_GEOM_CENT);
    }

    //add fine grid level
    mead_data_t * int_minmax = mmpbsa_utils::interaction_minmax(receptorCrds,ligandCrds);//order {min,max}
    mead_data_t intSize[3];
    for(size_t i = 0;i<3;i++)
        intSize[i] = int_minmax[3+i] - int_minmax[i];
    mead_data_t maxIntSize = max(max(intSize[0],intSize[1]),intSize[2]);
    mead_data_t intCenter[3];
    for(size_t i = 0;i<3;i++)
        intCenter[i] = (int_minmax[3+i]+int_minmax[i])/2;
    mead_data_t fine_dim = maxIntSize/mead_data_t(fine_grid_spacing) + mead_data_t(fine_grid_spacing);
    int fine_grid_dim = int(floor(fine_dim) + 1);
    if(fine_grid_dim % 2 == 0)
        fine_grid_dim++;

    fdm.add_level(fine_grid_dim,fine_grid_spacing,ON_CENT_OF_INTR);
    fdm.resolve(Coord(geoCenter[0],geoCenter[1],geoCenter[2]),Coord(intCenter[0],intCenter[1],intCenter[2]));

    printf("Int min max: (%f, %f, %f) (%f, %f, %f)",int_minmax[0],int_minmax[1],int_minmax[2],int_minmax[3],int_minmax[4],int_minmax[5]);

    delete [] int_minmax;
    return fdm;
}

mmpbsa::EMap mmpbsa::MeadInterface::full_EMap(const mmpbsa::EmpEnerFun& efun, const std::valarray<mmpbsa::Vector>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mead_data_t>& radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (mmpbsa::MeadException)
{
	throw mmpbsa::MMPBSAException("mmpbsa::MeadInterface::full_EMap: Deprecated!");
    mmpbsa::EMap returnMe(&efun,crds);
    mmpbsa_t * pbsa_values = pbsa_solvation(efun,crds,fdm,radii,residueMap,interactionStrength);
    returnMe.set_elstat_solv(pbsa_values[0]);
    returnMe.set_area(pbsa_values[1]);
    returnMe.set_sasol(pbsa_values[1]*surfTension+surfOffset);
    delete [] pbsa_values;
    return returnMe;
}

mmpbsa::EMap mmpbsa::MeadInterface::full_EMap(const std::vector<mmpbsa::atom_t>& atoms, const mmpbsa::forcefield_t& ff, const std::valarray<mmpbsa::Vector>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mead_data_t>& radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (mmpbsa::MeadException)
{
    mmpbsa::EMap returnMe(atoms,ff,crds);
    mmpbsa_t * pbsa_values = pbsa_solvation(atoms,ff,crds,fdm,radii,residueMap,interactionStrength);
    returnMe.set_elstat_solv(pbsa_values[0]);
    returnMe.set_area(pbsa_values[1]);
    returnMe.set_sasol(pbsa_values[1]*surfTension+surfOffset);
    delete [] pbsa_values;
    return returnMe;
}

mmpbsa_t* mmpbsa::MeadInterface::pbsa_solvation(const mmpbsa::EmpEnerFun& efun, const std::valarray<mmpbsa::Vector>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mead_data_t>& radii,
        const std::map<std::string,std::string>& residueMap,
        const mmpbsa_t& interactionStrength, const mmpbsa_t& exclusionRadius) throw (mmpbsa::MeadException)
{
	mmpbsa_t* returnMe;
	std::vector<mmpbsa::atom_t> atoms;
	mmpbsa::forcefield_t ff;
	efun.extract_force_field(ff);
	efun.extract_atom_structs(atoms);

	returnMe = pbsa_solvation(atoms, ff,crds,fdm,radii,residueMap,interactionStrength,exclusionRadius);

	destroy(&ff);
	return returnMe;
}

mmpbsa_t* mmpbsa::MeadInterface::pbsa_solvation(const std::vector<mmpbsa::atom_t>& atoms, const mmpbsa::forcefield_t& ff,
		const std::valarray<mmpbsa::Vector>& crds,
		const FinDiffMethod& fdm, const std::map<std::string,mead_data_t>& radii,
		const std::map<std::string,std::string>& residueMap,
		const mmpbsa_t& interactionStrength, const mmpbsa_t& exclusionRadius) throw (mmpbsa::MeadException)
{
    mmpbsa_t * returnMe = new mmpbsa_t[2];enum{esol = 0, area};

    mmpbsa::MeadInterface MI;//stores bond info for atoms not in radii.
    //PB
    AtomSet atmSet;
    std::vector<mmpbsa::atom_t>::const_iterator atom;

    size_t i = 0;
    for(atom = atoms.begin();atom != atoms.end();atom++,i++)
    {
        Atom currAtom;
        currAtom.atname = atom->name;
        currAtom.resname = "FOO";//FIX ME???
        currAtom.resnum = i+1;//FIX ME???
        currAtom.coord = ToCoord(crds[i]);
        currAtom.charge = atom->charge;
        currAtom.rad = mmpbsa_utils::lookup_radius(currAtom.atname,radii);

        if(currAtom.rad < 0.1 || currAtom.rad > 3.0)
            std::cerr << "WARNING: strange radius, " << currAtom.rad << ", for atom "
                    << currAtom.atname << " (index = " << i << ")";
        atmSet.insert(currAtom);
    }
    //Solvent energy
    ChargeDist rho(new AtomChargeSet(atmSet));
    DielectricEnvironment eps(new TwoValueDielectricByAtoms(atmSet,1.0,80.0,1.4));
    ElectrolyteEnvironment ely( new ElectrolyteByAtoms(atmSet,interactionStrength,exclusionRadius));
    ElstatPot phi_solv(fdm, eps, rho, ely);
    phi_solv.solve();
    mmpbsa_t prod_sol = mmpbsa_t(phi_solv * rho);
    
    //Electrolyte
    ElectrolyteEnvironment ely_ref( new UniformElectrolyte(0.0));
    DielectricEnvironment  eps_ref( new UniformDielectric(1.0));
    ElstatPot phi_ref(fdm, eps_ref, rho, ely_ref);
    phi_ref.solve();
    mmpbsa_t prod_ref = mmpbsa_t(phi_ref * rho);

    returnMe[esol] = (prod_sol - prod_ref) / 2.0;

    atmSet.clear();//Free a little memory before gobbling a lot for molsurf (which hopefully will be replaced)
    size_t numCoords = crds.size();
    REAL_T xs[numCoords],ys[numCoords],zs[numCoords];
    REAL_T rads[numCoords];
    
    //Surface Area
    for(size_t i = 0;i<numCoords;i++)
    {
        xs[i] = crds[i].x();
        ys[i] = crds[i].y();
        zs[i] = crds[i].z();
        rads[i] = mmpbsa_utils::lookup_radius(atoms.at(i).name,MI.brad) + 1.4;//SA radii are not necessarily the same as PB radii
    }
    returnMe[area] = molsurf(xs,ys,zs,rads,numCoords,0);//replace with molsurf stuff42;//

    return returnMe;
}

Coord ToCoord(const mmpbsa::Vector& v)
{
	return Coord(v.x(),v.y(),v.z());
}


