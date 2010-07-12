#include "MeadInterface.h"

MeadInterface::MeadInterface() {
}

MeadInterface::MeadInterface(const MeadInterface& orig) {
}

MeadInterface::~MeadInterface() {
}

FinDiffMethod MeadInterface::createFDM(const std::valarray<mmpbsa_t>& complexCrds,
        const std::valarray<mmpbsa_t>& receptorCrds, const std::valarray<mmpbsa_t>& ligandCrds,
        const int& outbox_grid_dim, const mmpbsa_t& fine_grid_spacing) throw (MMPBSAException)
{
    using std::max;
    
    if(complexCrds.size() == 0 || receptorCrds.size() == 0 || ligandCrds.size() == 0)
        throw MMPBSAException("Trivial coordinates supplied to createFDM must be 3-D.",DATA_FORMAT_ERROR);
    
    if(complexCrds.size() % 3 != 0 || receptorCrds.size() % 3 != 0 || ligandCrds.size() % 3 != 0)
        throw MMPBSAException("Coordinates supplied to createFDM must be 3-D.",DATA_FORMAT_ERROR);

    //Obtain Complex dimensions and size
    mmpbsa_t comMax[3] = {complexCrds[0],complexCrds[1],complexCrds[2]};
    mmpbsa_t comMin[3] = {complexCrds[0],complexCrds[1],complexCrds[2]};
    mmpbsa_t comSize[3];
    mmpbsa_t maxComSize;
    mmpbsa_t geoCenter[3];
    for(size_t i = 3;i<complexCrds.size();i+=3)
    {
        for(size_t j = 0;j<3;j++)
        {
            if(complexCrds[i+j] > comMax[j])
                comMax[j] = complexCrds[i+j];
            if(complexCrds[i+j] < comMin[j])
                comMin[j] = complexCrds[i+j];
        }
    }

    for(size_t i = 0;i<3;i++)
    {
        comSize[i] = comMax[i] - comMin[i];
        geoCenter[i] = (comMax[i]+comMin[i])/2;
    }
    
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
    Coord * int_minmax = mmpbsa_utils::interaction_minmax(receptorCrds,ligandCrds);//order {min,max}
    Coord intSize = int_minmax[1] - int_minmax[0];
    float maxIntSize = max(max(intSize.x,intSize.y),intSize.z);
    Coord intCenter = (int_minmax[1]+int_minmax[0])/2;
    float fine_dim = maxIntSize/float(fine_grid_spacing) + float(fine_grid_spacing);
    int fine_grid_dim = int(floor(fine_dim) + 1);
    if(fine_grid_dim % 2 == 0)
        fine_grid_dim++;

    fdm.add_level(fine_grid_dim,fine_grid_spacing,ON_CENT_OF_INTR);

//    /*unused in pyamber*/
//    Coord gCenter;
//    for(size_t i = 0;i<ligandCrds.size();i+=3)
//    {
//        gCenter.x += ligandCrds[i];
//        gCenter.y += ligandCrds[i+1];
//        gCenter.z += ligandCrds[i+2];
//    }

    fdm.resolve(Coord(geoCenter[0],geoCenter[1],geoCenter[2]),intCenter);
    return fdm;
}

EMap MeadInterface::full_EMap(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mmpbsa_t>* radii,
        const std::map<std::string,std::string>& residueMap,const mmpbsa_t& interactionStrength,
        const mmpbsa_t& surfTension, const mmpbsa_t& surfOffset) throw (MMPBSAException)
{
    EMap returnMe(&efun,crds);
    mmpbsa_t * pbsa_values = pbsa_solvation(efun,crds,fdm,radii,residueMap,interactionStrength);
    returnMe.set_elstat_solv(pbsa_values[0] * surfTension + surfOffset);
    returnMe.set_area(pbsa_values[1]);

    delete [] pbsa_values;
    return returnMe;
}

mmpbsa_t* MeadInterface::pbsa_solvation(const EmpEnerFun& efun, const std::valarray<mmpbsa_t>& crds,
        const FinDiffMethod& fdm, const std::map<std::string,mmpbsa_t>* radii,
        const std::map<std::string,std::string>& residueMap,
        const mmpbsa_t& interactionStrength, const mmpbsa_t& exclusionRadius) throw (MMPBSAException)
{
    if(crds.size() % 3 != 0)
        throw MMPBSAException("Coordinates supplied to pbsa_solvation must be "
                "3-D coordinates.", DATA_FORMAT_ERROR);

    mmpbsa_t * returnMe = new mmpbsa_t[2];size_t esol = 0;size_t area = 1;

    AtomSet atmSet;
    const mmpbsa_io::SanderParm * parminfo = efun.parminfo;
    for(size_t i = 0;i<parminfo->natom;i++)
    {
        Atom currAtom;
        currAtom.atname = parminfo->atom_names[i];
        currAtom.resname = "DRC";//FIX ME???
        currAtom.resnum = i+1;//FIX ME???
        currAtom.coord = Coord(crds[3*i],crds[3*i+1],crds[3*i+2]);
        currAtom.charge = parminfo->charges[i];
        if(radii)

            currAtom.rad = mmpbsa_utils::lookup_radius(currAtom.atname,(*radii));
        else
            currAtom.rad = bondi_lookup(currAtom.atname);

        if(currAtom.rad < 0.1 || currAtom.rad > 3.0)
            fprintf(stderr,"WARNING: strange radius, %f, for atom %s (index = %d)", currAtom.rad, currAtom.atname.c_str(),i);
        atmSet.insert(currAtom);
    }
    AtomChargeSet rho(atmSet);
    TwoValueDielectricByAtoms eps(atmSet,1.0,80.0,1.4);
    ElectrolyteByAtoms ely(atmSet,interactionStrength,exclusionRadius);
    FinDiffElstatPot phi_solv(fdm, &eps, &rho, &ely);
    phi_solv.solve();
    mmpbsa_t prod_sol = mmpbsa_t(rho * phi_solv);
    UniformElectrolyte ely_ref(0.0);
    UniformDielectric eps_ref(1.0);
    FinDiffElstatPot phi_ref(fdm, &eps_ref, &rho, &ely_ref);
    phi_ref.solve();
    mmpbsa_t prod_ref = mmpbsa_t(rho * phi_ref);
    //# esol will already be in kcal/mole because of Amber's charge units
    returnMe[esol] = (prod_sol - prod_ref) / 2.0;

    returnMe[area] = MMPBSA_PI;

    return returnMe;
}

mmpbsa_t MeadInterface::bondi_lookup(const std::string& atomName)
{
    //map used for clarity in reading the code.
    std::map<std::string,mmpbsa_t> brad;
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

    if(brad.find(mmpbsa_utils::trimString(atomName)) != brad.end())
        return brad[atomName];
    if(brad.find(atomName.substr(0,1)) != brad.end())
        return brad[atomName.substr(0,1)];

    fprintf(stderr,"WARNING: Could not find a bondi radius for \"%s\". "
            "Using zero instead.",atomName.c_str());

    return 0;
}


