#include "SanderIO.h"

sanderio::SanderParm::SanderParm() {
    natom = 0;// total number of atoms
    ntypes = 0;//< total number of distinct atom types
    nbonh = 0;// number of bonds containing hydrogen
    mbona = 0;// number of bonds not containing hydrogen
    ntheth = 0;// number of angles containing hydrogen
    mtheta = 0;// number of angles not containing hydrogen
    nphih = 0;// number of dihedrals containing hydrogen
    mphia = 0;// number of dihedrals not containing hydrogen
    nhparm = 0;// currently not used
    nparm = 0;// set to 1 if LES is used
    nnb = 0;// number of excluded atoms (=NEXT)
    nres = 0;// number of residues
    nbona = 0;// MBONA + number of constraint bonds
    ntheta = 0;// MTHETA + number of constraint angles
    nphia = 0;// MPHIA + number of constraint dihedrals
    numbnd = 0;// number of unique bond types
    numang = 0;// number of unique angle types
    nptra = 0;// number of unique dihedral types
    natyp = 0;// number of atom types in parameter file, see SOLTY below
    nphb = 0;// number of distinct 10-12 hydrogen bond pair types
    ifpert = 0;// set to 1 if perturbation info is to be read in
    nbper = 0;// number of bonds to be perturbed
    ngper = 0;// number of angles to be perturbed
    ndper = 0;//number of dihedrals to be perturbed
    mbper = 0;// number of bonds with atoms completely in perturbed group
    mgper = 0;// number of angles with atoms completely in perturbed group
    mdper = 0;// number of dihedrals with atoms completely in perturbed groups
    ifbox = 0;// set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = 0;// number of atoms in the largest residue
    ifcap = 0;// set to 1 if the CAP option from edit was specified
    nextra = 0;// number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    iptres = 0;//   last residue that is considered part of solute (base 1 index)
    nspm = 0;//     total number of molecules
    nspsol = 0;//   the first solvent "molecule" (base 1 index)


    //pointers are stored here in case the vectors are swapped around later
    //between difference classes. These arrays are LARGE; the goal here is to
    //save memory/time if that happens.
    titles = "";
    atom_names;
    charges;
    masses;
    atom_type_indices;
    number_excluded_atoms;
    nonbonded_parm_indices;
    residue_labels;
    residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants;
    bond_equil_values;
    angle_force_constants;
    dihedral_force_constants;
    angle_equil_values;
    dihedral_periodicities;
    dihedral_phases;
    soltys;//solubility?
    lennard_jones_acoefs;
    lennard_jones_bcoefs;
    bonds_inc_hydrogen;
    bonds_without_hydrogen;
    angles_inc_hydrogen;
    angles_without_hydrogen;
    dihedrals_inc_hydrogen;
    dihedrals_without_hydrogen;
    excluded_atoms_list;
    hbond_acoefs;
    hbond_bcoefs;
    hbcuts;
    amber_atom_types;
    tree_chain_classifications;
    join_array;
    irotats;
    radius_sets;
    radii;
    screen;/**/
}

sanderio::SanderParm::SanderParm(const SanderParm& orig) {
    natom = orig.natom;// total number of atoms
    ntypes = orig.ntypes;// total number of distinct atom types
    nbonh = orig.nbonh;// number of bonds containing hydrogen
    mbona = orig.mbona;// number of bonds not containing hydrogen
    ntheth = orig.ntheth;// number of angles containing hydrogen
    mtheta = orig.mtheta;// number of angles not containing hydrogen
    nphih = orig.nphih;// number of dihedrals containing hydrogen
    mphia = orig.mphia;// number of dihedrals not containing hydrogen
    nhparm = orig.nhparm;// currently not used
    nparm = orig.nparm;// set to 1 if LES is used
    nnb = orig.nnb;// number of excluded atoms (=NEXT)
    nres = orig.nres;// number of residues
    nbona = orig.nbona;// MBONA + number of constraint bonds
    ntheta = orig.ntheta;// MTHETA + number of constraint angles
    nphia = orig.nphia;// MPHIA + number of constraint dihedrals
    numbnd = orig.numbnd;// number of unique bond types
    numang = orig.numang;// number of unique angle types
    nptra = orig.nptra;// number of unique dihedral types
    natyp = orig.natyp;// number of atom types in parameter file, see SOLTY below
    nphb = orig.nphb;// number of distinct 10-12 hydrogen bond pair types
    ifpert = orig.ifpert;// set to 1 if perturbation info is to be read in
    nbper = orig.nbper;// number of bonds to be perturbed
    ngper = orig.ngper;// number of angles to be perturbed
    ndper = orig.ndper;//number of dihedrals to be perturbed
    mbper = orig.mbper;// number of bonds with atoms completely in perturbed group
    mgper = orig.mgper;// number of angles with atoms completely in perturbed group
    mdper = orig.mdper;// number of dihedrals with atoms completely in perturbed groups
    ifbox = orig.ifbox;// set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = orig.nmxrs;// number of atoms in the largest residue
    ifcap = orig.ifcap;// set to 1 if the CAP option from edit was specified
    nextra = orig.nextra;// number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    iptres = orig.iptres;//   last residue that is considered part of solute (base 1 index)
    nspm = orig.nspm;//     total number of molecules
    nspsol = orig.nspsol;//   the first solvent "molecule" (base 1 index)

    //pointers are stored here in case the vectors are swapped around later
    //between difference classes. These arrays are LARGE; the goal here is to
    //save memory/time if that happens.
    titles = orig.titles;
    atom_names = orig.atom_names;
    charges = orig.charges;
    masses = orig.masses;
    atom_type_indices = orig.atom_type_indices;
    number_excluded_atoms = orig.number_excluded_atoms;
    nonbonded_parm_indices = orig.nonbonded_parm_indices;
    residue_labels = orig.residue_labels;
    residue_pointers = orig.residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = orig.bond_force_constants;
    bond_equil_values = orig.bond_equil_values;
    angle_force_constants = orig.angle_force_constants;
    angle_equil_values = orig.angle_equil_values;
    dihedral_force_constants = orig.dihedral_force_constants;
    dihedral_periodicities = orig.dihedral_periodicities;
    dihedral_phases = orig.dihedral_phases;
    soltys = orig.soltys;//solubility?
    lennard_jones_acoefs = orig.lennard_jones_acoefs;
    lennard_jones_bcoefs = orig.lennard_jones_bcoefs;
    bonds_inc_hydrogen = orig.bonds_inc_hydrogen;
    bonds_without_hydrogen = orig.bonds_without_hydrogen;
    angles_inc_hydrogen = orig.angles_inc_hydrogen;
    angles_without_hydrogen = orig.angles_without_hydrogen;
    dihedrals_inc_hydrogen = orig.dihedrals_inc_hydrogen;
    dihedrals_without_hydrogen = orig.dihedrals_without_hydrogen;
    excluded_atoms_list = orig.excluded_atoms_list;
    hbond_acoefs = orig.hbond_acoefs;
    hbond_bcoefs = orig.hbond_bcoefs;
    hbcuts = orig.hbcuts;
    amber_atom_types = orig.amber_atom_types;
    tree_chain_classifications = orig.tree_chain_classifications;
    join_array = orig.join_array;
    irotats = orig.irotats;
    radius_sets = orig.radius_sets;
    radii = orig.radii;
    screen = orig.screen;
    /*
    atom_names = new std::valarray<std::string>(*(orig.atom_names));
    charges = new std::valarray<mmpbsa_t>(*(orig.charges));
    masses = new std::valarray<mmpbsa_t>(*(orig.masses));
    atom_type_indices = new std::valarray<int>(*(orig.atom_type_indices));
    number_excluded_atoms = new std::valarray<int>(*(orig.number_excluded_atoms));
    nonbonded_parm_indices = new std::valarray<int>(*(orig.nonbonded_parm_indices));
    residue_labels = new std::valarray<std::string>(*(orig.residue_labels));
    residue_pointers = new std::valarray<int>(*(orig.residue_pointers));//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = new std::valarray<mmpbsa_t>(*(orig.bond_force_constants));
    bond_equil_values = new std::valarray<mmpbsa_t>(*(orig.bond_equil_values));
    angle_force_constants = new std::valarray<mmpbsa_t>(*(orig.angle_force_constants));
    angle_equil_values = new std::valarray<mmpbsa_t>(*(orig.angle_equil_values));
    dihedral_force_constants = new std::valarray<mmpbsa_t>(*(orig.dihedral_force_constants));
    dihedral_periodicities = new std::valarray<mmpbsa_t>(*(orig.dihedral_periodicities));
    dihedral_phases = new std::valarray<mmpbsa_t>(*(orig.dihedral_phases));
    soltys = new std::valarray<mmpbsa_t>(*(orig.soltys));//solubility?
    lennard_jones_acoefs = new std::valarray<mmpbsa_t>(*(orig.lennard_jones_acoefs));
    lennard_jones_bcoefs = new std::valarray<mmpbsa_t>(*(orig.lennard_jones_bcoefs));
    bonds_inc_hydrogen = new std::valarray<int>(*(orig.bonds_inc_hydrogen));
    bonds_without_hydrogen = new std::valarray<int>(*(orig.bonds_without_hydrogen));
    angles_inc_hydrogen = new std::valarray<int>(*(orig.angles_inc_hydrogen));
    angles_without_hydrogen = new std::valarray<int>(*(orig.angles_without_hydrogen));
    dihedrals_inc_hydrogen = new std::valarray<int>(*(orig.dihedrals_inc_hydrogen));
    dihedrals_without_hydrogen = new std::valarray<int>(*(orig.dihedrals_without_hydrogen));
    excluded_atoms_list = new std::valarray<int>(*(orig.excluded_atoms_list));
    hbond_acoefs = new std::valarray<mmpbsa_t>(*(orig.hbond_acoefs));
    hbond_bcoefs = new std::valarray<mmpbsa_t>(*(orig.hbond_bcoefs));
    hbcuts = new std::valarray<mmpbsa_t>(*(orig.hbcuts));
    amber_atom_types = new std::valarray<std::string>(*(orig.amber_atom_types));
    tree_chain_classifications = new std::valarray<std::string>(*(orig.tree_chain_classifications));
    join_array = new std::valarray<int>(*(orig.join_array));
    irotats = new std::valarray<int>(*(orig.irotats));
    radius_sets = new std::valarray<std::string>(*(orig.radius_sets));
    radii = new std::valarray<mmpbsa_t>(*(orig.radii));
    screen = new std::valarray<mmpbsa_t>(*(orig.screen));*/

}

sanderio::SanderParm& sanderio::SanderParm::operator=(const sanderio::SanderParm& orig)
{
    if(this == &orig)
        return *this;

    natom = orig.natom;// total number of atoms
    ntypes = orig.ntypes;// total number of distinct atom types
    nbonh = orig.nbonh;// number of bonds containing hydrogen
    mbona = orig.mbona;// number of bonds not containing hydrogen
    ntheth = orig.ntheth;// number of angles containing hydrogen
    mtheta = orig.mtheta;// number of angles not containing hydrogen
    nphih = orig.nphih;// number of dihedrals containing hydrogen
    mphia = orig.mphia;// number of dihedrals not containing hydrogen
    nhparm = orig.nhparm;// currently not used
    nparm = orig.nparm;// set to 1 if LES is used
    nnb = orig.nnb;// number of excluded atoms (=NEXT)
    nres = orig.nres;// number of residues
    nbona = orig.nbona;// MBONA + number of constraint bonds
    ntheta = orig.ntheta;// MTHETA + number of constraint angles
    nphia = orig.nphia;// MPHIA + number of constraint dihedrals
    numbnd = orig.numbnd;// number of unique bond types
    numang = orig.numang;// number of unique angle types
    nptra = orig.nptra;// number of unique dihedral types
    natyp = orig.natyp;// number of atom types in parameter file, see SOLTY below
    nphb = orig.nphb;// number of distinct 10-12 hydrogen bond pair types
    ifpert = orig.ifpert;// set to 1 if perturbation info is to be read in
    nbper = orig.nbper;// number of bonds to be perturbed
    ngper = orig.ngper;// number of angles to be perturbed
    ndper = orig.ndper;//number of dihedrals to be perturbed
    mbper = orig.mbper;// number of bonds with atoms completely in perturbed group
    mgper = orig.mgper;// number of angles with atoms completely in perturbed group
    mdper = orig.mdper;// number of dihedrals with atoms completely in perturbed groups
    ifbox = orig.ifbox;// set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = orig.nmxrs;// number of atoms in the largest residue
    ifcap = orig.ifcap;// set to 1 if the CAP option from edit was specified
    nextra = orig.nextra;// number of "extra points" (atom type of EP)

    //data only used with solvents
    iptres = orig.iptres;//   last residue that is considered part of solute (base 1 index)
    nspm = orig.nspm;//     total number of molecules
    nspsol = orig.nspsol;//   the first solvent "molecule" (base 1 index)

    //pointers are stored here in case the vectors are swapped around later
    //between difference classes. These arrays are LARGE; the goal here is to
    //save memory/time if that happens.
    titles = orig.titles;
    atom_names = orig.atom_names;
    charges = orig.charges;
    masses = orig.masses;
    atom_type_indices = orig.atom_type_indices;
    number_excluded_atoms = orig.number_excluded_atoms;
    nonbonded_parm_indices = orig.nonbonded_parm_indices;
    residue_labels = orig.residue_labels;
    residue_pointers = orig.residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = orig.bond_force_constants;
    bond_equil_values = orig.bond_equil_values;
    angle_force_constants = orig.angle_force_constants;
    angle_equil_values = orig.angle_equil_values;
    dihedral_force_constants = orig.dihedral_force_constants;
    dihedral_periodicities = orig.dihedral_periodicities;
    dihedral_phases = orig.dihedral_phases;
    soltys = orig.soltys;//solubility?
    lennard_jones_acoefs = orig.lennard_jones_acoefs;
    lennard_jones_bcoefs = orig.lennard_jones_bcoefs;
    bonds_inc_hydrogen = orig.bonds_inc_hydrogen;
    bonds_without_hydrogen = orig.bonds_without_hydrogen;
    angles_inc_hydrogen = orig.angles_inc_hydrogen;
    angles_without_hydrogen = orig.angles_without_hydrogen;
    dihedrals_inc_hydrogen = orig.dihedrals_inc_hydrogen;
    dihedrals_without_hydrogen = orig.dihedrals_without_hydrogen;
    excluded_atoms_list = orig.excluded_atoms_list;
    hbond_acoefs = orig.hbond_acoefs;
    hbond_bcoefs = orig.hbond_bcoefs;
    hbcuts = orig.hbcuts;
    amber_atom_types = orig.amber_atom_types;
    tree_chain_classifications = orig.tree_chain_classifications;
    join_array = orig.join_array;
    irotats = orig.irotats;
    radius_sets = orig.radius_sets;
    radii = orig.radii;
    screen = orig.screen;

    return *this;
}

void sanderio::SanderParm::raw_read_amber_parm(std::string file)
{

    using std::fstream;
    using std::string;

    fstream prmtopFile(file.c_str(),fstream::in);

    if(!prmtopFile.is_open())
        throw SanderIOException(strcat("Could not open: ",file.c_str()));

    string currentLine = getNextLine(prmtopFile);
    if(!strcmp(currentLine.substr(0,9).c_str(),"%VERSION"))
        throw SanderIOException(file.append(" is malformed. %VERSION is missing.").c_str());

    string flag;
    string format;

    //Fill Arrays
    int flagCounter = 0;
    while(prmtopFile.good())
    {
        flagCounter++;
        currentLine = getNextLine(prmtopFile);//should be FLAG
        if(currentLine.substr(0,5) != "%FLAG")
        {
            if(prmtopFile.eof())//there was whitespace before EOF
                return;

            char* flagLocation;
            sprintf(flagLocation,"%i",flagCounter);
            throw SanderIOException(file.append(" is malformed. %FLAG is "
                    "missing. Flag #").append(flagLocation).c_str(),BROKEN_PRMTOP_FILE);
        }
        else
        {
            flag = currentLine.substr(5);
            trimString(flag);
        }

        currentLine = getNextLine(prmtopFile);//should be FORMAT
        if(currentLine.substr(0,7) != "%FORMAT")
        {
            char* flagLocation;
            sprintf(flagLocation,"%i",flagCounter);
            throw SanderIOException(file.append(" is malformed. %FORMAT is "
                    "missing.").append(flagLocation).c_str());
        }
        else
        {
            format = currentLine.substr(7);
            trimString(format);
        }

        parseParmtopFile(prmtopFile,flag,format);
    }



}

bool sanderio::SanderParm::sanityCheck() throw (SanderIOException)
{
    using std::string;

    bool thereIsAProblem = false;

    char* error;
    int minAtomType = 1;
    if(!rangeCheck(atom_type_indices,minAtomType,ntypes))
    {
        sprintf(error,"Incorrect number of atom types. There should be %d types but "
                "the data ranges from %d to %d .",ntypes,atom_type_indices.min(),atom_type_indices.max());
        throw SanderIOException(error,INVALID_PRMTOP_DATA);
    }


    int icobot = 1;
    if(nphb != 0)
        icobot = -nphb;

    int maxParmIndex = int(0.5*ntypes*(ntypes+1));
    if(!rangeCheck(nonbonded_parm_indices,icobot,maxParmIndex))
    {
        sprintf(error,"Incorrect number of Non bonded parameter indices. "
                "They should be between %d and %d",icobot,int(0.5*ntypes*(ntypes+1)));
        throw SanderIOException(error,INVALID_PRMTOP_DATA);
    }

    if(nres > 0 && residue_pointers[0] != 1)
    {
        sprintf(error,"The first residue pointer must equal one. Instead, it equals %d",residue_pointers[0]);
        throw SanderIOException(error,INVALID_PRMTOP_DATA);
    }
    //check the order of the residue_pointers. Each element must be greater
    //      than the previous one.
    int lastResiduePointer = 0;
    for(int i = 0;i<residue_pointers.size();i++)
        if(lastResiduePointer < residue_pointers[i])
            lastResiduePointer = residue_pointers[i];
        else
            throw SanderIOException("RESIDUE_POINTER is out of order.");

    //Check Bonds
    int atomsPerBond = 2;
    try{bondCheck(bonds_inc_hydrogen,natom,nbonh,numbnd,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("bonds_inc_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(bonds_without_hydrogen,natom,mbona,numbnd,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("bonds_without_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    atomsPerBond = 3;
    try{bondCheck(angles_inc_hydrogen,natom,ntheth,numang,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("angles_inc_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(angles_without_hydrogen,natom,mtheta,numang,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("angles_without_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    atomsPerBond = 4;
    try{bondCheck(dihedrals_inc_hydrogen,natom,nphih,nptra,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("dihedrals_inc_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(dihedrals_without_hydrogen,natom,mphia,nptra,atomsPerBond);
    }catch(SanderIOException sioe){
        std::string message("dihedrals_wihtout_hydrogen data error: ");
        throw SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    //Check atom exclusions.
    int exclusionIndex = 0;
    for(int i = 0;i<natom;i++)
    {
        if(number_excluded_atoms[i] < 0 || number_excluded_atoms[i] > natom-1)
        {
            sprintf(error,"Number of excluded atoms must be between 0 and %d",natom-1);
            throw SanderIOException(error,INVALID_PRMTOP_DATA);
        }
        for(int j = exclusionIndex;j<exclusionIndex+number_excluded_atoms[i];j++)
        {
            int currentExcludedAtom = excluded_atoms_list[j];
            if(currentExcludedAtom != 0 &&
                    !(currentExcludedAtom > i && currentExcludedAtom <= natom))
            {
                sprintf(error,"For atom %d, excluded_atom_list element %d is "
                        "not in the range %d to %d",i+1,currentExcludedAtom,i+2,natom);
                throw SanderIOException(error,INVALID_PRMTOP_DATA);
            }
        }
        exclusionIndex += number_excluded_atoms[i];
    }


    //Check solvant parameters
    if(ifbox > 0)
    {
        if(iptres == 0 || nspm ==0 || nspsol == 0)
            throw SanderIOException("This is a solvent. Thus, 3 solvent pointers are needed."
                    ,BROKEN_PRMTOP_FILE);
        if(atoms_per_molecule.size() == 0 || atoms_per_molecule.sum() != natom)
        {
            sprintf(error,"A solvent is being used. Therefore, one must"
                    "provide the number of atoms per molecule and it cannot exceed"
                    "the total number of atoms (%d)",natom);
            throw SanderIOException(error,INVALID_PRMTOP_DATA);
        }

        if(box_dimensions.size() != 4)
            throw SanderIOException("This is a solvent. Thus, 4 BOX_DIMENSIONS are needed.",
                    INVALID_PRMTOP_DATA);

        //Check and warn about unsupported parameters.
        if(ifcap > 0)
        {
            std::cerr << "Warning: ifcap > 0, but water cap not supported" << std::endl;
            thereIsAProblem = true;
        }
        if(ifpert > 0)
        {
            std::cerr << "Warning: ifpert > 0, but perturbation not supported" << std::endl;
            thereIsAProblem = true;
        }
        if(nparm == 1)
        {
            std::cerr << "Warning: nparm == 1, but Locally Enhanced Sampling  not supported" << std::endl;
            thereIsAProblem = true;
        }

    }
    return !thereIsAProblem;//sanity check returns true if everything checks out.

}//end sanitycheck

void sanderio::SanderParm::parseParmtopFile(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
    //is the file still able to be read?
    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file.");

    //determine which flag is being used.
    //Sorry, but this simply has to be a long if..else chain.
    //If the flag is not known, the function returns, because some prmtop files
    //have parameters that are not used here, such as solvents.
    if(flag == "POINTERS")
        return loadPointers(prmtopFile,flag,format);
    else if(flag == "ATOM_NAME")
        loadPrmtopData(prmtopFile,atom_names,natom,format);
    else if(flag == "CHARGE")
        loadPrmtopData(prmtopFile,charges,natom,format);
    else if(flag == "MASS")
        loadPrmtopData(prmtopFile,masses,natom,format);
    else if(flag == "ATOM_TYPE_INDEX")
        loadPrmtopData(prmtopFile,atom_type_indices,natom,format);
    else if(flag == "NUMBER_EXCLUDED_ATOMS")
        loadPrmtopData(prmtopFile,number_excluded_atoms,natom,format);
    else if(flag == "NONBONDED_PARM_INDEX")
        loadPrmtopData(prmtopFile,nonbonded_parm_indices,ntypes*ntypes,format);
    else if(flag == "RESIDUE_LABEL")
        loadPrmtopData(prmtopFile,residue_labels,nres,format);
    else if(flag == "RESIDUE_POINTER")
        loadPrmtopData(prmtopFile,residue_pointers,nres,format);
    else if(flag == "BOND_FORCE_CONSTANT")
        loadPrmtopData(prmtopFile,bond_force_constants,numbnd,format);
    else if(flag == "BOND_EQUIL_VALUE")
        loadPrmtopData(prmtopFile,bond_equil_values,numbnd,format);
    else if(flag == "ANGLE_FORCE_CONSTANT")
        loadPrmtopData(prmtopFile,angle_force_constants,numang,format);
    else if(flag == "ANGLE_EQUIL_VALUE")
        loadPrmtopData(prmtopFile,angle_equil_values,numang,format);
    else if(flag == "DIHEDRAL_FORCE_CONSTANT")
        loadPrmtopData(prmtopFile,dihedral_force_constants,nptra,format);
    else if(flag == "DIHEDRAL_PERIODICITY")
        loadPrmtopData(prmtopFile,dihedral_periodicities,nptra,format);
    else if(flag == "DIHEDRAL_PHASE")
        loadPrmtopData(prmtopFile,dihedral_phases,nptra,format);
    else if(flag == "SOLTY")
        loadPrmtopData(prmtopFile,soltys,natyp,format);
    else if(flag == "LENNARD_JONES_ACOEF")
        loadPrmtopData(prmtopFile,lennard_jones_acoefs,int(0.5*ntypes*(ntypes+1)),format);
    else if(flag == "LENNARD_JONES_BCOEF")
        loadPrmtopData(prmtopFile,lennard_jones_bcoefs,int(0.5*ntypes*(ntypes+1)),format);
    else if(flag == "BONDS_INC_HYDROGEN")
        loadPrmtopData(prmtopFile,bonds_inc_hydrogen,nbonh*3,format);
    else if(flag == "BONDS_WITHOUT_HYDROGEN")
        loadPrmtopData(prmtopFile,bonds_without_hydrogen,nbona*3,format);
    else if(flag == "ANGLES_INC_HYDROGEN")
        loadPrmtopData(prmtopFile,angles_inc_hydrogen,4*ntheth,format);
    else if(flag == "ANGLES_WITHOUT_HYDROGEN")
        loadPrmtopData(prmtopFile,angles_without_hydrogen,4*ntheta,format);
    else if(flag == "DIHEDRALS_INC_HYDROGEN")
        loadPrmtopData(prmtopFile,dihedrals_inc_hydrogen,5*nphih,format);
    else if(flag == "DIHEDRALS_WITHOUT_HYDROGEN")
        loadPrmtopData(prmtopFile,dihedrals_without_hydrogen,5*nphia,format);
    else if(flag == "EXCLUDED_ATOMS_LIST")
        loadPrmtopData(prmtopFile,excluded_atoms_list,nnb,format);
    else if(flag == "HBOND_ACOEF")
        loadPrmtopData(prmtopFile,hbond_acoefs,nphb,format);
    else if(flag == "HBOND_BCOEF")
        loadPrmtopData(prmtopFile,hbond_bcoefs,nphb,format);
    else if(flag == "HBCUT")
        loadPrmtopData(prmtopFile,hbcuts,nphb,format);
    else if(flag == "AMBER_ATOM_TYPE")
        loadPrmtopData(prmtopFile,amber_atom_types,natom,format);
    else if(flag == "TREE_CHAIN_CLASSIFICATION")
        loadPrmtopData(prmtopFile,tree_chain_classifications,natom,format);
    else if(flag == "JOIN_ARRAY")
        loadPrmtopData(prmtopFile,join_array,natom,format);
    else if(flag == "IROTAT")
        loadPrmtopData(prmtopFile,irotats,natom,format);
    else if(flag == "RADIUS_SET")
    {
        radius_sets = getNextLine(prmtopFile);
        trimString(radius_sets);
    }
    else if(flag == "RADII")
        loadPrmtopData(prmtopFile,radii,natom,format);
    else if(flag == "SCREEN")
        loadPrmtopData(prmtopFile,screen,natom,format);
    else if(flag == "TITLE")
    {
        titles = getNextLine(prmtopFile);
        trimString(titles);
    }
    else if(flag == "SOLVENT_POINTERS")
        loadSolventPointers(prmtopFile,flag,format);
    else if(flag == "ATOMS_PER_MOLECULE")
        loadPrmtopData(prmtopFile,atoms_per_molecule,nspm,format);
    else if(flag == "BOX_DIMENSIONS")
        loadPrmtopData(prmtopFile,box_dimensions,4,format);
    else
        std::cout << flag << " was not parsed." << std::endl;


}

void sanderio::SanderParm::loadPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
    using std::string;

    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file. loadPointers");

    std::valarray<int> pointers(0,31);
    loadPrmtopData(prmtopFile,pointers,31,format);


    //populate parameters
    int i=0;
    natom = pointers[i++];// total number of atoms
    ntypes = pointers[i++];// total number of distinct atom types
    nbonh = pointers[i++];// number of bonds containing hydrogen
    mbona = pointers[i++];// number of bonds not containing hydrogen
    ntheth = pointers[i++];// number of angles containing hydrogen
    mtheta = pointers[i++];// number of angles not containing hydrogen
    nphih = pointers[i++];// number of dihedrals containing hydrogen
    mphia = pointers[i++];// number of dihedrals not containing hydrogen
    nhparm = pointers[i++];// currently not used
    nparm = pointers[i++];// set to 1 if LES is used
    nnb = pointers[i++];// number of excluded atoms (=NEXT)
    nres = pointers[i++];// number of residues
    nbona = pointers[i++];// MBONA + number of constraint bonds
    ntheta = pointers[i++];// MTHETA + number of constraint angles
    nphia = pointers[i++];// MPHIA + number of constraint dihedrals
    numbnd = pointers[i++];// number of unique bond types
    numang = pointers[i++];// number of unique angle types
    nptra = pointers[i++];// number of unique dihedral types
    natyp = pointers[i++];// number of atom types in parameter file, see SOLTY below
    nphb = pointers[i++];// number of distinct 10-12 hydrogen bond pair types
    ifpert = pointers[i++];// set to 1 if perturbation info is to be read in
    nbper = pointers[i++];// number of bonds to be perturbed
    ngper = pointers[i++];// number of angles to be perturbed
    ndper = pointers[i++];//number of dihedrals to be perturbed
    mbper = pointers[i++];// number of bonds with atoms completely in perturbed group
    mgper = pointers[i++];// number of angles with atoms completely in perturbed group
    mdper = pointers[i++];// number of dihedrals with atoms completely in perturbed groups
    ifbox = pointers[i++];// set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = pointers[i++];// number of atoms in the largest residue
    ifcap = pointers[i++];// set to 1 if the CAP option from edit was specified
    nextra = pointers[i++];// number of "extra points" (atom type of EP)


}

void sanderio::SanderParm::loadSolventPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
    using std::string;

    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file. loadPointers");

    std::valarray<int> pointers(0,3);
    loadPrmtopData(prmtopFile,pointers,3,format);

    int i = 0;
    iptres = pointers[i++];///<   last residue that is considered part of solute (base 1 index)
    nspm = pointers[i++];///<     total number of molecules
    nspsol = pointers[i++];
}

template <class T> bool sanderio::SanderParm::rangeCheck(const std::valarray<T>& array,
    const T& min, const T& max)
{
    if(array.max() > max)
        return false;

    if(array.min() < min)
        return false;

    return true;
}

template <class T> bool sanderio::SanderParm::bondCheck(const std::valarray<T>& array,
        const int& natoms, const int& nbonds, const int& ntypes,
        const int& atomsPerBond)
{
    char * error;
    int maxi = (natoms-1)*3;
    if(array.size() != nbonds*(atomsPerBond+1))
    {
        std::string message("Incorrect number of bonds. Expected ");
        sprintf(error,"%d",nbonds*(atomsPerBond+1));
        throw SanderIOException(message.append(error),
            INVALID_PRMTOP_DATA);
    }

    //check bond code range in respect to their specific bonds.
    for(int i = 0;i<nbonds*(atomsPerBond+1);i+=atomsPerBond+1)
    {
        for(int j = i;j<i+atomsPerBond;j++)
        {
            int absbnd = abs(array[j]);
            if(absbnd > maxi)
            {
                std::string message("Bond code exceed max of ");
                sprintf(error,"%d",maxi);
                message.append(error);
                throw SanderIOException(message,INVALID_PRMTOP_DATA);
            }
        }

        if(array[i+atomsPerBond] > ntypes)
        {
            std::string message("Bond code exceed max of ");
            sprintf(error,"%d",ntypes);
            message.append(error);
            throw SanderIOException(message,INVALID_PRMTOP_DATA);
        }
    }

    return true;
}


template <class T> void sanderio::SanderParm::loadPrmtopData(std::fstream& prmtopFile,
        std::valarray<T>& array,int size,const std::string& format)
{
    using std::valarray;
    using std::string;

    //ensure array is of the correct size (or exists).
    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file. loadPrmtopData");

    //use the format to obtain the array dimensions (in the 2-D sense).
    int numberOfColumns = 0;
    int columnWidth = 0;
    sscanf(format.c_str(),"%*c%d%*c%d",&numberOfColumns,&columnWidth);

    if(!loadValarray(prmtopFile,array,size,columnWidth,numberOfColumns))
        throw SanderIOException("Could not load Parameter Data",BROKEN_PRMTOP_FILE);

}


std::string sanderio::read_crds(std::fstream& crdFile, std::valarray<mmpbsa_t>& crds)
{
    using std::string;

    if(!crdFile.good())
        throw SanderIOException("Cannot open coordinate file",FILE_READ_ERROR);

    string title = getNextLine(crdFile);
    string strNatoms = getNextLine(crdFile);
    trimString(strNatoms);
    int natoms = 0;
    sscanf(strNatoms.c_str(),"%d",&natoms);

    if(!loadValarray(crdFile,crds,natoms*3,12,8))
        throw SanderIOException("Coordinate file is too short.",FILE_READ_ERROR);

    return title;
}

void sanderio::write_crds(const char* fileName,const std::valarray<mmpbsa_t>& crds,
    const char* title)
{
    using std::valarray;
    using std::slice;

    if(crds.size() % 3 != 0)
        throw SanderIOException("The number of elements in the coordinate array "
                "must be a multiple of 3, ie 3-dimensions.",DATA_FORMAT_ERROR);

    int natoms = int(crds.size()/3);

    std::fstream outFile(fileName,std::ios::out);
    char* format = "%12.7f";//format of the coordinate data

    if(!outFile.good())
        throw SanderIOException(std::string("Could not open: ").append(fileName),FILE_READ_ERROR);

    char strOutput[12];//used for outputting Fortran formatted strings with sprintf.

    outFile << title << std::endl;

    sprintf(strOutput,"%5d",natoms);//number of atoms
    outFile << strOutput << std::endl;

    int m;
    mmpbsa_t dblOutput;
    //save data in rows of 6
    for(m = 0;m<crds.size() - 6;m+=6)
    {
        valarray<mmpbsa_t> row = crds[slice(m,6,1)];//m-th row

        for(int i = 0;i<6;i++)
        {
            dblOutput = row[i];
            sprintf(strOutput,format,dblOutput);
            outFile << strOutput;
        }
        outFile << std::endl;
    }

    //save the possibly incomplete last row.
    if(m<crds.size())
    {
        for(m;m<crds.size();m++)
        {
            dblOutput = crds[m];
            sprintf(strOutput,format,dblOutput);
            outFile << strOutput;
        }
        outFile << std::endl;
    }

    outFile.close();
}

std::string sanderio::get_traj_title(std::fstream& trajFile)
{
    trajFile.seekg(0,std::ios::beg);
    return getNextLine(trajFile);
}

bool sanderio::get_next_snap(std::fstream& trajFile, std::valarray<mmpbsa_t>& snapshot,
    const int& natoms,bool isPeriodic)
{
    bool returnMe = loadValarray(trajFile,snapshot,natoms*3,8,10);
    if(isPeriodic)
        getNextLine(trajFile);//ignoring periodic box information
    return returnMe;
}

void sanderio::skip_next_snap(std::fstream& trajFile, const int& natoms, bool isPeriodic)
{
    int numlines = int(natoms*3/10);
    if(natoms*3 % 10)
        numlines++;
    if(isPeriodic)
        numlines++;
    trajFile.seekg(numlines,trajFile.cur);
}

template <class T> bool sanderio::loadValarray(std::fstream& dataFile,
        std::valarray<T>& dataArray, const int& arrayLength, const int& width,
        const int& numberOfColumns)
{
    using std::string;

    if(dataFile.eof())
        return false;

    if(dataArray.size() != arrayLength)
        dataArray.resize(arrayLength,0);

    int lineIndex = 0;
    float fltCurrentData = 0;
    int dataIndex = 0;

     for(int i = 0;i<arrayLength;i++)
    {
        if(dataFile.eof())
            throw SanderIOException("Data file ended in the middle of the "
                    "data.",BROKEN_TRAJECTORY_FILE);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
            char* error;
            sprintf(error,"Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #%d is only"
                    "%d characters long.",lineIndex+1,currentLine.size());
            std::cerr << error << std::endl;
        }

        //tokenize line into data. put data into valarray.
        int currentColumn = 0;
        while(currentLine.size() > 0 && currentColumn < numberOfColumns)
        {
            string currentData = currentLine.substr(0,width);
            sscanf(currentData.c_str(),"%f",&fltCurrentData);
            dataArray[dataIndex++] = fltCurrentData;
            currentLine.erase(0,width);
            currentColumn++;
        }

        lineIndex++;
        if(dataIndex == arrayLength)
            return true;//in case the file has too much data, ie periodic box
    }


    return true;

}

template <> bool sanderio::loadValarray<std::string>(std::fstream& dataFile,
        std::valarray<std::string>& dataArray, const int& arrayLength, const int& width,
        const int& numberOfColumns)
{
    using std::string;

    if(dataFile.eof())
        return false;

    if(dataArray.size() != arrayLength)
        dataArray.resize(arrayLength,"");

    int lineIndex = 0;
    int dataIndex = 0;

    for(int i = 0;i<arrayLength;i++)
    {
        if(dataFile.eof())
            throw SanderIOException("Data file ended in the middle of the "
                    "data.",BROKEN_TRAJECTORY_FILE);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
            char* error;
            sprintf(error,"Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #%d is only"
                    "%d characters long.",lineIndex+1,currentLine.size());
            std::cerr << error << std::endl;
        }

        //tokenize line into data. put data into valarray.
        int currentColumn = 0;
        while(currentLine.size() > 0 && currentColumn < numberOfColumns)
        {
            dataArray[dataIndex++] = currentLine.substr(0,width);
            currentLine.erase(0,width);
            currentColumn++;
        }

        lineIndex++;
        if(dataIndex == arrayLength)
            return true;//in case the file has too much data, ie periodic box
    }


    return true;
}

std::string getNextLine(std::fstream& file)
{
    if(!file.good())
        throw SanderIOException("Could not read from file");

    std::string returnMe;
    getline(file,returnMe);
    return returnMe;
}

void sanderio::trimString(std::string& bean)
{
    using std::string;
    int pos = bean.find_last_not_of(' ');
    if(pos != string::npos)
    {
        bean.erase(pos + 1);
        pos = bean.find_first_not_of(' ');
        if(pos != string::npos) bean.erase(0, pos);
    }
    else
        bean.erase(bean.begin(), bean.end());
}
