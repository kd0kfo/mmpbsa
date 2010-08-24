#include "SanderParm.h"


mmpbsa::SanderParm::SanderParm() {
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
    initializeArrays();
}

mmpbsa::SanderParm::SanderParm(const SanderParm& orig) {
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
    atom_names.resize(orig.atom_names.size());//gnu's c++ compiler does not copy elements if valarrays' sizes differ :-(
    atom_names = orig.atom_names;
    charges.resize(orig.charges.size());
    charges = orig.charges;
    masses.resize(orig.masses.size());
    masses = orig.masses;
    atom_type_indices.resize(orig.atom_type_indices.size());
    atom_type_indices = orig.atom_type_indices;
    number_excluded_atoms.resize(orig.number_excluded_atoms.size());
    number_excluded_atoms = orig.number_excluded_atoms;
    nonbonded_parm_indices.resize(orig.nonbonded_parm_indices.size());
    nonbonded_parm_indices = orig.nonbonded_parm_indices;
    nonbonded_parm_mask.resize(orig.nonbonded_parm_mask.size());
    nonbonded_parm_mask = orig.nonbonded_parm_mask;
    residue_labels.resize(orig.residue_labels.size());
    residue_labels = orig.residue_labels;
    residue_pointers.resize(orig.residue_pointers.size());
    residue_pointers = orig.residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants.resize(orig.bond_force_constants.size());
    bond_force_constants = orig.bond_force_constants;
    bond_equil_values.resize(orig.bond_equil_values.size());
    bond_equil_values = orig.bond_equil_values;
    angle_force_constants.resize(orig.angle_force_constants.size());
    angle_force_constants = orig.angle_force_constants;
    angle_equil_values.resize(orig.angle_equil_values.size());
    angle_equil_values = orig.angle_equil_values;
    dihedral_force_constants.resize(orig.dihedral_force_constants.size());
    dihedral_force_constants = orig.dihedral_force_constants;
    dihedral_periodicities.resize(orig.dihedral_periodicities.size());
    dihedral_periodicities = orig.dihedral_periodicities;
    dihedral_phases.resize(orig.dihedral_phases.size());
    dihedral_phases = orig.dihedral_phases;
    soltys.resize(orig.soltys.size());
    soltys = orig.soltys;//solubility?
    lennard_jones_acoefs.resize(orig.lennard_jones_acoefs.size());
    lennard_jones_acoefs = orig.lennard_jones_acoefs;
    lennard_jones_bcoefs.resize(orig.lennard_jones_bcoefs.size());
    lennard_jones_bcoefs = orig.lennard_jones_bcoefs;
    bonds_inc_hydrogen.resize(orig.bonds_inc_hydrogen.size());
    bonds_inc_hydrogen = orig.bonds_inc_hydrogen;
    bonds_without_hydrogen.resize(orig.bonds_without_hydrogen.size());
    bonds_without_hydrogen = orig.bonds_without_hydrogen;
    angles_inc_hydrogen.resize(orig.angles_inc_hydrogen.size());
    angles_inc_hydrogen = orig.angles_inc_hydrogen;
    angles_without_hydrogen.resize(orig.angles_without_hydrogen.size());
    angles_without_hydrogen = orig.angles_without_hydrogen;
    dihedrals_inc_hydrogen.resize(orig.dihedrals_inc_hydrogen.size());
    dihedrals_inc_hydrogen = orig.dihedrals_inc_hydrogen;
    dihedrals_without_hydrogen.resize(orig.dihedrals_without_hydrogen.size());
    dihedrals_without_hydrogen = orig.dihedrals_without_hydrogen;
    dihedral_h_mask.resize(orig.dihedral_h_mask.size());
    dihedral_h_mask = orig.dihedral_h_mask;
    dihedral_mask.resize(orig.dihedral_mask.size());
    dihedral_mask = orig.dihedral_mask;
    excluded_atoms_list.resize(orig.excluded_atoms_list.size());
    excluded_atoms_list = orig.excluded_atoms_list;
    hbond_acoefs.resize(orig.hbond_acoefs.size());
    hbond_acoefs = orig.hbond_acoefs;
    hbond_bcoefs.resize(orig.hbond_bcoefs.size());
    hbond_bcoefs = orig.hbond_bcoefs;
    hbcuts.resize(orig.hbcuts.size());
    hbcuts = orig.hbcuts;
    amber_atom_types.resize(orig.amber_atom_types.size());
    amber_atom_types = orig.amber_atom_types;
    tree_chain_classifications.resize(orig.tree_chain_classifications.size());
    tree_chain_classifications = orig.tree_chain_classifications;
    join_array.resize(orig.join_array.size());
    join_array = orig.join_array;
    irotats.resize(orig.irotats.size());
    irotats = orig.irotats;
    radius_sets.resize(orig.radius_sets.size());
    radius_sets = orig.radius_sets;
    radii.resize(orig.radii.size());
    radii = orig.radii;
    screen.resize(orig.screen.size());
    screen = orig.screen;
    atoms_per_molecule.resize(orig.atoms_per_molecule.size());
    atoms_per_molecule = orig.atoms_per_molecule;
    box_dimensions.resize(orig.box_dimensions.size());
    box_dimensions = orig.box_dimensions;

}

mmpbsa::SanderParm& mmpbsa::SanderParm::operator=(const mmpbsa::SanderParm& orig)
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
    atom_names.resize(orig.atom_names.size());//gnu's c++ compiler does not copy elements if valarrays' sizes differ :-(
    atom_names = orig.atom_names;
    charges.resize(orig.charges.size());
    charges = orig.charges;
    masses.resize(orig.masses.size());
    masses = orig.masses;
    atom_type_indices.resize(orig.atom_type_indices.size());
    atom_type_indices = orig.atom_type_indices;
    number_excluded_atoms.resize(orig.number_excluded_atoms.size());
    number_excluded_atoms = orig.number_excluded_atoms;
    nonbonded_parm_indices.resize(orig.nonbonded_parm_indices.size());
    nonbonded_parm_indices = orig.nonbonded_parm_indices;
    nonbonded_parm_mask.resize(orig.nonbonded_parm_mask.size());
    nonbonded_parm_mask = orig.nonbonded_parm_mask;
    residue_labels.resize(orig.residue_labels.size());
    residue_labels = orig.residue_labels;
    residue_pointers.resize(orig.residue_pointers.size());
    residue_pointers = orig.residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants.resize(orig.bond_force_constants.size());
    bond_force_constants = orig.bond_force_constants;
    bond_equil_values.resize(orig.bond_equil_values.size());
    bond_equil_values = orig.bond_equil_values;
    angle_force_constants.resize(orig.angle_force_constants.size());
    angle_force_constants = orig.angle_force_constants;
    angle_equil_values.resize(orig.angle_equil_values.size());
    angle_equil_values = orig.angle_equil_values;
    dihedral_force_constants.resize(orig.dihedral_force_constants.size());
    dihedral_force_constants = orig.dihedral_force_constants;
    dihedral_periodicities.resize(orig.dihedral_periodicities.size());
    dihedral_periodicities = orig.dihedral_periodicities;
    dihedral_phases.resize(orig.dihedral_phases.size());
    dihedral_phases = orig.dihedral_phases;
    soltys.resize(orig.soltys.size());
    soltys = orig.soltys;//solubility?
    lennard_jones_acoefs.resize(orig.lennard_jones_acoefs.size());
    lennard_jones_acoefs = orig.lennard_jones_acoefs;
    lennard_jones_bcoefs.resize(orig.lennard_jones_bcoefs.size());
    lennard_jones_bcoefs = orig.lennard_jones_bcoefs;
    bonds_inc_hydrogen.resize(orig.bonds_inc_hydrogen.size());
    bonds_inc_hydrogen = orig.bonds_inc_hydrogen;
    bonds_without_hydrogen.resize(orig.bonds_without_hydrogen.size());
    bonds_without_hydrogen = orig.bonds_without_hydrogen;
    angles_inc_hydrogen.resize(orig.angles_inc_hydrogen.size());
    angles_inc_hydrogen = orig.angles_inc_hydrogen;
    angles_without_hydrogen.resize(orig.angles_without_hydrogen.size());
    angles_without_hydrogen = orig.angles_without_hydrogen;
    dihedrals_inc_hydrogen.resize(orig.dihedrals_inc_hydrogen.size());
    dihedrals_inc_hydrogen = orig.dihedrals_inc_hydrogen;
    dihedrals_without_hydrogen.resize(orig.dihedrals_without_hydrogen.size());
    dihedrals_without_hydrogen = orig.dihedrals_without_hydrogen;
    dihedral_h_mask.resize(orig.dihedral_h_mask.size());
    dihedral_h_mask = orig.dihedral_h_mask;
    dihedral_mask.resize(orig.dihedral_mask.size());
    dihedral_mask = orig.dihedral_mask;
    excluded_atoms_list.resize(orig.excluded_atoms_list.size());
    excluded_atoms_list = orig.excluded_atoms_list;
    hbond_acoefs.resize(orig.hbond_acoefs.size());
    hbond_acoefs = orig.hbond_acoefs;
    hbond_bcoefs.resize(orig.hbond_bcoefs.size());
    hbond_bcoefs = orig.hbond_bcoefs;
    hbcuts.resize(orig.hbcuts.size());
    hbcuts = orig.hbcuts;
    amber_atom_types.resize(orig.amber_atom_types.size());
    amber_atom_types = orig.amber_atom_types;
    tree_chain_classifications.resize(orig.tree_chain_classifications.size());
    tree_chain_classifications = orig.tree_chain_classifications;
    join_array.resize(orig.join_array.size());
    join_array = orig.join_array;
    irotats.resize(orig.irotats.size());
    irotats = orig.irotats;
    radius_sets.resize(orig.radius_sets.size());
    radius_sets = orig.radius_sets;
    radii.resize(orig.radii.size());
    radii = orig.radii;
    screen.resize(orig.screen.size());
    screen = orig.screen;
    atoms_per_molecule.resize(orig.atoms_per_molecule.size());
    atoms_per_molecule = orig.atoms_per_molecule;
    box_dimensions.resize(orig.box_dimensions.size());
    box_dimensions = orig.box_dimensions;

    return *this;
}

void mmpbsa::SanderParm::initializeArrays()
{
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
    dihedral_h_mask;
    dihedral_mask;
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
    screen;
    atoms_per_molecule;
    box_dimensions;
    }

void mmpbsa::SanderParm::raw_read_amber_parm(const std::string& file) throw (mmpbsa::SanderIOException)
{
    using std::fstream;
    using std::string;
    using namespace mmpbsa_utils;

    fstream prmtopFile;
    mmpbsa_io::fileopen(file.c_str(),fstream::in,prmtopFile);

    try
    {
        raw_read_amber_parm(prmtopFile);
    }
    catch(mmpbsa::SanderIOException sioe)
    {
        std::ostringstream buff;
        buff << sioe.what() << std::endl << "Parmtop File = " << file;
        throw mmpbsa::SanderIOException(buff,sioe.getErrType());
    }
    prmtopFile.close();
}

void mmpbsa::SanderParm::raw_read_amber_parm(std::fstream& prmtopFile) throw (mmpbsa::SanderIOException)
{
    using std::fstream;
    using std::string;
    using namespace mmpbsa_utils;
    using mmpbsa_io::getNextLine;

    if(!prmtopFile.is_open())
        throw mmpbsa::SanderIOException("Could not open parmtop file.",mmpbsa::BROKEN_PRMTOP_FILE);

    string currentLine = getNextLine(prmtopFile);
    if(!strcmp(currentLine.substr(0,9).c_str(),"%VERSION"))
        throw mmpbsa::SanderIOException("Parmtop file is malformed. %%VERSION is missing.",mmpbsa::BROKEN_PRMTOP_FILE);

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

            std::ostringstream flagLocation;
            flagLocation << "Parmtop file is malformed. %FLAG is "
                    "missing. Flag #" << flagCounter;
            throw mmpbsa::SanderIOException(flagLocation,mmpbsa::BROKEN_PRMTOP_FILE);
        }
        else
        {
            flag = currentLine.substr(5);
            flag = trimString(flag);
        }

        currentLine = getNextLine(prmtopFile);//should be FORMAT
        if(currentLine.substr(0,7) != "%FORMAT")
        {
            std::ostringstream flagLocation;
            flagLocation << "Parmtop file is malformed. %FORMAT is "
                    "missing. Format #" << flagCounter;
            throw mmpbsa::SanderIOException(flagLocation,mmpbsa::BROKEN_PRMTOP_FILE);
        }
        else
        {
            format = trimString(currentLine.substr(7));
        }

        parseParmtopFile(prmtopFile,flag,format);
    }
}

bool mmpbsa::SanderParm::sanityCheck() throw (mmpbsa::SanderIOException)
{
    using std::string;

    bool thereIsAProblem = false;

    std::ostringstream error;

    if(natom == 0)
    {
        std::cerr << "WARNING: SanderParm had natom = 0, i.e. trivial "
                "SanderParm. Check to see if this was intentional, as the "
                "program is continuing, albeit with the SanderParm object failing"
                "the sanity check." << std::endl;
        return false;
    }

    if(!rangeCheck(atom_type_indices,size_t(1),ntypes))
    {
        error << "Incorrect number of atom types. There should be "
                << ntypes << " types but the data ranges from "
                << atom_type_indices.min() 
                << " to " << atom_type_indices.max() << ".";
        throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
    }


    size_t maxParmIndex = size_t(0.5*ntypes*(ntypes+1));
    if(!rangeCheck(nonbonded_parm_indices,size_t(1),maxParmIndex))//nonbonded_parm_indices are positive here. The values that were negative in the prmtop file are indicated in the nonbonded_parm_mask
    {
        error << "Incorrect number of Non bonded parameter indices. "
                "They should be between 1 and " << maxParmIndex;
        throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
    }

    if(nres > 0 && residue_pointers[0] != 1)
    {
        error << "The first residue pointer must equal one. Instead, it equals "
                << residue_pointers[0];
        throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
    }

    //check the order of the residue_pointers. Each element must be greater
    //      than the previous one.
    size_t lastResiduePointer = 0;
    for(size_t i = 0;i<residue_pointers.size();i++)
        if(lastResiduePointer < residue_pointers[i])
            lastResiduePointer = residue_pointers[i];
        else
            throw mmpbsa::SanderIOException("RESIDUE_POINTER is out of order.");

    //Check Bonds
    int atomsPerBond = 2;
    try{bondCheck(bonds_inc_hydrogen,natom,nbonh,numbnd,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("bonds_inc_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(bonds_without_hydrogen,natom,mbona,numbnd,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("bonds_without_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    atomsPerBond = 3;
    try{bondCheck(angles_inc_hydrogen,natom,ntheth,numang,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("angles_inc_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(angles_without_hydrogen,natom,mtheta,numang,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("angles_without_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    atomsPerBond = 4;
    try{bondCheck(dihedrals_inc_hydrogen,natom,nphih,nptra,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("dihedrals_inc_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    try{bondCheck(dihedrals_without_hydrogen,natom,mphia,nptra,atomsPerBond);
    }catch(mmpbsa::SanderIOException sioe){
        std::string message("dihedrals_wihtout_hydrogen data error: ");
        throw mmpbsa::SanderIOException(message.append(sioe.what()),sioe.getErrType());
    }

    //Check atom exclusions.
    size_t exclusionIndex = 0;
    for(size_t i = 0;i<natom;i++)
    {
        if(number_excluded_atoms[i] < 0 || number_excluded_atoms[i] > natom-1)
        {
            error << "Number of excluded atoms must be between 0 and " << natom-1;
            throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
        }
        for(size_t j = exclusionIndex;j<exclusionIndex+number_excluded_atoms[i];j++)
        {
            size_t currentExcludedAtom = excluded_atoms_list[j];
            if(currentExcludedAtom != 0 &&
                    !(currentExcludedAtom > i && currentExcludedAtom <= natom))
            {
                error << "For atom " << i+1 << ", excluded_atom_list element "
                        << currentExcludedAtom << " is not in the range "
                        << i+2 << " to " << natom;
                throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
            }
        }
        exclusionIndex += number_excluded_atoms[i];
    }


    //Check solvent parameters
    if(ifbox > 0)
    {
        if(iptres == 0 || nspm ==0 || nspsol == 0)
            throw mmpbsa::SanderIOException("This is a solvent. Thus, 3 solvent pointers are needed."
                    ,mmpbsa::BROKEN_PRMTOP_FILE);
        if(atoms_per_molecule.size() == 0 || atoms_per_molecule.sum() != natom)
        {
            error << "A solvent is being used. Therefore, one must "
                    "provide the number of atoms per molecule and it cannot exceed "
                    "the total number of atoms (" << natom << ")";
            throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
        }

        if(box_dimensions.size() != 4)
            throw mmpbsa::SanderIOException("This is a solvent. Thus, 4 BOX_DIMENSIONS are needed.",
                    mmpbsa::INVALID_PRMTOP_DATA);

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

void mmpbsa::SanderParm::parseParmtopFile(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
  using std::slice;
  using std::valarray;
  using namespace mmpbsa_utils;
  using mmpbsa_io::getNextLine;

    //is the file still able to be read?
    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file.");

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
        loadPrmtopMaskedData(prmtopFile,nonbonded_parm_indices,nonbonded_parm_mask,ntypes*ntypes,format);
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
      {
        loadPrmtopData(prmtopFile,bonds_inc_hydrogen,nbonh*3,format);
	valarray<size_t> one(1,nbonh);
	bonds_inc_hydrogen[slice(2,nbonh,3)] -= one;
	one *= 3;
	bonds_inc_hydrogen[slice(0,nbonh,3)] /= one;
	bonds_inc_hydrogen[slice(1,nbonh,3)] /= one;

      }
    else if(flag == "BONDS_WITHOUT_HYDROGEN")
      {
        loadPrmtopData(prmtopFile,bonds_without_hydrogen,nbona*3,format);
	valarray<size_t> one(1,nbona);
	bonds_without_hydrogen[slice(2,nbona,3)] -= one;
	one *= 3;
	bonds_without_hydrogen[slice(0,nbona,3)] /= one;
	bonds_without_hydrogen[slice(1,nbona,3)] /= one;
      }
    else if(flag == "ANGLES_INC_HYDROGEN")
      {
        loadPrmtopData(prmtopFile,angles_inc_hydrogen,4*ntheth,format);
	valarray<size_t> one(1,ntheth);
	angles_inc_hydrogen[slice(3,ntheth,4)] -= one;
	one *= 3;
	angles_inc_hydrogen[slice(0,ntheth,4)] /= one;
	angles_inc_hydrogen[slice(1,ntheth,4)] /= one;
	angles_inc_hydrogen[slice(2,ntheth,4)] /= one;
      }
    else if(flag == "ANGLES_WITHOUT_HYDROGEN")
      {
        loadPrmtopData(prmtopFile,angles_without_hydrogen,4*ntheta,format);
	valarray<size_t> one(1,ntheta);
	angles_without_hydrogen[slice(3,ntheta,4)] -= one;
	one *= 3;
	angles_without_hydrogen[slice(0,ntheta,4)] /= one;
	angles_without_hydrogen[slice(1,ntheta,4)] /= one;
	angles_without_hydrogen[slice(2,ntheta,4)] /= one;
      }
    else if(flag == "DIHEDRALS_INC_HYDROGEN")
      {
        loadPrmtopMaskedData(prmtopFile,dihedrals_inc_hydrogen,dihedral_h_mask,5*nphih,format);
	valarray<size_t> one(1,nphih);
	dihedrals_inc_hydrogen[slice(4,nphih,5)] -= one;
	one *= 3;
	dihedrals_inc_hydrogen[slice(0,nphih,5)] /= one;
	dihedrals_inc_hydrogen[slice(1,nphih,5)] /= one;
	dihedrals_inc_hydrogen[slice(2,nphih,5)] /= one;
	dihedrals_inc_hydrogen[slice(3,nphih,5)] /= one;
      }
    else if(flag == "DIHEDRALS_WITHOUT_HYDROGEN")
      {
	loadPrmtopMaskedData(prmtopFile,dihedrals_without_hydrogen,dihedral_mask,5*nphia,format);
	valarray<size_t> one(1,nphia);
	dihedrals_without_hydrogen[slice(4,nphia,5)] -= one;
	one *= 3;
	dihedrals_without_hydrogen[slice(0,nphia,5)] /= one;
	dihedrals_without_hydrogen[slice(1,nphia,5)] /= one;
	dihedrals_without_hydrogen[slice(2,nphia,5)] /= one;
	dihedrals_without_hydrogen[slice(3,nphia,5)] /= one;
      }
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
        radius_sets = trimString(radius_sets);
    }
    else if(flag == "RADII")
        loadPrmtopData(prmtopFile,radii,natom,format);
    else if(flag == "SCREEN")
        loadPrmtopData(prmtopFile,screen,natom,format);
    else if(flag == "TITLE")
    {
        titles = getNextLine(prmtopFile);
        titles = trimString(titles);
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

void mmpbsa::SanderParm::loadPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
    using std::string;

    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file. loadPointers");

    std::valarray<size_t> pointers(size_t(0),31);
    loadPrmtopData(prmtopFile,pointers,31,format);


    //populate parameters
    size_t i=0;
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

void mmpbsa::SanderParm::loadSolventPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format)
{
    using std::string;

    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file. loadPointers");

    std::valarray<size_t> pointers(size_t(0),3);
    loadPrmtopData(prmtopFile,pointers,3,format);

    size_t i = 0;
    iptres = pointers[i++];///<   last residue that is considered part of solute (base 1 index)
    nspm = pointers[i++];///<     total number of molecules
    nspsol = pointers[i++];
}

template <class T> bool mmpbsa::SanderParm::rangeCheck(const std::valarray<T>& array,
    const T& min, const T& max)
{
    if(array.size() == 0)
        return false;

    if(array.max() > max)
        return false;

    if(array.min() < min)
        return false;

    return true;
}

template <class T> bool mmpbsa::SanderParm::bondCheck(const std::valarray<T>& array,
        const size_t& natoms, const size_t& nbonds, const size_t& ntypes,
        const size_t& atomsPerBond)
{
    std::ostringstream error;
    size_t maxi = (natoms-1)*3;
    if(array.size() != nbonds*(atomsPerBond+1))
    {
        error << "Incorrect number of bonds. Expected "
            << nbonds*(atomsPerBond+1);
        throw mmpbsa::SanderIOException(error,
            mmpbsa::INVALID_PRMTOP_DATA);
    }

    //check bond code range in respect to their specific bonds.
    for(size_t i = 0;i<nbonds*(atomsPerBond+1);i+=atomsPerBond+1)
    {
        for(size_t j = i;j<i+atomsPerBond;j++)
        {
            size_t absbnd = abs(array[j]);
            if(absbnd > maxi)
            {
                error << "Bond code " << absbnd << " exceeded max of "
                        << maxi << ". (i = " << i << ", j = " << j << ")";
                throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
            }
        }

        if(array[i+atomsPerBond] > ntypes)
        {
            error << "Bond code of " << array[i+atomsPerBond] <<  " exceeded number of types "
                    << "(" << ntypes << "). (i = " << i << ")";
            throw mmpbsa::SanderIOException(error,mmpbsa::INVALID_PRMTOP_DATA);
        }
    }

    return true;
}

void mmpbsa::SanderParm::loadPrmtopData(std::fstream& prmtopFile,
        std::valarray<std::string>& array,size_t size,const std::string& format)
{
    using std::valarray;
    using std::string;

    //ensure array is of the correct size (or exists).
    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file. loadPrmtopData");

    //use the format to obtain the array dimensions (in the 2-D sense).
    size_t numberOfColumns = 0;
    size_t columnWidth = 0;
    
    //sanity check for the following sscanf. no buffer overflows here.
    if(format.size() > 10 || *(format.begin()) != '(' || *(format.end()-1) != ')')
    {
        throw mmpbsa::SanderIOException("Improper format: "+format,mmpbsa::BROKEN_PRMTOP_FILE);
    }
    sscanf(format.c_str(),"%*c%d%*c%d",&numberOfColumns,&columnWidth);

    if(!mmpbsa_io::loadValarray(prmtopFile,array,size,columnWidth,numberOfColumns))
        throw mmpbsa::SanderIOException("Could not load Parameter Data",mmpbsa::BROKEN_PRMTOP_FILE);

}

template <class T> void mmpbsa::SanderParm::loadPrmtopData(std::fstream& prmtopFile,
        std::valarray<T>& array,size_t size,const std::string& format)
{
    using std::valarray;
    using std::string;

    //ensure array is of the correct size (or exists).
    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file. loadPrmtopData");

    //use the format to obtain the array dimensions (in the 2-D sense).
    size_t numberOfColumns = 0;
    size_t columnWidth = 0;
    
    //sanity check for the following sscanf. no buffer overflows here.
    if(format.size() > 10 || *(format.begin()) != '(' || *(format.end()-1) != ')')
    {
        throw mmpbsa::SanderIOException("Improper format: "+format,mmpbsa::BROKEN_PRMTOP_FILE);
    }
    sscanf(format.c_str(),"%*c%d%*c%d",&numberOfColumns,&columnWidth);

    if(!mmpbsa_io::loadValarray(prmtopFile,array,size,columnWidth,numberOfColumns))
        throw mmpbsa::SanderIOException("Could not load Parameter Data",mmpbsa::BROKEN_PRMTOP_FILE);

}



template <class T> void mmpbsa::SanderParm::loadPrmtopMaskedData(std::fstream& prmtopFile,
        std::valarray<T>& array,std::valarray<bool>& maskArray,size_t size,const std::string& format)
{
    using std::valarray;
    using std::string;

    //ensure array is of the correct size (or exists).
    if(!prmtopFile.good())
        throw mmpbsa::SanderIOException("Cannot read from file. loadPrmtopData");

    //use the format to obtain the array dimensions (in the 2-D sense).
    size_t numberOfColumns = 0;
    size_t columnWidth = 0;
    //sanity check for the following sscanf. no buffer overflows here.
    if(format.size() > 10 || *(format.begin()) != '(' || *(format.end()-1) != ')')
    {
        throw mmpbsa::SanderIOException("Improper format: "+format,mmpbsa::BROKEN_PRMTOP_FILE);
    }
    sscanf(format.c_str(),"%*c%d%*c%d",&numberOfColumns,&columnWidth);

    valarray<int> tmpArray(size);
    if(!mmpbsa_io::loadValarray(prmtopFile,tmpArray,size,columnWidth,numberOfColumns))
        throw mmpbsa::SanderIOException("Could not load Parameter Data",mmpbsa::BROKEN_PRMTOP_FILE);

    if(maskArray.size() != size)
        maskArray.resize(size);
    maskArray = tmpArray < 0;//the mask indicates whether or not the value is negative, which is a flag used in BondWalker.walk(...)
    tmpArray = abs(tmpArray);
    if(array.size() != size)
        array.resize(size);

    for(size_t i = 0;i<size;i++)
        array[i] = tmpArray[i];

}
