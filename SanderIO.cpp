#include "SanderIO.h"

SanderParm::SanderParm() {
    natom = 0;///< total number of atoms
    ntypes = 0;///< total number of distinct atom types
    nbonh = 0;///< number of bonds containing hydrogen
    mbona = 0;///< number of bonds not containing hydrogen
    ntheth = 0;///< number of angles containing hydrogen
    mtheta = 0;///< number of angles not containing hydrogen
    nphih = 0;///< number of dihedrals containing hydrogen
    mphia = 0;///< number of dihedrals not containing hydrogen
    nhparm = 0;///< currently not used
    nparm = 0;///< set to 1 if LES is used
    nnb = 0;///< number of excluded atoms (=NEXT)
    nres = 0;///< number of residues
    nbona = 0;///< MBONA + number of constraint bonds
    ntheta = 0;///< MTHETA + number of constraint angles
    nphia = 0;///< MPHIA + number of constraint dihedrals
    numbnd = 0;///< number of unique bond types
    numang = 0;///< number of unique angle types
    nptra = 0;///< number of unique dihedral types
    natyp = 0;///< number of atom types in parameter file, see SOLTY below
    nphb = 0;///< number of distinct 10-12 hydrogen bond pair types
    ifpert = 0;///< set to 1 if perturbation info is to be read in
    nbper = 0;///< number of bonds to be perturbed
    ngper = 0;///< number of angles to be perturbed
    ndper = 0;///<number of dihedrals to be perturbed
    mbper = 0;///< number of bonds with atoms completely in perturbed group
    mgper = 0;///< number of angles with atoms completely in perturbed group
    mdper = 0;///< number of dihedrals with atoms completely in perturbed groups
    ifbox = 0;///< set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = 0;///< number of atoms in the largest residue
    ifcap = 0;///< set to 1 if the CAP option from edit was specified
    nextra = 0;///< number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    iptres = 0;///<   last residue that is considered part of solute (base 1 index)
    nspm = 0;///<     total number of molecules
    nspsol = 0;///<   the first solvent "molecule" (base 1 index)

    /*
    //pointers are stored here in case the vectors are swapped around later
    //between difference classes. These arrays are LARGE; the goal here is to
    //save memory/time if that happens.
    titles = "";
    atom_names = new std::valarray<std::string>;
    charges = new std::valarray<double>;
    masses = new std::valarray<double>;
    atom_type_indices = new std::valarray<int>;
    number_excluded_atoms = new std::valarray<int>;
    nonbonded_parm_indices = new std::valarray<int>;
    residue_labels = new std::valarray<std::string>;
    residue_pointers = new std::valarray<int>;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = new std::valarray<double>;
    bond_equil_values = new std::valarray<double>;
    angle_force_constants = new std::valarray<double>;
    dihedral_force_constants = new std::valarray<double>;
    angle_equil_values = new std::valarray<double>;
    dihedral_periodicities = new std::valarray<double>;
    dihedral_phases = new std::valarray<double>;
    soltys = new std::valarray<double>;//solubility?
    lennard_jones_acoefs = new std::valarray<double>;
    lennard_jones_bcoefs = new std::valarray<double>;
    bonds_inc_hydrogen = new std::valarray<int>;
    bonds_without_hydrogen = new std::valarray<int>;
    angles_inc_hydrogen = new std::valarray<int>;
    angles_without_hydrogen = new std::valarray<int>;
    dihedrals_inc_hydrogen = new std::valarray<int>;
    dihedrals_without_hydrogen = new std::valarray<int>;
    excluded_atoms_list = new std::valarray<int>;
    hbond_acoefs = new std::valarray<double>;
    hbond_bcoefs = new std::valarray<double>;
    hbcuts = new std::valarray<double>;
    amber_atom_types = new std::valarray<std::string>;
    tree_chain_classifications = new std::valarray<std::string>;
    join_array = new std::valarray<int>;
    irotats = new std::valarray<int>;
    radius_sets = new std::valarray<std::string>;
    radii = new std::valarray<double>;
    screen = new std::valarray<double>;/**/
}

SanderParm::SanderParm(const SanderParm& orig) {
    natom = orig.natom;///< total number of atoms
    ntypes = orig.ntypes;///< total number of distinct atom types
    nbonh = orig.nbonh;///< number of bonds containing hydrogen
    mbona = orig.mbona;///< number of bonds not containing hydrogen
    ntheth = orig.ntheth;///< number of angles containing hydrogen
    mtheta = orig.mtheta;///< number of angles not containing hydrogen
    nphih = orig.nphih;///< number of dihedrals containing hydrogen
    mphia = orig.mphia;///< number of dihedrals not containing hydrogen
    nhparm = orig.nhparm;///< currently not used
    nparm = orig.nparm;///< set to 1 if LES is used
    nnb = orig.nnb;///< number of excluded atoms (=NEXT)
    nres = orig.nres;///< number of residues
    nbona = orig.nbona;///< MBONA + number of constraint bonds
    ntheta = orig.ntheta;///< MTHETA + number of constraint angles
    nphia = orig.nphia;///< MPHIA + number of constraint dihedrals
    numbnd = orig.numbnd;///< number of unique bond types
    numang = orig.numang;///< number of unique angle types
    nptra = orig.nptra;///< number of unique dihedral types
    natyp = orig.natyp;///< number of atom types in parameter file, see SOLTY below
    nphb = orig.nphb;///< number of distinct 10-12 hydrogen bond pair types
    ifpert = orig.ifpert;///< set to 1 if perturbation info is to be read in
    nbper = orig.nbper;///< number of bonds to be perturbed
    ngper = orig.ngper;///< number of angles to be perturbed
    ndper = orig.ndper;///<number of dihedrals to be perturbed
    mbper = orig.mbper;///< number of bonds with atoms completely in perturbed group
    mgper = orig.mgper;///< number of angles with atoms completely in perturbed group
    mdper = orig.mdper;///< number of dihedrals with atoms completely in perturbed groups
    ifbox = orig.ifbox;///< set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = orig.nmxrs;///< number of atoms in the largest residue
    ifcap = orig.ifcap;///< set to 1 if the CAP option from edit was specified
    nextra = orig.nextra;///< number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    iptres = orig.iptres;///<   last residue that is considered part of solute (base 1 index)
    nspm = orig.nspm;///<     total number of molecules
    nspsol = orig.nspsol;///<   the first solvent "molecule" (base 1 index)

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
    charges = new std::valarray<double>(*(orig.charges));
    masses = new std::valarray<double>(*(orig.masses));
    atom_type_indices = new std::valarray<int>(*(orig.atom_type_indices));
    number_excluded_atoms = new std::valarray<int>(*(orig.number_excluded_atoms));
    nonbonded_parm_indices = new std::valarray<int>(*(orig.nonbonded_parm_indices));
    residue_labels = new std::valarray<std::string>(*(orig.residue_labels));
    residue_pointers = new std::valarray<int>(*(orig.residue_pointers));//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = new std::valarray<double>(*(orig.bond_force_constants));
    bond_equil_values = new std::valarray<double>(*(orig.bond_equil_values));
    angle_force_constants = new std::valarray<double>(*(orig.angle_force_constants));
    angle_equil_values = new std::valarray<double>(*(orig.angle_equil_values));
    dihedral_force_constants = new std::valarray<double>(*(orig.dihedral_force_constants));
    dihedral_periodicities = new std::valarray<double>(*(orig.dihedral_periodicities));
    dihedral_phases = new std::valarray<double>(*(orig.dihedral_phases));
    soltys = new std::valarray<double>(*(orig.soltys));//solubility?
    lennard_jones_acoefs = new std::valarray<double>(*(orig.lennard_jones_acoefs));
    lennard_jones_bcoefs = new std::valarray<double>(*(orig.lennard_jones_bcoefs));
    bonds_inc_hydrogen = new std::valarray<int>(*(orig.bonds_inc_hydrogen));
    bonds_without_hydrogen = new std::valarray<int>(*(orig.bonds_without_hydrogen));
    angles_inc_hydrogen = new std::valarray<int>(*(orig.angles_inc_hydrogen));
    angles_without_hydrogen = new std::valarray<int>(*(orig.angles_without_hydrogen));
    dihedrals_inc_hydrogen = new std::valarray<int>(*(orig.dihedrals_inc_hydrogen));
    dihedrals_without_hydrogen = new std::valarray<int>(*(orig.dihedrals_without_hydrogen));
    excluded_atoms_list = new std::valarray<int>(*(orig.excluded_atoms_list));
    hbond_acoefs = new std::valarray<double>(*(orig.hbond_acoefs));
    hbond_bcoefs = new std::valarray<double>(*(orig.hbond_bcoefs));
    hbcuts = new std::valarray<double>(*(orig.hbcuts));
    amber_atom_types = new std::valarray<std::string>(*(orig.amber_atom_types));
    tree_chain_classifications = new std::valarray<std::string>(*(orig.tree_chain_classifications));
    join_array = new std::valarray<int>(*(orig.join_array));
    irotats = new std::valarray<int>(*(orig.irotats));
    radius_sets = new std::valarray<std::string>(*(orig.radius_sets));
    radii = new std::valarray<double>(*(orig.radii));
    screen = new std::valarray<double>(*(orig.screen));*/

}

SanderParm & SanderParm::operator=(const SanderParm& orig)
{
    if(this == &orig)
        return *this;

    natom = orig.natom;///< total number of atoms
    ntypes = orig.ntypes;///< total number of distinct atom types
    nbonh = orig.nbonh;///< number of bonds containing hydrogen
    mbona = orig.mbona;///< number of bonds not containing hydrogen
    ntheth = orig.ntheth;///< number of angles containing hydrogen
    mtheta = orig.mtheta;///< number of angles not containing hydrogen
    nphih = orig.nphih;///< number of dihedrals containing hydrogen
    mphia = orig.mphia;///< number of dihedrals not containing hydrogen
    nhparm = orig.nhparm;///< currently not used
    nparm = orig.nparm;///< set to 1 if LES is used
    nnb = orig.nnb;///< number of excluded atoms (=NEXT)
    nres = orig.nres;///< number of residues
    nbona = orig.nbona;///< MBONA + number of constraint bonds
    ntheta = orig.ntheta;///< MTHETA + number of constraint angles
    nphia = orig.nphia;///< MPHIA + number of constraint dihedrals
    numbnd = orig.numbnd;///< number of unique bond types
    numang = orig.numang;///< number of unique angle types
    nptra = orig.nptra;///< number of unique dihedral types
    natyp = orig.natyp;///< number of atom types in parameter file, see SOLTY below
    nphb = orig.nphb;///< number of distinct 10-12 hydrogen bond pair types
    ifpert = orig.ifpert;///< set to 1 if perturbation info is to be read in
    nbper = orig.nbper;///< number of bonds to be perturbed
    ngper = orig.ngper;///< number of angles to be perturbed
    ndper = orig.ndper;///<number of dihedrals to be perturbed
    mbper = orig.mbper;///< number of bonds with atoms completely in perturbed group
    mgper = orig.mgper;///< number of angles with atoms completely in perturbed group
    mdper = orig.mdper;///< number of dihedrals with atoms completely in perturbed groups
    ifbox = orig.ifbox;///< set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = orig.nmxrs;///< number of atoms in the largest residue
    ifcap = orig.ifcap;///< set to 1 if the CAP option from edit was specified
    nextra = orig.nextra;///< number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    iptres = orig.iptres;///<   last residue that is considered part of solute (base 1 index)
    nspm = orig.nspm;///<     total number of molecules
    nspsol = orig.nspsol;///<   the first solvent "molecule" (base 1 index)

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
    charges = new std::valarray<double>(*(orig.charges));
    masses = new std::valarray<double>(*(orig.masses));
    atom_type_indices = new std::valarray<int>(*(orig.atom_type_indices));
    number_excluded_atoms = new std::valarray<int>(*(orig.number_excluded_atoms));
    nonbonded_parm_indices = new std::valarray<int>(*(orig.nonbonded_parm_indices));
    residue_labels = new std::valarray<std::string>(*(orig.residue_labels));
    residue_pointers = new std::valarray<int>(*(orig.residue_pointers));//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    bond_force_constants = new std::valarray<double>(*(orig.bond_force_constants));
    bond_equil_values = new std::valarray<double>(*(orig.bond_equil_values));
    angle_force_constants = new std::valarray<double>(*(orig.angle_force_constants));
    angle_equil_values = new std::valarray<double>(*(orig.angle_equil_values));
    dihedral_force_constants = new std::valarray<double>(*(orig.dihedral_force_constants));
    dihedral_periodicities = new std::valarray<double>(*(orig.dihedral_periodicities));
    dihedral_phases = new std::valarray<double>(*(orig.dihedral_phases));
    soltys = new std::valarray<double>(*(orig.soltys));//solubility?
    lennard_jones_acoefs = new std::valarray<double>(*(orig.lennard_jones_acoefs));
    lennard_jones_bcoefs = new std::valarray<double>(*(orig.lennard_jones_bcoefs));
    bonds_inc_hydrogen = new std::valarray<int>(*(orig.bonds_inc_hydrogen));
    bonds_without_hydrogen = new std::valarray<int>(*(orig.bonds_without_hydrogen));
    angles_inc_hydrogen = new std::valarray<int>(*(orig.angles_inc_hydrogen));
    angles_without_hydrogen = new std::valarray<int>(*(orig.angles_without_hydrogen));
    dihedrals_inc_hydrogen = new std::valarray<int>(*(orig.dihedrals_inc_hydrogen));
    dihedrals_without_hydrogen = new std::valarray<int>(*(orig.dihedrals_without_hydrogen));
    excluded_atoms_list = new std::valarray<int>(*(orig.excluded_atoms_list));
    hbond_acoefs = new std::valarray<double>(*(orig.hbond_acoefs));
    hbond_bcoefs = new std::valarray<double>(*(orig.hbond_bcoefs));
    hbcuts = new std::valarray<double>(*(orig.hbcuts));
    amber_atom_types = new std::valarray<std::string>(*(orig.amber_atom_types));
    tree_chain_classifications = new std::valarray<std::string>(*(orig.tree_chain_classifications));
    join_array = new std::valarray<int>(*(orig.join_array));
    irotats = new std::valarray<int>(*(orig.irotats));
    radius_sets = new std::valarray<std::string>(*(orig.radius_sets));
    radii = new std::valarray<double>(*(orig.radii));
    screen = new std::valarray<double>(*(orig.screen));*/

    return *this;
}


SanderParm SanderParm::raw_read_amber_parm(std::string file)
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
    SanderParm returnMe;//this will store the data read from file. This takes
        //alot of memory. So it is create only *after* the file is correctly
        //opened.

    //Fill Arrays
    while(prmtopFile.good())
    {
        currentLine = getNextLine(prmtopFile);//should be FLAG
        if(!strcmp(currentLine.substr(0,5).c_str(),"%FLAG"))
            throw SanderIOException(file.append(" is malformed. %FLAG is missing.").c_str());
        else
        {
            flag = currentLine.substr(5);
            trimString(flag);
        }

        currentLine = getNextLine(prmtopFile);//should be FORMAT
        if(!strcmp(currentLine.substr(0,7).c_str(),"%FORMAT"))
            throw SanderIOException(file.append(" is malformed. %FORMAT is missing.").c_str());
        else
        {
            format = currentLine.substr(7);
            trimString(format);
        }

        returnMe.parseValarray(prmtopFile,flag.c_str(),format.c_str());
    }



}

void SanderParm::parseValarray(std::fstream& prmtopFile,const char* flag,
            const char* format)
{
    //is the file still able to be read?
    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file.");

    //determine which flag is being used.
    //Sorry, but this simply has to be a long if..else chain.
    //If the flag is not known, the function returns, because some prmtop files
    //have parameters that are not used here, such as solvents.
    if(strcmp(flag,"POINTERS"))
        return loadPointers(prmtopFile,flag,format);
    else if(strcmp(flag,"ATOM_NAME"))
        loadArray(prmtopFile,atom_names,natom);
    else if(strcmp(flag,"CHARGE"))
        loadArray(prmtopFile,charges,natom);
    else if(strcmp(flag,"MASS"))
        loadArray(prmtopFile,masses,natom);
    else if(strcmp(flag,"ATOM_TYPE_INDEX"))
        loadArray(prmtopFile,atom_type_indices,natom);
    else if(strcmp(flag,"NUMBER_EXCLUDED_ATOMS"))
        loadArray(prmtopFile,number_excluded_atoms,natom);
    else if(strcmp(flag,"NONBONDED_PARM_INDEX"))
        loadArray(prmtopFile,nonbonded_parm_indices,ntypes*ntypes);
    else if(strcmp(flag,"RESIDUE_LABEL"))
        loadArray(prmtopFile,residue_labels,nres);
    else if(strcmp(flag,"RESIDUE_POINTER"))
        loadArray(prmtopFile,residue_pointers,nres);
    else if(strcmp(flag,"BOND_FORCE_CONSTANT"))
        loadArray(prmtopFile,bond_force_constants,numbnd);
    else if(strcmp(flag,"BOND_EQUIL_VALUE"))
        loadArray(prmtopFile,bond_equil_values,numbnd);
    else if(strcmp(flag,"ANGLE_FORCE_CONSTANT"))
        loadArray(prmtopFile,angle_force_constants,numang);
    else if(strcmp(flag,"ANGLE_EQUIL_VALUE"))
        loadArray(prmtopFile,angle_equil_values,numang);
    else if(strcmp(flag,"DIHEDRAL_FORCE_CONSTANT"))
        loadArray(prmtopFile,dihedral_force_constants,nptra);
    else if(strcmp(flag,"DIHEDRAL_PERIODICITY"))
        loadArray(prmtopFile,dihedral_periodicities,nptra);
    else if(strcmp(flag,"DIHEDRAL_PHASE"))
        loadArray(prmtopFile,dihedral_phases,nptra);
    else if(strcmp(flag,"SOLTY"))
        loadArray(prmtopFile,soltys,natyp);
    else if(strcmp(flag,"LENNARD_JONES_ACOEF"))
        loadArray(prmtopFile,lennard_jones_acoefs,(int) 0.5*ntypes*(ntypes+1));
    else if(strcmp(flag,"LENNARD_JONES_BCOEF"))
        loadArray(prmtopFile,lennard_jones_bcoefs,(int) 0.5*ntypes*(ntypes+1));
    else if(strcmp(flag,"BONDS_INC_HYDROGEN"))
        loadArray(prmtopFile,bonds_inc_hydrogen,nbonh);
    else if(strcmp(flag,"BONDS_WITHOUT_HYDROGEN"))
        loadArray(prmtopFile,bonds_without_hydrogen,nbona);
    else if(strcmp(flag,"ANGLES_INC_HYDROGEN"))
        loadArray(prmtopFile,angles_inc_hydrogen,ntheth);
    else if(strcmp(flag,"ANGLES_WITHOUT_HYDROGEN"))
        loadArray(prmtopFile,angles_without_hydrogen,ntheta);
    else if(strcmp(flag,"DIHEDRALS_INC_HYDROGEN"))
        loadArray(prmtopFile,dihedrals_inc_hydrogen,nphih);
    else if(strcmp(flag,"DIHEDRALS_WITHOUT_HYDROGEN"))
        loadArray(prmtopFile,dihedrals_without_hydrogen,nphia);
    else if(strcmp(flag,"EXCLUDED_ATOMS_LIST"))
        loadArray(prmtopFile,excluded_atoms_list,nnb);
    else if(strcmp(flag,"HBOND_ACOEF"))
        loadArray(prmtopFile,hbond_acoefs,nphb);
    else if(strcmp(flag,"HBOND_BCOEF"))
        loadArray(prmtopFile,hbond_bcoefs,nphb);
    else if(strcmp(flag,"HBCUT"))
        loadArray(prmtopFile,hbcuts,nphb);
    else if(strcmp(flag,"AMBER_ATOM_TYPE"))
        loadArray(prmtopFile,amber_atom_types,natom);
    else if(strcmp(flag,"TREE_CHAIN_CLASSIFICATION"))
        loadArray(prmtopFile,tree_chain_classifications,natom);
    else if(strcmp(flag,"JOIN_ARRAY"))
        loadArray(prmtopFile,join_array,natom);
    else if(strcmp(flag,"IROTAT"))
        loadArray(prmtopFile,irotats,natom);
    else if(strcmp(flag,"RADIUS_SET"))
        loadArray(prmtopFile,radius_sets,natom);
    else if(strcmp(flag,"RADII"))
        loadArray(prmtopFile,radii,natom);
    else if(strcmp(flag,"SCREEN"))
        loadArray(prmtopFile,screen,natom);
    else
        std::cout << flag << " was not parsed." << std::endl;


}

void SanderParm::loadPointers(std::fstream& prmtopFile,const char* flag,
            const char* format)
{
    using std::string;
    
    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file. loadPointers");

    std::valarray<int> pointers(0,31);
    char next;//first character of next line
    int i = 0;
    do
    {
        //"peek" to see if are finished parsing, ie next flag begins.
        prmtopFile.get(next);
        if(next == '%')
        {
            prmtopFile.unget();

            //check to see if the correct number of data values have been read.
            if(i != pointers.size() - 1)
                throw SanderIOException("SanderParm has read to many or too "
                    "few data values.");

            return;
        }

        prmtopFile.unget();//undo the "peek"
        string currentLine = getNextLine(prmtopFile);
        
        std::stringstream ss(currentLine);//will be used to tokenize the elements.
        string tempString;//Non-whitespace strings will be stored here
        while (ss >> tempString)
        {
            double numericalValue = atoi(tempString.c_str());

            if(numericalValue == 0 && tempString[0] != '0')
                throw SanderIOException(tempString.append(" is not a int.").c_str());
            pointers[i++] = numericalValue;
        }

    }while(prmtopFile.good());


    //populate parameters
    i=0;
    natom = pointers[i++];///< total number of atoms
    ntypes = pointers[i++];///< total number of distinct atom types
    nbonh = pointers[i++];///< number of bonds containing hydrogen
    mbona = pointers[i++];///< number of bonds not containing hydrogen
    ntheth = pointers[i++];///< number of angles containing hydrogen
    mtheta = pointers[i++];///< number of angles not containing hydrogen
    nphih = pointers[i++];///< number of dihedrals containing hydrogen
    mphia = pointers[i++];///< number of dihedrals not containing hydrogen
    nhparm = pointers[i++];///< currently not used
    nparm = pointers[i++];///< set to 1 if LES is used
    nnb = pointers[i++];///< number of excluded atoms (=NEXT)
    nres = pointers[i++];///< number of residues
    nbona = pointers[i++];///< MBONA + number of constraint bonds
    ntheta = pointers[i++];///< MTHETA + number of constraint angles
    nphia = pointers[i++];///< MPHIA + number of constraint dihedrals
    numbnd = pointers[i++];///< number of unique bond types
    numang = pointers[i++];///< number of unique angle types
    nptra = pointers[i++];///< number of unique dihedral types
    natyp = pointers[i++];///< number of atom types in parameter file, see SOLTY below
    nphb = pointers[i++];///< number of distinct 10-12 hydrogen bond pair types
    ifpert = pointers[i++];///< set to 1 if perturbation info is to be read in
    nbper = pointers[i++];///< number of bonds to be perturbed
    ngper = pointers[i++];///< number of angles to be perturbed
    ndper = pointers[i++];///<number of dihedrals to be perturbed
    mbper = pointers[i++];///< number of bonds with atoms completely in perturbed group
    mgper = pointers[i++];///< number of angles with atoms completely in perturbed group
    mdper = pointers[i++];///< number of dihedrals with atoms completely in perturbed groups
    ifbox = pointers[i++];///< set to 1 if standard periodic box, 2 when truncated octahedral
    nmxrs = pointers[i++];///< number of atoms in the largest residue
    ifcap = pointers[i++];///< set to 1 if the CAP option from edit was specified
    nextra = pointers[i++];///< number of "extra points" (atom type of EP)


}

template <class T> void SanderParm::loadArray(std::fstream& prmtopFile,
        std::valarray<T> array,int size)
{
    using std::valarray;
    using std::string;

    //ensure array is of the correct size (or exists).
    if(!prmtopFile.good())
        throw SanderIOException("Cannot read from file. loadArray");

    if(array.size() != size)
        array.resize(size);

    char next;//first character of next line
    int currentIndex = 0;
    do
    {
        //"peek" to see if are finished parsing, ie next flag begins.
        prmtopFile.get(next);
        if(next == '%')
        {
            prmtopFile.unget();

            //check to see if the correct number of data values have been read.
            if(currentIndex != array.size() - 1)
                throw SanderIOException("SanderParm has read to many or too "
                    "few data values.");

            return;
        }

        prmtopFile.unget();//undo the "peek"
        string currentLine = getNextLine(prmtopFile);
        
        std::stringstream ss(currentLine);//will be used to tokenize the elements.
        string tempString;//Non-whitespace strings will be stored here
        while (ss >> tempString)
        {
            double numericalValue = atof(tempString.c_str());

            if(numericalValue == 0 && tempString[0] != '0')
                throw SanderIOException(tempString.append(" is not a double.").c_str());
            array[currentIndex++] = numericalValue;
        }

    }while(prmtopFile.good());

    //check to see if the correct number of data values have been read.
    if(currentIndex != array.size() - 1)
        throw SanderIOException("SanderParm has read to many or too "
            "few data values.");
}

std::string SanderParm::getNextLine(std::fstream& file)
{
    if(!file.good())
        throw SanderIOException("Could not read from file");

    std::string returnMe;
    getline(file,returnMe);
    return returnMe;
}

void trimString(std::string& bean)
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


