/**
 * SanderIO
 *
 * Handles IO for sander parameter files. Parameters are stored in the
 *      "SanderParm" class.
 *
 * Created by David Coss on June 16, 2010, 11:32 AM
 */


#ifndef SANDERIO_H
#define	SANDERIO_H

//Standard Includes
#include <string>
#include <fstream>
#include <valarray>
#include <streambuf>
#include <sstream>
#include <iostream>

//project specific stuff
#include "mmpbsa_exceptions.h"

class SanderParm {
public:
    /**
     * Default Constructor.
     *
     */
    SanderParm();

    /**
     * Copy Constructor.
     *
     * @param orig SanderParm to be copied
     */
    SanderParm(const SanderParm& orig);

    //Destructor
    ~SanderParm(){}

    /**
     * Reads the AMBER Topology Parameter files (prmtop) and
     * @param file
     * @return
     */
    void raw_read_amber_parm(std::string file);

    /**
     * Determines if the data loaded by loadSolventPointers makes sense.
     * Returns true if it passes and false if it fails but for a reason that is
     * not crash worthy, for example when warnings are sent to cerr.
     * 
     * @return bool
     */
    bool sanityCheck() throw (SanderIOException);

    SanderParm & operator=(const SanderParm& orig);


    //Here begin the many data types (public)

    //public to be more easily accessed. change later(?)
    //Data read from parameter file
    int natom;///< total number of atoms
    int ntypes;///< total number of distinct atom types
    int nbonh;///< number of bonds containing hydrogen
    int mbona;///< number of bonds not containing hydrogen
    int ntheth;///< number of angles containing hydrogen
    int mtheta;///< number of angles not containing hydrogen
    int nphih;///< number of dihedrals containing hydrogen
    int mphia;///< number of dihedrals not containing hydrogen
    int nhparm;///< currently not used
    int nparm;///< set to 1 if LES is used
    int nnb;///< number of excluded atoms (=NEXT)
    int nres;///< number of residues
    int nbona;///< MBONA + number of constraint bonds
    int ntheta;///< MTHETA + number of constraint angles
    int nphia;///< MPHIA + number of constraint dihedrals
    int numbnd;///< number of unique bond types
    int numang;///< number of unique angle types
    int nptra;///< number of unique dihedral types
    int natyp;///< number of atom types in parameter file, see SOLTY below
    int nphb;///< number of distinct 10-12 hydrogen bond pair types
    int ifpert;///< set to 1 if perturbation info is to be read in
    int nbper;///< number of bonds to be perturbed
    int ngper;///< number of angles to be perturbed
    int ndper;///<number of dihedrals to be perturbed
    int mbper;///< number of bonds with atoms completely in perturbed group
    int mgper;///< number of angles with atoms completely in perturbed group
    int mdper;///< number of dihedrals with atoms completely in perturbed groups
    int ifbox;///< set to 1 if standard periodic box, 2 when truncated octahedral
    int nmxrs;///< number of atoms in the largest residue
    int ifcap;///< set to 1 if the CAP option from edit was specified
    int nextra;///< number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    int iptres;///<   last residue that is considered part of solute (base 1 index)
    int nspm;///<     total number of molecules
    int nspsol;///<   the first solvent "molecule" (base 1 index)

    //Pointers to arrays that are read from the prmtop file
    //pointers are stored here in case the vectors are swapped around later
    //between difference classes. These arrays are LARGE; the goal here is to
    //save memory/time if that happens.
    std::string titles;
    std::valarray<std::string> atom_names;
    std::valarray<double> charges;
    std::valarray<double> masses;
    std::valarray<int> atom_type_indices;
    std::valarray<int> number_excluded_atoms;
    std::valarray<int> nonbonded_parm_indices;
    std::valarray<std::string> residue_labels;
    std::valarray<int> residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    std::valarray<double> bond_force_constants;
    std::valarray<double> bond_equil_values;
    std::valarray<double> angle_force_constants;
    std::valarray<double> angle_equil_values;
    std::valarray<double> dihedral_force_constants;
    std::valarray<double> dihedral_periodicities;
    std::valarray<double> dihedral_phases;
    std::valarray<double> soltys;//solubility?
    std::valarray<double> lennard_jones_acoefs;
    std::valarray<double> lennard_jones_bcoefs;
    std::valarray<int> bonds_inc_hydrogen;
    std::valarray<int> bonds_without_hydrogen;
    std::valarray<int> angles_inc_hydrogen;
    std::valarray<int> angles_without_hydrogen;
    std::valarray<int> dihedrals_inc_hydrogen;
    std::valarray<int> dihedrals_without_hydrogen;
    std::valarray<int> excluded_atoms_list;
    std::valarray<double> hbond_acoefs;
    std::valarray<double> hbond_bcoefs;
    std::valarray<double> hbcuts;
    std::valarray<std::string> amber_atom_types;
    std::valarray<std::string> tree_chain_classifications;
    std::valarray<int> join_array;
    std::valarray<int> irotats;
    std::string radius_sets;
    std::valarray<double> radii;
    std::valarray<double> screen;

    //solvent arrays
    std::valarray<double> atoms_per_molecule;
    std::valarray<double> box_dimensions;

private:
    /**
     * Fills the appropriate valarray with data, based on the provided flag and
     *      format. Format is used to get the number of columns and character
     *      width of the column. Columns are *not* whitespace delimited.
     *
     * @param
     * @param
     * @param
     * @return
     */
    void parseValarray(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Loads the POINTERS from the file to the necessary variable in SanderParm
     * For SOLVENT_POINTERS, see loadSolventPointers
     * 
     * @param prmtopFile
     * @param flag
     * @param format
     */
    void loadPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Loads the SOLVENT_POINTERS from the file to the necessary variable in SanderParm
     * 
     * @param prmtopFile
     * @param flag
     * @param format
     */
    void loadSolventPointers(std::fstream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Reads the data
     * 
     * @param prmtopFile
     * @param array pointer to new array. Contents of array will be replaced
     *      with new values. If there is a disagreement between array.size() and
     *      the parameter "size", array is deleted and replaced with new valarray(size)
     * @param size expected size of array.
     */
    template <class T> static void loadArray(std::fstream& prmtopFile,
        std::valarray<T>& array, int size,const std::string& format);

    /**
     * Loads data from the file to the given array. Data in array will be replaced
     *      by data read from file.
     * 
     * @param prmtopFile
     * @param array
     * @param size
     * @param format
     */
    static void loadArray(std::fstream& prmtopFile,
        std::valarray<std::string>& array, int size,const std::string& format);


    /**
     * Determines if all of the values in array fall between max and min.
     * array's template class must have operator> and operator<
     *
     * Returns true if everything is between max and min, false otherwise
     * 
     * @param array
     * @param min
     * @param max
     * @return
     */
    template <class T> static bool rangeCheck(const std::valarray<T>& array, const T& min,
        const T& max);

    template <class T> static bool bondCheck(const std::valarray<T>& array,
        const int& natoms, const int& nbonds, const int& ntype,
        const int& atomsPerBond);

};

std::string read_crds(std::fstream& crdFile, std::valarray<double>& crds);

void write_crds(const char* fileName,const std::valarray<double>& crds,
    const char* title = "");

/**
 * Removes white space at the beginning and end of a string.
 * @param str
 */
void trimString(std::string& bean);//why doesn't stl have this??

/**
     * Gets the next line with data, ie empty, whitespace lines are ignored.
     * @param file
     * @return
     */
static std::string getNextLine(std::fstream& file);

#endif	//SANDERIO_H

