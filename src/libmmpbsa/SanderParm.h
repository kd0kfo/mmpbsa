/**
 * @class mmpbsa::SanderParm
 * @brief Sander parameter container.
 *
 * SanderParm encapsulates data contained in a parmtop file.
 * The entire dataset of the parmtop file is stored in 
 * valarrays based on the catagory in the parmtop file.
 * The use of valarrays facilitates mathematical functions
 * which opearate on the topology data through valarray
 * operator overloads.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */

#ifndef SANDERPARM_H
#define	SANDERPARM_H

#include <string>
#include <sstream>
#include <valarray>
#include <iostream>

#include "mmpbsa_exceptions.h"
#include "globals.h"

namespace mmpbsa{
//forward declarations


/**
 * @class mmpbsa::SanderIOException
 * @brief Exception to be used when there is a problem reading Sander files.
 *
 * Created by David Coss, 2010
 */
class SanderIOException : public MMPBSAException {
public:
    SanderIOException(const std::string& error) : MMPBSAException( error){}
    
    SanderIOException(const std::string& error, const MMPBSAErrorTypes& errorType)
        : MMPBSAException(error,errorType){}
    SanderIOException(const std::ostringstream& error,const MMPBSAErrorTypes& errorType)
        : MMPBSAException(error,errorType){}

    const char* identifier(){return "SanderIO Error";}
};


class SanderParm {
public:
    /**
     * Default Constructor.
     * SanderParm stores all of the parameter data from the Sander Parmtop file.
     * This object may be passed along within the mmpbsa program to centralize
     * the storage of data, preventing redundancy.
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
     * Reads the AMBER Topology Parameter files (prmtop)
     *
     * @param file
     * @return
     */
    void raw_read_amber_parm(const std::string& file) throw (mmpbsa::SanderIOException);

    /**
     * Reads the AMBER Topology Parameter files (prmtop)
     *
     * @param prmtopFile
     */
    void raw_read_amber_parm(std::iostream& prmtopFile) throw (mmpbsa::SanderIOException);

    /**
     * Determines if the data loaded by loadSolventPointers makes sense.
     * Returns true if it passes and false if it fails but for a reason that is
     * not crash worthy, for example when warnings are sent to cerr.
     *
     * @return bool
     */
    bool sanityCheck() throw (SanderIOException);

    void initialize_arrays();

    SanderParm & operator=(const SanderParm& orig);

    //Here begin the many data types (public)

    //public to be more easily accessed. change later(?)
    //Data read from parameter file

    //size_t is used because the data from the prmtop file *actually* represents
    //an index in an array. Those values will actually be used as indices of
    //vectors and valarrays here. Therefore size_t is the logical, portable choice.
    size_t natom;///< total number of atoms
    size_t ntypes;///< total number of distinct atom types
    size_t nbonh;///< number of bonds containing hydrogen
    size_t mbona;///< number of bonds not containing hydrogen
    size_t ntheth;///< number of angles containing hydrogen
    size_t mtheta;///< number of angles not containing hydrogen
    size_t nphih;///< number of dihedrals containing hydrogen
    size_t mphia;///< number of dihedrals not containing hydrogen
    size_t nhparm;///< currently not used
    size_t nparm;///< set to 1 if LES is used
    size_t nnb;///< number of excluded atoms (=NEXT)
    size_t nres;///< number of residues
    size_t nbona;///< MBONA + number of constraint bonds
    size_t ntheta;///< MTHETA + number of constraint angles
    size_t nphia;///< MPHIA + number of constraint dihedrals
    size_t numbnd;///< number of unique bond types
    size_t numang;///< number of unique angle types
    size_t nptra;///< number of unique dihedral types
    size_t natyp;///< number of atom types in parameter file, see SOLTY below
    size_t nphb;///< number of distinct 10-12 hydrogen bond pair types
    size_t ifpert;///< set to 1 if perturbation info is to be read in
    size_t nbper;///< number of bonds to be perturbed
    size_t ngper;///< number of angles to be perturbed
    size_t ndper;///<number of dihedrals to be perturbed
    size_t mbper;///< number of bonds with atoms completely in perturbed group
    size_t mgper;///< number of angles with atoms completely in perturbed group
    size_t mdper;///< number of dihedrals with atoms completely in perturbed groups
    size_t ifbox;///< set to 1 if standard periodic box, 2 when truncated octahedral
    size_t nmxrs;///< number of atoms in the largest residue
    size_t ifcap;///< set to 1 if the CAP option from edit was specified
    size_t nextra;///< number of "extra points" (atom type of EP)

    //Data calculated after parsing the parameter.
    size_t iptres;///<   last residue that is considered part of solute (base 1 index)
    size_t nspm;///<     total number of molecules
    size_t nspsol;///<   the first solvent "molecule" (base 1 index)

    std::string titles;
    std::valarray<std::string> atom_names;
    std::valarray<mmpbsa_t> charges;
    std::valarray<mmpbsa_t> masses;
    std::valarray<size_t> atom_type_indices;//pointers are size_t
    std::valarray<size_t> number_excluded_atoms;
    std::valarray<size_t> nonbonded_parm_indices;//negative means 10-12 parameter arrays, positive means 6-12
    std::valarray<bool> nonbonded_parm_mask;//for use with finding the correct nonbonded parameter array.
    std::valarray<std::string> residue_labels;
    std::valarray<size_t> residue_pointers;//pointer means location in the array
       //not c++ pointer. This is an amber name from the prmtop file
       //(cf %FLAG RESIDUE_POINTER)
    std::valarray<mmpbsa_t> bond_force_constants;
    std::valarray<mmpbsa_t> bond_equil_values;
    std::valarray<mmpbsa_t> angle_force_constants;
    std::valarray<mmpbsa_t> angle_equil_values;
    std::valarray<mmpbsa_t> dihedral_force_constants;
    std::valarray<mmpbsa_t> dihedral_periodicities;
    std::valarray<mmpbsa_t> dihedral_phases;
    std::valarray<mmpbsa_t> soltys;//solubility?
    std::valarray<mmpbsa_t> lennard_jones_acoefs;
    std::valarray<mmpbsa_t> lennard_jones_bcoefs;
    std::valarray<size_t> bonds_inc_hydrogen;//some of the bond codes can be negative.
    std::valarray<size_t> bonds_without_hydrogen;
    std::valarray<size_t> angles_inc_hydrogen;
    std::valarray<size_t> angles_without_hydrogen;
    std::valarray<size_t> dihedrals_inc_hydrogen;
    std::valarray<size_t> dihedrals_without_hydrogen;
    std::valarray<bool> dihedral_h_mask;//some dihedral values can be negative to indicate types of bond. Values are true if the sander data is negative. false otherwise.
    std::valarray<bool> dihedral_mask;//some dihedral values can be negative to indicate types of bond. Values are true if the sander data is negative. false otherwise.
    std::valarray<size_t> excluded_atoms_list;
    std::valarray<mmpbsa_t> hbond_acoefs;
    std::valarray<mmpbsa_t> hbond_bcoefs;
    std::valarray<mmpbsa_t> hbcuts;
    std::valarray<std::string> amber_atom_types;
    std::valarray<std::string> tree_chain_classifications;
    std::valarray<size_t> join_array;
    std::valarray<size_t> irotats;
    std::string radius_sets;
    std::valarray<mmpbsa_t> radii;
    std::valarray<mmpbsa_t> screen;

    //solvent arrays
    std::valarray<size_t> atoms_per_molecule;
    std::valarray<size_t> box_dimensions;

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
    void parseParmtopFile(std::iostream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Loads the POINTERS from the file to the necessary variable in SanderParm
     * For SOLVENT_POINTERS, see loadSolventPointers
     *
     * @param prmtopFile
     * @param flag
     * @param format
     */
    void loadPointers(std::iostream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Loads the SOLVENT_POINTERS from the file to the necessary variable in SanderParm
     *
     * @param prmtopFile
     * @param flag
     * @param format
     */
    void loadSolventPointers(std::iostream& prmtopFile,const std::string& flag,
            const std::string& format);

    /**
     * Reads the data. An optional offset is given to shift data value, for example
     * with 1-indexed atom indices.
     *
     * @param prmtopFile
     * @param array pointer to new array. Contents of array will be replaced
     *      with new values. If there is a disagreement between array.size() and
     *      the parameter "size", array is deleted and replaced with new valarray(size)
     * @param size expected size of array.
     * @param offset If this value is non-zero, this offset will be added to the data.
     */
    template <class T> static void loadPrmtopData(std::iostream& prmtopFile,
        std::valarray<T>& array, size_t size,const std::string& format);

    static void loadPrmtopData(std::iostream& prmtopFile,std::valarray<std::string>& array,
    size_t size,const std::string& format);

    /**
     * Reads parmtop bond data that is specific to Dihedral Bonds. Unlike other
     * parmtop bond data, dihedral values can be negative. The negative sign is
     * a flag used for bond types in BondWalker. Therefore, there are two arrays
     * for dihedral bond data. One array, called array in this method, will hold
     * the absolute value of the dihedral bond data. The other array, called
     * maskArray here, is a boolean array, where true indicates the prmtop data
     * was negative. Programmatically, maskArray is useful because it could be used
     * as a mask_array, in the valarray sense of the word, to retrieve the respective
     * data quickly.
     *
     * @param prmtopFile
     * @param array
     * @param maskArray
     * @param size
     * @param format
     */
    template <class T> void loadPrmtopMaskedData(std::iostream& prmtopFile,
        std::valarray<T>& array,std::valarray<bool>& maskArray,size_t size,const std::string& format);

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
        const size_t& natoms, const size_t& nbonds, const size_t& ntype,
        const size_t& atomsPerBond);


};

/**
 * Looks up the residue of the atom and returns the label.
 * If a non-zero pointer is given to resIndex, its data is
 * replaced with the residue index (which is one-indexed!).
 */
std::string getResidueLabel(const size_t& index, const std::valarray<size_t>& res_ranges, const mmpbsa::SanderParm* parm, size_t* residueIndex);

}//end namespace mmpbsa



#endif	/* SANDERPARM_H */

