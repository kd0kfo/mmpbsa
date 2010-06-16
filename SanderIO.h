/**
 * SanderIO
 *
 * Handles IO for sander parameter files. Parameters are stored in the
 *      "SanderParm" class.
 *
 * Created by David Coss on June 16, 2010, 11:32 AM
 */

#include <exception>

#ifndef SANDERIO_H
#define	SANDERIO_H

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
    virtual ~SanderParm();


private:
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
};


class SanderIOException : public std::exception {
    SanderIOException();
    virtual ~SanderIOException();

};
#endif	/* SANDERIO_H */

