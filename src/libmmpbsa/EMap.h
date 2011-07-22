/**
 * @class mmpbsa::EMap
 * @brief Energy Abstraction class.
 *
 * EMap abstracts the energy data into one object, allowing calculations
 * to be performed on the energy data, including gathering statistics
 * and finding free energy. Using Emap, all of the data can be treated
 * as one mathematical object on which arithmetic can be performed, including
 * square root.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */


#ifndef EMAP_H
#define EMAP_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include <valarray>
#include "structs.h"
#include "globals.h"

//Forward declaration in namespace mmpbsa_utils
namespace mmpbsa_utils{ class XMLNode;}

namespace mmpbsa{

//forward declarations
class Vector;
class EmpEnerFun;

class EMap{
public:

    static const char DEFAULT_XML_TAG[];

    /**
     * Default Constructor for Energy Map Object
     * @see EMap(const EmpEnerFun* efun,const std::valarray<mmpbsa_t>& crds)
     *
     */
    EMap();

    /**
     * Copy constructor
     *
     * @param orig
     */
    EMap(const EMap& orig);

    /**
     * Creates an Energy Map object for the given Energy Function.
     * Energy maps are used to provide arithmetic operations to Energy data.
     * Operator overloads provided in the EMap class allow one to do arithmetic
     * on correspoding energy values in the Energy object.
     *
     * @param efun
     * @param crds
     */
    EMap(const EmpEnerFun* efun,const std::valarray<Vector>& crds);

    EMap(const std::vector<atom_t>& atoms, const mmpbsa::forcefield_t& ff, const std::valarray<Vector>& crds);

    ~EMap(){}

    /**
     * Sets all data equal to zero
     */
    void clear();

    const mmpbsa_t& set_elstat_solv(const mmpbsa_t& val){return elstat_solv = val;}
    const mmpbsa_t& set_area(const mmpbsa_t& val){return area = val;}
    const mmpbsa_t& set_sasol(const mmpbsa_t& val){return sasol = val;}
    /**
     * ostream output for convienent writing of energy data.
     *
     * Format:
     * BOND
     * ANGLE
     * DIHED
     * VDW14
     * ELE14
     * VDWAALS
     * VACELE
     *
     * @param rhs
     * @return
     */

    //operators
    EMap& operator=(const EMap& rhs);
    EMap operator+(const EMap& rhs)const;
    EMap& operator+=(const EMap& rhs);
    EMap operator-(const EMap& rhs)const;
    EMap& operator-=(const EMap& rhs);
    EMap operator*(const EMap& rhs)const;
    EMap& operator*=(const EMap& rhs);
    EMap operator/(const mmpbsa_t& rhs)const;
    EMap& operator/=(const mmpbsa_t& rhs);

    bool operator>(const EMap& rhs)const;
    bool operator==(const EMap& rhs)const;
    bool operator<(const EMap& rhs)const{return !(*this > rhs || *this == rhs);}
    bool operator>=(const EMap& rhs)const{return *this > rhs || *this == rhs;}
    bool operator<=(const EMap& rhs)const{return !(*this > rhs);}
    bool operator!=(const EMap& rhs)const{return !(*this == rhs);}

    friend std::ostream& operator<<(std::ostream& theStream, const mmpbsa::EMap& toWrite);

    /**
     * Creates an XML document node using the energy data. The name of the data
     * in captial letters is the name of each energy data tag.
     * 
     * @param name -- optional std::string title for the XMLNode (default = "energy")
     * @return
     */
    mmpbsa_utils::XMLNode* toXML(const std::string& name = DEFAULT_XML_TAG)const;

    /**
     * Creates an EMap object based on the given XML Node. Assumes that
     * the data tags are the children of xmlEnergy.
     * 
     * @param xmlEnergy
     * @return
     */
    static EMap loadXML(const mmpbsa_utils::XMLNode* xmlEnergy);

    /**
     * Gives the net electrostatic energy, which is the sum of ele14 and vacele.
     */
    mmpbsa_t total_elec_energy()const{return ele14 + vacele;}

    /**
     * Gives the net Van der Waals energy, which is the sum: vdw14 + vdwaals.
     */
    mmpbsa_t total_vdw_energy()const{return vdw14 + vdwaals;}

    /**
     * Gives the net internal energy, which is the sum: angle + bond + dihed
     */
    mmpbsa_t total_internal_energy()const{return angle + bond + dihed;}

    /**
     * Gives the total gas energy, which is the sum: total_elec_energy() + total_vdw_energy() + total_internal_energy()
     */
    mmpbsa_t total_gas_energy()const{return total_elec_energy() + total_vdw_energy() + total_internal_energy();}

    mmpbsa_t bond;///<Bond Energy
    mmpbsa_t angle;///<Angle Energy
    mmpbsa_t dihed;///<Dihedral Energy
    mmpbsa_t vdw14;///<Van der Waals Energy between 1st and 4th atom in dihedral
    mmpbsa_t ele14;///<Electrostatic Energy between 1st and 4th atom in dihedral
    mmpbsa_t vdwaals;///<Van der Walls Energy
    mmpbsa_t vacele;///<Electrostatic Energy
    mmpbsa_t elstat_solv;///<Electrostatic contribution to solvation energy.
    mmpbsa_t area;///<Solvant Accessible Surface Area of Molecule
    mmpbsa_t sasol;///<Surface area contribution to solvation energy.

    bool molsurf_failed;

};

}//end namespace mmpbsa
mmpbsa::EMap abs(const mmpbsa::EMap& obj);
mmpbsa::EMap sqrt(const mmpbsa::EMap& obj);


#endif
