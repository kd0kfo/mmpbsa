#ifndef EMAP_H
#define EMAP_H

#include <valarray>
#include "mmpbsa_utils.h"
#include "EmpEnerFun.h"

namespace mmpbsa{

class EMap{
public:
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
    EMap(const EmpEnerFun* efun,const std::valarray<mmpbsa_t>& crds);

    ~EMap(){}

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
    friend EMap abs(const mmpbsa::EMap& obj);
    friend EMap sqrt(const mmpbsa::EMap& obj);

    /**
     * Creates an XML document node using the energy data. The name of the data
     * in captial letters is the name of each energy data tag.
     * 
     * @param name -- optional std::string title for the XMLNode (default = "energy")
     * @return
     */
    mmpbsa_utils::XMLNode* toXML(const std::string& name = "energy")const;

    /**
     * Creates an EMap object based on the given XML Node.
     * 
     * @param xmlEnergy
     * @return
     */
    static EMap loadXML(const mmpbsa_utils::XMLNode* xmlEnergy);

protected:
    mmpbsa_t bond;
    mmpbsa_t angle;
    mmpbsa_t dihed;
    mmpbsa_t vdw14;
    mmpbsa_t ele14;
    mmpbsa_t vdwaals;
    mmpbsa_t vacele;
    mmpbsa_t elstat_solv;
    mmpbsa_t area;
    mmpbsa_t sasol;

};

}//end namespace mmpbsa


#endif
