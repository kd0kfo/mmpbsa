#include "EMap.h"

mmpbsa::EMap::EMap()
{
    bond = 0;
    angle = 0;
    dihed = 0;
    vdw14 = 0;
    ele14 = 0;
    vdwaals = 0;
    vacele = 0;
    elstat_solv = 0;
    area = 0;
    sasol = 0;
}


mmpbsa::EMap::EMap(const mmpbsa::EMap& orig)
{
    bond = orig.bond;
    angle = orig.angle;
    dihed = orig.dihed;
    vdw14 = orig.vdw14;
    ele14 = orig.ele14;
    vdwaals = orig.vdwaals;
    vacele = orig.vacele;
    elstat_solv = orig.elstat_solv;
    area = orig.area;
    sasol = orig.sasol;
}

mmpbsa::EMap::EMap(const mmpbsa::EmpEnerFun* efun, const std::valarray<mmpbsa_t>& crds)
{
    if(efun == 0)
        throw mmpbsa::MMPBSAException("An attempt was made to create an EMap with a null"
                " EmpEnerFun pointer.",mmpbsa::UNKNOWN_ERROR);

    bond = efun->total_bond_energy(crds);
    angle = efun->total_angle_energy(crds);
    dihed = efun->total_dihedral_energy(crds);
    vdw14 = efun->total_vdw14_energy(crds);
    ele14 = efun->total_elstat14_energy(crds);
    vdwaals = efun->total_vdwaals_energy(crds);
    vacele = efun->total_elstat_energy(crds);
    elstat_solv = 0;
    area = 0;
    sasol = 0;
}

namespace mmpbsa{
std::ostream& operator<<(std::ostream& theStream, const mmpbsa::EMap& toWrite)
{
    std::streamsize prevPrecision = theStream.precision();
    std::ios::fmtflags prevFloatfield = theStream.floatfield;
    theStream.precision(8);
    theStream.setf(std::ios::fixed,std::ios::floatfield);
    theStream << "BOND " << toWrite.bond << std::endl;
    theStream << "ANGLE " << toWrite.angle << std::endl;
    theStream << "DIHED " << toWrite.dihed << std::endl;
    theStream << "VDW14 " << toWrite.vdw14 << std::endl;
    theStream << "ELE14 " << toWrite.ele14 << std::endl;
    theStream << "VACELE " << toWrite.vacele << std::endl;
    theStream << "VDWAALS " << toWrite.vdwaals << std::endl;
    theStream << "PBSOL " << toWrite.elstat_solv << std::endl;
    theStream << "SASOL " << toWrite.sasol << std::endl;
    theStream << "area " << toWrite.area;
    theStream.precision(prevPrecision);
    theStream.setf(prevFloatfield,std::ios::floatfield);
    return theStream;
}

mmpbsa::EMap abs(const mmpbsa::EMap& obj)
{
    mmpbsa::EMap returnMe = obj;
    returnMe.angle = std::abs(returnMe.angle);
    returnMe.area = std::abs(returnMe.area);
    returnMe.bond = std::abs(returnMe.bond);
    returnMe.dihed = std::abs(returnMe.dihed);
    returnMe.ele14 = std::abs(returnMe.ele14);
    returnMe.elstat_solv = std::abs(returnMe.elstat_solv);
    returnMe.sasol = std::abs(returnMe.sasol);
    returnMe.vacele = std::abs(returnMe.vacele);
    returnMe.vdw14 = std::abs(returnMe.vdw14);
    returnMe.vdwaals = std::abs(returnMe.vdwaals);

    return returnMe;
}

mmpbsa::EMap sqrt(const mmpbsa::EMap& obj)
{
    mmpbsa::EMap returnMe = obj;
    returnMe.angle = std::sqrt(returnMe.angle);
    returnMe.area = std::sqrt(returnMe.area);
    returnMe.bond = std::sqrt(returnMe.bond);
    returnMe.dihed = std::sqrt(returnMe.dihed);
    returnMe.ele14 = std::sqrt(returnMe.ele14);
    returnMe.elstat_solv = std::sqrt(returnMe.elstat_solv);
    returnMe.sasol = std::sqrt(returnMe.sasol);
    returnMe.vacele = std::sqrt(returnMe.vacele);
    returnMe.vdw14 = std::sqrt(returnMe.vdw14);
    returnMe.vdwaals = std::sqrt(returnMe.vdwaals);

    return returnMe;
}

}

mmpbsa::EMap& mmpbsa::EMap::operator=(const mmpbsa::EMap& rhs)
{
    if(this == &rhs)
        return *this;

    bond = rhs.bond;
    angle = rhs.angle;
    dihed = rhs.dihed;
    vdw14 = rhs.vdw14;
    ele14 = rhs.ele14;
    vdwaals = rhs.vdwaals;
    vacele = rhs.vacele;
    elstat_solv = rhs.elstat_solv;
    area = rhs.area;
    sasol = rhs.sasol;
    return *this;
}

mmpbsa::EMap mmpbsa::EMap::operator+(const mmpbsa::EMap& rhs)const
{
    mmpbsa::EMap returnMe;
    returnMe.bond = rhs.bond+bond;//i know. rhs. But addition is commutative...
    returnMe.angle = rhs.angle+angle;
    returnMe.dihed = rhs.dihed+dihed;
    returnMe.vdw14 = rhs.vdw14+vdw14;
    returnMe.ele14 = rhs.ele14+ele14;
    returnMe.vdwaals = rhs.vdwaals+vdwaals;
    returnMe.vacele = rhs.vacele+vacele;
    returnMe.elstat_solv = rhs.elstat_solv+elstat_solv;
    returnMe.area = rhs.area+area;
    returnMe.sasol = rhs.sasol+sasol;
    return returnMe;
}

mmpbsa::EMap& mmpbsa::EMap::operator+=(const mmpbsa::EMap& rhs)
{
    bond += rhs.bond;
    angle += rhs.angle;
    dihed += rhs.dihed;
    vdw14 += rhs.vdw14;
    ele14 += rhs.ele14;
    vdwaals += rhs.vdwaals;
    vacele += rhs.vacele;
    elstat_solv += rhs.elstat_solv;
    area += rhs.area;
    sasol += rhs.sasol;
    return *this;
}

mmpbsa::EMap mmpbsa::EMap::operator-(const mmpbsa::EMap& rhs)const
{
    EMap returnMe;
    returnMe.bond = bond-rhs.bond;
    returnMe.angle = angle-rhs.angle;
    returnMe.dihed = dihed-rhs.dihed;
    returnMe.vdw14 = vdw14-rhs.vdw14;
    returnMe.ele14 = ele14-rhs.ele14;
    returnMe.vdwaals = vdwaals-rhs.vdwaals;
    returnMe.vacele = vacele-rhs.vacele;
    returnMe.elstat_solv = elstat_solv-rhs.elstat_solv;
    returnMe.area = area-rhs.area;
    returnMe.sasol = sasol-rhs.sasol;
    return returnMe;
}


mmpbsa::EMap& mmpbsa::EMap::operator-=(const mmpbsa::EMap& rhs)
{
    bond -= rhs.bond;
    angle -= rhs.angle;
    dihed -= rhs.dihed;
    vdw14 -= rhs.vdw14;
    ele14 -= rhs.ele14;
    vdwaals -= rhs.vdwaals;
    vacele -= rhs.vacele;
    elstat_solv -= rhs.elstat_solv;
    area -= rhs.area;
    sasol -= rhs.sasol;

    return *this;
}


mmpbsa::EMap mmpbsa::EMap::operator*(const mmpbsa::EMap& rhs)const
{
    mmpbsa::EMap returnMe;
    returnMe.bond = bond*rhs.bond;
    returnMe.angle = angle*rhs.angle;
    returnMe.dihed = dihed*rhs.dihed;
    returnMe.vdw14 = vdw14*rhs.vdw14;
    returnMe.ele14 = ele14*rhs.ele14;
    returnMe.vdwaals = vdwaals*rhs.vdwaals;
    returnMe.vacele = vacele*rhs.vacele;
    returnMe.elstat_solv = elstat_solv*rhs.elstat_solv;
    returnMe.area = area*rhs.area;
    returnMe.sasol = sasol*rhs.sasol;
    return returnMe;
}

mmpbsa::EMap& mmpbsa::EMap::operator*=(const mmpbsa::EMap& rhs)
{
    bond *= rhs.bond;
    angle *= rhs.angle;
    dihed *= rhs.dihed;
    vdw14 *= rhs.vdw14;
    ele14 *= rhs.ele14;
    vdwaals *= rhs.vdwaals;
    vacele *= rhs.vacele;
    elstat_solv *= rhs.elstat_solv;
    area *= rhs.area;
    sasol *= rhs.sasol;

    return *this;
}

mmpbsa::EMap mmpbsa::EMap::operator/(const mmpbsa_t& rhs)const
{
    mmpbsa::EMap returnMe;
    returnMe.bond = bond/rhs;
    returnMe.angle = angle/rhs;
    returnMe.dihed = dihed/rhs;
    returnMe.vdw14 = vdw14/rhs;
    returnMe.ele14 = ele14/rhs;
    returnMe.vdwaals = vdwaals/rhs;
    returnMe.vacele = vacele/rhs;
    returnMe.elstat_solv = elstat_solv/rhs;
    returnMe.area = area/rhs;
    returnMe.sasol = sasol/rhs;
    return returnMe;
}

mmpbsa::EMap& mmpbsa::EMap::operator/=(const mmpbsa_t& rhs)
{
    bond /= rhs;
    angle /= rhs;
    dihed /= rhs;
    vdw14 /= rhs;
    ele14 /= rhs;
    vdwaals /= rhs;
    vacele /= rhs;
    elstat_solv /= rhs;
    area /= rhs;
    sasol /= rhs;

    return *this;
}

bool mmpbsa::EMap::operator>(const mmpbsa::EMap& rhs)const
{
    return bond > rhs.bond ||
            angle > rhs.angle ||
            dihed > rhs.dihed ||
            vdw14 > rhs.vdw14 ||
            ele14 > rhs.ele14 ||
            vdwaals > rhs.vdwaals ||
            vacele > rhs.vacele ||
            elstat_solv > rhs.elstat_solv ||
            area > rhs.area ||
            sasol > rhs.sasol;
}

bool mmpbsa::EMap::operator==(const mmpbsa::EMap& rhs)const
{
    return bond == rhs.bond &&
            angle == rhs.angle &&
            dihed == rhs.dihed &&
            vdw14 == rhs.vdw14 &&
            ele14 == rhs.ele14 &&
            vdwaals == rhs.vdwaals &&
            vacele == rhs.vacele &&
            elstat_solv == rhs.elstat_solv &&
            area == rhs.area &&
            sasol == rhs.sasol;
}


mmpbsa_utils::XMLNode* mmpbsa::EMap::toXML(const std::string& name)const
{
    mmpbsa_utils::XMLNode* theNode = new mmpbsa_utils::XMLNode(name);

    std::ostringstream value;
    value << MMPBSA_FORMAT << bond;
    theNode->insertChild("BOND",value.str());value.str("");
    value << MMPBSA_FORMAT << angle;
    theNode->insertChild("ANGLE",value.str());value.str("");
    value << MMPBSA_FORMAT << dihed;
    theNode->insertChild("DIHED",value.str());value.str("");
    value << MMPBSA_FORMAT << vdw14;
    theNode->insertChild("VDW14",value.str());value.str("");
    value << MMPBSA_FORMAT << ele14;
    theNode->insertChild("ELE14",value.str());value.str("");
    value << MMPBSA_FORMAT << vacele;
    theNode->insertChild("VACELE",value.str());value.str("");
    value << MMPBSA_FORMAT << vdwaals;
    theNode->insertChild("VDWAALS",value.str());value.str("");
    value << MMPBSA_FORMAT << elstat_solv;
    theNode->insertChild("PBSOL",value.str());value.str("");
    value << MMPBSA_FORMAT << sasol;
    theNode->insertChild("SASOL",value.str());value.str("");
    value << MMPBSA_FORMAT << area;
    theNode->insertChild("AREA",value.str());

    return theNode;
}

mmpbsa::EMap mmpbsa::EMap::loadXML(const mmpbsa_utils::XMLNode* xmlEnergy)
{

    mmpbsa::EMap returnMe;
    if(xmlEnergy == 0)
        return returnMe;
    std::string data_type,data;
    for(const mmpbsa_utils::XMLNode* it = xmlEnergy->children;it != 0;it = it->siblings)
    {
        data_type = it->getName();
        std::istringstream data(it->getText());
        if(data_type == "BOND")
        {
            data >> MMPBSA_FORMAT >> returnMe.bond;
            continue;
        }

        if(data_type == "ANGLE")
        {
            data >> MMPBSA_FORMAT >> returnMe.angle;
            continue;
        }

        if(data_type == "DIHED")
        {
            data >> MMPBSA_FORMAT >> returnMe.dihed;
            continue;
        }

        if(data_type == "VDW14")
        {
            data >> MMPBSA_FORMAT >> returnMe.vdw14;
            continue;
        }

        if(data_type == "ELE14")
        {
            data >> MMPBSA_FORMAT >> returnMe.ele14;
            continue;
        }

        if(data_type == "VACELE")
        {
            data >> MMPBSA_FORMAT >> returnMe.vacele;
            continue;
        }

        if(data_type == "VDWAALS")
        {
            data >> MMPBSA_FORMAT >> returnMe.vdwaals;
            continue;
        }

        if(data_type == "PBSOL")
        {
            data >> MMPBSA_FORMAT >> returnMe.elstat_solv;
            continue;
        }

        if(data_type == "SASOL")
        {
            data >> MMPBSA_FORMAT >> returnMe.sasol;
            continue;
        }

        if(data_type == "AREA")
        {
            data >> MMPBSA_FORMAT >> returnMe.area;
            continue;
        }

        std::cerr << "Emap::loadXML was given an unknown data type (" <<
                data_type << "), which will be ignored.";
    }

    return returnMe;
}


