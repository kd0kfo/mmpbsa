/* 
 * Molecular dynamics wrapper
 *
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */

#ifndef ENERGYINFO_H
#define	ENERGYINFO_H

#include <cmath>
#include <valarray>
#include <string>
#include <sstream>
#include <fstream>

#include "StringTokenizer.h"
#include "mmpbsa_io.h"
#include "mmpbsa_utils.h"

namespace mmpbsa{

class EnergyInfo : public std::valarray<mmpbsa_t> {
public:
    /**
     * MD Energy Class.
     * Wraps the energy data, providing arithmetic over the set of energy data.
     * 
     */
    EnergyInfo() : std::valarray<mmpbsa_t>(0.0,total_parameters){}
    virtual ~EnergyInfo(){}

    /**
     * Loads first energy info from the provided mdout file.
     * 
     * @param fileName
     * @return 
     */
    void get_first_energyinfo(const char* fileName);

    /**
     * Loads next energy info from the provided mdout file.
     * 
     * @param mdoutFile
     * @return
     */
    void get_next_energyinfo(std::fstream& mdoutFile);

    /**
     * Stores the provided energy value in the correct place in the EnergyInfo class
     * Returns true if the identifier is a correct energy type. Returns false otherwise.
     * 
     * @param identifier
     * @param value
     * @return
     */
    bool loadEnergyValue(const std::string& identifier,const mmpbsa_t& value);
    bool loadEnergyValue(const std::string& identifier,const std::string& value);


    /**
     * Returns true if at least one element is greater than the corresponding
     * element of rhs. Otherwise, false is returned.
     * 
     * @param rhs
     * @return
     */
    bool operator>(const EnergyInfo& rhs)const;

    friend std::ostream& operator<<(std::ostream& theStream, const mmpbsa::EnergyInfo& data);
    friend std::fstream& operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data);

    /**
     * energy information from sander/dynlib.f:prntmd
     *
     */
    enum EnergyInfoIndices {
        nstep, /*step number*/
        time,
        temp,
        press,
        etot,
        ektot,
        eptot,
        bond,
        angle,
        dihed,
        nb14,
        eel14,
        vdwaals,
        eelec,
        ehbond,
        egb,
        restraint,
        esurf,
        nonconst_pot,
        dvdl,
        rms_pbs,
        ekcmt,
        virial,
        volume,
        eksolt,
        eksolv,
        epol,
        e3bod,
        diprms,
        dipitr,
        dipole_temp,
        density,
        ewalderr,
        rmsdvalue,
        tgtrmsd,
        total_parameters/*gives the total number of energy data types in the EnergyInfo class. For use with loops internally.*/
    };

    
};

class AveRmsEnerInfo
{
public:
    AveRmsEnerInfo();
    AveRmsEnerInfo(const AveRmsEnerInfo& orig);
    AveRmsEnerInfo(const EnergyInfo& avgs, const EnergyInfo& rmses);

    virtual ~AveRmsEnerInfo(){}

    /**
     * Loads next energy info from the provided mdout file.
     *
     * @param mdoutFile
     * @return
     */
    void get_next_energyinfo(std::fstream& mdoutFile);
    void clear();


    AveRmsEnerInfo& operator=(const AveRmsEnerInfo& rhs);

private:
    EnergyInfo avg;//averages
    EnergyInfo rms;//Root Mean Squares
    EnergyInfo relrms;//rms/abs(avg)

};

/**
 * Since this file provides the classes to handle Energy Info, the definitions and
 * declarations are provided here as well, though placed in the expected mmpbsa_io
 * namespace.
 *
 * @param mdout
 * @param energyinfos
 * @param avginfos
 */
void mdout2enerinfos(std::fstream& mdout, std::valarray<EnergyInfo>& energyinfos,
        std::valarray<AveRmsEnerInfo>& avginfos);

/**
 * Return final energy from the mdout file of a minization run
 * 
 * @param mdout
 * @return
 */
float get_minimized_energy(std::fstream& mdout) throw (SanderIOException);

};//end namespace mmpbsa

mmpbsa::EnergyInfo sqrt(const mmpbsa::EnergyInfo& rhs);

#endif	/* ENERGYINFO_H */

