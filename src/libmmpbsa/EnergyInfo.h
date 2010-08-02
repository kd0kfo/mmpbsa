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

class EnergyInfo {
public:
    /**
     * MD Energy Class.
     * Wraps the energy data, providing arithmetic over the set of energy data.
     * 
     */
    EnergyInfo();
    EnergyInfo(const EnergyInfo& orig);
    virtual ~EnergyInfo();

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

    const std::valarray<mmpbsa_t>& getEnergyData(){return energydata;}

    /**
     * Copies the provided data array to the energy data object.
     * Note: it is the user's resposibility to ensure the data is in the correct
     * place in the array.
     * 
     * @param newData
     */
    void setEnergyData(const std::valarray<mmpbsa_t>& newData);

    /**
     * Clears the energy data
     *
     */
    void clear();

    EnergyInfo& operator=(const EnergyInfo& rhs);
    EnergyInfo operator+(const EnergyInfo& rhs)const;
    void operator+=(const EnergyInfo& rhs);
    EnergyInfo operator-(const EnergyInfo& rhs)const;
    void operator-=(const EnergyInfo& rhs);
    EnergyInfo operator/(const EnergyInfo& rhs)const;
    void operator/=(const EnergyInfo& rhs);



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


private:
    std::valarray<mmpbsa_t> energydata;//keys correspond to above listed variables.
};

class AveRmsEnerInfo
{
public:
    AveRmsEnerInfo();
    AveRmsEnerInfo(const AveRmsEnerInfo& orig);
    AveRmsEnerInfo(const EnergyInfo& avgs, const EnergyInfo& rmses);

    virtual ~AveRmsEnerInfo(){}

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

#endif	/* ENERGYINFO_H */

