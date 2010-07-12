/* 
 * File:   EnergyInfo.h
 * Author: dcoss
 *
 * Created on June 21, 2010, 1:49 PM
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


class EnergyInfo {
public:
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
    void setEnergyData(const std::valarray<mmpbsa_t>& newData){energydata = newData;}
    void clear();

    EnergyInfo& operator=(const EnergyInfo& rhs);
    EnergyInfo operator+(const EnergyInfo& rhs)const;
    void operator+=(const EnergyInfo& rhs);
    EnergyInfo operator-(const EnergyInfo& rhs)const;
    void operator-=(const EnergyInfo& rhs);
    EnergyInfo operator/(const EnergyInfo& rhs)const;
    void operator/=(const EnergyInfo& rhs);



    static const int nstep = 0;//step number
    //energy information from sander/dynlib.f:prntmd
    static const int time = 1;
    static const int temp = 2;
    static const int press = 3;
    static const int etot = 4;
    static const int ektot = 5;
    static const int eptot = 6;
    static const int bond = 7;
    static const int angle = 8;
    static const int dihed = 9;
    static const int nb14 = 10;
    static const int eel14 = 11;
    static const int vdwaals = 12;
    static const int eelec = 13;
    static const int ehbond = 14;
    static const int egb = 15;
    static const int restraint = 16;
    static const int esurf = 17;
    static const int nonconst_pot = 18;
    static const int dvdl = 19;
    static const int rms_pbs = 20;
    static const int ekcmt  = 21;
    static const int virial = 22;
    static const int volume = 23;
    static const int eksolt = 24;
    static const int eksolv = 25;
    static const int epol = 26;
    static const int e3bod = 27;
    static const int diprms = 28;
    static const int dipitr = 29;
    static const int dipole_temp = 30;
    static const int density = 31;
    static const int ewalderr = 32;
    static const int rmsdvalue = 33;
    static const int tgtrmsd = 34;

    //gives the total number of energy data types in the EnergyInfo class. For use with loops internally.
    static const int total_parameters = 35;


protected:
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

//Since this file provides the classes to handle Energy Info, the definitions and
//declarations are provided here as well, though placed in the expected mmpbsa_io
//namespace.
void mdout2enerinfos(std::fstream& mdout, std::valarray<EnergyInfo>& energyinfos,
        std::valarray<AveRmsEnerInfo>& avginfos);

/**
 * Return final energy from the mdout file of a minization run
 * 
 * @param mdout
 * @return
 */
float get_minimized_energy(std::fstream& mdout) throw (SanderIOException);

#endif	/* ENERGYINFO_H */

