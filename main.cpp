/* 
 * File:   main.cpp
 * Author: dcoss
 *
 * Created on June 16, 2010, 10:55 AM
 */

#include <cstdlib>
#include <iostream>
#include <valarray>
#include <fstream>

#include "mmpbsa_exceptions.h"
#include "SanderIO.h"
#include "EnergyInfo.h"
#include "Energy.h"

int testsubmain(int argc, char** argv);
int testenergy(int argc, char** argv);

std::fstream myOutput;

/*
 * 
 */
int main(int argc, char** argv)
{
    try 
    {
        myOutput;
        myOutput.open("/ibis/users_linux/dcoss/working_dir/mmpbsa-test/myoutput.txt",std::ios::out);
        int returnMe = testenergy(argc, argv);
        myOutput.close();
    }    
    catch (MMPBSAException e)
    {
        std::cerr << e.identifier() << ": " << e.what() << std::endl;
        return e.getErrType();
    }
}

int testenergy(int argc, char** argv)
{
    using namespace std;
    using namespace sanderio;
    string dir = "/ibis/users_linux/dcoss/working_dir/mmpbsa-test/";
    string testdir = dir;
    string prmTopFilename = dir + "mm668036_022W.top";
    //string crdFilename = "";
    string trajFilename = dir + "NPTMDprod.mdcrd"; //dir + "polyAT_vac_md1_12Acut.mdcrd";
    string mdoutFilename = dir + "polyAT_vac_md1_12Acut.out";

    SanderParm * sp = new SanderParm;
    sp->raw_read_amber_parm(prmTopFilename);
    
    if (sp->sanityCheck())
        cout << prmTopFilename << " loaded correctly and is a valid prmtop file" << std::endl;
    else
        throw MMPBSAException("Error in loading ParmInfo",INVALID_PRMTOP_DATA);

    //EmpEnerFun * efun = new EmpEnerFun(sp);
    EmpEnerFun efun(sp);
    
    fstream trajFile(trajFilename.c_str(),ios::in);
    string trajTitle = get_traj_title(trajFile);
    printf("Opened trajectory file <%s> with title \"%s\"\n",trajFilename.c_str(),trajTitle.c_str());

    std::valarray<mmpbsa_t> snapshot(sp->natom*3);
    if(get_next_snap(trajFile, snapshot, sp->natom,true))
        cout << prmTopFilename << " loaded first snapshot" << endl;
    else
        throw MMPBSAException("Error in loading snapshot",BROKEN_TRAJECTORY_FILE);

    valarray<bool> keepers(false, sp->natom);
    keepers[slice(0, 1415, 1)] = valarray<bool>(true,1415);
    EmpEnerFun stripped = efun.stripEnerFun(keepers, true);
    cout << stripped.ereport(snapshot) << endl;

    valarray<mmpbsa_t> reccrds(1348*3);
    reccrds = (snapshot[slice(0,1348*3,1)]);
    valarray<bool> reckeepers(false,1415);
    for(size_t i = 0;i<1348;i++)
        reckeepers[i] = true;
    EmpEnerFun recstripped = stripped.stripEnerFun(reckeepers,true);
    cout << recstripped.ereport(reccrds) << endl;

    valarray<bool> ligkeepers(false,1415);
    for(size_t i = 1348;i<1415;i++)
        ligkeepers[i] = true;
    valarray<mmpbsa_t> ligcrds(67*3);
    ligcrds = (snapshot[slice(1348*3,67*3,1)]);
    EmpEnerFun ligstripped = stripped.stripEnerFun(ligkeepers,true);
    cout << ligstripped.ereport(ligcrds) << endl;

    cout.flush();
    cerr.flush();
    trajFile.close();
    //delete efun;
    delete sp;
    return 0;
}

int testsubmain(int argc, char** argv)
{
    using namespace std;
    using namespace sanderio;
    string dir = "/ibis/users_linux/dcoss/working_dir/mmpbsa-test/";
    string testdir = dir;
    string prmTopFilename = dir + "mm668036_022W.top";
    //string crdFilename = "";
    string trajFilename = dir + "NPTMDprod.mdcrd"; //dir + "polyAT_vac_md1_12Acut.mdcrd";
    string mdoutFilename = dir + "polyAT_vac_md1_12Acut.out";

    SanderParm * sp = new SanderParm;
    sp->raw_read_amber_parm(prmTopFilename);
    cout << "Loaded Parameters" << endl;

    if (sp->sanityCheck())
        cout << prmTopFilename << " loaded correctly and is a valid prmtop file" << std::endl;
    else
        cout << prmTopFilename << " loaded but is not a valid prmtop file" << std::endl;

    //        valarray<double> crds;
    //        fstream crdFile(crdFilename.c_str(), ios::in);
    //        string title = read_crds(crdFile, crds);
    //        printf("Read %d coordinates.\n",int(crds.size()/3));
    //
    //        printf("First coordinate is (%f,%f,%f)\n", crds[0], crds[1], crds[2]);std::cout.flush();
    //        write_crds(testdir.append("testWriteCrds.inpcrd").c_str(), crds, title.c_str());
    //
    std::fstream trajFile(trajFilename.c_str(), std::ios::in);
    std::valarray<mmpbsa_t> trajArray(0.0, sp->natom * 3);
    string trajTitle = get_traj_title(trajFile);
    bool isPeriodic = true;
    std::cout << "Title (" << trajTitle << ")" << std::endl;
    if (get_next_snap(trajFile, trajArray, sp->natom, isPeriodic))
        printf("Loaded 1st snapshot\n");
    else
        printf("Failed to loadsnapshot\n");
    std::cout.flush();

    if (get_next_snap(trajFile, trajArray, sp->natom, isPeriodic))
        printf("Loaded 2nd snapshot\n");
    else
        printf("Failed to loadsnapshot\n");

    cout.flush();

    //        EnergyInfo ei;
    //        ei.get_first_energyinfo(mdoutFilename.c_str());
    //
    //        printf("Loaded energyinfo\n");
    //        cout.flush();

    EmpEnerFun * efun = new EmpEnerFun(sp);
    

    delete sp;
    delete efun;
    cout << "Everything went well." << endl;

    return 0;
}



