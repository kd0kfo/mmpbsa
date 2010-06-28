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

using namespace std;

/*
 * 
 */
int main(int argc, char** argv)
{

    try {

        using namespace std;
        using namespace sanderio;


        string dir = "/net/sbarray1/share/ibis/users_linux/dcoss/working_dir/tut1/";
        string testdir = "/ibis/users_linux/dcoss/working_dir/test_mmpbsa/";
        string prmTopFilename = dir + "polyAT_wat.prmtop";
        string crdFilename = dir + "polyAT_wat.inpcrd";
        string trajFilename = dir + "polyAT_vac_md1_12Acut.mdcrd";
        string mdoutFilename = dir + "polyAT_vac_md1_12Acut.out";

        SanderParm* sp = new SanderParm();
        sp->raw_read_amber_parm(prmTopFilename);
        cout << "Loaded Parameters" << endl;
        
        if (sp->sanityCheck())
            cout << prmTopFilename << " loaded correctly and is a valid prmtop file" << std::cout;
        else
            cout << prmTopFilename << " loaded but is not a valid prmtop file" << std::cout;

        valarray<double> crds;
        fstream crdFile(crdFilename.c_str(), ios::in);
        string title = read_crds(crdFile, crds);

        printf("First coordinate is (%f,%f,%f)\n", crds[0], crds[1], crds[2]);std::cout.flush();
        write_crds(testdir.append("testWriteCrds.inpcrd").c_str(), crds, title.c_str());

        std::fstream trajFile(trajFilename.c_str(), std::ios::in);
        std::valarray<double> trajArray(0.0, 638);
        string trajTitle = get_traj_title(trajFile);
        std::cout << "Title (" << trajTitle << ")" << std::endl;
        if (get_next_snap(trajFile, trajArray, 638))
            printf("Loaded 1st snapshot\n");
        else
            printf("Failed to loadsnapshot\n");

        if (get_next_snap(trajFile, trajArray, 638))
            printf("Loaded 2nd snapshot\n");
        else
            printf("Failed to loadsnapshot\n");

        cout.flush();

        EnergyInfo ei;
        ei.get_first_energyinfo(mdoutFilename.c_str());

        printf("Loaded energyinfo\n");
        cout.flush();

        EmpEnerFun * efun = new EmpEnerFun(sp);

        delete sp;
        cout << "Everything went well." << endl;
    }
  catch(MMPBSAException e)
    {
      std::cerr << e.identifier() << ": " << e.what() << std::endl;
      return e.getErrType();
    }
    
  return 0;
}



