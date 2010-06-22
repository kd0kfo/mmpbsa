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

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

  try
    {
      using std::string;
      using namespace sanderio;
      
      string dir = "/net/sbarray1/share/ibis/users_linux/dcoss/working_dir/tut1/";
      string prmTopFilename = dir + "polyAT_wat.prmtop";
      string crdFilename = dir + "polyAT_wat.inpcrd";
      string trajFilename = dir + "polyAT_vac_md1_12Acut.mdcrd";
      string mdoutFilename = dir + "polyAT_vac_md1_12Acut.out";

      SanderParm sp;
      sp.raw_read_amber_parm(prmTopFilename);
      std::cout << "Loaded Parameters" << std::endl;
      bool test = sp.sanityCheck();
      
      if(test)
          std::cout << prmTopFilename << " loaded correctly and is a valid prmtop file" << std::cout;
      else
          std::cout << prmTopFilename  << " loaded but is not a valid prmtop file" << std::cout;

      std::valarray<double> crds;
      std::fstream crdFile(crdFilename.c_str(),std::ios::in);
      std::string title = read_crds(crdFile, crds);

      printf("First coordinate is (%f,%f,%f)",crds[0],crds[1],crds[2]);
      write_crds("/ibis/users_linux/dcoss/tmp/testWriteCrds.inpcrd",crds,title.c_str());

      std::fstream trajFile(trajFilename.c_str(),std::ios::in);
      std::valarray<double> trajArray(0.0,638);
      string trajTitle = get_traj_title(trajFile);
      std::cout << "Title (" << trajTitle << ")" << std::endl;
      if(get_next_snap(trajFile, trajArray,638))
          printf("Loaded 1st snapshot\n");
      else
          printf("Failed to loadsnapshot\n");

      if(get_next_snap(trajFile, trajArray,638))
          printf("Loaded 2nd snapshot\n");
      else
          printf("Failed to loadsnapshot\n");

      cout.flush();

      EnergyInfo ei;
      ei.get_first_energyinfo(mdoutFilename.c_str());

      printf("Loaded energyinfo\n");
      cout.flush();

    }
  catch(SanderIOException sioe)
    {
      std::cerr << "SanderIO Error: " << sioe.what() <<std::endl;
      return sioe.getErrType();
    }

  return 0;
}



