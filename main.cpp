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

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

  try
    {
      using std::string;
      string dir = "/net/sbarray1/share/ibis/users_linux/dcoss/working_dir/tut1/";
      string prmTopFilename = dir + "polyAT_wat.prmtop";
      string crdFilename = dir + "polyAT_wat.inpcrd";

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

      printf("First coordinate is (%e,%e,%e)",crds[0],crds[1],crds[2]);
      write_crds("/ibis/users_linux/dcoss/tmp/testWriteCrds.inpcrd",crds,title.c_str());
    }
  catch(SanderIOException sioe)
    {
      std::cerr << "SanderIO Error: " << sioe.what() <<std::endl;
      return sioe.getErrType();
    }

  return 0;
}



