/* 
 * File:   main.cpp
 * Author: dcoss
 *
 * Created on June 16, 2010, 10:55 AM
 */

#include <cstdlib>
#include <iostream>
#include "mmpbsa_exceptions.h"
#include "SanderIO.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

  try
    {
      SanderParm sp;
      sp.raw_read_amber_parm("polyAT_cio.prmtop");
      std::cout << "Done" << std::endl;
    }
  catch(SanderIOException sioe)
    {
      std::cerr << sioe.what() <<std::endl;
      return sioe.getErrType();
    }

  return 0;
}



