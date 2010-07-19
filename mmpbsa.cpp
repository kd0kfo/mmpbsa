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
#include "mmpbsa_io.h"
#include "EnergyInfo.h"
#include "Energy.h"
#include "MeadInterface.h"

#include "MEAD/FinDiffMethod.h"

int realDeal(int argc, char** argv);
void printSnapshot(const EMap& complexEMap, const EMap& receptorEMap,
        const EMap& ligandEMap, std::fstream& outFile);
int testsubmain(int argc, char** argv);
int testenergy(int argc, char** argv);
int testPolyAT(int argc, char** argv);
void parseArgs(int argc, char** argv, MeadInterface mi);
void parseParameter(std::string arg, MeadInterface mi);
void parseFlag(std::string flag, MeadInterface mi);


std::fstream myOutput;
std::fstream prmtopFile;
std::fstream trajFile;
std::fstream radiiFile;
std::fstream outputFile;

/*
 * 
 */
int main(int argc, char** argv)
{
    try 
    {
        myOutput;
        myOutput.open("/ibis/users_linux/dcoss/working_dir/tut1/polyAT-mmpbsa/myoutput.txt",std::ios::out);
        int returnMe = realDeal(argc, argv);
        myOutput.close();
        return returnMe;
    }    
    catch (MMPBSAException e)
    {
        std::cerr << e.identifier() << ": " << e.what() << std::endl;
        return e.getErrType();
    }
}

int realDeal(int argc, char** argv)
{
    using std::valarray;
    using std::slice;
    using std::map;
    
    MeadInterface mi;
    mi.istrength = 0;
    mi.surf_tension =  0.00542;// kcal/mol/Ang^2
    mi.surf_offset = 0.92;// kcal/mol

    parseArgs(argc,argv,mi);

    //load and check the parmtop file.
    mmpbsa_io::SanderParm * sp = new mmpbsa_io::SanderParm;
    sp->raw_read_amber_parm(::prmtopFile);
    if(!sp->sanityCheck())
        throw MMPBSAException("Parmtop file is insane.",INVALID_PRMTOP_DATA);

    //Create energy function with the parmtop data. This energy function will
    //have everything in it. Receptor and ligand (and of course the whole complex)
    //will be stripped out.
    EmpEnerFun entireEFun(sp);

    valarray<bool> complexKeepers(false,sp->natom);
    valarray<bool> receptorKeepers(false,sp->natom);
    valarray<bool> ligandKeepers(false,sp->natom);
    valarray<size_t> receptorStartPos(size_t(0),1);//would be better as a command line argument, but I don't know yet what this will look like. To add later.
    valarray<size_t> ligandStartPos(size_t(1),1);//array containing the starting positions of ligands and receiptors
    size_t bottom,top;
    size_t receptorSize = 0;
    size_t ligandSize = 0;

    for(size_t i = 0;i<receptorStartPos.size();i++)
    {
        size_t currPos = receptorStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currReceptor(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorSize += top-bottom;
    }
    for(size_t i = 0;i<ligandStartPos.size();i++)
    {
        size_t currPos = ligandStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currLigand(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currLigand;
        ligandKeepers[slice(bottom,top-bottom,1)] = currLigand;
        ligandSize += top-bottom;
    }
    size_t complexSize = receptorSize+ligandSize;

    EmpEnerFun complexEFun = entireEFun.stripEnerFun(complexKeepers,true);
    EmpEnerFun receptorEFun = entireEFun.stripEnerFun(receptorKeepers,true);
    EmpEnerFun ligandEFun = entireEFun.stripEnerFun(ligandKeepers,true);

    //load radii data, if available
    map<std::string,mmpbsa_t> radii;//later, check to see if radii.size() > 0 before calling full_EMap(...)
    map<std::string,std::string> residues;
    if(radiiFile.is_open())
        mmpbsa_io::read_siz_file(radiiFile,radii, residues);

    if(!trajFile.good())
        throw MMPBSAException("Unable to read from trajectory file",BROKEN_TRAJECTORY_FILE);

    using namespace mmpbsa_io;
    get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
    int snapcounter = 0;
    valarray<mmpbsa_t> snapshot(sp->natom*3);
    valarray<mmpbsa_t> complexSnap(complexSize*3);
    valarray<mmpbsa_t> receptorSnap(receptorSize*3);
    valarray<mmpbsa_t> ligandSnap(ligandSize*3);
    while(!trajFile.eof())
    {
        try{
            if(get_next_snap(trajFile, snapshot, sp->natom,true))
                printf("Snapshot #%d has been loaded.\n",++snapcounter);
            else
                throw MMPBSAException("Error in loading snapshot",BROKEN_TRAJECTORY_FILE);
        }
        catch(MMPBSAException e)
        {
            if(e.getErrType() == UNEXPECTED_EOF)
            {
              std::cerr << "End of Snapshots Reached" << std::endl;
              return 0;
            }
        }
        //process snapshot.
        size_t complexCoordIndex = 0;
        size_t receptorCoordIndex = 0;
        size_t ligandCoordIndex = 0;
        for(size_t i = 0;i<sp->natom;i++)
        {
            std::slice_array<mmpbsa_t> currCoord = snapshot[slice(3*i,3,1)];
            if(complexKeepers[i])
                complexSnap[slice(3*complexCoordIndex++,3,1)] = currCoord;
            if(receptorKeepers[i])
                receptorSnap[slice(3*receptorCoordIndex++,3,1)] = currCoord;
            if(ligandKeepers[i])
                ligandSnap[slice(3*ligandCoordIndex++,3,1)] = currCoord;
        }

        FinDiffMethod fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
        map<std::string,mmpbsa_t> * pradii = 0;//don't delete!!!
        if(radii.size())//if the radii map is empty, use MeadInterface's lookup table.
            pradii = &radii;

        EMap complexEMap = MeadInterface::full_EMap(complexEFun,complexSnap,fdm,
                pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
        EMap receptorEMap = MeadInterface::full_EMap(receptorEFun,receptorSnap,fdm,
                pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
        EMap ligandEMap = MeadInterface::full_EMap(ligandEFun,ligandSnap,fdm,
                pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);

        printSnapshot(complexEMap,receptorEMap,ligandEMap,::outputFile);
    }


    delete sp;

    ::outputFile.close();
    ::prmtopFile.close();
    ::radiiFile.close();
    ::trajFile.close();

    return 0;
}

void printSnapshot(const EMap& complexEMap, const EMap& receptorEMap, const EMap& ligandEMap, std::fstream& outFile)
{
    outFile << "COMPLEX" << std::endl;
    outFile << complexEMap  << std::endl;
    outFile << "RECEPTOR" << std::endl;
    outFile << receptorEMap  << std::endl;
    outFile << "LIGAND" << std::endl;
    outFile << ligandEMap  << std::endl << std::endl;
}

void parseArgs(int argc, char** argv, MeadInterface mi)
{
    using std::string;
    for(int i = 1;i<argc;i++)
    {
        string currArg = argv[i];
        
        if(currArg.substr(0,2) == "--")
            currArg.erase(currArg.begin(),currArg.begin()+1);

        if(currArg.find("=") != string::npos)
            parseParameter(currArg,mi);
        else
            parseFlag(currArg,mi);
    }
}

void parseParameter(std::string arg, MeadInterface mi)
{
    if(arg.find("=") == string::npos)
        return parseFlag(arg,mi);//oops that's a flag.
    
    if(arg.find("=") != arg.rfind("="))
    {
        char error[256];
        sprintf(error,"Multiple occurance of \"=\" in parameter: %s\n",arg.c_str());
        throw MMPBSAException(error,COMMAND_LINE_ERROR);
    }
    
    using std::string;
    string name = arg.substr(0,arg.find("="));
    string value = arg.substr(arg.find("=")+1);

    if(name == "prmtop" || name == "parmtop")
    {
        ::prmtopFile.open(value.c_str(),std::ios::in);
        return;
    }
    else if(name == "coordinates" || name == "traj" || name == "mdcrd")
    {
        ::trajFile.open(value.c_str(),std::ios::in);
        return;
    }
    else if(name == "out" || name == "output")
    {
        ::outputFile.open(value.c_str(),std::ios::out);
        return;
    }
    else if(name == "radii")
    {
        ::radiiFile.open(value.c_str(),std::ios::in);
        return;
    }
    else if(name == "istrength")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.istrength = mmpbsa_t(fValue);
        return;
    }
    else if(name == "surf_offset")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.surf_offset = mmpbsa_t(fValue);
        return;
    }
    else if(name == "surf_tension")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.surf_tension = mmpbsa_t(fValue);
        return;
    }

    fprintf(stderr,"I don't know what to do with the parameter %s (=%s)\n",name.c_str(),value.c_str());
}

void parseFlag(std::string flag, MeadInterface mi)
{
    if(flag.find("=") != flag.npos)
    {
        parseParameter(flag,mi);//Oops that's a flag.
        return;
    }
    //no flags yet :-/
    fprintf(stderr,"I don't know what to do with the flag \"%s\"\n",flag.c_str());
}

int testPolyAT(int argc, char** argv)
{
    using namespace std;

  try{
    mmpbsa_io::SanderParm *sp = new mmpbsa_io::SanderParm;
    sp->raw_read_amber_parm("/ibis/users_linux/dcoss/working_dir/tut1/polyAT_vac.prmtop");

    if(!sp->sanityCheck()){
      std::cout << "Sander parmtop file is not sane." << std::endl;
      return BROKEN_PRMTOP_FILE;
    }

    EmpEnerFun efun(sp);
    std::valarray<mmpbsa_t> crds(sp->natom * 3);
    std::fstream crdFile("/ibis/users_linux/dcoss/working_dir/tut1/results_12Acut/polyAT_vac_md1_12Acut.mdcrd");
    mmpbsa_io::read_crds(crdFile,crds);

    EMap eMap(&efun,crds);

    cout << eMap << endl;

    delete sp;
  }
  catch(MMPBSAException e)
    {
      std::cout << e.identifier() << " " << e.what() << std::endl;
      return e.getErrType();
    }

  return 0;
}

int testenergy(int argc, char** argv)
{
    using namespace std;
    using namespace mmpbsa_io;
    string dir = "/ibis/users_linux/dcoss/working_dir/mmpbsa-test/";
    string testdir = dir;
    string prmTopFilename = dir + "mm668036_022W.top";
    //string crdFilename = "";
    string trajFilename = dir + "NPTMDprod.mdcrd"; //dir + "polyAT_vac_md1_12Acut.mdcrd";
    string mdoutFilename = dir + "polyAT_vac_md1_12Acut.out";
    string radiiFilename = dir + "my_parse_delphi.siz";
    string resultFilename = dir + "mmpbsa-results.out";
    fstream resultFile(resultFilename.c_str(),ios::out);

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

    int snapNumber = 0;
    while(!trajFile.eof())
    {
        printf("Snapshot #%d\n",++snapNumber);
        std::valarray<mmpbsa_t> snapshot(sp->natom*3);
        if(get_next_snap(trajFile, snapshot, sp->natom,true))
            printf("Snapshot loaded\n");
        else
            throw MMPBSAException("Error in loading snapshot",BROKEN_TRAJECTORY_FILE);

        valarray<bool> keepers(false, sp->natom);
        keepers[slice(0, 1415, 1)] = valarray<bool>(true,1415);
        EmpEnerFun stripped = efun.stripEnerFun(keepers, true);

        valarray<mmpbsa_t> reccrds(1348*3);
        reccrds = (snapshot[slice(0,1348*3,1)]);
        valarray<bool> reckeepers(false,1415);
        for(size_t i = 0;i<1348;i++)
            reckeepers[i] = true;
        EmpEnerFun recstripped = stripped.stripEnerFun(reckeepers,true);

        valarray<bool> ligkeepers(false,1415);
        for(size_t i = 1348;i<1415;i++)
            ligkeepers[i] = true;
        valarray<mmpbsa_t> ligcrds(67*3);
        ligcrds = (snapshot[slice(1348*3,67*3,1)]);
        EmpEnerFun ligstripped = stripped.stripEnerFun(ligkeepers,true);

        EMap complexMap(&stripped,snapshot);
        EMap receptorMap(&recstripped,reccrds);
        EMap ligandMap(&ligstripped,ligcrds);

        fstream radiiFile(radiiFilename.c_str(),ios::in);
        map<string,mmpbsa_t> radii;
        map<string,string> radiiResidues;
        mmpbsa_io::read_siz_file(radiiFile,radii,radiiResidues);

        FinDiffMethod fdm = MeadInterface::createFDM(snapshot[slice(0,1415*3,1)],reccrds,ligcrds);

        //#PB params
        mmpbsa_t interactionStrength = 0.0;

        //#SA params
        mmpbsa_t surf_tension = 0.00542; //# kcal/mol/Ang^2
        mmpbsa_t surf_offset = 0.92; //# kcal/mol

        resultFile << "COMPLEX" << std::endl;
        EMap com_emap = MeadInterface::full_EMap(stripped,snapshot[slice(0,1415*3,1)],
                fdm,&radii,radiiResidues,interactionStrength,surf_tension,surf_offset);
        resultFile << com_emap << std::endl;

        resultFile << "RECEPTOR" << endl;
        EMap rec_emap = MeadInterface::full_EMap(recstripped,reccrds,fdm,&radii,
                radiiResidues,interactionStrength,surf_tension,surf_offset);
        resultFile << rec_emap << std::endl;

        resultFile << "LIGAND" << endl;
        EMap lig_emap = MeadInterface::full_EMap(ligstripped,ligcrds,fdm,&radii,
                radiiResidues,interactionStrength,surf_tension,surf_offset);
        resultFile << lig_emap << std::endl;

        resultFile << std::endl;//separate snapshots with a blank line
    }//end interating through snapshots

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
    using namespace mmpbsa_io;
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



