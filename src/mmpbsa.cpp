#include "mmpbsa.h"
#include <iomanip>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{
  using mmpbsa::MMPBSAState;

#ifdef USE_PTHREADS
  pthread_mutex_init(&mmpbsa_mutex,NULL);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
#endif

#ifdef USE_MPI
  std::string mmpbsa_output_filename;
  mmpbsa_utils::mpi_init_hosts(&argc,&argv,mpi_rank,mpi_size);
  if(mpi_rank == 0)
    data_list = new mmpbsa_utils::XMLNode(MMPBSA_XML_TITLE);
  else
    data_list = 0;
#endif

  std::cerr << PACKAGE_STRING;

#ifdef USE_MPI
  if(mpi_rank == 0)
    mpi_processes_running = mpi_size;
  std::cerr << " with MPI on host number " << mpi_rank << " (" << getpid() << ")";
#endif

  std::cerr <<" started on " << mmpbsa_utils::get_human_time() << std::endl;

  try
    {
      ::timeAtPreviousCheckpoint = 0;
      mmpbsa_boinc_init();//must be called before any other BOINC routines. If BOINC is not used, nothing will happen.
      int retval = 0;

      //Try to create a queue based on a provided queue file.
      getQueueFile(processQueue,argc,argv);

      //If no queue was found, run based on arguments in argv. These should
      //correspond to what is placed in the queue XML file.
      if(processQueue.size() == 0)
	::processQueue = parseArgs(argc,argv);
      if(processQueue.size() == 0)
        {
	  std::cout << helpString() << std::endl;
        }
      else
        {
	  int queuePosition= 0;
	  for(std::vector<MMPBSAState>::iterator job = processQueue.begin();
	      job != processQueue.end();job++,queuePosition++)
            {
	      switch(job->currentProcess)
                {
		case MMPBSAState::SANDER:
		  restart_sander(*job,job->currentSI);
		  if(job->placeInQueue > queuePosition)
		    {
		      queuePosition++;
		      continue;
		    }
		  else if(job->placeInQueue < queuePosition)
		    {
		      job->currentSI.completed = false;
		      job->fractionDone = 0;
		      job->placeInQueue = queuePosition;
		    }
		  retval = sander_run(*job,job->currentSI);
		  break;
		case MMPBSAState::MOLSURF:
		  if(job->placeInQueue > queuePosition)
		    {
		      queuePosition++;
		      continue;
		    }
		  else if(job->placeInQueue < queuePosition)
		    {
		      job->placeInQueue = queuePosition;
		    }
		  molsurf_run(*job);
		  break;
		case MMPBSAState::MMPBSA:
		  restart_mmpbsa(*job);
		  if(job->placeInQueue > queuePosition)
		    {
		      queuePosition++;
		      continue;
		    }
		  else if(job->placeInQueue < queuePosition)
		    {
		      job->fractionDone = 0;
		      job->placeInQueue = queuePosition;
		    }
		  retval = mmpbsa_run(*job,job->currentMI);
#ifdef USE_MPI
		  if(has_filename(MMPBSA_OUT_TYPE,*job))
		    mmpbsa_output_filename = get_filename(MMPBSA_OUT_TYPE,*job);
		  else if(has_filename(SANDER_MDOUT_TYPE,*job))
		    mmpbsa_output_filename = get_filename(SANDER_MDOUT_TYPE,*job);
#endif
		  break;
                }
	      if(retval)
		break;
            }
        }
#ifdef USE_BOINC
      if(retval)
	std::cerr << "BOINC Error: " << boincerror(retval) << std::endl;
      boinc_finish(retval);
#endif

#ifdef USE_MPI
      mmpbsa_utils::mpi_dump_data(data_list,mmpbsa_output_filename);
      delete data_list;
      MPI_Finalize();
      std::cerr << mpi_rank << " finalized." << std::endl;
#endif
      std::cerr << PACKAGE_STRING <<" finished (" << retval << ") on " << mmpbsa_utils::get_human_time() << "\n" << std::endl;
      return retval;
    }    
  catch (mmpbsa::MMPBSAException e)
    {
      std::cerr << e.identifier() << ": " << e.what() << std::endl;
#ifdef USE_BOINC
      boinc_finish(e.getErrType());
#endif
#ifdef USE_MPI
      MPI_Finalize();
      std::cerr << mpi_rank << " finalized." << std::endl;
#endif
      std::cerr << PACKAGE_STRING <<" finished (" << e.getErrType() << ") on " << mmpbsa_utils::get_human_time() << "\n" << std::endl;
      return e.getErrType();

    }
}

void writePDB(const std::vector<mmpbsa::atom_t>& atoms,const mmpbsa::forcefield_t& ff,
	      const std::valarray<mmpbsa::Vector>crds,const mmpbsa::MMPBSAState& currState,const std::string& molecule)
{
  std::fstream pdbFile;
  std::string filename = "default";
  if(has_filename(MMPBSA_OUT_TYPE,currState))
    filename = get_filename(MMPBSA_OUT_TYPE,currState);
  else if(has_filename(SANDER_MDOUT_TYPE,currState))
    filename = get_filename(SANDER_MDOUT_TYPE,currState);
  size_t dotLoc = filename.find_last_of('.');
  if(dotLoc != std::string::npos)
    filename.erase(dotLoc);
  filename += "_" + molecule;
  std::ostringstream buff;
  buff << currState.currentSnap;
  filename += + "_snap-" + buff.str() + ".pdb";

  pdbFile.open(filename.c_str(),std::ios::out);
  if(!pdbFile.is_open())
    throw mmpbsa::MMPBSAException("writePDB could not open for writing: " + filename,mmpbsa::FILE_IO_ERROR);

  streamPDB(pdbFile,atoms,ff,crds);
  pdbFile.close();
}

#include <unistd.h>
void study_cpu_time()
{
#if 0//def USE_BOINC
  double debug_cpu_time;
  pid_t the_pid = getpid();
  boinc_wu_cpu_time(debug_cpu_time);
  std::cout << "MMPBSA time: " << linux_cpu_time(the_pid) << " BOINC time: " << debug_cpu_time << std::endl;
#endif
}

int write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState)
{
#ifdef USE_MPI
  return mmpbsa_utils::mpi_write_mmpbsa_data(energy_data,currState,mpi_rank,data_list);
#else

  string filename;
  if(has_filename(MMPBSA_OUT_TYPE,currState))
    filename = get_filename(MMPBSA_OUT_TYPE,currState);
  else if(has_filename(SANDER_MDOUT_TYPE,currState))
    filename = get_filename(SANDER_MDOUT_TYPE,currState);
  else
    filename = "mmpbsa-output.xml";

  std::string data = energy_data.toString();
  std::ios::openmode the_mode = std::ios::out;
  if(filename.find(".gz") != std::string::npos || filename.find(".tar") != std::string::npos || filename.find(".tgz") != std::string::npos)
    the_mode |= std::ios::binary;
  std::fstream out_file(filename.c_str(),the_mode);
  if(!out_file.good())
    throw mmpbsa::MMPBSAException("write_mmpbsa_data: unable to open " + filename + " for writing.",mmpbsa::FILE_IO_ERROR);
  mmpbsa_io::smart_write(out_file,data.c_str(),data.size(),&filename);
  out_file.close();
  return 0;

#endif
}

mmpbsa_utils::XMLNode* read_mmpbsa_data(const mmpbsa::MMPBSAState& currState)
{
  using namespace mmpbsa;
  using std::string;
  
  string filename;
  std::stringstream data;
  std::ios::openmode the_mode = std::ios::in;
  
  if(has_filename(MMPBSA_OUT_TYPE,currState))
    filename = get_filename(MMPBSA_OUT_TYPE,currState);
  else if(has_filename(SANDER_MDOUT_TYPE,currState))
    filename = get_filename(SANDER_MDOUT_TYPE,currState);
  else
    throw MMPBSAException("read_mmpbsa_data: Not provided a data filename",COMMAND_LINE_ERROR);
	
  if(filename.find(".gz") != string::npos || filename.find(".tar") != string::npos || filename.find(".tgz") != string::npos)
    the_mode |= std::ios::binary;
  std::fstream in_file(filename.c_str(),the_mode);
  if(!in_file.good())
    throw MMPBSAException("read_mmpbsa_data: unable to open " + filename + " for writing.",FILE_IO_ERROR);
  mmpbsa_io::smart_read(data,in_file,&filename);
  in_file.close();
  return mmpbsa_utils::XMLParser::parse(data);
}

void get_sander_forcefield(mmpbsa::MMPBSAState& currState,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list,mmpbsa_io::trajectory_t& trajfile)
{
  using namespace mmpbsa;
  using std::valarray;
  using std::slice;

  //load and check the parmtop file.
  mmpbsa::SanderParm * sp = new mmpbsa::SanderParm;
  if(!has_filename(MMPBSA_TOPOLOGY_TYPE,currState))
    throw mmpbsa::MMPBSAException("get_sander_forcefield: no parmtop file.",BROKEN_PRMTOP_FILE);
  sp->raw_read_amber_parm(get_filename(MMPBSA_TOPOLOGY_TYPE,currState));
  if(!currState.trustPrmtop)
    if(!sp->sanityCheck())
      throw MMPBSAException("get_sander_forcefield: Parmtop file, " + get_filename(SANDER_PRMTOP_TYPE,currState)
			    + " is insane.",INVALID_PRMTOP_DATA);

  //Create energy function with the parmtop data. This energy function will
  //have everything in it. Receptor and ligand will be stripped out.
  EmpEnerFun entireEFun(sp);

  *split_ff = new mmpbsa::forcefield_t[MMPBSAState::END_OF_MOLECULES];
  *atom_lists = new std::vector<atom_t>[MMPBSAState::END_OF_MOLECULES];

  valarray<bool> complexKeepers(false,sp->natom);//array of atoms to keep.
  valarray<bool> receptorKeepers(false,sp->natom);
  valarray<bool> ligandKeepers(false,sp->natom);
  size_t bottom,top;
  size_t receptorSize = 0;
  size_t ligandSize = 0;

  //Prepare a list of the beginnings and ends of receptors and ligands
  for(std::set<size_t>::const_iterator currPos = currState.receptorStartPos.begin();
      currPos != currState.receptorStartPos.end();currPos++)
    {
      bottom = entireEFun.mol_ranges[2* *currPos];
      top = entireEFun.mol_ranges[2* *currPos+1];
      valarray<bool> currReceptor(true,top-bottom);
      complexKeepers[slice(bottom,top-bottom,1)] = currReceptor;
      receptorKeepers[slice(bottom,top-bottom,1)] = currReceptor;
      receptorSize += top-bottom;
    }
  for(std::set<size_t>::const_iterator currPos = currState.ligandStartPos.begin();
      currPos != currState.ligandStartPos.end();currPos++)
    {
      bottom = entireEFun.mol_ranges[2* *currPos];
      top = entireEFun.mol_ranges[2* *currPos+1];
      valarray<bool> currLigand(true,top-bottom);
      complexKeepers[slice(bottom,top-bottom,1)] = currLigand;
      ligandKeepers[slice(bottom,top-bottom,1)] = currLigand;
      ligandSize += top-bottom;
    }
  size_t complexSize = receptorSize+ligandSize;

  //Separate the molecules
  EmpEnerFun complexEFun = entireEFun.stripEnerFun(complexKeepers,true);
  complexEFun.extract_force_field((*split_ff)[MMPBSAState::COMPLEX]);
  complexEFun.extract_atom_structs((*atom_lists)[MMPBSAState::COMPLEX]);
  EmpEnerFun receptorEFun = entireEFun.stripEnerFun(receptorKeepers,true);
  receptorEFun.extract_force_field((*split_ff)[MMPBSAState::RECEPTOR]);
  receptorEFun.extract_atom_structs((*atom_lists)[MMPBSAState::RECEPTOR]);
  EmpEnerFun ligandEFun = entireEFun.stripEnerFun(ligandKeepers,true);
  ligandEFun.extract_force_field((*split_ff)[MMPBSAState::LIGAND]);
  ligandEFun.extract_atom_structs((*atom_lists)[MMPBSAState::LIGAND]);

  mol_list.resize(sp->natom,MMPBSAState::END_OF_MOLECULES);//In this case, it is solvent
  mol_list[receptorKeepers] = MMPBSAState::RECEPTOR;
  mol_list[ligandKeepers] = MMPBSAState::LIGAND;

  trajfile.natoms = mol_list.size();
  trajfile.ifbox = sp->ifbox;

  // These Energy Functions loose scope and their SanderParameters should be deleted.
  delete complexEFun.parminfo;complexEFun.parminfo = 0;
  delete receptorEFun.parminfo;receptorEFun.parminfo = 0;
  delete ligandEFun.parminfo;ligandEFun.parminfo = 0;
  delete sp;
}

bool should_calculate_snapshot(const size_t& currentSnap, const std::vector<size_t>& snapList)
{
#ifndef USE_MPI
  return true;
#else
  if(snapList.size() == 0)
    return (currentSnap-1) % mpi_size == mpi_rank;

  for(size_t i = 0;i<snapList.size();i++)
    {
      //Is the requested snap shot in the list?
      if(snapList.at(i) == currentSnap)
	if(i % mpi_size == mpi_rank)//If so should this host do it?
	  return true;
	else
	  return false;
    }

  return false;//If the snapshot is not in the list, don't do it.

#endif

}


void get_forcefield(mmpbsa::MMPBSAState& currState,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list,mmpbsa_io::trajectory_t& trajfile)
{
#ifdef USE_GROMACS
  if(!has_filename(MMPBSA_TOPOLOGY_TYPE,currState))
    throw mmpbsa::MMPBSAException("get_forcefield: no parmtop file.",mmpbsa::BROKEN_PRMTOP_FILE);
  std::string filename = get_filename(MMPBSA_TOPOLOGY_TYPE,currState);
  if(filename.find(".tpr") != std::string::npos)
    {
      std::set<size_t> *receptor_start,*ligand_start;
      receptor_start = (currState.receptorStartPos.size()) ? &currState.receptorStartPos : 0;
      ligand_start = (currState.ligandStartPos.size()) ? &currState.ligandStartPos : 0;
      mmpbsa_io::get_gromacs_forcefield(filename.c_str(),split_ff,atom_lists,mol_list,receptor_start,ligand_start);
      return;
    }
#endif
  get_sander_forcefield(currState,split_ff,atom_lists,mol_list,trajfile);
}

void dump_crds(const std::vector<mmpbsa::atom_t>& atoms,
	       const std::valarray<mmpbsa::Vector>& crds,
	       std::map<std::string,float> radii)
{
  size_t size = crds.size();
  std::cout << "# Start Coordinates" << std::endl;
  for(size_t i = 0;i<size;i++)
    std::cout << crds[i].x() << " " << crds[i].y() << " " << crds[i].z() << " " << mmpbsa_utils::lookup_radius(atoms[i].name,radii) << std::endl;
  std::cout << "# End Coordinates" << std::endl;
}


int molsurf_run(mmpbsa::MMPBSAState& currState)
{
  using std::valarray;
  using std::map;
  using namespace mmpbsa;
  mmpbsa_t area_value,msms_value;
  std::vector<size_t>::const_iterator curr_snap;
  //Setup Trajectory Structure
  if(!has_filename(MMPBSA_TRAJECTORY_TYPE,currState))
    throw mmpbsa::MMPBSAException("molsurf_run: no trajectory file was given.",BROKEN_TRAJECTORY_FILE);
  mmpbsa_io::trajectory_t trajFile = mmpbsa_io::open_trajectory(get_filename(MMPBSA_TRAJECTORY_TYPE,currState),currState.keep_traj_in_mem);



  //Get forcefield data
  mmpbsa::forcefield_t* split_ff = 0;
  std::vector<atom_t>* atom_lists = 0;
  std::valarray<MMPBSAState::MOLECULE> mol_list;
  get_forcefield(currState,&split_ff,&atom_lists,mol_list,trajFile);


  //load radii data, if available
  map<std::string,float> radii;//later, check to see if radii.size() > 0 before calling full_EMap(...)
  map<std::string,std::string> residues;
  if(has_filename(RADII_TYPE,currState))
    {
      std::string radiiFilename = get_filename(RADII_TYPE,currState);
      std::fstream radiiFile(radiiFilename.c_str(),std::ios::in);
      std::stringstream radiiData;
      mmpbsa_io::smart_read(radiiData,radiiFile,&radiiFilename);
      if(currState.verbose)
	std::cout << "Reading radii from " << radiiFilename  << std::endl;
      radiiFile.close();
      mmpbsa_io::read_siz_file(radiiData,radii, residues);
    }
  else
    {
      mmpbsa::MeadInterface default_radii;
      radii = default_radii.brad;
    }

  //setup trajectory storage
  get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
  if(currState.snapList.size() == 0)
    {
      std::fstream trajstream(trajFile.sander_filename->c_str(),std::ios::in);
      size_t num_snaps = mmpbsa_io::count_snapshots(trajstream,trajFile.natoms,trajFile.ifbox > 0);
      trajstream.close();
      for(size_t filler = 0;filler < num_snaps;filler++)
	currState.snapList.push_back(filler + 1);
    }

  curr_snap = currState.snapList.begin();
  for(;curr_snap != currState.snapList.end();curr_snap++)
    {
      valarray<mmpbsa::Vector> snapshot(mol_list.size());
      mmpbsa_io::seek(trajFile,*curr_snap);
      mmpbsa_io::get_next_snap(trajFile,snapshot);
      size_t complexSize,receptorSize,ligandSize;
      receptorSize = ligandSize = 0;
      for(size_t i = 0;i<mol_list.size();i++)
	{
	  if(mol_list[i] == MMPBSAState::RECEPTOR)
	    receptorSize++;
	  else if(mol_list[i] == MMPBSAState::LIGAND)
	    ligandSize++;
	}
      complexSize = receptorSize + ligandSize;
      valarray<mmpbsa::Vector> complexSnap(complexSize);
      valarray<mmpbsa::Vector> receptorSnap(receptorSize);
      valarray<mmpbsa::Vector> ligandSnap(ligandSize);

      //separate coordinates
      size_t complexCoordIndex = 0;
      size_t receptorCoordIndex = 0;
      size_t ligandCoordIndex = 0;
      for(size_t i = 0;i<mol_list.size();i++)
	{
	  Vector& currCoord = snapshot[i];
	  if(mol_list[i] == MMPBSAState::RECEPTOR)
	    {
	      complexSnap[complexCoordIndex++] = currCoord;
	      receptorSnap[receptorCoordIndex++] = currCoord;
	    }
	  else if(mol_list[i] == MMPBSAState::LIGAND)
	    {
	      complexSnap[complexCoordIndex++] = currCoord;
	      ligandSnap[ligandCoordIndex++] = currCoord;
	    }
	}

      std::vector<atom_t>* atoms;
      valarray<mmpbsa::Vector> *snap;
      switch(currState.currentMolecule)
	{
	case MMPBSAState::RECEPTOR:
	  atoms = &atom_lists[MMPBSAState::RECEPTOR];
	  snap = &receptorSnap;
	  break;
	case MMPBSAState::COMPLEX:
	  atoms = &atom_lists[MMPBSAState::COMPLEX];
	  snap = &complexSnap;
	  break;
	case MMPBSAState::LIGAND:
	  atoms = &atom_lists[MMPBSAState::LIGAND];
	  snap = &ligandSnap;
	  break;
	default:
	  throw mmpbsa::MMPBSAException("molsurf_run: invalid molecule type");
	}
      std::cout << MeadInterface::molsurf_area(*atoms,*snap,radii);//molsurf
      //		std::cout << MeadInterface::msms_area(*atoms,*snap,radii);
      std::cout << std::endl;
    }//end of iteration through snap list
  return 0;
}

int mmpbsa_run(mmpbsa::MMPBSAState& currState, mmpbsa::MeadInterface& mi)
{
  using std::valarray;
  using std::vector;
  using std::slice;
  using std::map;
  using namespace mmpbsa;
  using mmpbsa::MMPBSAState;
  using namespace mmpbsa_io;

  std::cout << "Starting MMPBSA calculation ";
#ifdef USE_MPI
  if(mpi_rank == 0)
    mpi_processes_running = mpi_size;
  std::cout << " with MPI on host number " << mpi_rank << " (" << getpid() << ")";
#endif

	std::cout << std::endl;

	if(!has_filename(MMPBSA_OUT_TYPE,currState) && !has_filename(SANDER_MDOUT_TYPE,currState))
	  currState.filename_map[SANDER_MDOUT_TYPE] = "mmpbsa-output.xml";

    //Upon restart, MMPBSA needs to reload energy that was calculated previously and then
    //append new data to it.
    mmpbsa_utils::XMLParser previousEnergyData(new mmpbsa_utils::XMLNode("mmpbsa_energy"));
    if(!currState.overwrite)
      {
	try{
	  mmpbsa_utils::XMLNode* old_data = read_mmpbsa_data(currState);
	  if(old_data != NULL)
	    previousEnergyData.setHead(old_data);
	}
	catch(mmpbsa::MMPBSAException xmlpe)
	  {
	    if(xmlpe.getErrType() != mmpbsa::FILE_IO_ERROR)
	      {
		std::cerr << "Previous energy data is corrupt. Overwriting." << std::endl;
	      }
	    previousEnergyData.setHead(new mmpbsa_utils::XMLNode("mmpbsa_energy"));
	  }
      }

  //Setup Trajectory Structure
  if(!has_filename(MMPBSA_TRAJECTORY_TYPE,currState))
    throw mmpbsa::MMPBSAException("mmpbsa_run: no trajectory file was given.",BROKEN_TRAJECTORY_FILE);
  mmpbsa_io::trajectory_t trajFile = mmpbsa_io::open_trajectory(get_filename(MMPBSA_TRAJECTORY_TYPE,currState),currState.keep_traj_in_mem);



  //Get forcefield data
  mmpbsa::forcefield_t* split_ff = 0;
  std::vector<atom_t>* atom_lists = 0;
  std::valarray<MMPBSAState::MOLECULE> mol_list;
  get_forcefield(currState,&split_ff,&atom_lists,mol_list,trajFile);


  //load radii data, if available
  map<std::string,mead_data_t> radii = mi.brad;// if there is no radii file, use MeadInterface's default.
  map<std::string,std::string> residues;
  if(has_filename(RADII_TYPE,currState))
    {
      std::string radiiFilename = get_filename(RADII_TYPE,currState);
      std::fstream radiiFile(radiiFilename.c_str(),std::ios::in);
      std::stringstream radiiData;
      mmpbsa_io::smart_read(radiiData,radiiFile,&radiiFilename);
      std::cout << "Reading radii from " << radiiFilename  << std::endl;
      radiiFile.close();
      mmpbsa_io::read_siz_file(radiiData,radii, residues);
    }

  //setup trajectory storage
  get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
  valarray<mmpbsa::Vector> snapshot(mol_list.size());
  size_t complexSize,receptorSize,ligandSize;
  receptorSize = ligandSize = 0;
  for(size_t i = 0;i<mol_list.size();i++)
    {
      if(mol_list[i] == MMPBSAState::RECEPTOR)
	receptorSize++;
      else if(mol_list[i] == MMPBSAState::LIGAND)
	ligandSize++;
    }
  complexSize = receptorSize + ligandSize;
  valarray<mmpbsa::Vector> complexSnap(complexSize);
  valarray<mmpbsa::Vector> receptorSnap(receptorSize);
  valarray<mmpbsa::Vector> ligandSnap(ligandSize);


  //if the program is resuming a previously started calculation, advance to the
  //last snapshot.
  if(currState.currentSnap)//zero = one = start from beginning.
    {
      mmpbsa_io::seek(trajFile,currState.currentSnap-1);
    }
  else
    {
      trajFile.curr_snap = currState.currentSnap = 1;//start from beginning if currentSnap was originally zero.
      currState.currentMolecule = MMPBSAState::COMPLEX;
    }

  // Report number of snapshots to be calculated.
#ifdef USE_MPI
  if(mpi_rank == MMPBSA_MASTER)
    {
      printf("Calculating %u Snapshots\n",currState.snapList.size());fflush(stdout);
      fflush(stdout);
    }
#else
  printf("Calculating %Lu Snapshots\n",currState.snapList.size());fflush(stdout);
  fflush(stdout);
#endif

  mmpbsa_utils::XMLNode* outputXML = previousEnergyData.getHead();

  //Walk through the snapshots. This is where MMPBSA is actually done.
  while(!mmpbsa_io::eof(trajFile))
    {
      try{
	//if a list of snaps to be run is provided, check to see if this snapshot
	//should be used. Remember: snapcounter is 1-indexed.
	//
	//Additionally, check to see if the snapshot should be run by this node
	//if multiple nodes are used, e.g. MPI
	if(currState.snapList.size())//check if the current snapshot should be skipped
	  {
	    if(!mmpbsa_utils::contains(currState.snapList,currState.currentSnap) || !should_calculate_snapshot(currState.currentSnap, currState.snapList))
	      {
		try
		  {
		    mmpbsa_io::seek(trajFile,trajFile.curr_snap+1);
		    std::cout << "Skipping Snapshot #" << currState.currentSnap << std::endl;
		    currState.currentSnap += 1;
		    currState.currentMolecule = MMPBSAState::COMPLEX;
		    continue;
		  }
		catch(mmpbsa::MMPBSAException mmpbsae)
		  {
		    if(mmpbsae.getErrType() == mmpbsa::FILE_IO_ERROR && mmpbsa_io::eof(trajFile))
		      break;// Done with snapshots.
		    printf("Throwing exception type: %d as opposed to %d  EOF = %d Message: %s\n",mmpbsae.getErrType(),mmpbsa::FILE_IO_ERROR,((mmpbsa_io::eof(trajFile))?1:0),mmpbsae.what());fflush(stdout);
		    
		    throw mmpbsae;
		  }
	      }
	  }
        
	//Can the snapshot be loaded, and if so, should it be run by this node, if multiple
	//nodes are used.
	if(get_next_snap(trajFile, snapshot) && should_calculate_snapshot(currState.currentSnap, currState.snapList))
	  std::cout << "Running Snapshot #" << currState.currentSnap << std::endl;
	else if(!should_calculate_snapshot(currState.currentSnap,currState.snapList))
	  {
	    std::cout << "Skipping snapshot #" << currState.currentSnap++ << std::endl;
	    continue;
	  }
	else
	  {
	    std::ostringstream error;
	    error << "mmpbsa_run: Error in loading snapshot #"  << ++(currState.currentSnap) << std::endl;
	    throw MMPBSAException(error,BROKEN_TRAJECTORY_FILE);
	  }
      }
      catch(MMPBSAException e)
        {
	  if(e.getErrType() == UNEXPECTED_EOF)
	    break;//return 0;
	  throw e;
        }

      study_cpu_time();

      //retrieve snapshot ID and start a snapshot XML dataset.
      std::ostringstream strSnapNumber;
      strSnapNumber << currState.currentSnap;//Node used to output energy data
      mmpbsa_utils::XMLNode* snapshotXML = new mmpbsa_utils::XMLNode("snapshot");
      snapshotXML->insertChild("ID",strSnapNumber.str());
      if(mi.snap_list_offset != 0)
        {
	  strSnapNumber.str("");
	  strSnapNumber << mi.snap_list_offset;
	  snapshotXML->insertChild("snap_list_offset",strSnapNumber.str());
        }

      //separate coordinates
      size_t complexCoordIndex = 0;
      size_t receptorCoordIndex = 0;
      size_t ligandCoordIndex = 0;
      for(size_t i = 0;i<mol_list.size();i++)
        {
	  Vector& currCoord = snapshot[i];
	  if(mol_list[i] == MMPBSAState::RECEPTOR)
            {
	      complexSnap[complexCoordIndex++] = currCoord;
	      receptorSnap[receptorCoordIndex++] = currCoord;
            }
	  else if(mol_list[i] == MMPBSAState::LIGAND)
            {
	      complexSnap[complexCoordIndex++] = currCoord;
	      ligandSnap[ligandCoordIndex++] = currCoord;
            }
        }

      //write PDB information, if requested.
      if(currState.savePDB)
        {
	  writePDB(atom_lists[MMPBSAState::COMPLEX],split_ff[MMPBSAState::COMPLEX],complexSnap,currState,"complex");
	  writePDB(atom_lists[MMPBSAState::RECEPTOR],split_ff[MMPBSAState::RECEPTOR],receptorSnap,currState,"receptor");
	  writePDB(atom_lists[MMPBSAState::LIGAND],split_ff[MMPBSAState::LIGAND],ligandSnap,currState,"ligand");
        }

      //Check to see if we just ended one snap shot and need to start the
      //next snap shot at the beginning (which is the whole COMPLEX. If
      //so, reset the current Molecule to Complex and increment the snap count.
      if(currState.currentMolecule == MMPBSAState::END_OF_MOLECULES)
        {
	  currState.currentMolecule = MMPBSAState::COMPLEX;
	  ::updateMMPBSAProgress(currState,1);
	  currState.currentSnap += 1;
	  report_boinc_progress();
	  continue;//Restarted program at the end of a snapshot. So, move on.
        }
        
      // Iterate through the three parts of the complex and calculate energies
      for(;currState.currentMolecule < MMPBSAState::END_OF_MOLECULES;++currState.currentMolecule)
        {
	  FinDiffMethod fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
	  const std::valarray<Vector> *curr_crds;
	  mmpbsa_t energy;
	  std::string mol_name;
	  int molsurf_error_flag = 0;

	  switch(currState.currentMolecule)
	    {
	    case MMPBSAState::COMPLEX:
	      curr_crds = &complexSnap;
	      mol_name = "COMPLEX";
	      break;
	    case MMPBSAState::RECEPTOR:
	      curr_crds = &receptorSnap;
	      mol_name = "RECEPTOR";
	      break;
	    case MMPBSAState::LIGAND:
	      curr_crds = &ligandSnap;
	      mol_name = "LIGAND";
	      break;
	    default:
	      throw MMPBSAException("mmpbsa_run: invalid molecule in switch.",mmpbsa::DATA_FORMAT_ERROR);
	    }

	  std::cout << "Calculating " << mol_name << std::endl;
	  // Constructor performs MM
	  EMap results(atom_lists[currState.currentMolecule],split_ff[currState.currentMolecule],*curr_crds);

	  // PB
	  energy = MeadInterface::pb_solvation(atom_lists[currState.currentMolecule],*curr_crds,fdm,radii,residues,mi.istrength);
	  results.set_elstat_solv(energy);

	  // SA
#ifdef _WIN32
	  energy = MeadInterface::molsurf_windows32(atom_lists[currState.currentMolecule],*curr_crds,radii,&molsurf_error_flag);
#else //posix
	  energy = MeadInterface::molsurf_posix(atom_lists[currState.currentMolecule],*curr_crds,radii,&molsurf_error_flag);
#endif
	  if(molsurf_error_flag != 0)
	    {
	      results.molsurf_failed = true;
	      results.set_area(0.0);
	      results.set_sasol(0.0);
	    }
	  else
	    {
	      results.molsurf_failed = false;
	      results.set_area(energy);
	      results.set_sasol(energy*mi.surf_tension+mi.surf_offset);
	    }

	  // Current molecule calculation has finished. Checkpoint.
	  thread_safe_checkpoint(mol_name.c_str(),results, currState,snapshotXML, NULL);

        }

      //MMPBSA is complete. Save the state and update the process on the
      //status, if monitoring is being done, e.g. BOINC.
      checkpoint_mmpbsa(currState);
      outputXML->insertChild(snapshotXML);
      currState.currentMolecule = MMPBSAState::COMPLEX;//Reset current molecule
      currState.currentSnap += 1;

#ifndef USE_MPI
      write_mmpbsa_data(previousEnergyData,currState);
#endif

      //If there was a snapshot list, has the list been completed?
      if(currState.snapList.size())
	if(*(currState.snapList.end()-1) == currState.currentSnap - 1)//decrement currentSnap because it was increment above.
	  break;


    }//end of snapshot loop
  study_cpu_time();

  currState.fractionDone = 1.0;
  checkpoint_mmpbsa(currState);

  for(size_t i = 0;i<MMPBSAState::END_OF_MOLECULES;i++)
    destroy(&split_ff[i]);
  delete [] split_ff;
  delete [] atom_lists;

#ifdef USE_MPI
  mpi_processes_running--;
  //For master node, wait for other processes.
  MPI_Status status;
  if(mpi_rank == MMPBSA_MASTER)
    {
      int reporter[2];//{sender,status}
      while(mpi_processes_running != 0)
    	{
	  reporter[1] = 1;
	  MPI_Recv(reporter,2, MPI_INT,MPI_ANY_SOURCE, mmpbsa_utils::STATUS, MPI_COMM_WORLD, &status);
	  while(reporter[1] != 0)
	    {
	      mmpbsa_utils::mpi_recv_mmpbsa_data(mpi_rank,reporter[0],mpi_size,currState,data_list,data_fragments);
	      reporter[1]--;
	    }
	  mpi_processes_running--;
    	}
    }
  else
    {
      int my_status[2];
      my_status[0] = mpi_rank;
      my_status[1] = 0;
      if(previousEnergyData.getHead() != 0 && previousEnergyData.getHead()->children != 0)
	my_status[1] = (int)ceil(float(previousEnergyData.toString().size())/MMPBSA_MPI_MAX_BUFFER);

      MPI_Send(my_status,2,MPI_INT,MMPBSA_MASTER,mmpbsa_utils::STATUS,MPI_COMM_WORLD);
    }
#endif

  write_mmpbsa_data(previousEnergyData,currState);

  return 0;
}

int sander_run(mmpbsa::MMPBSAState& currState,mmpbsa::SanderInterface& si)
{
  using namespace mmpbsa;
  if(si.completed || currState.currentProcess == MMPBSAState::MMPBSA)
    return 0;

  std::cout << "Starting molecular dynamics" << std::endl;
  int retval = si.start(currState.filename_map);
  if(retval)
    return retval;

  int status;
  double debug_cpu_time;
  while(true)//si has a fork running. Monitor it.
    {
      if(si.poll(status))
        {
	  if(status)
            {
	      std::cerr << "Sander(" << si.getPID()
                        << ") had a problem: ";
	      std::cerr.setf(std::ios::hex);
	      std::cerr << "0x" << status;//show status in hex.
	      std::cerr.setf(std::ios::dec);
	      std::cerr << std::endl;
	      return EXIT_CHILD_FAILED;
            }
	  break;
        }
#ifdef USE_BOINC
      ::poll_boinc_messages(si);
      if(::boinc_time_to_checkpoint())
        {

	  ::netFractionDone = overallFractionDone();
	  checkpoint_sander(currState,si);
	  std::cout << "Checkpointed." << std::endl;
        }
      std::cout << "Sander time: " << si.cpu_time();
      boinc_report_app_status(si.start_time()+si.cpu_time(),::timeAtPreviousCheckpoint,overallFractionDone());
      boinc_wu_cpu_time(debug_cpu_time);
      std::cout << " Boinc time: " << debug_cpu_time << std::endl;
      boinc_sleep(SanderInterface::pollPeriod);
#endif
    }
  si.completed = true;
  currState.fractionDone = 1.0;
  checkpoint_sander(currState,si);

  return 0;
}

std::vector<mmpbsa::MMPBSAState> parseArgs(int argc, char** argv)
{
  using std::string;
  using mmpbsa::MMPBSAState;
  std::map<string,string> pseudo_queue;
  MMPBSAState job_equivalent;
  std::vector<MMPBSAState>returnMe;
  string name,value;
  int retval;
  for(int i = 1;i<argc;i++)
    {
      string currArg = argv[i];
      if(currArg.substr(0,2) == "--")
	currArg.erase(currArg.begin(),currArg.begin()+2);

      if(currArg.find("=") != string::npos)
        {
	  name = currArg.substr(0,currArg.find("="));
	  value = currArg.substr(currArg.find("=")+1);
	  pseudo_queue[name] = value;
        }
      else
        {
	  pseudo_queue[currArg] = "";
        }
    }

  retval = parseParameter(pseudo_queue,job_equivalent);
  if(retval == 0)
    returnMe.push_back(job_equivalent);
  return returnMe;
}

int parseParameter(std::map<std::string,std::string> args, mmpbsa::MMPBSAState& currState)
{
  using mmpbsa_utils::loadListArg;
  using namespace mmpbsa;
  std::string resolved_filename;
  int returnMe = 0;
  MeadInterface& mi = currState.currentMI;
  SanderInterface& si = currState.currentSI;
  currState.receptorStartPos.insert(0);//in case these are not set manually by the use. This is the default.
  currState.ligandStartPos.insert(1);
  for(std::map<std::string,std::string>::const_iterator it = args.begin();it != args.end();it++)
    {
      if (it->second.find("=") != string::npos)
    	{
	  std::ostringstream error;
	  error << "Multiple occurrence of \"=\" in parameter: " <<
	    it->first << " = " << it->second << std::endl;
	  throw mmpbsa::MMPBSAException(error, mmpbsa::COMMAND_LINE_ERROR);
    	}

      std::istringstream buff(it->second);
      if(it->first == "rec_list")
    	{
	  currState.receptorStartPos.clear();
	  loadListArg(it->second, currState.receptorStartPos,1);
    	}
      else if(it->first == "lig_list"){
	currState.ligandStartPos.clear();
	loadListArg(it->second, currState.ligandStartPos,1);
      }
      else if(it->first == "snap_list" || it->first == "snaplist")
    	{
	  currState.snapList.clear();
	  loadListArg(it->second, currState.snapList);
    	}
      else if(it->first == "help" || it->first == "h")
    	{
	  std::cout << helpString() << std::endl;
	  return 1;
    	}
      else if(it->first == "version")
    	{
	  std::cout << PACKAGE_STRING << std::endl;
	  return 1;
    	}
      else if(it->first == "verbose")
    	{
	  std::istringstream vbuff(it->second);
	  vbuff >> mmpbsa_verbosity;
	  if(vbuff.fail())
	    throw mmpbsa::MMPBSAException("parse_parameters: \"" + it->second + "\" is an invalid verbosity level.",
					  mmpbsa::COMMAND_LINE_ERROR);
    	}
      else if(it->first == "sample_queue")
    	{
	  if(!it->second.size())
	    {
	      std::cout << "Parameter \"sample_queue\" requires a file name\"" << std::endl;
	      return 2;
	    }
	  std::string queueFilename;
	  mmpbsa_io::resolve_filename(it->second,queueFilename);
	  ::sampleQueue(queueFilename);
	  return 1;
    	}
      else if (it->first == "trust_prmtop") {
	currState.trustPrmtop = true;
      }
      else if (it->first == "mmpbsa_only") {
	si.completed = true;
	currState.currentProcess = mmpbsa::MMPBSAState::MMPBSA;
      }
      else if(it->first == "surface_area"){
	if(it->second.size() == 0)
	  {
	    std::cerr << "surface_area flag needs a molecule type (COMPLEX, RECEPTOR, LIGAND)" << std::endl;
	    return 2;
	  }
	if(it->second == "COMPLEX")
	  currState.currentMolecule = mmpbsa::MMPBSAState::COMPLEX;
	else if(it->second == "RECEPTOR")
	  currState.currentMolecule = mmpbsa::MMPBSAState::RECEPTOR;
	else if(it->second == "LIGAND")
	  currState.currentMolecule = mmpbsa::MMPBSAState::LIGAND;
	else
	  {
	    std::cerr << "surface_area flag needs a molecule type (COMPLEX, RECEPTOR, LIGAND)" << std::endl;
	    return 2;
	  }
	currState.currentProcess = mmpbsa::MMPBSAState::MOLSURF;
      }
      else if (it->first == "md_only")
    	{
	  currState.MDOnly = true;
	  return 0;
    	}
      else if(it->first == "weight")
    	{
	  std::istringstream buff(it->second);
	  buff >> currState.weight;
    	}
      else if(it->first == "id" || it->first == "prereq")//this is used by the queue system only. Not needed for calculations.
	continue;
      else if (it->first == "istrength")
    	{
	  buff >> MMPBSA_FORMAT >> mi.istrength;
    	}
      else if (it->first == "surf_offset")
    	{
	  buff >> MMPBSA_FORMAT >> mi.surf_offset;
    	}
      else if (it->first == "surf_tension")
    	{
	  buff >> MMPBSA_FORMAT >> mi.surf_tension;
    	}
      else if (it->first == "trust_prmtop")
    	{
	  currState.trustPrmtop = true;
    	}
      else if(it->first == "multithread" || it->first == "nthreads")
    	{
#ifndef USE_PTHREADS
	  std::cerr << "Warning: not compiled with threads. Ignoring multithread tag." << std::endl;
	  continue;
#endif
	  buff >> mi.multithread;
	  if(buff.fail())
	    {
	      std::cerr << "Warning: '" << it->second << "' is not a valid value for the 'multithread' flag. Not using multithreading." << std::endl;
	      mi.multithread = 0;
	    }
    	}
      else if(it->first == "traj_in_memory")
    	{
	  int bool_buff;
	  buff >> bool_buff;
	  if(buff.fail())
	    {
	      std::cerr << "Warning: could not determine whether or not to keep trajectory in memory based on:  " << it->first << " = " << it->second << " Using default: " << ((currState.keep_traj_in_mem) ? "true" : "false") << std::endl;
	    }
	  currState.keep_traj_in_mem = (bool_buff != 0);
    	}
      else if(it->first == "save_pdb")
    	{
	  currState.savePDB = (it->second != "0");
    	}
      else if(it->first == "snap_list_offset")
    	{
	  buff >> mi.snap_list_offset;
	  if(buff.fail())
	    throw mmpbsa::MMPBSAException("parse_parameters: \"" + it->second + "\" is an invalid verbosity level.",
					  mmpbsa::COMMAND_LINE_ERROR);
    	}
	else if(it->first == "overwrite")
	  {
	    if(it->second.size() > 0)
	      {
		if(it->second == "n" || it->second == "N" || it->second == "0")
		  {
		    currState.overwrite = false;
		    continue;
		  }
	      }
	    currState.overwrite = true;
	  }
    	else//assuming if the argument is not one of the parameters listed above, it's a filename
    	{
	  mmpbsa_io::resolve_filename(it->second,resolved_filename);
	  currState.filename_map[it->first] = resolved_filename;
    	}

    }
  return returnMe;
}

std::string helpString()
{
  std::string returnMe = "MMPBSA Calculations\n";
#ifdef USE_BOINC
  returnMe += "Compiled with BOINC\n";
#endif

  return returnMe + "Usage: ./mmpbsa [parameters]\n"
    "\n"
    "\nParameters:\n"
    "\ntraj=<trajectory file>"
    "\ntop=<topology file>"
    "\nradii=<radii file>"
    "\n\tSIZ radii file. If no file is values are used "
    "\n\tfrom a lookup table built into mmpbsa"
    "\nistrength=<strength value>"
    "\n\t(default = 0)"
    "\nsurf_offset=<surface offset>"
    "\n\t(default = 0.92 kcal/mol)"
    "\nsurf_tension=<surface tension value>"
    "\n\t(default = 0.00542 kcal/mol/Ang^2)"
    "\nrec_list=<comma separated list>"
    "\n\tList of beginning atoms of residues, where the"
    "\n\tlength is deduced from the parmtop file. "
    "\n\tNote: this is a one-indexed list."
    "\n\tExample: rec_list=0,13,25 would use recptors"
    "\n\tbeginning with the first, fourteenth and"
    "\n\ttwenty-sixth atoms"
    "\nlig_list=<comma separated list>"
    "\n\tList of beginning atoms of ligand."
    "\n\tSee also, rec_list."
    "\nsnap_list=<comma separated list>"
    "\n\t1-indexed list of snapshots to be included."
    "\n\tIf this option is not used, all snapshots"
    "\n\tare calculated."
    "\ntrust_prmtop"
    "\n\tOverride the Parmtop sanity check."
    "\n\tUse with caution!"
    "\nsample_queue=<filename>"
    "\n\tCreates a sample queue XML file.";
}

bool restart_sander(mmpbsa::MMPBSAState& restartState, mmpbsa::SanderInterface& si)
{
  using mmpbsa_utils::XMLParser;
  using mmpbsa::MMPBSAState;
  XMLParser xmlDoc;

  if(!has_filename(CHECKPOINT_FILE_TYPE,restartState))
    return false;

  try{
    xmlDoc.parse(get_filename(CHECKPOINT_FILE_TYPE,restartState));
    if(xmlDoc.mainTag() == MMPBSASTATE_TAG)
      restartState.currentProcess = MMPBSAState::MMPBSA;
  }
  catch(mmpbsa::XMLParserException xpe)
    {
      std::cerr << "restart_sander: Did not open " << get_filename(CHECKPOINT_FILE_TYPE,restartState) << std::endl;
      if(xpe.getErrType() == mmpbsa::FILE_IO_ERROR)
	return false;
      else
	throw xpe;
    }
  std::map<std::string,std::string> checkMap = XMLParser::mapNode(xmlDoc.getHead());

  bool usedAllParameters = true;
  std::string tag = "";
  double runtime,cpu_time,start_time;
  runtime = cpu_time = start_time = 0;
  for(std::map<std::string,std::string>::iterator it = checkMap.begin();
      it != checkMap.end();it++)
    {
      tag = it->first;
      std::istringstream buff(it->second);
      //chain to load parameters into correct variables.
      if(tag == "pid")
        {
	  int newPID = 0;
	  buff >> newPID;
	  si.setPID(newPID);
        }
      else if(tag == "cpu_time")
        {
	  buff >> std::scientific >> cpu_time;
        }
      else if(tag == "start_time")
        {
	  buff >> std::scientific >> start_time;
        }
      else if(tag == "runtime")
        {
	  buff >> std::scientific >> runtime;
        }
      else if(tag == "finished")
        {
	  buff >> si.completed;
        }
      else if(tag == "queue_position")
        {
	  buff >> restartState.placeInQueue;
        }
      else
        {
	  usedAllParameters = false;
        }
    }
  if(si.completed)
    si.set_start_time(runtime);
  else
    {
      si.set_start_time(cpu_time+start_time);
    }
  ::netCPUTime = si.start_time();
  return usedAllParameters;
}

bool restart_mmpbsa(mmpbsa::MMPBSAState& restartState)
{
  using mmpbsa::MMPBSAState;

  if(!has_filename(CHECKPOINT_FILE_TYPE,restartState))
    return false;

  using mmpbsa_utils::XMLParser;
  XMLParser xmlDoc;
  try{
    xmlDoc.parse(get_filename(CHECKPOINT_FILE_TYPE,restartState));
  }
  catch(mmpbsa::XMLParserException xpe)
    {
      if(xpe.getErrType() == mmpbsa::FILE_IO_ERROR)
	return false;
      else
	throw xpe;
    }
  std::map<std::string,std::string> checkMap = XMLParser::mapNode(xmlDoc.getHead());

  bool usedAllParameters = true;
  std::string tag = "";
  for(std::map<std::string,std::string>::iterator it = checkMap.begin();
      it != checkMap.end();it++)
    {
      tag = it->first;
      std::istringstream buff(it->second);
      //chain to load parameters into correct variables.
      if(tag == "current_molecule")
        {
	  int currMol = 0;
	  buff >> currMol;
	  switch(currMol)
            {
	    case MMPBSAState::COMPLEX:
	      restartState.currentMolecule = MMPBSAState::COMPLEX;
	      break;
	    case MMPBSAState::LIGAND:
	      restartState.currentMolecule = MMPBSAState::LIGAND;
	      break;
	    case MMPBSAState::RECEPTOR:
	      restartState.currentMolecule = MMPBSAState::RECEPTOR;
	      break;
	    default:
	      restartState.currentMolecule = MMPBSAState::END_OF_MOLECULES;
	      usedAllParameters = false;
	      break;
            }
        }//end "current_molecule" case
      else if(tag == "current_snap")
        {
	  buff >> restartState.currentSnap;
        }//end "current_snap" case
      else if(tag == "checkpoint_counter")
        {
	  buff >> restartState.checkpointCounter;
        }
      else if(tag == "queue_position")
        {
	  buff >> restartState.placeInQueue;
        }
      else if(tag == "save_pdb" && it->second != "0")
        {
	  restartState.savePDB == true;
        }
      else
        {
	  usedAllParameters = false;
        }
    }
  return usedAllParameters;
}

void checkpoint_sander(mmpbsa::MMPBSAState& saveState, mmpbsa::SanderInterface& si)
{
  std::fstream sanderProgress(HEARTBEAT_FILENAME,std::ios::in);
  if(sanderProgress.good())
    {
      double completed,remaining;
      sanderProgress >> completed >> remaining;
      saveState.fractionDone = completed/(completed+remaining);
      sanderProgress.close();
    }
  std::map<std::string,std::string> checkMap;
  std::ostringstream buff;
  buff << si.getPID();
  checkMap["pid"] = buff.str();
  buff.str("");buff << si.start_time();
  checkMap["start_time"] = buff.str();
  buff.str("");buff << si.cpu_time();
  checkMap["cpu_time"] = buff.str();
  buff.str("");buff << si.netRuntime();
  checkMap["runtime"] = buff.str();
  buff.str("");buff << si.completed;
  checkMap["finished"] = buff.str();
  buff.str("");buff << saveState.placeInQueue;
  checkMap["queue_position"] = buff.str();
  buff.str("");buff << saveState.fractionDone;
  checkMap["stage_fraction_done"] = buff.str();
  mmpbsa_utils::XMLParser xmlDoc(MDSTATE_TAG,checkMap);

  send_status_message(si, saveState.fractionDone,::timeAtPreviousCheckpoint);
  checkpoint_out(saveState,xmlDoc);

}

void checkpoint_mmpbsa(mmpbsa::MMPBSAState& saveState)
{
    
  std::map<std::string,std::string> checkMap;
  std::ostringstream buff;
  buff << saveState.currentMolecule;
  checkMap["current_molecule"] = buff.str();
  buff.str("");buff << saveState.currentSnap;
  checkMap["current_snap"] = buff.str();
  buff.str("");buff << saveState.checkpointCounter;
  checkMap["checkpoint_counter"] = buff.str();
  buff.str("");buff << saveState.placeInQueue;
  checkMap["queue_position"] = buff.str();

  if(saveState.savePDB)
    checkMap["save_pdb"] = "1";
  mmpbsa_utils::XMLParser xmlDoc(MMPBSASTATE_TAG,checkMap);
  report_boinc_progress();
  checkpoint_out(saveState,xmlDoc);
}

void checkpoint_out(mmpbsa::MMPBSAState& saveState,mmpbsa_utils::XMLParser& xmlDoc)
{
  saveState.checkpointCounter++;
  if(has_filename(CHECKPOINT_FILE_TYPE,saveState))
    xmlDoc.write(get_filename(CHECKPOINT_FILE_TYPE,saveState));
#ifdef USE_BOINC
  boinc_checkpoint_completed();
#endif
}


int mmpbsa_boinc_init()
{
#ifndef USE_BOINC
  return 0;
#else
  BOINC_OPTIONS options;
  memset(&options, 0, sizeof(options));
  options.main_program = true;
  options.check_heartbeat = true;
  options.handle_process_control = true;

#ifdef USE_GRAPHICS
  options.backwards_compatible_graphics = true;
  gshmem =(MMPBSA_SHMEM*)boinc_graphics_make_shmem("mmpbsa",sizeof(MMPBSA_SHMEM));
  if(!gshmem)
    std::cerr << "Could not create shared memory for mmpbsa." << std::endl;
  update_gshmem();
  boinc_register_timer_callback(update_gshmem);
#endif

  return boinc_init_options(&options);
#endif
}

void send_status_message(mmpbsa::SanderInterface& si, double frac_done, 
			 double checkpoint_cpu_time)
{
#ifdef USE_BOINC
  double current_cpu_time =  si.start_time() + si.cpu_time();
  boinc_report_app_status(current_cpu_time,checkpoint_cpu_time,frac_done);
#else
  return;
#endif
}

void poll_boinc_messages(mmpbsa::SanderInterface& si)
{
#ifdef USE_BOINC
  BOINC_STATUS status;
  boinc_get_status(&status);
  if (status.no_heartbeat) {
    si.kill();
    exit(0);
  }
  if (status.quit_request) {
    si.kill();
    exit(0);
  }
  if (status.abort_request) {
    si.kill();
    exit(0);
  }
  if (status.suspended) {
    if (!si.isSuspended()) {
      si.stop();
    }
  } else {
    if (si.isSuspended()) {
      si.resume();
    }
  }
#endif
  int a = 3;
  a *= 42;
}

std::vector<mmpbsa::MMPBSAState>& getQueueFile(std::vector<mmpbsa::MMPBSAState>& queue_vector, int argc,char** argv)
{
  using mmpbsa_utils::XMLParser;
  using mmpbsa::MMPBSAState;
  std::vector<MMPBSAState>& returnMe = queue_vector;

  XMLParser queueXML;
  std::string xmlFilename = "";
  std::string arg,name,value;
  for(int i = 1;i<argc;i++)
    {
      arg = argv[i];
      if(arg.find("=") == std::string::npos)
	continue;
      name = arg.substr(0,arg.find("="));
      value = arg.substr(arg.find("=")+1);
      if(name == "queue")
        {
	  mmpbsa_io::resolve_filename(value,xmlFilename);
	  break;
        }
    }

  if(xmlFilename.size() == 0)
    return returnMe;

  try
    {
      queueXML.parse(xmlFilename);
    }
  catch(mmpbsa::XMLParserException xpe)
    {
      std::cerr << "getQueueFile: Did not open "<< xmlFilename << std::endl;
      if(xpe.getErrType() == mmpbsa::FILE_IO_ERROR)
	return returnMe;
      else
	throw xpe;
    }

  const mmpbsa_utils::XMLNode * head = queueXML.getHead();
  if(head == 0)
    return returnMe;

  if(head->getName() == MMPBSA_QUEUE_TITLE)
    head = head->children;

  int queuePosition = 0;
  for(const mmpbsa_utils::XMLNode* sibling = head;sibling;sibling = sibling->siblings)
    {
      MMPBSAState nodeState;
      mmpbsa::SanderInterface si;
      std::map<std::string,std::string> tags = XMLParser::mapNode(sibling);
      if(sibling->getName() == "mmpbsa")
	nodeState.currentProcess = MMPBSAState::MMPBSA;
      else if(sibling->getName() == "molecular_dynamics" || sibling->getName() == "moledyn")
	nodeState.currentProcess = MMPBSAState::SANDER;
      else
	throw mmpbsa::MMPBSAException("getQueueFile: " + sibling->getName() + " is an unknown process type.",mmpbsa::BAD_XML_TAG);
      parseParameter(tags,nodeState);
      nodeState.placeInQueue = queuePosition++;
      returnMe.push_back(nodeState);
    }

  return returnMe;
}

void sampleQueue(const std::string& filename)
{
  using std::map;
  using mmpbsa_utils::XMLParser;
  using mmpbsa_utils::XMLNode;

  XMLNode* mmpbsaXML = new XMLNode("mmpbsa");
  mmpbsaXML->insertChild(MMPBSA_TOPOLOGY_TYPE,"sander_prmtop_file.prmtop");
  mmpbsaXML->insertChild(MMPBSA_TRAJECTORY_TYPE,"sander_snapshot_file.mdcrd");
  mmpbsaXML->insertChild("radii","DelPhi_radii_file.siz");
  mmpbsaXML->insertChild("mmpbsa_out","mmpbsa-result-output.out");
  mmpbsaXML->insertChild("snap_list","1,3");

#ifndef USE_BOINC
  XMLParser writeMe(mmpbsaXML);
  writeMe.write(filename);
#else
  XMLNode* sanderXML = new XMLNode("molecular_dynamics");
  sanderXML->insertChild("mdin","sander_input.in");
  sanderXML->insertChild(SANDER_MDOUT_TYPE,"sander_output.out");
  sanderXML->insertChild("restart","sander_restart.rst");
  sanderXML->insertChild(SANDER_INPCRD_TYPE,"sander_input_coordinates.inpcrd");
  sanderXML->insertChild(SANDER_PRMTOP_TYPE,"sander_prmtop_file.prmtop");
  sanderXML->insertChild("mdcrd","sander_snapshot_file.mdcrd");
  sanderXML->insertChild("checkpoint","checkpoint_file_name.xml");

  XMLNode* theDoc = new XMLNode(MMPBSA_QUEUE_TITLE);
  theDoc->insertChild(sanderXML);
  theDoc->insertChild(mmpbsaXML);

  XMLParser writeMe(theDoc);
  writeMe.write(filename);
#endif

}

void updateMMPBSAProgress(mmpbsa::MMPBSAState& currState,const double& increment)
{
  if(currState.snapList.size())
    {
      currState.fractionDone += increment/currState.snapList.size();
    }
  else
    {
      currState.fractionDone = 0;//if we do not know how many snapshots we're reading, do we know the fraction done?
    }
}

double overallFractionDone()
{
  double totalWeight = 0;
  double completed = 0;
  for(size_t i = 0;i< ::processQueue.size();i++)
    {
      totalWeight += processQueue[i].weight;
      completed += processQueue[i].fractionDone*processQueue[i].weight;
    }
  return completed/totalWeight;
}

void report_boinc_progress()
{
#ifdef USE_BOINC
  ::netFractionDone = overallFractionDone();
  //boinc_fraction_done(completed/totalWeight);
  double cpu_time;
  ::boinc_wu_cpu_time(cpu_time);
  ::boinc_report_app_status(cpu_time,::timeAtPreviousCheckpoint,::netFractionDone);
#else
  return;
#endif
}

void update_gshmem()
{
#if defined(USE_BOINC) && defined(USE_GRAPHICS)
  if(!gshmem)
    return;
  gshmem->update_time = dtime();//without this, the graphics app will think we've died.
  if(gshmem->countdown > 0)
    gshmem->countdown--;
  else
    return;
  gshmem->fraction_done = ::netFractionDone;
  gshmem->cpu_time = ::netCPUTime + boinc_elapsed_time();
  boinc_get_status(&gshmem->status);
#else
  return;
#endif
}

#if 0//why use threads if it's not supported in mead or molsurf
void *do_mmpbsa_calculation_thread(void* args)
{
  struct mmpbsa_thread_arg * thread_args = (mmpbsa_thread_arg*) args;
  if(thread_args == 0)
    throw mmpbsa::MMPBSAException("do_mmpbsa_calculation_thread: given null pointer for argument structure");
  mmpbsa::EMap theEMap = mmpbsa::MeadInterface::full_EMap(*thread_args->atoms,*thread_args->ff,*thread_args->snap,
							  *thread_args->fdm,
							  *thread_args->pradii,*thread_args->residues,thread_args->mi->istrength,
							  thread_args->mi->surf_tension,thread_args->mi->surf_offset);
  thread_safe_checkpoint(thread_args->next_mole,thread_args->mole_name,
			 theEMap,*thread_args->currState,thread_args->snapshotXML,thread_args->mmpbsa_mutex);
#ifdef USE_PTHREADS
  std::cout << "Finished " << thread_args->mole_name << std::endl;
  delete thread_args;
#endif
}

int do_mmpbsa_calculation(void* thread_object,int useMultithread,
			  const std::vector<mmpbsa::atom_t>& atoms, const mmpbsa::forcefield_t& ff,
			  const std::valarray<mmpbsa::Vector>& Snap,
			  const FinDiffMethod& fdm,const std::map<std::string,float>& radii,
			  const std::map<std::string,std::string>& residues,
			  const mmpbsa::MeadInterface& mi,
			  mmpbsa::MMPBSAState& currState,mmpbsa_utils::XMLNode* snapshotXML,
			  mmpbsa::MMPBSAState::MOLECULE next_mole, const char* mole_name, void * mmpbsa_mutex)
{
  struct mmpbsa_thread_arg * pass_these = new struct mmpbsa_thread_arg;
  pass_these->atoms = &atoms;
  pass_these->ff = &ff;
  pass_these->snap = &Snap;
  pass_these->fdm = &fdm;
  pass_these->pradii = &radii;
  pass_these->residues = &residues;
  pass_these->mi = &mi;
  pass_these->currState = &currState;
  pass_these->snapshotXML = snapshotXML;
  pass_these->next_mole = next_mole;
  pass_these->mmpbsa_mutex = mmpbsa_mutex;
  pass_these->mole_name = mole_name;


#ifdef USE_PTHREADS
  if(useMultithread)
    {
      std::cout << "Running thread " << mole_name << "(" << &pass_these << ", " << pass_these->snap <<  ")" << std::endl;
      return pthread_create((pthread_t*)thread_object,&attr,do_mmpbsa_calculation_thread,(void*) pass_these);
    }
#endif
  std::cout << "Calculating " << mole_name << std::endl;
  (*do_mmpbsa_calculation_thread)((void*) pass_these);
  return 0;
}
#endif

void thread_safe_checkpoint(const char* mole_name,
			    const mmpbsa::EMap& EMap, mmpbsa::MMPBSAState& currState,
			    mmpbsa_utils::XMLNode* snapshotXML, void * mmpbsa_mutex)
{
  if(snapshotXML == 0)
    throw mmpbsa::MMPBSAException("mmpbsa_update: given a null pointer for the snapshot output file.");

#ifdef USE_PTHREADS
  pthread_mutex_t* pMutex = (pthread_mutex_t*)mmpbsa_mutex;
  if(pMutex != 0)
    pthread_mutex_lock(pMutex);
#endif

  snapshotXML->insertChild(EMap.toXML(((mole_name) ? mole_name : "UNKNOWN_MOLECULE")));
  ::updateMMPBSAProgress(currState,0.33333333);
  checkpoint_mmpbsa(currState);
  study_cpu_time();
#ifdef USE_PTHREADS
  if(pMutex != 0)
    pthread_mutex_unlock(pMutex);
#endif
}


