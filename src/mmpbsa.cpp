#include "mmpbsa.h"

int main(int argc, char** argv)
{
	using mmpbsa::MMPBSAState;

#ifdef USE_PTHREADS
	pthread_mutex_init(&mmpbsa_mutex,NULL);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
#endif

	std::cerr << PACKAGE_STRING <<" started on " << mmpbsa_utils::get_human_time() << std::endl;

	try
    {
        ::timeAtPreviousCheckpoint = 0;
        mmpbsa_boinc_init();//must be called before any other BOINC routines. If BOINC is not used, nothing will happen.
        int retval = 0;

        ::processQueue = getQueueFile(argc,argv);
        //If no queue was found, run based on arguments in argv. These should
        //correspond to what is placed in the queue XML file.
        if(processQueue.size() == 0)
        {
            MMPBSAState currState;
            std::map<std::string,std::string> argMap = parseArgs(argc,argv);
            mmpbsa::SanderInterface si;
            retval = parseParameter(argMap,currState,si);
            if(retval == 0)
            {
            	sander_run(currState,si);
                if(!retval && !currState.MDOnly)
                {
                    mmpbsa::MeadInterface mi;
                    retval = parseParameter(argMap,currState,mi);
                    if(!retval)
                        retval = mmpbsa_run(currState,mi);
                }
            }
            else
                retval--;
        }
        else
        {
            int queuePosition= 0;
            for(std::vector<MMPBSAState>::iterator it = processQueue.begin();
                    it != processQueue.end();it++)
            {
                switch(it->currentProcess)
                {
                    case MMPBSAState::SANDER:
                        restart_sander(*it,it->currentSI);
                        if(it->placeInQueue > queuePosition)
                        {
                            queuePosition++;
                            continue;
                        }
                        else if(it->placeInQueue < queuePosition)
                        {
                            it->currentSI.completed = false;
                            it->fractionDone = 0;
                            it->placeInQueue = queuePosition;
                        }
                        retval = sander_run(*it,it->currentSI);
                        break;
                    case MMPBSAState::MMPBSA:
                        restart_mmpbsa(*it);
                        if(it->placeInQueue > queuePosition)
                        {
                            queuePosition++;
                            continue;
                        }
                        else if(it->placeInQueue < queuePosition)
                        {
                            it->fractionDone = 0;
                            it->placeInQueue = queuePosition;
                        }
                        retval = mmpbsa_run(*it,it->currentMI);
                        break;
                }
                if(retval)
                    break;
                queuePosition++;
            }
        }
#ifdef USE_BOINC
        if(retval)
            std::cerr << "BOINC Error: " << boincerror(retval) << std::endl;
        boinc_finish(retval);
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
        std::cerr << PACKAGE_STRING <<" finished (" << e.getErrType() << ") on " << mmpbsa_utils::get_human_time() << "\n" << std::endl;
        return e.getErrType();

    }
}

void writePDB(const std::vector<mmpbsa::atom_t>& atoms,const mmpbsa::forcefield_t& ff,
		const std::valarray<mmpbsa_t>crds,const mmpbsa::MMPBSAState& currState,const std::string& molecule)
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

void write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState)
{
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
}

mmpbsa_utils::XMLNode* read_mmpbsa_data(const mmpbsa::MMPBSAState& currState)
{
	const string& filename = get_filename(SANDER_MDOUT_TYPE,currState);
	std::stringstream data;
	std::ios::openmode the_mode = std::ios::in;
	if(filename.find(".gz") != std::string::npos || filename.find(".tar") != std::string::npos || filename.find(".tgz") != std::string::npos)
		the_mode |= std::ios::binary;
	std::fstream in_file(filename.c_str(),the_mode);
	if(!in_file.good())
		throw mmpbsa::MMPBSAException("read_mmpbsa_data: unable to open " + filename + " for writing.",mmpbsa::FILE_IO_ERROR);
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
    if(!has_filename(SANDER_PRMTOP_TYPE,currState))
    	throw mmpbsa::MMPBSAException("get_sander_forcefield: no parmtop file.",BROKEN_PRMTOP_FILE);
    sp->raw_read_amber_parm(get_filename(SANDER_PRMTOP_TYPE,currState));
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
    for(size_t i = 0;i<currState.receptorStartPos.size();i++)
    {
        size_t currPos = currState.receptorStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currReceptor(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorSize += top-bottom;
    }
    for(size_t i = 0;i<currState.ligandStartPos.size();i++)
    {
        size_t currPos = currState.ligandStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
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

    delete sp;
}

void get_forcefield(mmpbsa::MMPBSAState& currState,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list,mmpbsa_io::trajectory_t& trajfile)
{
#ifdef USE_GROMACS
	if(!has_filename(SANDER_PRMTOP_TYPE,currState))
		throw mmpbsa::MMPBSAException("get_forcefield: no parmtop file.",mmpbsa::BROKEN_PRMTOP_FILE);
	std::string filename = get_filename(SANDER_PRMTOP_TYPE,currState);
	if(filename.find(".tpr") != std::string::npos)
	{
		mmpbsa_io::get_gromacs_forcefield(filename.c_str(),split_ff,atom_lists,mol_list);
		return;
	}
#endif
	get_sander_forcefield(currState,split_ff,atom_lists,mol_list,trajfile);
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

    std::cout << "Starting MMPBSA calculation " << std::endl;


    if(!has_filename(SANDER_MDOUT_TYPE,currState))
        currState.filename_map[SANDER_MDOUT_TYPE] = "mmpbsa-output.xml";

    //Upon restart, MMPBSA needs to reload energy that was calculated previously and then
    //append new data to it.
    mmpbsa_utils::XMLParser previousEnergyData;
    try{
    	mmpbsa_utils::XMLNode* old_data = read_mmpbsa_data(currState);
    	if(old_data == 0)
    		previousEnergyData.setHead(new mmpbsa_utils::XMLNode("mmpbsa_energy"));
    	else
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

    //Setup Trajectory Structure
    if(!has_filename(SANDER_INPCRD_TYPE,currState))
    	throw mmpbsa::MMPBSAException("mmpbsa_run: no trajectory file was given.",BROKEN_TRAJECTORY_FILE);
    mmpbsa_io::trajectory_t trajFile = mmpbsa_io::open_trajectory(get_filename(SANDER_INPCRD_TYPE,currState));



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
        std::cout << "Reading radii from " << radiiFilename  << std::endl;
        radiiFile.close();
        mmpbsa_io::read_siz_file(radiiData,radii, residues);
    }

    //setup trajectory storage
    get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
    valarray<mmpbsa_t> snapshot(mol_list.size()*3);
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
    valarray<mmpbsa_t> complexSnap(complexSize*3);
    valarray<mmpbsa_t> receptorSnap(receptorSize*3);
    valarray<mmpbsa_t> ligandSnap(ligandSize*3);


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

    mmpbsa_utils::XMLNode* outputXML = previousEnergyData.getHead();

    //Walk through the snapshots. This is where MMPBSA is actually done.
    while(!mmpbsa_io::eof(trajFile))
    {
        try{
            //if a list of snaps to be run is provided, check to see if this snapshot
            //should be used. Remember: snapcounter is 1-indexed.
            if(currState.snapList.size())//check if the current snapshot should be skipped
                if(!mmpbsa_utils::contains(currState.snapList,currState.currentSnap))
                {
                    mmpbsa_io::seek(trajFile,trajFile.curr_snap+1);
                    std::cout << "Skipping Snapshot #" << currState.currentSnap << std::endl;
                    currState.currentSnap += 1;
                    currState.currentMolecule = MMPBSAState::COMPLEX;
                    continue;
                }
        
            if(get_next_snap(trajFile, snapshot))
                std::cout << "Running Snapshot #" << currState.currentSnap << std::endl;
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
            {
            	write_mmpbsa_data(previousEnergyData,currState);
            	return 0;
            }
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
            std::slice_array<mmpbsa_t> currCoord = snapshot[slice(3*i,3,1)];
            if(mol_list[i] == MMPBSAState::RECEPTOR)
            {
            	complexSnap[slice(3*complexCoordIndex++,3,1)] = currCoord;
            	receptorSnap[slice(3*receptorCoordIndex++,3,1)] = currCoord;
            }
            else if(mol_list[i] == MMPBSAState::LIGAND)
            {
            	complexSnap[slice(3*complexCoordIndex++,3,1)] = currCoord;
                ligandSnap[slice(3*ligandCoordIndex++,3,1)] = currCoord;
            }
        }

                //write PDB information, if requested.
        if(currState.savePDB)
        {
        	writePDB(atom_lists[MMPBSAState::COMPLEX],split_ff[MMPBSAState::COMPLEX],complexSnap,currState,"complex");
        	writePDB(atom_lists[MMPBSAState::RECEPTOR],split_ff[MMPBSAState::RECEPTOR],receptorSnap,currState,"receptor");
        	writePDB(atom_lists[MMPBSAState::LIGAND],split_ff[MMPBSAState::LIGAND],ligandSnap,currState,"ligand");
        }

        //This is the section where we actually do the MMPBSA calculations.
        //For the sake of multithreading, one FDM will be made and
        //two pointers will be created for ligand and receptor. If
        //multithreading is used, eg on a multicore processor,
        //those pointers will point to unique FDM objects. If
        //multithreading is not used, they will point to comp_fdm.
        //Number of threads is specified in the MeadInterface object,
        //to allow the number of threads to be dynamically set.
        FinDiffMethod *ligand_fdm,*recpt_fdm,comp_fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
    	ligand_fdm = recpt_fdm = &comp_fdm;
    	if(mi.multithread)
        {
        	ligand_fdm = new FinDiffMethod;
        	*ligand_fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
        	recpt_fdm = new FinDiffMethod;
        	*recpt_fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
        }

        //Radii information
        map<std::string,float>* pradii = &(mi.brad);//don't delete!!!
        if(radii.size())//if the radius map is empty, use MeadInterface's lookup table.
            pradii = &radii;

        //output-ing is broken up by section, in case the program needs to be
        //monitored or paused.

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
        
        //MMPBSA on Complex
        int comp_thread_retval, recpt_thread_retval, ligand_retval;
    	//map<std::string,float> recpt_radii = *pradii;

#ifdef USE_PTHREADS
      pthread_t comp_thread,recpt_thread,ligand_thread;
#else
      unsigned short comp_thread, recpt_thread, ligand_thread,mmpbsa_mutex;//will be ignored. Dummy variable name
#endif
        if(currState.currentMolecule == MMPBSAState::COMPLEX)
        {
        	do_mmpbsa_calculation((void*) &comp_thread, mi.multithread, atom_lists[MMPBSAState::COMPLEX], split_ff[MMPBSAState::COMPLEX], complexSnap,comp_fdm,*pradii,residues,mi,currState,snapshotXML,MMPBSAState::RECEPTOR,"COMPLEX",(void*)&mmpbsa_mutex);
        }

        //MMPBSA on Receptor
        if(currState.currentMolecule == MMPBSAState::RECEPTOR || mi.multithread)
        {
        	do_mmpbsa_calculation((void*) &recpt_thread, mi.multithread,  atom_lists[MMPBSAState::RECEPTOR], split_ff[MMPBSAState::RECEPTOR], receptorSnap,*recpt_fdm,*pradii,residues,mi,currState,snapshotXML,MMPBSAState::LIGAND,"RECEPTOR",(void*)&mmpbsa_mutex);
        }

        //MMPBSA on Ligand
        if(currState.currentMolecule == MMPBSAState::LIGAND || mi.multithread)
        {
        	do_mmpbsa_calculation((void*) &ligand_thread,mi.multithread,  atom_lists[MMPBSAState::LIGAND], split_ff[MMPBSAState::LIGAND], ligandSnap,*ligand_fdm,*pradii,residues,mi,currState,snapshotXML,MMPBSAState::END_OF_MOLECULES,"LIGAND",(void*)&mmpbsa_mutex);
        }

#ifdef USE_PTHREADS
        if(mi.multithread)
        {
			std::cout << "Waiting to join threads.... ";
			pthread_join(ligand_thread,NULL);
			pthread_join(comp_thread,NULL);
			pthread_join(recpt_thread,NULL);
			std::cout << "Finished."<< std::endl;
        }
#endif
        //MMPBSA is complete. Save the state and update the process on the
        //status, if monitoring is being done, e.g. BOINC.
        checkpoint_mmpbsa(currState);
        outputXML->insertChild(snapshotXML);
        currState.currentMolecule = MMPBSAState::COMPLEX;//Reset current molecule
        currState.currentSnap += 1;

        write_mmpbsa_data(previousEnergyData,currState);

        //clean up FinDiffMethods
        if(ligand_fdm != 0 && ligand_fdm != &comp_fdm)
        	delete ligand_fdm;
        if(recpt_fdm != 0 && recpt_fdm != ligand_fdm && recpt_fdm != &comp_fdm)
        	delete recpt_fdm;

    }//end of snapshot loop
    study_cpu_time();

    currState.fractionDone = 1.0;
    checkpoint_mmpbsa(currState);

    write_mmpbsa_data(previousEnergyData,currState);


    for(size_t i = 0;i<MMPBSAState::END_OF_MOLECULES;i++)
    	destroy(&split_ff[i]);
    delete [] split_ff;
    delete [] atom_lists;
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

std::map<std::string,std::string> parseArgs(int argc, char** argv)
{
    using std::string;
    std::map<std::string,std::string> returnMe;

    string name,value;
    for(int i = 1;i<argc;i++)
    {
        string currArg = argv[i];
        if(currArg.substr(0,2) == "--")
            currArg.erase(currArg.begin(),currArg.begin()+2);

        if(currArg.find("=") != string::npos)
        {
            name = currArg.substr(0,currArg.find("="));
            value = currArg.substr(currArg.find("=")+1);
            returnMe[name] = value;
        }
        else if(currArg == "verbose")
        {
        	mmpbsa_verbosity = 1;
        }
        else
        {
            returnMe[currArg] = "";
        }
    }
    return returnMe;
}

int parseParameter(std::map<std::string,std::string> args, mmpbsa::MMPBSAState& currState, mmpbsa::MeadInterface& mi)
{
    for(std::map<std::string,std::string>::const_iterator it = args.begin();it != args.end();it++)
    {
        if (it->second.find("=") != std::string::npos)
        {
            std::ostringstream error;
            error << "Multiple occurrence of \"=\" in parameter: "
                    << it->first << "=" << it->second << std::endl;
            throw mmpbsa::MMPBSAException(error, mmpbsa::COMMAND_LINE_ERROR);
        }

        std::istringstream buff(it->second);
        if (it->first == "istrength")
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
        else if(it->first == "multithread")
        {
#ifndef USE_PTHREADS
        	std::cerr << "Warning: not compiled with threads. multithread tag." << std::endl;
        	continue;
#endif
        	buff >> mi.multithread;
        	if(buff.fail())
        	{
        		std::cerr << "Warning: could not determine number of threads from " << it->first << " = " << it->second << "  Not using multithreading." << std::endl;
        		mi.multithread = 0;
        	}
        }
        else if (it->first == "help" || it->first == "h")
		{
			std::cout << helpString() << std::endl;
			return 1;
		}
		else if(it->first == "version")
		{
			std::cout << PACKAGE_STRING << std::endl;
			return 1;
		}
		else if(it->first == "save_pdb")
		{
			currState.savePDB = (it->second != "0");
		}
		else if(it->first == "verbose")
		{
			buff >> mmpbsa_verbosity;
			if(buff.fail())
				throw mmpbsa::MMPBSAException("parse_parameters: \"" + it->second + "\" is an invalid verbosity level.",
						mmpbsa::COMMAND_LINE_ERROR);
		}
		else if(it->first == "snap_list_offset")
		{
			buff >> mi.snap_list_offset;
			if(buff.fail())
				throw mmpbsa::MMPBSAException("parse_parameters: \"" + it->second + "\" is an invalid verbosity level.",
						mmpbsa::COMMAND_LINE_ERROR);
		}

    }
    return 0;
}

int parseParameter(std::map<std::string,std::string> args, mmpbsa::MMPBSAState& currState, mmpbsa::SanderInterface& si)
{
    using mmpbsa_utils::loadListArg;
    std::string resolved_filename;
    int returnMe = 0;
    currState.receptorStartPos.push_back(0);//in case these are not set manually by the use. This is the default.
    currState.ligandStartPos.push_back(1);
    for(std::map<std::string,std::string>::const_iterator it = args.begin();it != args.end();it++)
    {
        if (it->second.find("=") != string::npos)
        {
            std::ostringstream error;
            error << "Multiple occurrence of \"=\" in parameter: " <<
                    it->first << " = " << it->second << std::endl;
            throw mmpbsa::MMPBSAException(error, mmpbsa::COMMAND_LINE_ERROR);
        }

        if (it->first == "rec_list") {
            currState.receptorStartPos.clear();
            loadListArg(it->second, currState.receptorStartPos,1);
        } else if (it->first == "lig_list") {
            currState.ligandStartPos.clear();
            loadListArg(it->second, currState.ligandStartPos,1);
        } else if (it->first == "snap_list" || it->first == "snaplist") {
            currState.snapList.clear();
            loadListArg(it->second, currState.snapList);
        } else if (it->first == "help" || it->first == "h") {
            std::cout << helpString() << std::endl;
            return 1;
        } else if(it->first == "version") {
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
            return 0;
        }
        else if (it->first == "mmpbsa_only") {
            si.completed = true;
            currState.currentProcess = mmpbsa::MMPBSAState::MMPBSA;
            return 0;
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
        else if(it->first == "id" || it->first == "prereq" || it->first == "snap_list_offset")//this is used by the queue system only. Not needed for calculations.
        	continue;
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
    "\ncreate_sample=<filename>"
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

std::vector<mmpbsa::MMPBSAState> getQueueFile(int argc,char** argv)
{
    using mmpbsa_utils::XMLParser;
    using mmpbsa::MMPBSAState;
    std::vector<MMPBSAState> returnMe;

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

    if(head->getName() == MMPBSA_XML_TITLE)
        head = head->children;

    int queuePosition = 0;
    for(const mmpbsa_utils::XMLNode* sibling = head;sibling;sibling = sibling->siblings)
    {
        MMPBSAState nodeState;
        mmpbsa::SanderInterface si;
        std::map<std::string,std::string> tags = XMLParser::mapNode(sibling);
        if(sibling->getName() == "mmpbsa")
        {
            mmpbsa::MeadInterface mi;
            parseParameter(tags,nodeState,si);
            parseParameter(tags,nodeState,mi);
            nodeState.currentMI = mi;
            nodeState.currentProcess = MMPBSAState::MMPBSA;
            nodeState.placeInQueue = queuePosition++;
            returnMe.push_back(nodeState);
        }
        else if(sibling->getName() == "molecular_dynamics" || sibling->getName() == "moledyn")
        {
            parseParameter(tags,nodeState,si);
            nodeState.currentSI = si;
            nodeState.currentProcess = MMPBSAState::SANDER;
            nodeState.placeInQueue = queuePosition++;
            returnMe.push_back(nodeState);
        }
        else
        	throw mmpbsa::MMPBSAException("getQueueFile: " + sibling->getName() + " is an unknown process type.",mmpbsa::BAD_XML_TAG);
    }

    return returnMe;
}

void sampleQueue(const std::string& filename)
{
    using std::map;
    using mmpbsa_utils::XMLParser;
    using mmpbsa_utils::XMLNode;

    XMLNode* sanderXML = new XMLNode("molecular_dynamics");
    sanderXML->insertChild("mdin","sander_input.in");
    sanderXML->insertChild("mdout","sander_output.out");
    sanderXML->insertChild("restart","sander_restart.rst");
    sanderXML->insertChild("inpcrd","sander_input_coordinates.inpcrd");
    sanderXML->insertChild("prmtop","sander_prmtop_file.prmtop");
    sanderXML->insertChild("mdcrd","sander_snapshot_file.mdcrd");
    sanderXML->insertChild("checkpoint","checkpoint_file_name.xml");

    XMLNode* mmpbsaXML = new XMLNode("mmpbsa");
    mmpbsaXML->insertChild("prmtop","sander_prmtop_file.prmtop");
    mmpbsaXML->insertChild("inpcrd","sander_snapshot_file.mdcrd");
    mmpbsaXML->insertChild("radii","DelPhi_radii_file.siz");
    mmpbsaXML->insertChild("mmpbsa_out","mmpbsa-result-output.out");
    mmpbsaXML->insertChild("snap_list","1,3");
    mmpbsaXML->insertChild("checkpoint","checkpoint_file_name.xml");
    
    XMLNode* theDoc = new XMLNode(MMPBSA_XML_TITLE);
    theDoc->insertChild(sanderXML);
    theDoc->insertChild(mmpbsaXML);
    
    XMLParser writeMe(theDoc);
    writeMe.write(filename);
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

bool has_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state)
{
	return the_state.filename_map.find(filetype) != the_state.filename_map.end();
}

const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state) throw (mmpbsa::MMPBSAException)
{
	if(the_state.filename_map.find(filetype) == the_state.filename_map.end())
		throw mmpbsa::MMPBSAException("get_filename: No filename provided for the file type " + filetype,mmpbsa::COMMAND_LINE_ERROR);
	return the_state.filename_map.find(filetype)->second;
}

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
		const std::valarray<mmpbsa_t>& Snap,
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

void thread_safe_checkpoint(const mmpbsa::MMPBSAState::MOLECULE& next_mole,const char* mole_name,
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
	currState.currentMolecule = next_mole;
	::updateMMPBSAProgress(currState,0.33333333);
	checkpoint_mmpbsa(currState);
	study_cpu_time();
#ifdef USE_PTHREADS
	if(pMutex != 0)
		pthread_mutex_unlock(pMutex);
#endif
}


