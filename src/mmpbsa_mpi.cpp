#include "mmpbsa_mpi.h"

void mmpbsa_utils::mpi_init_hosts(int* argc, char*** argv, int& mpi_rank,int& mpi_size)
{
	//Initialize
	MPI_Init(argc,argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
}

void mmpbsa_utils::mpi_store_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list)
{
	mmpbsa_utils::XMLNode *id_node,*it;
	int index = 0;

	if(data_list == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store energy data into an empty list (which is a null pointer).",mmpbsa::NULL_POINTER);

	it = energy_data.getHead();
	if(it == 0 || it->children == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store empty energy data." ,mmpbsa::DATA_FORMAT_ERROR);

	data_list->insertChild(it->children);
	energy_data.getHead()->children = 0;
}

void mmpbsa_utils::mpi_store_mmpbsa_data(char* data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list)
{
	mmpbsa_utils::XMLNode* head;
	std::stringstream buff;

	if(data_list == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store in an empty data list, which is a null pointer.",mmpbsa::NULL_POINTER);

	buff << data;
	head = mmpbsa_utils::XMLParser::parse(buff);
	mmpbsa_utils::XMLParser parser(head);
	mpi_store_mmpbsa_data(parser,data_index,data_list);
	parser.detachHead();
}


int mmpbsa_utils::mpi_send_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const int& mpi_rank)
{
	std::string str_data = energy_data.toString();
	size_t data_length = str_data.size()+1;// + MMPBSA_MPI_ID_WORD;
	char* data = (char*)calloc(data_length,sizeof(char));//add 1 to str_data.size() for the null character at the end and then add to the beginning for mpi_rank of the sender.
	strcpy(data,str_data.c_str());

	int returnMe = MPI_Send(data, data_length*sizeof(char), MPI_CHAR,MMPBSA_MASTER,mmpbsa_utils::DATA, MPI_COMM_WORLD);
	free(data);
	return returnMe;
}

int mmpbsa_utils::mpi_recv_mmpbsa_data(const int& my_rank, const int& source_rank,
		const int& mpi_size, const mmpbsa::MMPBSAState& currState,
		mmpbsa_utils::XMLNode* data_list)
{
	if(mpi_size == 0 || my_rank != 0)
		return 0;

	char *mpi_data = (char*)calloc(MMPBSA_MPI_MAX_BUFFER,sizeof(char));
	MPI_Status status;

	//Receive Data from node
	int returnMe = MPI_Recv(mpi_data, MMPBSA_MPI_MAX_BUFFER, MPI_CHAR,
			source_rank, mmpbsa_utils::DATA, MPI_COMM_WORLD, &status);

	//Is data good?
	if(returnMe != 0)
	{
		std::ostringstream error;
		int errlen;
		MPI_Error_string(returnMe, mpi_data, &errlen);

		error << "mmpbsa_utils::mpi_recv_mmpbsa_data: MPI_Recv had a problem" << std::endl;
		error << "Message (" << returnMe << "): " << mpi_data;
		free(mpi_data);
		throw mmpbsa::MMPBSAException(error,mmpbsa::MPI_ERROR);
	}

	//Store data and write if needed
	mpi_store_mmpbsa_data(mpi_data,42,data_list);
	free(mpi_data);
	return 0;

}

int mmpbsa_utils::mpi_write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState,
		const int& mpi_rank, mmpbsa_utils::XMLNode* data_list)
{
	if(mpi_rank == MMPBSA_MASTER)
	{
		mpi_store_mmpbsa_data(energy_data,mpi_rank,data_list);
		//while(mpi_update_data_file(currState,mpi_rank,data_list,next_node))1;
		return 0;
	}
	else
	{
		if(energy_data.getHead() == 0 || energy_data.getHead()->children == 0)
			return 0;
		return mpi_send_mmpbsa_data(energy_data,mpi_rank);
	}
}

void mmpbsa_utils::mpi_finish_output(const int& mpi_rank,const int& mpi_size,const std::string& mmpbsa_output_filename)
{
	if(mpi_rank != 0 || mpi_size == 1)
		return;
	std::fstream output(mmpbsa_output_filename.c_str(),std::ios::out | std::ios::app);

	if(!output.good())
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_finish_output: Could not open file for writing.",mmpbsa::FILE_IO_ERROR);

	output << "</" << MMPBSA_XML_TITLE << ">" << std::endl;
	output.close();
}

int mpi_mmpbsa_sorter(mmpbsa_utils::XMLNode  *test_element, mmpbsa_utils::XMLNode  *partition_element)
{
	std::string te_ID,pe_ID;
	mmpbsa_utils::XMLNode* it;
	for(it = test_element->children;it != 0;it = it->siblings)
		if(it->getName() == "ID")
		{
			te_ID = it->getText();
			break;
		}
	for(it = partition_element->children;it != 0;it = it->siblings)
		if(it->getName() == "ID")
		{
			pe_ID = it->getText();
			break;
		}
	if(te_ID.size() == 0 || pe_ID.size() == 0)
	{
		std::cerr << "mpi_mmpbsa_sorter: could not get id from one of the nodes" << std::endl;
		exit(1);
	}

	std::istringstream buff;
	int peid,teid;
	buff.str(pe_ID);
	buff >> peid;
	if(buff.fail())
	{
		std::cerr << "mpi_mmpbsa_sorter: could not get id from partition element" << std::endl;
		exit(1);
	}
	buff.clear();buff.str(te_ID);
	buff >> teid;
	if(buff.fail())
	{
		std::cerr << "mpi_mmpbsa_sorter: could not get id from test element" << std::endl;
		exit(1);
	}

	int returnMe = 0;
	if(teid == peid)
		returnMe += 1;
	if(teid > peid)
		returnMe += 2;
	if(teid < peid)
		returnMe += 4;

	return returnMe;
}

void mmpbsa_utils::mpi_dump_data(mmpbsa_utils::XMLNode* data_list,const std::string& filename)
{
	if(data_list == 0)
		return;
	std::fstream output(filename.c_str(),std::ios::out | std::ios::app);
	if(!output.good())
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_dump_data: Could not open file for writing.",mmpbsa::FILE_IO_ERROR);

	std::vector<mmpbsa_utils::XMLNode *> snaps = quick_sort(data_list, mpi_mmpbsa_sorter);
	//std::vector<mmpbsa_utils::XMLNode *>::const_iterator snap = snaps.begin();
	//for(;snap != snaps.end();snap++)
	mmpbsa_utils::XMLNode *temp;
	for(size_t i = 0;i<snaps.size();i++)
		if(snaps.at(i) != 0)
		{
			temp = snaps.at(i)->siblings;
			snaps.at(i)->siblings = 0;
			output << snaps.at(i)->toString() << std::endl;
			snaps.at(i)->siblings = temp;
		}
	output.close();
}

#if 0//probably don't want this
void insert(mmpbsa_utils::mpi_node* new_node,mmpbsa_utils::mpi_node** node_list_addr)
{
	size_t next_index;
	mmpbsa_utils::mpi_node *node_list, *next_node;

	if(node_list_addr == 0)
		throw mmpbsa::MMPBSAException("insert: Cannot insert an mpi_node with a list into which it should be inserted.",mmpbsa::NULL_POINTER);

	node_list = node_list_addr[0];
	if(node_list == 0)
	{
		*node_list_addr = new_node;
		return;
	}

	//If the new node should go at the head of the queue.
	if(node_list->index >= new_node->index)
	{
		next_node = new_node;
		while(next_node->next != 0)
			next_node = next_node->next;
		next_node->next = node_list;
		*node_list_addr = new_node;
		return;
	}

	//If the new node should not go at the head of the queue.
	while(node_list != 0)
	{
		next_node = node_list->next;
		if(next_node == 0)
		{
			node_list->next = new_node;
			return;
		}
	}

}

size_t mmpbsa_utils::mpi_snap_per_host(const size_t& total_snapshots, const int& rank, const int& size)
{
	if(size == 1)
			return total_snapshots;

	if(size == 0)
		throw mmpbsa::MMPBSAExceptions("mmpbsa_utils::mpi_snap_per_host: Requested number of snap shots for zero hosts. Cannot divide by zero.",mmpbsa::DATA_FORMAT_ERROR);

	//Master node (rank = 0) should do the fewest snapshots since it has more disk overhead.
	if(rank == 0)
		if(snapshots % size == 0)
			return mpi_snap_per_host(total_snapshots,rank,size-1);
		else
			return snapshots % size;

	return ceil(total_snapshots/size);

}

#endif
