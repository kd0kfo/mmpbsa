#include "mmpbsa_mpi.h"

void init(mmpbsa_utils::mpi_node* node)
{
	if(node == 0)
		throw mmpbsa::MMPBSAException("init: mpi_node to be initialized is a null pointer.",mmpbsa::NULL_POINTER);
	node->data = 0;
	node->index= 0;
}

void destroy(mmpbsa_utils::mpi_node* node)
{
	if(node == 0)
		return;
	if(node->data != 0)
	{
		free(node->data);
		node->data = 0;
	}
	node->index = 0;
}

void mmpbsa_utils::mpi_init_hosts(int* argc, char*** argv, int& mpi_rank,int& mpi_size)
{
	//Initialize
	MPI_Init(argc,argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
}

void mmpbsa_utils::mpi_store_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data,const int& data_index,
		mmpbsa_utils::mpi_data_list* data_list)
{
	mmpbsa_utils::XMLNode *id_node,*it;
	int index = 0;

	if(data_list == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store energy data into an empty list (which is a null pointer).",mmpbsa::NULL_POINTER);

	it = energy_data.getHead();
	if(it == 0 || it->children == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store empty energy data." ,mmpbsa::DATA_FORMAT_ERROR);

#if 0//probably don't want this
	it = it->children;//this should be the snapshot tag

	//iterate through the snapshot's children tags. Only 1st snapshot is stored.
	//Call this function for each snapshot
	id_node = 0;
	for(it = it->children;it != 0;it = it->siblings)
	{
		if(it->getName() == "ID" || it->getName() == "id")
		{
			id_node = it;
			break;
		}
	}

	//If an index is listed for the snapshot, set index equal to its value
	if(id_node != 0)
	{
		std::istringstream buff(id_node->getText());
		buff >> index;
		if(buff.fail())
			throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Invalid index in snapshot data: " + id_node->getText(),mmpbsa::DATA_FORMAT_ERROR);
	}
#endif

	//Create char* and store data
	it = energy_data.getHead()->children;
	std::string str_data = it->toString();
	char* data = (char*)calloc(str_data.size() + 1,sizeof(char));
	strcpy(data,str_data.c_str());
	mpi_store_mmpbsa_data(data,data_index,data_list);
}

void mmpbsa_utils::mpi_store_mmpbsa_data(char* data,const int& data_index,
		mmpbsa_utils::mpi_data_list* data_list)
{
	mmpbsa_utils::mpi_node new_data;

	if(data_list == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_store_mmpbsa_data: Cannot store in an empty data list, which is a null pointer.",mmpbsa::NULL_POINTER);
	new_data.data = data;
	new_data.index = data_index;

	data_list->insert(new_data);
}

/**
 * Returns 1 if the head of the data_list was written. Zero otherwise.
 */
int mmpbsa_utils::mpi_update_data_file(const mmpbsa::MMPBSAState& currState,
		const int& mpi_rank, mmpbsa_utils::mpi_data_list* data_list, mmpbsa_utils::mpi_node* next_node)
{
	//Only master node writes files.
	if(mpi_rank != MMPBSA_MASTER)
		return 0;

	printf("data_list size: %d \n next_node index: %d \n",data_list->size(),next_node->index);
	if(data_list->size() > 0)
		printf("first data_list index: %d\n",data_list->begin()->index);

	if(data_list == 0 || data_list->size() == 0)
		return 0;
	if(next_node == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_update_data_file: Must have a place to store the last node written. Null pointer was provided.",mmpbsa::NULL_POINTER);

	if(next_node->index != data_list->begin()->index)
		return 0;

#if 0//probably don't want this
	//Decide if the head of the data_list should be written
	//If not, leave.
	if(currState.snapList.size() == 0)
		if(data_list->begin()->index != next_node->index + 1)
			return 0;

	//If index of the head of the data_list is after the index last node node,
	//write the head of data_list
	for(size_t i = 0;i<currState.snapList.size();i++)
		if(currState.snapList.at(i) == next_node->index)
			if(currState.snapList.size() == i+1 || currState.snapList.at(i+1) != data_list->begin()->index)
			{
				printf("i: %d size: %d i+1: %d first index: %d\n",i,currState.snapList.size(),currState.snapList.at(i+1),data_list->begin()->index);
				return 0;
			}
#endif

	//Open file
	std::fstream output;
	if(has_filename(MMPBSA_OUT_TYPE,currState))
		output.open(get_filename(MMPBSA_OUT_TYPE,currState).c_str(),std::ios::out | std::ios::app);
	else if(has_filename(SANDER_MDOUT_TYPE,currState))
		output.open(get_filename(SANDER_MDOUT_TYPE,currState).c_str(),std::ios::out | std::ios::app);
	else
		output.open("mmpbsa_mpi_output.xml",std::ios::out | std::ios::app);

	if(!output.good())
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_update_data_file: Could not open file for writting.",mmpbsa::FILE_IO_ERROR);

	//Write data
	if(data_list->begin()->index == 0)
		output << "<" << MMPBSA_XML_TITLE << ">" << std::endl;
	if(data_list->begin()->data != 0)
		output << data_list->begin()->data << std::endl;
	output.close();

	//Update list
	next_node->index++;
	if(data_list->begin()->data != 0)
		free(data_list->begin()->data);
	data_list->erase(data_list->begin());

	return 1;
}

int mmpbsa_utils::mpi_send_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const int& mpi_rank)
{
	//Create char* and store data
	if(energy_data.getHead() == 0)
		return 0;

	const mmpbsa_utils::XMLNode *it = energy_data.getHead()->children;
	std::string str_data = it->toString();
	size_t data_length = str_data.size()+1 + MMPBSA_MPI_ID_WORD;
	char* data = (char*)calloc(data_length,sizeof(char));//add 1 to str_data.size() for the null character at the end and then add to the beginning for mpi_rank of the sender.
	strcpy(&data[MMPBSA_MPI_ID_WORD],str_data.c_str());

	//tag data with mpi rank as little endian char's.
	for(size_t i = 0;i<MMPBSA_MPI_ID_WORD;i++)
	{
		data[i] = (char)((mpi_rank & (0xff << 2*i)) >> 2*i);
	}

	printf("%d is sending data (%d) to %d.\n",mpi_rank,data,MMPBSA_MASTER);
	int returnMe = MPI_Send(data, data_length*sizeof(char), MPI_CHAR,MMPBSA_MASTER,0, MPI_COMM_WORLD);
	free(data);
	printf("%d finished sending data (%d) to %d.\n",mpi_rank,data,MMPBSA_MASTER);
	return returnMe;
}

int mmpbsa_utils::mpi_recv_mmpbsa_data(const int& my_rank, const int& source_rank,
		const int& mpi_size, const mmpbsa::MMPBSAState& currState,
		mmpbsa_utils::mpi_data_list* data_list, mmpbsa_utils::mpi_node* next_node)
{
	if(mpi_size == 0 || my_rank != 0)
		return 0;

	char *write_data,*mpi_data = (char*)calloc(MMPBSA_MPI_MAX_BUFFER,sizeof(char));
	int sender_rank = 0;
	size_t write_data_size;
	MPI_Status status;

	//Receive Data from node
	int returnMe = MPI_Recv(mpi_data, MMPBSA_MPI_MAX_BUFFER, MPI_CHAR,
			source_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

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

	//Extract source rank.
	for(size_t i = 0;i<MMPBSA_MPI_ID_WORD;i++)
	{
		sender_rank += ((int)mpi_data[i]) << 2*i;
	}

	//Extract text to be written
	write_data_size = strlen(&mpi_data[MMPBSA_MPI_ID_WORD]);
	write_data = (char*)calloc(write_data_size,sizeof(char));
	strcpy(write_data,&mpi_data[MMPBSA_MPI_ID_WORD]);
	free(mpi_data);

	//Store data and write if needed
	mpi_store_mmpbsa_data(write_data,sender_rank,data_list);
	return mmpbsa_utils::mpi_update_data_file(currState,my_rank,data_list,next_node);//will free(write_data)

}

int mmpbsa_utils::mpi_write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState,
		const int& mpi_rank, mmpbsa_utils::mpi_data_list* data_list, mmpbsa_utils::mpi_node* next_node)
{
	if(mpi_rank == MMPBSA_MASTER)
	{
		mpi_store_mmpbsa_data(energy_data,mpi_rank,data_list);
		while(mpi_update_data_file(currState,mpi_rank,data_list,next_node))1;
		return 0;
	}
	else
		mpi_send_mmpbsa_data(energy_data,mpi_rank);
}

void mmpbsa_utils::mpi_finish_output(const int& mpi_rank,const int& mpi_size,const std::string& mmpbsa_output_filename)
{
	if(mpi_rank != 0 || mpi_size == 1)
		return;
	std::fstream output(mmpbsa_output_filename.c_str(),std::ios::out | std::ios::app);

	if(!output.good())
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_update_data_file: Could not open file for writting.",mmpbsa::FILE_IO_ERROR);

	output << "</" << MMPBSA_XML_TITLE << ">" << std::endl;
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
