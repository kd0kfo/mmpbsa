#include "mmpbsa_mpi.h"

#include <mpi.h>

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

void mmpbsa_utils::mpi_store_mmpbsa_data(const char* data,const int& data_index,
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
	//parser.detachHead();
}


int mmpbsa_utils::mpi_send_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const int& mpi_rank)
{
	std::string str_data = energy_data.toString(),substring;
	bool is_EOF;
	int returnMe = 0;
	size_t data_length;
	while(str_data.size() != 0)
	{
		if(str_data.size() >= MMPBSA_MPI_MAX_BUFFER)
		{
			substring = str_data.substr(0,MMPBSA_MPI_MAX_BUFFER);
			str_data.erase(0,MMPBSA_MPI_MAX_BUFFER);
			data_length = MMPBSA_MPI_MAX_BUFFER;
		}
		else
		  {
			substring = str_data;
			str_data = "";
			data_length = substring.size() + 1;//add 1 to str_data.size() for the null character at the end
		}
		char* data = (char*)calloc(data_length,sizeof(char));
		strcpy(data,substring.c_str());

		returnMe = MPI_Send(data, data_length*sizeof(char), MPI_CHAR,MMPBSA_MASTER,mmpbsa_utils::DATA, MPI_COMM_WORLD);
		free(data);
		if(returnMe)
			break;
	}
	return returnMe;
}

int mmpbsa_utils::mpi_recv_mmpbsa_data(const int& my_rank, const int& source_rank,
		const int& mpi_size, const mmpbsa::MMPBSAState& currState,
		mmpbsa_utils::XMLNode* data_list,std::map<int,std::string>& data_fragments)
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

	//If the data has a terminating null character, store data and write if needed.
	//Otherwise, store it for recombination later.
	if(mpi_data[MMPBSA_MPI_MAX_BUFFER-1] != 0)
	{
		if(data_fragments.find(source_rank) != data_fragments.end())
			data_fragments.find(source_rank)->second += mpi_data;
		else
		{
			std::pair<int,std::string> fragment;
			fragment.first = source_rank;
			fragment.second = std::string(mpi_data);
			data_fragments.insert(fragment);
		}
	}
	else
	{
		std::string whole_data = mpi_data;
		if(data_fragments.find(source_rank) != data_fragments.end())
		{
			whole_data = data_fragments.find(source_rank)->second + whole_data;
			data_fragments.erase(source_rank);
		}

		mpi_store_mmpbsa_data(whole_data.c_str(),42,data_list);
	}
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


#if 0//deprecated

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
#endif //deprecated

typedef struct  {
  bool operator() (const std::pair<int,mmpbsa_utils::XMLNode *>& lhs, const std::pair<int,mmpbsa_utils::XMLNode *>& rhs) const
  {return lhs.first<rhs.first;}
}xml_sorter;


void sort_this(std::set<std::pair<int,mmpbsa_utils::XMLNode *>,xml_sorter >& sorted_data,const mmpbsa_utils::XMLNode* data_list)
{
	std::istringstream buff;
	int id;
	std::pair<int,mmpbsa_utils::XMLNode *> new_pair;
	mmpbsa_utils::XMLNode *it, *snap_shot = data_list->children;
	for(;snap_shot != 0;snap_shot = snap_shot->siblings)
		for(it = snap_shot->children;it != 0;it = it->siblings)
			if(it->getName() == "ID")
			{
				buff.clear();
				buff.str(it->getText());
				buff >> id;
				if(!buff.fail())
				{
					new_pair.first = id;
					new_pair.second = snap_shot;
					sorted_data.insert(new_pair);
				}
				break;
			}
}

void mmpbsa_utils::mpi_dump_data(mmpbsa_utils::XMLNode* data_list,const std::string& filename)
{
	using namespace std;
	if(data_list == 0 || filename.size() == 0)
		return;
	fstream output(filename.c_str(),std::ios::out | std::ios::app);
	if(!output.good())
		throw mmpbsa::MMPBSAException("mmpbsa_utils::mpi_dump_data: Could not open file for writing.",mmpbsa::FILE_IO_ERROR);

	set<pair<int,mmpbsa_utils::XMLNode *>,xml_sorter > snaps;
	set<pair<int,mmpbsa_utils::XMLNode *>,xml_sorter >::const_iterator snap_shot;
	sort_this(snaps,data_list);
	//std::vector<mmpbsa_utils::XMLNode *>::const_iterator snap = snaps.begin();
	//for(;snap != snaps.end();snap++)
	mmpbsa_utils::XMLNode *temp;
	output << "<" << MMPBSA_XML_TITLE << ">" << std::endl;
	for(snap_shot = snaps.begin();snap_shot != snaps.end();snap_shot++)
		if(snap_shot->second != 0)
		{
			temp = snap_shot->second->siblings;
			snap_shot->second->siblings = 0;
			output << snap_shot->second->toString() << std::endl;
			snap_shot->second->siblings = temp;
		}
	output << "</" << MMPBSA_XML_TITLE << ">" << std::endl;
	output.close();
}

