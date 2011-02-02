#ifndef MMPBSA_MPI_H
#define MMPBSA_MPI_H

#include <stdlib.h>
#include <set>
#include <sstream>
#include <mpi.h>
#include <stdlib.h>

#include "libmmpbsa/mmpbsa_exceptions.h"
#include "libmmpbsa/MMPBSAState.h"

#define MMPBSA_MASTER 0
#define MMPBSA_MPI_MAX_BUFFER 4096
#define MMPBSA_MPI_ID_WORD 2 //Number of chars used to identify host. Thus, number of MPI ranks is equal to 2^(MMPBSA_MPI_ID_WORD*8)

namespace mmpbsa_utils
{

void mpi_init_hosts(int* argc, char*** argv, int& mpi_rank,int& mpi_size);
void mpi_store_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list);
void mpi_store_mmpbsa_data(char* data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list);

int mpi_send_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const int& mpi_rank);

int mpi_recv_mmpbsa_data(const int& my_rank, const int& source_rank,
		const int& mpi_size, const mmpbsa::MMPBSAState& currState,
		mmpbsa_utils::XMLNode* data_list);

int mpi_write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState,
		const int& mpi_rank, mmpbsa_utils::XMLNode* data_list);
void mpi_dump_data(mmpbsa_utils::XMLNode* data_list,const std::string& filename);
void mpi_finish_output(const int& mpi_rank,const int& mpi_size,const std::string& mmpbsa_output_filename);

}//end namespace mmpbsa_utils


#endif
