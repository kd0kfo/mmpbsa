#ifndef MMPBSA_MPI_H
#define MMPBSA_MPI_H

#include <stdlib.h>
#include <set>
#include <sstream>
#include <stdlib.h>

#include "libmmpbsa/XMLParser.h"
#include "libmmpbsa/mmpbsa_exceptions.h"
#include "libmmpbsa/MMPBSAState.h"

#define MMPBSA_MASTER 0
#define MMPBSA_MPI_MAX_BUFFER 1024//4096

namespace mmpbsa_utils
{

enum MMPBSA_MPI_TAG {DATA = 0,STATUS,NUM_OF_TAGS};

/**
 * Initializes hosts.
 *
 * Calls MPI_Init.
 */
void mpi_init_hosts(int* argc, char*** argv, int& mpi_rank,int& mpi_size);

/**
 * Stores energy data sent from the slave nodes to the master. The provided
 * energy data is inserted into the master data set, data_list
 */
void mpi_store_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list);

/**
 * Stores energy data sent from the slave nodes to the master. The provided
 * energy data is inserted into the master data set, data_list.
 *
 * data is loaded into an XMLParser and then a call is made to
 * void mpi_store_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list)
 */
void mpi_store_mmpbsa_data(const char* data,const int& data_index,
		mmpbsa_utils::XMLNode* data_list);

/**
 * Sends energy data to the master node, whose rank is defined above.
 * The rank provided as a argument is the rank of the slave node
 * sending the data.
 *
 * This is a blocking call, using MPI_Send
 */
int mpi_send_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const int& mpi_rank);

/**
 * This receives data from the slave node, source_rank, an stores it in
 * the XMLNode data_list. If this method is called by any node other
 * than the master, the function immediately returns 0 (and obviously does not block).
 *
 * This is a blocking call, using MPI_Recv
 */
int mpi_recv_mmpbsa_data(const int& my_rank, const int& source_rank,
		const int& mpi_size, const mmpbsa::MMPBSAState& currState,
		mmpbsa_utils::XMLNode* data_list,std::map<int,std::string>& data_fragments);

/**
 * 	Replaces the act of writing mmpbsa data in a manner
 *  dependent on the node's rank.
 *
 * If called by the master node, the energy data is stored
 * in data_list.
 *
 * If called by a slave node, the energy data is send to the master
 * node, using mpi_send_mmpbsa_data.
 *
 * This function should be used when a serial function would write the
 * mmpbsa data, because it abstracts away the fact that all of the
 * nodes are separated from each other. Then when the master node finishes,
 * it can actually perform the file write, using mpi_dump_data.
 */
int mpi_write_mmpbsa_data(mmpbsa_utils::XMLParser& energy_data, const mmpbsa::MMPBSAState& currState,
		const int& mpi_rank, mmpbsa_utils::XMLNode* data_list);

/**
 * Performs a file write using the data_list XMLNode. This should be called
 * by the master node when the calculation has completed, or by any node
 * if something goes wrong and data should be dumped for analysis.
 */
void mpi_dump_data(mmpbsa_utils::XMLNode* data_list,const std::string& filename);


}//end namespace mmpbsa_utils


#endif
