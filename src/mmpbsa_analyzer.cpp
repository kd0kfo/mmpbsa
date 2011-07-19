#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>

#include <cerrno>
#include <cstring>

#include "libmmpbsa/XMLParser.h"
#include "libmmpbsa/EMap.h"

static char mmpbsa_analyzer_doc[] = "mmpbsa_analyzer -- Program for analyzing and manipulating data produced by mmpbsa.";
static char mmpbsa_analyzer_usage[] = "Usage: mmpbsa_analyzer [options] [input file]\n\tIf no input file is supplied, standard input is used.\n\tOutput is sent to standard output.";
void args_usage();

enum NON_CHAR_OPTIONS{VERSION_OPT = 1};
enum ERR_CODES{INVALID_INPUT_PARAMETER,INVALID_CLI_ARGUMENT};
typedef struct{
	struct option getoptions;/* getopt.h option type. */
	std::string help_string,arg_type;
}mmpbsa_analyzer_options;

typedef struct {
	int verbosity;
	std::istream *input;
	std::ostream *output;
}mmpbsa_analyzer_arguments ;

enum MOLECULE_ENUM{COMPLEX=0,RECEPTOR,LIGAND,DELTA,NUM_MOLECULES};

typedef struct{
  mmpbsa_t ele, vdw, internal, gas, pbsur, pbsolv, area; 
}mmpbsa_analyzer_data;

typedef struct{
  size_t molsurf_dependent[NUM_MOLECULES], molsurf_independent[NUM_MOLECULES];
  mmpbsa_analyzer_data averages[NUM_MOLECULES],stddev[NUM_MOLECULES];
}mmpbsa_analyzer_average;


/**
 * Possible options provided to the command line.
 */
const mmpbsa_analyzer_options options[] = {
		{{"verbose",optional_argument,0,'v'},"Sets the verbosity of the program. Higher the level, the more verbose.","INTEGER"},
		{{"help",no_argument,0,'h'},"This help dialog",""},
		{{"version",no_argument,0,VERSION_OPT},"Display version",""},
		{{"usage",no_argument,0,0},"Same as --help",""},
		{{0,0,0,0},"",""}
};

void parse_opt(const int& key, mmpbsa_analyzer_arguments& args, char** argv)
{
	std::istringstream buff;
	switch(key)
	{
	case 'v':
	  if(optarg == 0)
	    {
	      args.verbosity = 1;
	      break;
	    }
		buff.str(optarg);
		buff >> args.verbosity;
		if(buff.fail())
		{
		  std::cerr << "Warning: " << argv[optind] << " is an invalid verbose level. Need integer, but got: " << optarg << std::endl;
			exit(INVALID_INPUT_PARAMETER);
		}
		break;
	case VERSION_OPT:
	  std::cout << "mmpbsa_analyzer " << PACKAGE_VERSION << std::endl;
	  exit(0);
	case 0:case '?':case 'h':
	  args_usage();
	  break;
	default:
		std::cerr << "Unknown flag: " << (char)key << std::endl;
		std::cerr << "run mmpbsa_analyzer --help for more information" << std::endl;
		exit(1);
	}
}

void mmpbsa_analyzer_defaults(mmpbsa_analyzer_arguments& args)
{
	args.input = &std::cin;
	args.output = &std::cout;
	args.verbosity = 0;
}

void args_usage()
{
	using namespace std;
	mmpbsa_analyzer_options curr_opts;
	size_t opt_counter;

	cout << mmpbsa_analyzer_doc << endl;
	cout << mmpbsa_analyzer_usage << endl;
	cout << endl;
	opt_counter = 0;
	curr_opts = options[opt_counter];
	while(curr_opts.getoptions.name != 0)
	{
		if(curr_opts.getoptions.val > 0x41)
			cout << "-" << (char)curr_opts.getoptions.val << ", ";
		else
			cout << "    ";
		cout << "--" << curr_opts.getoptions.name;
		if(curr_opts.arg_type.size() != 0)
		{
			cout << "=";
			if(curr_opts.getoptions.has_arg == optional_argument)
				cout << "[";
			cout << curr_opts.arg_type;
			if(curr_opts.getoptions.has_arg == optional_argument)
				cout << "]";
		}
		cout << "\t" << curr_opts.help_string << endl;
		curr_opts = options[++opt_counter];
	}
	exit(0);
}

void finalize_average(mmpbsa_analyzer_average& avg)
{
  for(size_t i = 0;i<NUM_MOLECULES;i++)
    {

    }
  
}

void average_data(const mmpbsa::EMap& new_emap, mmpbsa_analyzer_average& avg, mmpbsa::EMap& curr_delta,int molecule)
{
  mmpbsa_t curr_data;
  size_t i;//index into arrays 0 <= i < NUM_MOLECULES
  if(molecule < 0 || molecule >= NUM_MOLECULES)
    throw mmpbsa::MMPBSAException("average_data: invalid molecule enum value.",mmpbsa::DATA_FORMAT_ERROR);

  i = (size_t) molecule;
  curr_data = new_emap.total_elec_energy();
  avg.averages[i].ele += curr_data;
  avg.stddev[i].ele += curr_data*curr_data;
  
  curr_data = new_emap.total_vdw_energy();
  avg.averages[i].vdw += curr_data; 
  avg.stddev[i].vdw += curr_data*curr_data; 
  
  curr_data = new_emap.total_internal_energy();
  avg.averages[i].internal += curr_data;
  avg.stddev[i].internal += curr_data*curr_data;
  
  curr_data = new_emap.total_gas_energy();
  avg.averages[i].gas += curr_data;
  avg.stddev[i].gas += curr_data*curr_data;
  
  curr_data = new_emap.elstat_solv;
  avg.averages[i].pbsolv += curr_data;
  avg.stddev[i].pbsolv += curr_data*curr_data;

  avg.molsurf_independent[i]++;

  if(new_emap.molsurf_failed)
    {
      curr_data = new_emap.sasol;
      avg.averages[i].pbsur += curr_data;
      avg.stddev[i].pbsur += curr_data*curr_data;
      
      curr_data = new_emap.area;
      avg.averages[i].area += curr_data;
      avg.stddev[i].area += curr_data*curr_data;
      avg.molsurf_dependent[i]++;
    }

}

void sanity_check_average(mmpbsa_analyzer_average& avg)
{
	if(num_complex != num_receptor || num_receptor != num_ligand)
		cerr << "Warning: Number of components differs." << endl
			<< "Receptor: " << num_receptor << endl
			<< "Ligand: " << num_ligand << endl
			<< "Complex: " << num_complex << endl;

}

void init_average(mmpbsa_analyzer_average& avg)
{
  for(size_t i = 0;i<NUM_MOLECULES;i++)
    {
      molsurf_dependent[i] = molsurf_independent[i] = 0;
      
      avg.averages[i].ele = avg.stddev[i].ele = 0.0;
      avg.averages[i].vdw = avg.stddev[i].vdw = 0.0;
      avg.averages[i].internal = avg.stddev[i].internal = 0.0;
      avg.averages[i].gas = avg.stddev[i].gas = 0.0;
      avg.averages[i].pbsur = avg.stddev[i].pbsur = 0.0;
      avg.averages[i].pbsolv = avg.stddev[i].pbsolv = 0.0;
      avg.averages[i].area = avg.stddev[i].area = 0.0;
    }
}

void summarize_delta(std::ostream* output, mpbsa::EMap* delta, mmpbsa_t* gas_energies)
{
	using std::ios;
	if(output == 0)
		return;

	output->precision(2);
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << " ";output->flags(ios::internal|ios::fixed);
	output->width(18);
	*output << "DELTA" << std::endl;

	//ELE = vacele + ele14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "ELE";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].total_elec_energy();output->width(12);*output << delta[1].total_elec_energy() << std::endl;

	//VDW = vdwaals + vdw14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "VDW";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].total_vdw_energy();output->width(12);*output << delta[1].total_vdw_energy() << std::endl;

	//INT = vdwaals + vdw14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "INT";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].total_internal_energy();output->width(12);*output << delta[1].total_internal_energy() << std::endl;

	//Gas energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "GAS";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].total_gas_energy();output->width(12);*output << sqrt(*gas_energies - pow(delta[0].total_gas_energy(),2)) << std::endl;

	//Surface Area Solvation energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "PBSUR";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].sasol;output->width(12);*output << delta[1].sasol << std::endl;

	//Poisson-Boltzmann Solvation energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "PBSOLV";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].elstat_solv;output->width(12);*output << delta[1].elstat_solv << std::endl;

	//Surface Area
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "AREA";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << delta[0].area;output->width(12);*output << delta[1].area << std::endl;


}


void summarize_molecules(std::ostream* output,mmpbsa::EMap* complex,mmpbsa::EMap* receptor,mmpbsa::EMap* ligand,mmpbsa_t* gas_energies)
{
	using std::ios;
	if(output == 0)
		return;

	output->precision(2);
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << " ";output->flags(ios::internal|ios::fixed);
	output->width(18);
	*output << "COMPLEX";output->width(24);
	*output << "RECEPTOR";output->width(24);
	*output << "LIGAND" << std::endl;

	//ELE = vacele + ele14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "ELE";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].total_elec_energy();output->width(12);*output << complex[1].total_elec_energy();output->width(12);
	*output << receptor[0].total_elec_energy();output->width(12);*output << receptor[1].total_elec_energy();output->width(12);
	*output << ligand[0].total_elec_energy();output->width(12);*output << ligand[1].total_elec_energy() << std::endl;

	//VDW = vdwaals + vdw14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "VDW";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].total_vdw_energy();output->width(12);*output << complex[1].total_vdw_energy();output->width(12);
	*output << receptor[0].total_vdw_energy();output->width(12);*output << receptor[1].total_vdw_energy();output->width(12);
	*output << ligand[0].total_vdw_energy();output->width(12);*output << ligand[1].total_vdw_energy() << std::endl;

	//INT = vdwaals + vdw14
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "INT";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].total_internal_energy();output->width(12);*output << complex[1].total_internal_energy();output->width(12);
	*output << receptor[0].total_internal_energy();output->width(12);*output << receptor[1].total_internal_energy();output->width(12);
	*output << ligand[0].total_internal_energy();output->width(12);*output << ligand[1].total_internal_energy() << std::endl;

	//Gas energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "GAS";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].total_gas_energy();output->width(12);*output << sqrt(gas_energies[0]-pow(complex[0].total_gas_energy(),2));output->width(12);
	*output << receptor[0].total_gas_energy();output->width(12);*output << sqrt(gas_energies[1]-pow(receptor[0].total_gas_energy(),2));output->width(12);
	*output << ligand[0].total_gas_energy();output->width(12);*output << sqrt(gas_energies[2]-pow(ligand[0].total_gas_energy(),2)) << std::endl;

	//Surface Area Solvation energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "PBSUR";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].sasol;output->width(12);*output << complex[1].sasol;output->width(12);
	*output << receptor[0].sasol;output->width(12);*output << receptor[1].sasol;output->width(12);
	*output << ligand[0].sasol;output->width(12);*output << ligand[1].sasol << std::endl;

	//Poisson-Boltzmann Solvation energy
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "PBSOLV";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].elstat_solv;output->width(12);*output << complex[1].elstat_solv;output->width(12);
	*output << receptor[0].elstat_solv;output->width(12);*output << receptor[1].elstat_solv;output->width(12);
	*output << ligand[0].elstat_solv;output->width(12);*output << ligand[1].elstat_solv << std::endl;

	//Surface Area
	output->flags(ios::left|ios::fixed);output->width(12);
	*output << "AREA";output->flags(ios::right|ios::fixed);
	output->width(12);
	*output << complex[0].area;output->width(12);*output << complex[1].area;output->width(12);
	*output << receptor[0].area;output->width(12);*output << receptor[1].area;output->width(12);
	*output << ligand[0].area;output->width(12);*output << ligand[1].area << std::endl;


}

mmpbsa::EMap stddev_couple(const mmpbsa::EMap& energy, mmpbsa_t& gas)
{
	mmpbsa::EMap coupled_term = energy;
	coupled_term.ele14 += energy.vacele; coupled_term.vacele = 0;
	coupled_term.vdw14 += energy.vdwaals;coupled_term.vdwaals = 0;
    coupled_term.angle += energy.bond + energy.dihed;coupled_term.bond = coupled_term.dihed = 0;
    gas += pow(coupled_term.total_gas_energy(),2);
    return coupled_term*coupled_term;
}

void summarize(mmpbsa_analyzer_arguments& args)
{
	using namespace std;
	using mmpbsa_utils::XMLParser;
	using mmpbsa_utils::XMLNode;
	using mmpbsa::EMap;

	XMLNode *data = NULL;

	//some data does depend on molsurf results; others do not. Therefore the averaging is different if molsurf fails.
	mmpbsa_analyzer_average avg;
	size_t bad_area_com, bad_area_rec, bad_area_lig;
	
	EMap curr_delta;//place holder for calculating the DELTA energy for a given snapshot.

	if(args.input == NULL)
		throw mmpbsa::MMPBSAException("summarize: No input stream provided.",mmpbsa::FILE_IO_ERROR);
	if(args.output == NULL)
		throw mmpbsa::MMPBSAException("summarize: No output stream provided.",mmpbsa::FILE_IO_ERROR);

	bad_area_com = bad_area_rec = bad_area_lig = 0;
	
	*data = XMLParser::parse(*args.input);
	if(data == NULL)
		throw mmpbsa::MMPBSAException("summarize: Could not parse data from input stream.",mmpbsa::FILE_IO_ERROR);

	init_average(avg);

	if(data->children != 0)
	{
		XMLNode *molecule, *snap = data->children;
		for(;snap != 0;snap = snap->siblings)
		{
		  curr_delta.clear();
			for(molecule = snap->children;molecule != 0;molecule = molecule->siblings)
			{
				if(molecule->getName() == "COMPLEX")
				{
					curr = EMap::loadXML(molecule);
					average_data(curr,avg,curr_delta,COMPLEX);
					if(curr.molsurf_failed)
					  bad_area_com++;
				}
				else if(molecule->getName() == "RECEPTOR")
				{
					curr = EMap::loadXML(molecule);
					average_data(curr,avg,curr_delta,RECEPTOR);
					if(curr.molsurf_failed)
					  bad_area_rec++;
				}
				else if(molecule->getName() == "LIGAND")
				{
					curr = EMap::loadXML(molecule);
					average_data(curr,avg,curr_delta,LIGAND);
					if(curr.molsurf_failed)
						bad_area_lig++;
				}
				else if(args.verbosity)
				  {
				    *args.output << "Ignoring tag: " << molecule->getName() << std::endl;				
				  }
			}
			average_data(curr_delta,avg,curr_delta,DELTA);
		}
	}

	// Verify that the data collected makes sense
	sanity_check_average(avg);

	//Finish average and std. dev. calculations
	finalize_average(avg);
	
	// Display failed data information
	if(args.verbosity)
	  {
	    if(bad_area_com > 0)
	      {
		fprintf(stderr,"Data contains %d occurances where molsurf failed in complex.\n",bad_area_com);
	      }
	    if(bad_area_lig > 0)
	      {
		fprintf(stderr,"Data contains %d occurances where molsurf failed in ligand.\n",bad_area_lig);
	      }
	    if(bad_area_rec > 0)
	      {
		fprintf(stderr,"Data contains %d occurances where molsurf failed in receptor.\n",bad_area_rec);
	      }
	  }

	// Output results
	*args.output << "Summary of " << snapshot_counter << " snapshots." << std::endl;

	summarize_molecules(args.output,ecomplex,receptor,ligand,gas_energies);
	summarize_delta(args.output,delta,&gas_energies[3]);

	delete data;

}


int main(int argc, char** argv)
{
	mmpbsa_analyzer_arguments args;
	int getopt_retval,option_index;
	struct option * long_opts;
	mmpbsa_analyzer_options curr_ma_opt;
	size_t num_opts = 0;

	mmpbsa_analyzer_defaults(args);

	curr_ma_opt = options[num_opts];
	while(curr_ma_opt.getoptions.name != 0)
		curr_ma_opt = options[++num_opts];
	num_opts++;//for the null opt at the end.
	long_opts = new struct option[num_opts];
	for(size_t i = 0;i<num_opts;i++)
		long_opts[i] = options[i].getoptions;
	
	while((getopt_retval = getopt_long(argc,argv,"v::",long_opts,&option_index)) != -1)
		parse_opt(getopt_retval,args,argv);

	if(optind < argc)
	{
		std::fstream *input = new std::fstream;
		if(access(argv[optind],R_OK))
		{
			fprintf(stderr,"Cannot read (%d)%s\nReason: %s",optind,argv[optind],strerror(errno));
			return errno;
		}

		input->open(argv[optind],std::ios::in);
		if(!input->good())
		{
			std::ostringstream err;
			delete input;
			fprintf(stderr, "Could not open %s for reading.",argv[optind]);
			return 42;
		}
		args.input = input;
	}

	delete [] long_opts;

	summarize(args);
	
	if(args.input != &std::cin)
		delete args.input;
	return 0;

}
