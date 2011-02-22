#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>

#include "libmmpbsa/XMLParser.h"
#include "libmmpbsa/EMap.h"

static char mmpbsa_analyzer_doc[] = "mmpbsa_analyzer -- Program for analyzing and manipulating data produced by mmpbsa.";
static char mmpbsa_analyzer_usage[] = "Usage: mmpbsa_analyzer [options]";
void args_usage();

enum ANALYZER_COMMANDS{SUMMARIZE,NUM_COMMANDS};
enum ERR_CODES{INVALID_INPUT_PARAMETER,INVALID_CLI_ARGUMENT};
typedef struct{
	struct option getoptions;/* getopt.h option type. */
	std::string help_string,arg_type;
}mmpbsa_analyzer_options;

typedef struct {
	int verbosity;
	int command;
	std::string parameter,output_filename;
}mmpbsa_analyzer_arguments ;

/**
 * Possible options provided to the command line.
 */
const mmpbsa_analyzer_options options[] = {
		{{"summarize",required_argument,0,'s'},"Summarizes the data by averaging over all snapshots/frames.","FILE"},
		{{"verbose",optional_argument,0,'v'},"Sets the verbosity of the program. Higher the level, the more verbose.","INTEGER"},
		{{"output",required_argument,0,'o'},"Send output to the specified file. DEFAULT: standard output.","FILE"},
		{{"help",no_argument,0,'h'},"This help dialog",""},
		{{"usage",no_argument,0,0},"Same as --help",""},
		{{0,0,0,0},"",""}
};

void parse_opt(const int& key, mmpbsa_analyzer_arguments& args, char** argv)
{
	std::istringstream buff;
	switch(key)
	{
	case 's':
		args.command = SUMMARIZE;
		args.parameter = optarg;
		break;
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
	args.command = NUM_COMMANDS;
	args.parameter = args.output_filename = "";
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

void summarize_delta(std::ostream* output,mmpbsa::EMap* delta, mmpbsa_t* gas_energies)
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

	if(args.parameter.size() == 0)
		return;

	const string& input_filename = args.parameter;

	//open output. If a file is specified, use it. Otherwise, use stdout
	ostream* output = &std::cout;
	if(args.output_filename.size() != 0)
	{
		std::fstream* outfile = new std::fstream(args.output_filename.c_str(),std::ios::out);
		if(!outfile->good())
		{
			std::cerr << "Could not read from " << args.output_filename << std::endl;
			return;
		}
		output = outfile;
	}

	bool added_snapshot;
	size_t snapshot_counter = 0;
	XMLParser data;
	EMap curr_delta,curr,ecomplex[2],receptor[2],ligand[2],delta[2];//[0] = mean, [1] = stddev.
	size_t num_complex,num_receptor,num_ligand;
	mmpbsa_t gas_energies[4];//complex, receptor, ligand, delta
	num_complex = num_receptor = num_ligand = 0;
	data.parse(input_filename);

	for(size_t i = 0;i<4;i++)
		gas_energies[i] = 0;

	if(data.getHead() != 0 && data.getHead()->children != 0)
	{
		XMLNode *molecule, *snap = data.getHead()->children;
		for(;snap != 0;snap = snap->siblings)
		{
			EMap curr_delta;
			added_snapshot = false;
			for(molecule = snap->children;molecule != 0;molecule = molecule->siblings)
			{
				if(molecule->getName() == "COMPLEX")
				{
					curr = EMap::loadXML(molecule);
					ecomplex[0] += curr;
					ecomplex[1] += stddev_couple(curr,gas_energies[0]);
					curr_delta += curr;
					num_complex++;
				}
				else if(molecule->getName() == "RECEPTOR")
				{
					curr = EMap::loadXML(molecule);
					receptor[0] += curr;
					receptor[1] += stddev_couple(curr,gas_energies[1]);
					curr_delta -= curr;
					num_receptor++;
				}
				else if(molecule->getName() == "LIGAND")
				{
					curr = EMap::loadXML(molecule);
					ligand[0] += curr;
					ligand[1] += stddev_couple(curr,gas_energies[2]);
					curr_delta -= curr;
					num_ligand++;
				}
			}
			delta[0] += curr_delta;
			delta[1] += stddev_couple(curr_delta,gas_energies[3]);
			if(!added_snapshot)
			{
				snapshot_counter++;
				added_snapshot = true;
			}
		}
	}

	if(num_complex != num_receptor || num_receptor != num_ligand)
		cerr << "Warning: Number of components differs." << endl
			<< "Receptor: " << num_receptor << endl
			<< "Ligand: " << num_ligand << endl
			<< "Complex: " << num_complex << endl;

	for(size_t i = 0;i<4;i++)
		gas_energies[i] /= num_complex;

	mmpbsa_t useless;
	ecomplex[0] /= num_complex;ecomplex[1] /= num_complex;ecomplex[1] = sqrt(ecomplex[1] - stddev_couple(ecomplex[0],useless));
	receptor[0] /= num_complex;receptor[1] /= num_complex;receptor[1] = sqrt(receptor[1] - stddev_couple(receptor[0],useless));
	ligand[0] /= num_complex;ligand[1] /= num_complex;ligand[1] = sqrt(ligand[1] - stddev_couple(ligand[0],useless));
	delta[0] /= num_complex;delta[1] /= num_complex;delta[1] = sqrt(delta[1] - stddev_couple(delta[0],useless));

	*output << "Summary of " << snapshot_counter << " snapshots." << std::endl;

	summarize_molecules(output,ecomplex,receptor,ligand,gas_energies);
	summarize_delta(output,delta,&gas_energies[3]);

	if(output != &std::cout)
		delete output;

}


int main(int argc, char** argv)
{
	mmpbsa_analyzer_arguments args;
	int getopt_retval,option_index;
	struct option * long_opts;
	mmpbsa_analyzer_options curr_ma_opt;
	size_t num_opts = 0;

	if(argc < 2)
		args_usage();
	mmpbsa_analyzer_defaults(args);

	curr_ma_opt = options[num_opts];
	while(curr_ma_opt.getoptions.name != 0)
		curr_ma_opt = options[++num_opts];
	num_opts++;//for the null opt at the end.
	long_opts = new struct option[num_opts];
	for(size_t i = 0;i<num_opts;i++)
		long_opts[i] = options[i].getoptions;

	while(true)
	{
		getopt_retval = getopt_long(argc,argv,"s:v::",long_opts,&option_index);
		if(getopt_retval == -1)
			break;

		parse_opt(getopt_retval,args,argv);
	}

	delete [] long_opts;

	if(args.command == SUMMARIZE)
		summarize(args);
	return 0;

}
