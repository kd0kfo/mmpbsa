#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>


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


int main(int argc, char** argv)
{
	mmpbsa_analyzer_arguments args;
	int getopt_retval,option_index;
	struct option * long_opts;
	mmpbsa_analyzer_options curr_ma_opt;
	char* ch_command;
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

	std::cout << "command: " << args.command << " parameter: " << args.parameter << std::endl;
	return 0;

}
