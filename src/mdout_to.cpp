#include "libmmpbsa/XMLNode.h"
#include "libmmpbsa/EnergyInfo.h"
#include "libmmpbsa/EMap.h" // to use the same XMLNode tag as the mmpbsa energy data
#include "libmmpbsa/mmpbsa_exceptions.h"

#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <string>
#include <cerrno>
#include <sstream>

mmpbsa_utils::XMLNode* EnergyInfo2XML(const mmpbsa::EnergyInfo& data)
{
  using mmpbsa_utils::XMLNode;
  using mmpbsa::str_energy_type;
  
  std::ostringstream buff;
  size_t term_idx,energyinfo_size = data.size();
  XMLNode *energy_term = NULL, *returnMe;
  returnMe = new XMLNode(mmpbsa::EMap::DEFAULT_XML_TAG);
  
  for(term_idx = 0;term_idx < energyinfo_size;term_idx++)
    {
      if(str_energy_type(term_idx) == "UNKNOWN")
	continue;
      if(energy_term == NULL)
	{
	  energy_term = new XMLNode(str_energy_type(term_idx));
	  returnMe->children = energy_term;
	}
      else
	{
	  energy_term->siblings = new XMLNode(str_energy_type(term_idx));
	  energy_term = energy_term->siblings;
	}
      buff.clear();buff.str("");
      buff << data[term_idx];
      energy_term->setText(buff.str());
    }
  
  return returnMe;
}

void output_xml(std::vector<mmpbsa::EnergyInfo>& snapshots, mmpbsa::AveRmsEnerInfo& avg_rms)
{
  using mmpbsa_utils::XMLNode;
  std::vector<mmpbsa::EnergyInfo>::iterator curr_snap;
  std::vector<mmpbsa::EnergyInfo>::difference_type loc;
  std::ostringstream buff;
  XMLNode *head, *xmlsnap, *xmlavg;
 
  xmlsnap = NULL;
  xmlavg = NULL;
  head = new XMLNode("mdout");
  
  curr_snap = snapshots.begin();
  for(;curr_snap != snapshots.end();curr_snap++)
    {
      if(xmlsnap == NULL)
	{
	  xmlsnap = EnergyInfo2XML(*curr_snap);
	  head->children = xmlsnap;
	}
      else
	{
	  xmlsnap->siblings = EnergyInfo2XML(*curr_snap);
	  xmlsnap = xmlsnap->siblings;
	}
      buff.clear();buff.str("");
      loc = std::distance(snapshots.begin(),curr_snap);
      buff << (loc + 1);
      xmlsnap->insertChild("snapshot_number",buff.str());
    }

  xmlavg = EnergyInfo2XML(avg_rms.get_average());
  xmlavg->insertChild("snapshot_number","average");
  xmlavg->siblings = EnergyInfo2XML(avg_rms.get_rms());
  xmlavg->siblings->insertChild("snapshot_number","rms");

  // Attach data to mdout node to create the whole tree
  if(xmlavg != NULL)
    {
      XMLNode *children = head->children;
      while(children->siblings != NULL)
	children = children->siblings;
      children->siblings = xmlavg;
    }
  std::cout << head->toString("") << std::endl;

  delete head;// other nodes will be attached to head and thus deleted here too

}

std::string EnergyInfo2csv(const mmpbsa::EnergyInfo& energy)
{
  size_t idx,size = energy.size();
  std::ostringstream buffer;

  idx = 0;
  for(;idx <size;idx++)
    {
      if(mmpbsa::str_energy_type(idx) == "UNKNOWN")
	continue;
      if(idx != 0)
	buffer << ",";
      buffer << energy[idx];
    }
  return buffer.str();
}

void output_csv(std::vector<mmpbsa::EnergyInfo>& snapshots, mmpbsa::AveRmsEnerInfo& avg_rms)
{
  using mmpbsa_utils::XMLNode;
  using mmpbsa::str_energy_type;
  std::vector<mmpbsa::EnergyInfo>::iterator curr_snap;
  std::vector<mmpbsa::EnergyInfo>::difference_type loc;
  std::ostringstream buff;
  std::string csv_entry;
  size_t energy_term,energy_size = avg_rms.get_average().size();
  
  printf("\"index\",");
  energy_term = 0;
  for(;energy_term < energy_size;energy_term++)
    {
     if(mmpbsa::str_energy_type(energy_term) == "UNKNOWN")
	continue;
     if(energy_term != 0)
       printf(",");
     printf("\"%s\"",str_energy_type(energy_term).c_str());
    }
  printf("\n");
  fflush(stdout);
  
  curr_snap = snapshots.begin();
  for(;curr_snap != snapshots.end();curr_snap++)
    {
      csv_entry = EnergyInfo2csv(*curr_snap);
      loc = (std::distance(snapshots.begin(),curr_snap) + 1);
      printf("%ld,%s\n",(long)loc,csv_entry.c_str());
    }

  csv_entry = EnergyInfo2csv(avg_rms.get_average());
  printf("%ld,%s\n",std::distance(snapshots.begin(),curr_snap),csv_entry.c_str());
  csv_entry = EnergyInfo2csv(avg_rms.get_rms());
  printf("%ld,%s\n",std::distance(snapshots.begin(),curr_snap),csv_entry.c_str());
  fflush(stdout);
}

void load_input(std::vector<mmpbsa::EnergyInfo>& snapshots, mmpbsa::AveRmsEnerInfo& avg_rms,
		std::istream& input_stream)
{
  mmpbsa::EnergyInfo curr,avg,rms;
  int retval = 0,counter = 0;
  
  while(true)
    {
      curr.clear();
      retval = curr.get_next_energyinfo(input_stream);
      if(retval == 0)
	{
	  snapshots.push_back(curr);
	  if(input_stream.eof())
	    return;// There's no Average/RMS data. Hopefully this is a minimization mdout file...
	  counter++;
	  continue;
	}
      else if(retval == mmpbsa::UNEXPECTED_EOF && input_stream.eof())
	throw mmpbsa::MMPBSAException("Unexpected end of file reached for mdout file.",mmpbsa::UNEXPECTED_EOF);
      else if(retval == mmpbsa::UNEXPECTED_EOF)
	break;// in this case we reached the "A V E R A G E" line
      else
	{
	  std::ostringstream error;
	  error <<  "Problem loading snapshot of mdout file. Error code: " << retval;
	  throw mmpbsa::MMPBSAException(error);
	}
    }

  avg.get_next_energyinfo(input_stream);
  rms.get_next_energyinfo(input_stream);
  
  avg_rms = mmpbsa::AveRmsEnerInfo(avg,rms);
}

void print_help()
{
  printf("mdout_to -- Sander mdout file format converter.\n");
  printf("Reads an mdout file and writes the data using the specified format\n");
  printf("to standard output.\n");
  printf("Usage: mdout_to [options] <new format> <input file>\n\n");
  printf("Output formats: xml csv\n\n");
  printf("Options:\n");
  printf("--help, -h  \t This help message.\n");
  printf("--usage, -u \t Same as \"--help\".\n");
  printf("--version, v\t Prints the version information.\n");
	 
}

struct option long_opts[] = {
  {"help",0,NULL,'h'},
  {"usage",0,NULL,'h'},
  {"version",0,NULL,'v'},
  {NULL,0,NULL,0}
};
static const char short_opts[] = "h";


int main(int argc, char **argv)
{
  using std::vector;
  using mmpbsa::EnergyInfo;

  int opt_flag;
  const char *output_type = NULL;
  const char *input_filename = NULL;
  std::ifstream input_file;
  vector<EnergyInfo> snapshots;
  mmpbsa::AveRmsEnerInfo avg_rms;
  if(argc == 1)
    {
      print_help();
      return 0;
    }

  while((opt_flag = getopt_long(argc,argv,short_opts,long_opts,NULL)) != -1)
    {
      switch(opt_flag)
	{
	case 'h':
	  print_help();
	  exit(0);
	  break;
	default:
	  fprintf(stderr,"Unknown flag: %c",opt_flag);
	  if(optarg != NULL)
	    fprintf(stderr,"(%s)",optarg);
	  fprintf(stderr,"\nFor help, try %s --help\n",argv[0]);
	  exit(-1);
	  break;
	}
    }

  if(optind + 1 >= argc)
    {
      fprintf(stderr,"An output file type and input file are required.");
      fprintf(stderr,"\nFor help, try %s --help\n",argv[0]);
      exit(-1);
    }

  output_type = argv[optind++];
  input_filename = argv[optind++];

  if(access(input_filename,R_OK) == -1)
    {
      fprintf(stderr,"Could not read from %s\n",input_filename);
      fprintf(stderr,"Reason: %s\n",strerror(errno));
      exit(errno);
    }

  input_file.open(input_filename);

  try
    {
      load_input(snapshots,avg_rms,input_file);
  
      if(strncmp(output_type,"xml",3) == 0)
	output_xml(snapshots,avg_rms);
      else if(strncmp(output_type,"csv",3) == 0)
	output_csv(snapshots,avg_rms);
      else
	{
	  fprintf(stderr,"Unknown output type: %s\n",output_type);
	  exit(-1);
	}
    }
  catch(mmpbsa::MMPBSAException mmpbsae)
    {
      std::cerr << argv[0] << " encountered a problem. Error Message:" << std::endl;
      std::cerr << mmpbsae << std::endl;
      exit(mmpbsae.getErrType());
    }
  return 0;
}
