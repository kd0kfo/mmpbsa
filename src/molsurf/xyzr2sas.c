#include "molsurf.h"

#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>

int load_coordinates(REAL_T **xcrds, REAL_T **ycrds, REAL_T **zcrds, REAL_T **radii, long natoms, FILE *input)
{
  char buffer[4096],*tok;
  int test;
  long linecount = 0;
  if(natoms == 0)
    return 0;
  if(xcrds == NULL || ycrds == NULL || zcrds == NULL || radii == NULL || input == NULL)
    {
      fprintf(stderr,"load_coordinates was given a null pointer.\n");
      return -1;
    }

  *xcrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *ycrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *zcrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *radii = (REAL_T*)malloc(natoms*sizeof(REAL_T));

  while(!feof(input) && !ferror(input))
    {
      if(fgets(buffer,4096,input) == NULL)
	break;
      if(strlen(buffer) == 0 || buffer[0] == '#')//comment with # character
	continue;
      if(!isdigit(buffer[0]) && buffer[0] != '+' && buffer[0] != '-')
	continue;

      test = sscanf(buffer,"%lf %lf %lf %lf",&(*xcrds)[linecount],&(*ycrds)[linecount],&(*zcrds)[linecount],&(*radii)[linecount]);
      (*radii)[linecount] += 1.4;
      linecount++;
    }
  if(ferror(input))
    {
      fprintf(stderr,"File error. Reason: %s\n",strerror(errno));
      return -errno;
    }

  return linecount;
}

void print_help()
{
  printf("xyzr2sas -- Molecular Surface Area Calculator\n");
  printf("Usage: ./xyzr2sas [options]\n\n");
  printf("Options:\n");
  printf("--input,   -i <FILE>\tAtomic coordinates and radii (Default: standard input)\n");
  printf("--version, -v       \tDisplays the version\n");
  printf("--help,    -h       \tThis help message\n");
  printf("--usage             \tSame as \"--help\"\n");
}

const struct option long_opts[] = 
  {
    {"help",0,NULL,'h'},
    {"usage",0,NULL,'h'},
    {"version",0,NULL,'v'},
    {"input",1,NULL,'i'},
    {NULL,0,NULL,0}
  };
const char short_opts[] = "hi:v";

int main(int argc, char **argv)
{
  FILE *input = stdin;
  char num_atom_str[4096],curr_opt;
  long natoms;
  REAL_T *crds[3],*radii;
  REAL_T area;
  int retval;

  natoms = 0;
  for(;natoms < 3;natoms++)
    crds[natoms] = NULL;
  radii = NULL;

  while((curr_opt = getopt_long(argc,argv,short_opts,long_opts,NULL)) != -1)
    {
      switch(curr_opt)
	{
	case 'i':
	  {
	    input = fopen(argv[1],"r");
	    if(input == NULL)
	      {
		fprintf(stderr,"Could not open %s\nReason: %s\n",argv[1],strerror(errno));
		exit(errno);
	      }
	    break;
	  }
	case 'h':
	  print_help();
	  exit(0);
	  break;
	case 'v':
	  printf("xyz2sas %s\n",PACKAGE_VERSION);
	  exit(0);
	default:
	  fprintf(stderr,"Unknown flag: %c\n",curr_opt);
	  break;
	}
    }

  if(optind >= argc)
    {
      fprintf(stderr,"At least the number of atoms is needed.\n");
      return EINVAL;
    }
  else
    strncpy(num_atom_str,argv[optind],4096);

  natoms = atoi(num_atom_str);
  if(natoms == 0)
    {
      fprintf(stderr,"No atoms provided: %s\n",num_atom_str);
      return 0;
    }
  fprintf(stderr,"Loading %lu coordinates\n",natoms);
  retval = load_coordinates(&crds[0],&crds[1],&crds[2],&radii,natoms,input);
  if(retval != natoms)
    {
      fprintf(stderr,"Could not load all atoms (expected: %lu). Retval: %d\n",natoms,retval);
      if(retval < 0)
	return -retval;//errno
      return 1;
    }


  area = molsurf(crds[0],crds[1],crds[2],radii,natoms,0);

  natoms = 0;
  for(;natoms < 3;natoms++)
    free(crds[natoms]);
  free(radii);

  printf("%lf\n",area);

  return 0;

}
