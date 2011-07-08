#include "molsurf.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>

int load_coordinates(REAL_T **xcrds, REAL_T **ycrds, REAL_T **zcrds, REAL_T **radii, long natoms, FILE *input)
{
  long linecount = 0;
  if(natoms == 0)
    return 0;
  if(xcrds == NULL || ycrds == NULL || zcrds == NULL || radii == NULL || input == NULL)
    {
      fprintf(stderr,"given a null pointer.");
      return -1;
    }

  *xcrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *ycrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *zcrds = (REAL_T*)malloc(natoms*sizeof(REAL_T));
  *radii = (REAL_T*)malloc(natoms*sizeof(REAL_T));

  while(!feof(input) && !ferror(input))
    {
      char buffer[4096],*tok;
      int test;
      if(fgets(buffer,4096,input) == NULL)
	{
	  return linecount;
	}
      if(strlen(buffer) == 0 || buffer[0] == '#')//comment with # character
	continue;
      if(!isdigit(buffer[0]))
	continue;

      test = sscanf(buffer,"%lf %lf %lf %lf",&(*xcrds)[linecount],&(*ycrds)[linecount],&(*zcrds)[linecount],&(*radii)[linecount]);
      (*radii)[linecount] += 1.4;
      linecount++;
    }

  return linecount;
}

int main(int argc, char **argv)
{
  FILE *input = stdin;
  char num_atom_str[4096];
  long natoms;
  REAL_T *crds[3],*radii;
  REAL_T area;
  int retval;

  natoms = 0;
  for(;natoms < 3;natoms++)
    crds[natoms] = NULL;
  radii = NULL;

  if(argc > 2)
    {
      input = fopen(argv[1],"r");
      if(input == NULL)
	{
	  fprintf(stderr,"Could not open %s\nReason: %s\n",argv[1],strerror(errno));
	  return errno;
	}
      strncpy(num_atom_str,argv[2],4096);
    }
  else if(argc < 2)
    {
      fprintf(stderr,"At least the number of atoms is needed.\n");
      return EINVAL;
    }
  else
    strncpy(num_atom_str,argv[1],4096);

  natoms = atoi(num_atom_str);
  if(natoms == 0)
    {
      fprintf(stderr,"No atoms provided: %s\n",num_atom_str);
      return 0;
    }
  
  retval = load_coordinates(&crds[0],&crds[1],&crds[2],&radii,natoms,input);
  if(retval != natoms)
    {
      fprintf(stderr,"Could not load all atoms.\n");
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
