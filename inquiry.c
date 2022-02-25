/* FORTRAN interfaced-C function for extracting information from the .omf file*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

//void inquiryomf(int *xcount, int *ycount,
void inquiryomf_(int *xcount, int *ycount,
   int *zcount, int *startline, float *cellsize, float *zcellsize, 
   float *mscalefactor, char *omffilename) 
{
  FILE *omffp;
  char line[240], *p1, *p2;
  unsigned short int check=0;

  //  printf("size of xcount: %d\n",sizeof(signed int));

  //Strip trailing spaces
  int i = strlen(omffilename) - 1;
  while (isspace(omffilename[i]))
    omffilename[i--] = '\0';

  if ((omffp = fopen(omffilename, "r")) == NULL)
  {
    fprintf(stderr, "Inquiryomf: cannot open %s!\n", omffilename);
    exit(-1);
  }
  *startline=1;

  // Read lines until the start of data section
  while (fgets(line,240,omffp) != NULL && strstr(line,"Begin: Data") == NULL)
  {
    (*startline)++;

    // Strip spaces and convert to lowercase
    p1 = line;
    p2 = line;
    while(*p1 != 0) {
      if(isspace(*p1)) {
        ++p1;
      }
      else
      *p2++ = (char)tolower(*p1++); 
    }
    *p2 = 0; 

    if ((p1=strstr(line,"#xstepsize:")) != NULL)
    {
      p1+=strlen("#xstepsize:");
      if (!(*cellsize=atof(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x01;
    }
    if ((p1=strstr(line,"#zstepsize:")) != NULL)
    {
      p1+=strlen("#zstepsize:");
      if (!(*zcellsize=atof(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x02;
    }
    if ((p1=strstr(line,"#xnodes:")) != NULL)
    {
      p1+=strlen("#xnodes:");
      if (!(*xcount=atoi(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x04;
    }
    if ((p1=strstr(line,"#ynodes:")) != NULL)
    {
      p1+=strlen("#ynodes:");
      if (!(*ycount=atoi(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x08;
    }
    if ((p1=strstr(line,"#znodes:")) != NULL)
    {
      p1+=strlen("#znodes:");
      if (!(*zcount=atoi(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x10;
    }
    if ((p1=strstr(line,"#valuemultiplier:")) != NULL)
    {
      p1+=strlen("#valuemultiplier:");
      if (!(*mscalefactor=atof(p1)))
      {
        printf("Wrongly formatted preamble in %s file!\n", omffilename);
        exit(-1);
      }
      check += 0x20;
    }
  }
  
    if ((check ^ 0x3F))
    {
      printf("Wrongly formatted preamble in %s file!\n", omffilename);
      exit(-1);
    }
  if (feof(omffp))
  {
    printf("Missing data in %s file!\n", omffilename);
    exit(-1);
  }
}
