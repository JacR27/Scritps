#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
#include <ctype.h>
#include <unistd.h>
#include "deflate.h"
#include "bcl.h"

int lvalue, tvalue, wvalue, svalue, cvalue, nFiles;
int fflag, rflag, mflag, jflag, iflag, dflag;
char bclFormat[100];
char outExtention[100];
char printFile[50];
main( int argc, char *argv[])

    {
      extern int lvalue, tvalue, wvalue, svalue, cvalue;
      extern char bclFormat[];
      
      lvalue = 1;
      tvalue = 24;
      wvalue = 2;
      svalue = 2;
      cvalue = 302;
      strcpy(bclFormat, "bcl.filtered");
      strcpy(outExtention, "dbcl.2reads.filtered.gz");

      extern int fflag, rflag, mflag, jflag, iflag, dflag;
      fflag = 0; //filtering
      rflag = 0; // reduceingResolution
      mflag = 0; // remapping
      jflag = 0; // joining
      iflag = 0; // spliting
      dflag = 0; // demultiplexing

      
      int index;
      int c;

      opterr = 0;

      while ((c = getopt (argc, argv, "l:t:w:s:c:z:frmjid")) != -1)
	switch (c)
	  {
	  case 'l':
	    lvalue = atoi(optarg);
	    break;
	  case 't':
	    tvalue = atoi(optarg);
	    break;
	  case 'w':
	    wvalue = atoi(optarg);
	    break;
	  case 's':
	    svalue = atoi(optarg);
	    break;
	  case 'c':
	    cvalue = atoi(optarg);
	    break;
	  case 'f':
	    fflag = 1;
	    break;
	  case 'r':
	    rflag = 1;
	    break;
	  case 'm':
	    mflag = 1;
	    break;
	  case 'j':
	    jflag = 1;
	    break;
	  case 'i':
	    iflag = 1;
	    break;
	  case 'd':
	    dflag = 1;
	    break;
	  case '?':
	    if (optopt == 'l' || optopt == 'o' || optopt == 't' || optopt == 's' || optopt == 'f' || optopt == 'c')
	      fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	    else if (isprint (optopt))
	      fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	    else
	      fprintf (stderr,
		       "Unknown option character `\\x%x'.\n",
		       optopt);
	    return 1;
	  default:
	    abort ();
	  }

      printf ("lane = %d, cycles = %d, ntiles = %d, nswaths = %d, nsurfaces = %d, bclformat = %s, filter = %d, reduce resolution = %d, remap = %d, join = %d, split = %d, demultiplexed = %d\n",
	      lvalue ,cvalue, tvalue , wvalue , svalue ,bclFormat, fflag , rflag , mflag, jflag, iflag, dflag);

      for (index = optind; index < argc; index++)
	printf ("Non-option argument %s\n", argv[index]);
      
      	
      frr();
      
      return 0;
    }



void frr(void)
{
  #define MAXFOLDERNAMELEN 8
  printf("starting main\n");
  extern int cvalue;
  //extern char **folderNames;
  int nCycles;
  nCycles = (cvalue + 1) - 1;
  printf("%d\n", nCycles);
  char folderNames[nCycles][MAXFOLDERNAMELEN];
  int f;
  int t;
  for (f = 1, t = 0; f <= cvalue; ++f, ++t){
    printf("%d\n",f);
    sprintf(folderNames[t],"C%d.1/",f);
  }



#define MAXFILENAMELEN 50


  int i, j, k , n;
  nFiles= (svalue * 1 * 2);
  char bclNames[nFiles][MAXFILENAMELEN];
  extern int svalue, tvalue, wvalue, lvalue, nFiles;  
  n = 0;
  printf("nfiles %d\n",nFiles);
  for (i = 1; i <= svalue ; ++i){
    for (j = 2; j <= wvalue; ++j){
      for (k = 23; k <= tvalue; ++k){
	snprintf(bclNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.",lvalue,i,j,k);
	n++;
      }
    }
  }
  
  printf("got file names");

  

  #define MAXEXTENTIONLEN 100
    
  char fpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+MAXEXTENTIONLEN];
  fpath[0] = '\0';
  int u;
  int w;
  char inpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+MAXEXTENTIONLEN];
  char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+MAXEXTENTIONLEN];
  

  for( w=0;w<nFiles;w++){
    outpath[0] = '\0';
    strcat(outpath,bclNames[w]);
    strcat(outpath,outExtention);
    unsigned int nClusters;
    fpath[0] = '\0';
    strcpy(fpath, folderNames[0]);
    strcat(fpath, bclNames[w]);
    strcat(fpath,bclFormat);
    nClusters = getClusters(fpath);
    unsigned char *demultiplexed;
    unsigned char *baseCalls;
    baseCalls = malloc(nClusters);
    demultiplexed = malloc(nClusters*nCycles);
    printf("%d\n",nClusters);
    
    
    for (u=0;u<nCycles;u++){
      strcpy(inpath,folderNames[u]);
      strcat(inpath,bclNames[w]);
      strcat(inpath,bclFormat);
      readBaseCalls(nClusters, baseCalls, inpath);
      unsigned long z;
      for (z = 0; z < nClusters; z++){
	demultiplexed[z*nCycles+u] = baseCalls[z];
      }
    }
    
    printf("%d\n",nClusters*nCycles);
    printArray(demultiplexed,nClusters*nCycles,outpath,nClusters);
    free(demultiplexed);
    free(baseCalls);
  }
}



