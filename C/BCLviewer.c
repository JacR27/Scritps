#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bcl.h"



int rflag;
char printFile[50];

void printBCL(char filename[])
{
  unsigned long uniqueQ[256] = {0};
  unsigned long uniqueB[256] = {0};  
  int nClusters;
  int i;
  char cluster[50] = {'\0'};
  unsigned char * baseCalls, * qualities, *bases;
  nClusters = getClusters(filename);
  baseCalls = malloc(nClusters);
  bases = malloc(nClusters);
  qualities = malloc(nClusters);
  readBaseCalls(nClusters, baseCalls, filename);
  splitBaseCalls(baseCalls, nClusters, bases, qualities);
  free(baseCalls);
  charHist(uniqueQ, qualities,nClusters);
  charHist(uniqueB, bases,nClusters);
  printHist(uniqueQ);
  printHist(uniqueB);
  char con[] = {'A','C','G','T'};
  if (0){
  for (i = 0; i < nClusters; ++i){
    cluster[0] = '\0';
    if ((i+5) % 5 == 0)
	sprintf(cluster,"%010d\t",i);
      
    printf("%s%c-%d\t%s", cluster ,con[bases[i]], qualities[i], ((i+1) % 5) ? "" : "\n");
    //printf("%d%c-%d\t%s",(i % 5) ? i : i , con[bases[i]], qualities[i], (i % 5) ? "" : "\n");
    cluster[0] = '\0';
  }
  }
  free(qualities);
  free(bases);
}
     

void printcompressedBCL(char filename[])
{
  unsigned long uniqueQ[256] = {0};
  unsigned long uniqueB[256] = {0};  
  int nClusters;
  int i;
  char cluster[50] = {'\0'};
  unsigned char * baseCalls, * qualities, *bases;
  nClusters = getClusters(filename);
  baseCalls = malloc(nClusters);
  bases = malloc(nClusters);
  qualities = malloc(nClusters);

  readReducedBaseCalls(nClusters, baseCalls, filename);
  splitBaseCalls(baseCalls, nClusters, bases, qualities);
  free(baseCalls);
  charHist(uniqueQ, qualities,nClusters);
  charHist(uniqueB, bases,nClusters);
  printHist(uniqueQ);
  printHist(uniqueB);
  char con[] = {'A','C','G','T'};
  for (i = 0; i < nClusters; ++i){
    cluster[0] = '\0';
    if ((i+5) % 5 == 0)
  	sprintf(cluster,"%010d\t",i);
      
    printf("%s%c-%d\t%s", cluster ,con[bases[i]], qualities[i], ((i+1) % 5) ? "" : "\n");
    printf("%d%c-%d\t%s",(i % 5) ? i : i , con[bases[i]], qualities[i], (i % 5) ? "" : "\n");
    cluster[0] = '\0';
  }
  free(qualities);
  free(bases);
}


void printDenltiplexedBCL(char filename[],int readlen)
{
  unsigned long uniqueQ[256] = {0};
  unsigned long uniqueB[256] = {0};
  int nClusters;
  int i;
  char cluster[50] = {'\0'};
  unsigned char * baseCalls, * qualities, *bases;
  nClusters = getClusters(filename)*readlen;
  baseCalls = malloc(nClusters);
  bases = malloc(nClusters);
  qualities = malloc(nClusters);
  readBaseCalls(nClusters, baseCalls, filename);
  splitBaseCalls(baseCalls, nClusters, bases, qualities);
  free(baseCalls);
  charHist(uniqueQ, qualities,nClusters);
  charHist(uniqueB, bases,nClusters);
  printHist(uniqueQ);
  printHist(uniqueB);
  char con[] = {'A','C','G','T'};
  for (i = 0; i < nClusters; ++i){
    cluster[0] = '\0';
    if ((i+5) % 5 == 0)
      sprintf(cluster,"%010d\t",i);

    printf("%s%c-%d\t%s", cluster ,con[bases[i]], qualities[i], ((i+1) % 5) ? "" : "\n");
    //printf("%d%c-%d\t%s",(i % 5) ? i : i , con[bases[i]], qualities[i], (i % 5) ? "" : "\n");
    cluster[0] = '\0';
  }
  free(qualities);
  free(bases);
}




main( int argc, char *argv[])

    {
      extern int rflag;
      rflag = 0; // reduceingResolution
      int demultiplexed = 1;
      int index;
      int c;
      opterr = 0;
      while ((c = getopt (argc, argv, "p:rd:")) != -1)
	switch (c)
	  {
	  case 'd':
	    demultiplexed = atoi(optarg);
	  case 'r':
	    rflag = 1;
	    break;
          case 'p':
            strcpy(printFile, optarg);
	    break;
	  case '?':
	    if (optopt == 'p')
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

      for (index = optind; index < argc; index++)
	printf ("Non-option argument %s\n", argv[index]);
      
      
      if (!rflag)
	printBCL(printFile);
      else if (demultiplexed != 1)
	printDenltiplexedBCL(printFile, demultiplexed);
      else
      	printcompressedBCL(printFile);
      return 0;
    }

