#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
#include <ctype.h>
#include <unistd.h>
#include "bcl.h"

int mode, lvalue, tvalue, wvalue, svalue, cvalue, nFiles, cfvalue;
int fflag, rflag, mflag, jflag, iflag, dflag;

char bclFormat[30];

void frr(void)
{
  #define MAXFOLDERNAMELEN 8
  extern int cvalue;
  //extern char **folderNames;
  int nCycles;
  nCycles = (cvalue + 1) - cfvalue;
  printf("%d\n", nCycles);
  char folderNames[nCycles][MAXFOLDERNAMELEN];
  int f;
  int t;
  for (f = cfvalue, t = 0; f <= cvalue; ++f, ++t){
    printf("%d\n",f);
    sprintf(folderNames[t],"C%d.1/",f);
  }



#define MAXFILENAMELEN 50


  int i, j, k , n;
  nFiles= (svalue * wvalue * tvalue);
  char bclNames[nFiles][MAXFILENAMELEN];
  char filterNames[nFiles][MAXFILENAMELEN];

  extern int svalue, tvalue, wvalue, lvalue, nFiles;  
  extern int fflag; //filtering
  n = 0;
  printf("nfiles %d\n",nFiles);
  for (i = 1; i <= svalue ; ++i){
    for (j = 1; j <= wvalue; ++j){
      for (k = 1; k <= tvalue; ++k){
	if (fflag)
	  snprintf(filterNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.%s",lvalue,i,j,k,"filter");
	snprintf(bclNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.%s",lvalue,i,j,k,bclFormat);
	n++;
      }
    }
  }
  
  printf("got file names");



  

  int nQualites = 8;

  unsigned char qualityMap[] = {0,1,2,3,4,5,6,7};
  unsigned char origonalQualities[] = {0,14,14,14,34,34,34,42};
  printf("nQualities: %d\n",sizeof(origonalQualities));

  
  /*meta stats*/
  
  
  unsigned long uniqueQ[256] = {0};
  unsigned long uniqueB[256] = {0};
  unsigned long totalClusters;
  unsigned long totalPasses;
  totalClusters = 0;
  totalPasses = 0;
  char fpath[MAXFILENAMELEN+MAXFOLDERNAMELEN];
  extern int mode, rflag, mflag, jflag, iflag, dflag; //reducingResolution, remapping, joining, spliting, demultiplexing;
  fpath[0] = 0;
  int fi;
  printf("stuff %s\n",fpath);
  
  for (fi=0; fi <nFiles;++fi){ /* loop through tiles*/
  unsigned int nClusters;
  unsigned int nPasses;
  unsigned int nPostFiltering;
  unsigned char *filter;
  
  switch(mode) {
 case 1:
   printf("geting number of cluster passing filtering %s\n",filterNames[fi]);
   nClusters= getFilterClusters(filterNames[fi]);
   nPasses = getFilterPasses(filterNames[fi]);
   printf("%s: %d passed out of %d total clusters (%f%%)\n",filterNames[fi],nPasses,nClusters, ((float)nPasses/nClusters)*100);
   totalClusters += nClusters;
   totalPasses += nPasses;
   printf("%s: %d passed out of %d total clusters (%f%%)\n",filterNames[fi],totalPasses,totalClusters, ((float)totalPasses/totalClusters)*100);
   break;
 case 2:
   printf("generating file name\n");
   fpath[0] = '\0';
   strcpy(fpath, folderNames[0]);
   strcat(fpath, bclNames[fi]);
   printf("%s\n",fpath);
   nClusters = getClusters(fpath);
   printf("nclusters: %d\n", nClusters);  
   
   if (fflag){
      printf("getting filter file\n");
      filter = malloc(nClusters);
      nPasses = getFilterMask(filter, nClusters,filterNames[fi]);
      nPostFiltering = nPasses;
      printf("%d\n",nPasses);
   }
    else{
      nPostFiltering = nClusters;
    }

    int ci;
    printf("entering parrallel tile loop\n");
    #pragma omp parallel for
    for (ci = 0; ci <nCycles; ++ci){ /*loop through cycles*/

     
      unsigned char * baseCalls, *filteredBaseCalls, *bases, *qualities, *reducedQualities, *rejoined, *interleafed, *joined, *postFilteringBaseCalls;
      char path[MAXFILENAMELEN+MAXFOLDERNAMELEN];
      char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+10];
      path[0] = '\0';
      strcpy(path, folderNames[ci]);
      strcat(path, bclNames[fi]);

      baseCalls = malloc(nClusters);
      readBaseCalls(nClusters, baseCalls, path);

      if (fflag){
         outpath[0] = '\0';
	 strcpy(outpath, path);	
	 strcat(outpath, ".filtered");
	 filteredBaseCalls= malloc(nPasses);
	 printf("%d,%d, %d,%d,%d\n",nClusters, nPostFiltering,sizeof(filter),sizeof(filteredBaseCalls),nPasses);
         filterBaseCalls(nClusters, baseCalls, filter, filteredBaseCalls);
	 printArray(filteredBaseCalls, nPostFiltering, outpath ,nPostFiltering);
	 postFilteringBaseCalls = filteredBaseCalls;
	 free(baseCalls);
      }
      else{
         postFilteringBaseCalls = baseCalls;
      }

      if (mflag){
	bases = malloc(nPostFiltering);
        qualities = malloc(nPostFiltering);
        splitBaseCalls(postFilteringBaseCalls, nPostFiltering, bases, qualities);
        outpath[0] = '\0';
        strcpy(outpath, path);
        strcat(outpath, ".rm.gz");
        reducedQualities = malloc(nPostFiltering);
        reduceQualities(reducedQualities,origonalQualities,nPostFiltering,qualities,qualityMap,nQualites);
	interleafed = interleafe(reducedQualities, bases, nPostFiltering);
	joined = join(nPostFiltering * 2, interleafed, 4);
	printArray(joined, celingDev(nPostFiltering*2,4),outpath,nPostFiltering);
        free(joined);
	free(interleafed);
        free(bases);
        free(qualities);
        free(reducedQualities);
      }
     
      if (rflag){
	bases = malloc(nPostFiltering);
	qualities = malloc(nPostFiltering);
	splitBaseCalls(postFilteringBaseCalls, nPostFiltering, bases, qualities);    
        outpath[0] = '\0';
	strcpy(outpath, path);	
        strcat(outpath, ".rrr_42.gz");
	reducedQualities = malloc(nPostFiltering);
	reduceQualities(reducedQualities,origonalQualities, nPostFiltering, qualities, qualityMap,nQualites);
	rejoined = malloc(nPostFiltering);
	rejoin(rejoined, nPostFiltering, reducedQualities,bases);
	printArray(rejoined, nPostFiltering,outpath ,nPostFiltering);
	free(rejoined);
	free(bases);
	free(qualities);
	free(reducedQualities);
      }
      else
	;
      if (iflag){
	bases = malloc(nPostFiltering);
        qualities = malloc(nPostFiltering);
        splitBaseCalls(postFilteringBaseCalls, nPostFiltering, bases, qualities);
        outpath[0] = '\0';
        strcpy(outpath, path);
        strcat(outpath, ".q.gz");
        printArray(qualities, nPostFiltering,outpath ,nPostFiltering);
        outpath[0] = '\0';
        strcpy(outpath, path);
        strcat(outpath, ".b.gz");	
	printArray(bases, nPostFiltering,outpath ,nPostFiltering);
        free(bases);
        free(qualities);
        
      }
      free(postFilteringBaseCalls);
    }

    if (fflag)
      free(filter);
    break;
  case 3:
    printf("%s\n",filterNames[fi]);
    
    for (ci = 0; ci <nCycles; ++ci){ /*loop through cycles*/
      printf("nCycles %d %d\n",nCycles,ci);
      fpath[0] = '\0';
      strcpy(fpath, folderNames[ci]);
      strcat(fpath, bclNames[fi]);
      printf("%s\n",fpath);  
    }
    break;
  case 4:
    printf("generating file name\n");
    fpath[0] = '\0';
    strcpy(fpath, folderNames[0]);
    strcat(fpath, bclNames[fi]);
    printf("%s\n",fpath);
    nClusters = getClusters(fpath);
    printf("nclusters: %d\n", nClusters);
    #pragma omp parallel for
    for (ci = 0; ci <nCycles; ++ci){
      unsigned char * baseCalls, *bases, *qualities;
      char path[MAXFILENAMELEN+MAXFOLDERNAMELEN];
      char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+10];
      path[0] = '\0';
      strcpy(path, folderNames[ci]);
      strcat(path, bclNames[fi]);
      baseCalls = malloc(nClusters);
      readBaseCalls(nClusters, baseCalls, path);
      bases = malloc(nClusters);
      qualities = malloc(nClusters);
      splitBaseCalls(baseCalls, nClusters, bases, qualities);      	
      #pragma omp critical 
      {
	charHist(uniqueQ, qualities,nClusters);
	charHist(uniqueB, bases,nClusters);
      }
      free(bases);
      free(qualities);
      free(baseCalls);
    } 
  }
  }
  
  printHist(uniqueQ);
  printHist(uniqueB);
}


main( int argc, char *argv[])

    {
      extern int mode, lvalue, tvalue, wvalue, svalue, cvalue;
      extern char bclFormat[];
      
      lvalue = 1;
      tvalue = 24;
      wvalue = 2;
      svalue = 2;
      cvalue = 302;
      cfvalue = 1;
      mode = 0;
      strcpy(bclFormat, "bcl");
      

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

      while ((c = getopt (argc, argv, "o:l:t:w:s:c:z:frmjid")) != -1)
	switch (c)
	  {
	  case 'o':
	    mode = atoi(optarg);
	    if ((mode == 1) || (mode == 3))
	      fflag = 1;
	    break;
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
	  case 'z':
	    strcpy(bclFormat, optarg);
	    printf("%s\n",bclFormat);
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


