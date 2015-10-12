#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
#include <ctype.h>
#include <unistd.h>




/* parse, filter, demultiplex, reduse and split, bcl files */

int readHeader(FILE *input);

void printHist(unsigned long unique[]);

void charHist(unsigned long unique[] , unsigned char array[], int len);

void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[]);

void splitBaseCalls(unsigned char baseCalls[], unsigned int nBases, unsigned char bases[], unsigned char qualities[]);

void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[],unsigned char filter[],unsigned char filteredBaseCalls[]);

unsigned char * join(unsigned int len, unsigned char array[],int bytesToJoin);

void rejoin(unsigned char *reduced, unsigned int len, unsigned char qualities[],unsigned char bases[]);

void demultiplex(unsigned int len, unsigned char baseCalls[]);

void reduceQualities(unsigned char reducedQualities[], unsigned char origonalQualities[], unsigned int len, unsigned char qualities[], unsigned char qualityMap[], int nQualities);

void printArray(unsigned char array[], unsigned int len, char *fileName, unsigned int header);

unsigned char * interleafe(unsigned char array1[], unsigned char array2[], unsigned int len);

void Firstoutputs(void);

void analysis(void);

void printArrayStdout(unsigned char array[], unsigned int len, unsigned int header);

int getClusters(unsigned char filename[]);

int getFilterPasses(unsigned char filename[]);

int getFilterMask(unsigned char filter[], unsigned int len, unsigned char filename[]);

void stuff(void);

void printFileNames(char extention[]);

void frr(void);

int lvalue, tvalue, wvalue, svalue, cvalue;
int fflag, rflag, mflag, jflag, iflag, dflag;

main( int argc, char *argv[])

    {
      extern int lvalue, tvalue, wvalue, svalue, cvalue;
      
      lvalue = 1;
      tvalue = 8;
      wvalue = 3;
      svalue = 2;
      cvalue = 202;

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

      while ((c = getopt (argc, argv, "l:t:w:s:c:frmjid")) != -1)
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
	    if (optopt == 'l' || optopt == 't' || optopt == 's' || optopt == 'f' || optopt == 'c')
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

      printf ("lane = %d, cycles = %d, ntiles = %d, nswaths = %d, nsurfaces = %d, filter = %d, reduce resolution = %d, remap = %d, join = %d, split = %d, demultiplexed = %d\n",
	      lvalue ,cvalue, tvalue , wvalue , svalue , fflag , rflag , mflag, jflag, iflag, dflag);

      for (index = optind; index < argc; index++)
	printf ("Non-option argument %s\n", argv[index]);

      frr();
      
      return 0;
    }




void frr(void)
{
  #define MAXFILENAMELEN 20
  #define MAXFOLDERNAMELEN 8
  int nFiles;
  int i, j, k , n, f;  

  extern int lvalue, tvalue, wvalue, svalue, cvalue;
  
  extern int fflag, rflag, mflag, jflag, iflag, dflag; //filtering, reducingResolution, remapping, joining, spliting, demultiplexing;

  char folderNames[cvalue][MAXFOLDERNAMELEN];
  
  for (f = 1; f <= cvalue; ++f){
     sprintf(folderNames[f-1],"C%d.1/",f);
  }

  nFiles = (svalue * wvalue * tvalue);
  char filterNames[nFiles][MAXFILENAMELEN]; /*weast of memory when not filtering???*/
  char bclNames[nFiles][MAXFILENAMELEN];
  n = 0;
  for (i = 1; i <= svalue ; ++i){
    for (j = 1; j <= wvalue; ++j){
      for (k = 1; k <= tvalue; ++k){ 
	if (fflag){
	  sprintf(filterNames[n],"s_%d_%i%i%02d%s",lvalue,i,j,k,".filter");
	}
	sprintf(bclNames[n],"s_%d_%i%i%02d%s",lvalue,i,j,k,".bcl");
	n++;
      }
    }
  }


  int nQualites = 38;
  unsigned char origonalQualities[] = {0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41};
  
  printf("nQualities: %d\n",sizeof(origonalQualities));

  unsigned char qualityMap[] =        {0,7,7,7,7,7,11,11,11,11,11,11,22,22,22,22,22,22,22,27,27,27,27,27,32,32,32,32,32,37,37,37,37,37,42,42,42,42};
  

  /*meta stats*/
  unsigned long uniqueQ[256] = {0};
  unsigned long uniqueB[256] = {0};
  
  
  unsigned int nClusters;
  unsigned int nPasses;
  unsigned int nPostFiltering;
  int fi, ci;
  unsigned char * baseCalls, *filteredBaseCalls, *bases, *qualities, *filter, *postFilteringBaseCalls, *reducedQualities, *rejoined;
  
  char path[MAXFILENAMELEN+MAXFOLDERNAMELEN];
  char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+10];
  
  for (fi=0; fi <nFiles;++fi){ /* loop through tiles*/
  
    path[0] = '\0';
    strcpy(path, folderNames[0]);
    strcat(path, bclNames[fi]);
    nClusters = getClusters(path);
    
    
    if (fflag){
      filter = malloc(nClusters);
      nPasses = getFilterMask(filter, nClusters,filterNames[fi]);
      filteredBaseCalls= malloc(nPasses);
      nPostFiltering = nPasses;
    }
    
    else{
      nPostFiltering = nClusters;
    }
    
    baseCalls = malloc(nClusters);
    if (rflag){
       bases = malloc(nPostFiltering);
       qualities = malloc(nPostFiltering);
       reducedQualities = malloc(nPostFiltering);
    }

    for (ci = 0; ci <cvalue; ++ci){ /*loop through cycles*/
      path[0] = '\0';
      strcpy(path, folderNames[ci]);
      strcat(path, bclNames[fi]);
      readBaseCalls(nClusters, baseCalls, path);
      
      if (fflag){
         filterBaseCalls(nClusters, baseCalls, filter, filteredBaseCalls);
	 postFilteringBaseCalls = filteredBaseCalls;
      }
      else{
         postFilteringBaseCalls = baseCalls;
      }
    
      if (rflag){
        outpath[0] = '\0';
	strcpy(outpath, path);
	strcat(outpath, ".rr");
	splitBaseCalls(postFilteringBaseCalls, nPostFiltering, bases, qualities);
	reduceQualities(reducedQualities,origonalQualities, nPostFiltering, qualities, qualityMap,nQualites);
	
	charHist(uniqueQ, reducedQualities,nPostFiltering);
	charHist(uniqueB, bases,nPostFiltering);
	rejoined = malloc(nPostFiltering);
	rejoin(rejoined, nPostFiltering, reducedQualities,bases);
	printArray(rejoined, nPostFiltering,outpath ,nPostFiltering);
	free(rejoined);
	
	  
      }
      else{
         printf("check1");
         //reducedQualities = qualities;
      }
    }
    

    //printArray(qualities, nClusters,"qualities", nClusters);
    //charHist(qualities,nPostFiltering);
    
    
    free(postFilteringBaseCalls);
    if (fflag){
      free(filter);
      free(baseCalls);
    }
    
    if (rflag){
      free(bases);
      free(qualities);
      free(reducedQualities);
    }
    //
  }
  printHist(uniqueQ);
  printHist(uniqueB);
}



int celingDev(int dividend, int devisor)
{
  return ((dividend / devisor) + (dividend % devisor != 0));
}

/* read 32 bit little eden int from file */
int readHeader(FILE *input)
{
  #define BYTESININT32 4
  #define BITSINBYTE 8

  unsigned int nClusters;
  int c , i;
  
  nClusters = 0;

  for (i = 0; i<BYTESININT32; ++i){
    nClusters = nClusters + fgetc(input) * pow(2,i*BITSINBYTE);
  }
  return nClusters;
}

unsigned char * interleafe(unsigned char array1[], unsigned char array2[], unsigned int len)
{
  int j;
  unsigned char *interleved;
  
  interleved = malloc(len*2);
  for (j=0; j<len; ++j){
    interleved[j*2] = array1[j];
    interleved[(j*2)+1] = array2[j];
  }
  return interleved;
}
       
//interlevedQB = interleafe(qualities,bases,nClusters);

/* read n bytes from standard in to array */
void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[])
{
  #define BYTESININT32 4
  #define BITSINBYTE 8

  int i;
  int nClusters;
  unsigned char c;
  FILE * inputFile;
  inputFile = fopen(filename,"r");
  unsigned char rawNClusters[BYTESININT32];
  fread(rawNClusters,1, BYTESININT32,inputFile);
  nClusters =0;
  for (i = 0; i<BYTESININT32; ++i){
    nClusters = nClusters + rawNClusters[i] * pow(2,i*BITSINBYTE);
  }
  printf("filename: %s, nClusters %d\n",filename,nClusters);
  fread(baseCalls, 1, nBases,inputFile);
  fclose(inputFile);
}



void splitBaseCalls(unsigned char baseCalls[], unsigned int nBases, unsigned char bases[], unsigned char qualities[])
{
  int i;
  unsigned char c;
  for (i=0; i<nBases; ++i){
    c = baseCalls[i];
    bases[i] = (c % 4);
    qualities[i] = (c / 4);
  }
}

void charHist(unsigned long unique[] , unsigned char array[], int len)
{
  #define MAXUNIQUE 256
  int i;
  int nUnique;
  unsigned char c;
  //unsigned int unique[MAXUNIQUE] = {0};

  //nUnique = 0;
  for (i=0;i<len;++i){
    c = array[i];
    //if (!unique[c]){
    //   nUnique += 1;
    //}
    unique[c] +=1;
  }
  //unsigned char elements[nUnique];
  //unsigned int count[nUnique];
  //int n;
  //n = 0;
  //for (i=0; i<MAXUNIQUE;i++){
  //    if (unique[i]){
  //      elements[n] = i;
  //    count[n] = unique[i];
  //    n++;
  //   }
  //}
  //printHist(nUnique, elements, count);
}

void printHist(unsigned long unique[])
{
  #define MAXUNIQUE 256
  unsigned long i, c, tot;
  tot = 0;
  for (i=0;i<MAXUNIQUE;i++){
    if (c = unique[i]){
      printf("element: %lu, count %lu \n", i, c);
      tot +=c;  
    }
  }
  printf("total: %lu\n", tot);
}



/*Join bytes that do not use there maximum range of values */
unsigned char * join(unsigned int len, unsigned char array[], int bytesToJoin)
{
  #define BITSINBYTE 8 
  int i;
  int maxMod;
  int newLen;
  unsigned char *reduced;
  maxMod = bytesToJoin - 1;
  newLen = celingDev(len,bytesToJoin);
  reduced = calloc(newLen,1);
  //printf("%d",newLen);
  
  for (i = 0; i < len; ++i){
    reduced[i/bytesToJoin] = reduced[i/bytesToJoin] + array[i] * pow(2,(maxMod-(i%bytesToJoin))*(BITSINBYTE/bytesToJoin));
    //printf("base: %d index: %d newbase: %d\n",array[i],i/bytesToJoin,reduced[i/bytesToJoin]);
    }
  return reduced;
}

/*Join bytes that do not use there maximum range of values */
void rejoin(unsigned char *reduced, unsigned int len, unsigned char qualities[],unsigned char bases[])
{
  int i;
  for (i = 0; i < len; ++i)
    reduced[i] = (qualities[i] << 2) + bases[i];
}


 void reduceQualities(unsigned char reducedQualities[],unsigned char origonalQualities[], unsigned int len, unsigned char qualities[],unsigned char qualityMap[],int nQualities)
{
  int i , j;  
  for (i = 0; i < len; ++i){
    for (j = 1; j < nQualities; ++j){
      if (qualities[i] == origonalQualities[j]){
	reducedQualities[i] = qualityMap[j];
      }
    }
  }
}

 void printArray(unsigned char array[], unsigned int len, char *fileName,unsigned int header)
{
  int i;
  FILE *outHandle;
  outHandle = fopen(fileName,"w");
  fwrite(&header,sizeof(header),1,outHandle);
  fwrite(array, 1, len, outHandle);
  fclose(outHandle);
}


void printArrayStdout(unsigned char array[], unsigned int len, unsigned int header)
{
  write(1, &header,sizeof(header));
  write(1, array, len); 
}


/* read v3 filter file into array */
int getFilterMask(unsigned char filter[], unsigned int len, unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  int filterVersion; /* version, contained in header of filter file */
  int nPassed; /* number of passed bases */
  int i; /* loop variable */
  unsigned char c; /*temp var*/
  FILE *fp; 
  fp = fopen(filename,"r");
  readHeader(fp);
  filterVersion = readHeader(fp);
  nClusters = readHeader(fp);
  nPassed = 0;
  for (i = 0; i < nClusters; ++i){
    c = fgetc(fp);
    filter[i] = c;
    nPassed += c;
  }
  printf("Filter file:%s, version:%d number of clusters: %d, number of passed clusters %d (%f%%)\n",filename,filterVersion,nClusters,nPassed,(float)nPassed/nClusters *100);
  fclose(fp);
  
  return nPassed;
}

int getClusters(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  FILE *fp;
  fp = fopen(filename,"r");
  nClusters = readHeader(fp);
  return nClusters;
}
int getFilterPasses(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  int nPassed; /* number of passed bases */
  int i; /* loop variable */
  unsigned char c; /*temp var*/
  FILE *fp;
  fp = fopen(filename,"r");
  readHeader(fp);
  readHeader(fp);
  nClusters = readHeader(fp);
  nPassed = 0;
  for (i = 0; i < nClusters; ++i){
    c = fgetc(fp);
    nPassed += c;
  }
  fclose(fp);
  return nPassed;
}


void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[], unsigned char filter[], unsigned char filteredBaseCalls[])
{
  int j , i;
  j = 0;
  int nFilteredBases;
  for (i = 0; i < nBases; ++i){
    if (filter[i] == 1){
      filteredBaseCalls[j] = baseCalls[i];
      j += 1;
    }
  }
}

