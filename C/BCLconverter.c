#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
/* parse, filter, demultiplex, reduse and split, bcl files */

int readHeader(FILE *input);

void readBaseCalls(unsigned int nBases, unsigned char baseCalls[]);

void readAndSplitBaseCalls(unsigned int nBases, unsigned char bases[], unsigned char qualities[]);

void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[],unsigned char filter[]);

unsigned char * join(unsigned int len, unsigned char array[],int bytesToJoin);

void demultiplex(unsigned int len, unsigned char baseCalls[]);

unsigned char * reduceQualities(unsigned int len, unsigned char qualities[], unsigned char qualityMap[]);

void printArray(unsigned char array[], unsigned int len, char *fileName, unsigned int header);

unsigned char * interleafe(unsigned char array1[], unsigned char array2[], unsigned int len);

void Firstoutputs(void);

void sumfirst300(unsigned char array[]);

void analysis(void);

void printArrayStdout(unsigned char array[], unsigned int len, unsigned int header);

void printBCL(void);

int getFilterClusters(unsigned char filename[]);

int getFilterPasses(unsigned char filename[]);

int getFilterMask(unsigned char filter[], unsigned int len, unsigned char filename[]);

void stuff(void);

void printFileNames(char extention[]);



main( int argc, char *argv[])
{
  int n, m, x, l, ch;
  for (n = 1; n < argc; ++n){
    switch( (int)argv[n][0] ){
    case '/':
    case '-': 
      x = 0;
      l = strlen( argv[n] );
      for ( m=1; m <l; ++m ){
	ch = (int)argv[n][m];
	switch(ch){
	case 'p':
	  printBCL();
	  break;
	case 'n':
	  break;
	case 's':
	  stuff();
	  break;
	default:
	  printf("illegal option code = %c\n", ch);
	  break;
	}
      }
      break;
    default:
      printf( "test - %s\n", argv[n]);
      break;
    }
  }
}

void stuff(void)
{
  printFileNames(".filter");
}

void printFileNames(char extention[])
{
  #define MAXFILENAMELEN 20
  int nSurfaces, nSwaths, nTiles, nFiles;
  int i, j, k , n;
  char lane[] = "6";
  
  nSurfaces = 2;
  nSwaths = 3;
  nTiles = 8;
  
  nFiles = (nSurfaces * nSwaths * nTiles);
  char filenames[nFiles][MAXFILENAMELEN];
  n = 0;
  for (i = 1; i <= nSurfaces ; ++i){
    for (j = 1; j <= nSwaths; ++j){
      for (k = 1; k <= nTiles; ++k){         
	sprintf(filenames[n],"s_%s_%i%i%02d%s",lane,i,j,k,extention);
	n++;
      }
    }
  }
  int nClusters;
  nClusters = getFilterClusters(filenames[0]);
  unsigned char * filter;
  filter = malloc(nClusters);
  int nPasses;
  int fi;
  for (fi=0; fi <1;++fi){
    nPasses = getFilterMask(filter, nClusters,filenames[fi]);
  }
  int i;
  for (i = 0; i<nClusters; 
  free(filter);
  
}

void printBCL(void)
{
  int i;
  unsigned int nClusters;
  unsigned char *baseCalls;
  nClusters = readHeader(stdin);
  printf("Number of Clusters: %d\n", nClusters);
  readBaseCalls(nClusters, baseCalls);
  for (i = 0; i < nClusters; i++){
    printf("%d",baseCalls[i]);
    if ((i % 10) == 0){
      printf("\n");
    }
  }
}

void Firstoutputs(void)
{
  unsigned int nClusters;
  unsigned char *bases, *qualities, *jointBases, *jointQualities, *remapedQualities, *interlevedQB, *jointInterlevedQB;
  unsigned char qmap1[8] = {0,1,1,2,2,2,3,3};
  unsigned char qmap2[8] = {0,1,2,3,4,5,6,7};

  
  nClusters = readHeader(stdin);
  printf("Number of Clusters: %d\n", nClusters);
  
  bases = malloc(nClusters);
  qualities = malloc(nClusters);
  
  /* populate bases and qualities */
  readAndSplitBaseCalls(nClusters, bases, qualities);
  
  /* write qualities */
  printArray(qualities, nClusters, "qualities", nClusters);
  
  jointBases = join(nClusters, bases, 4);
  printArray(jointBases, celingDev(nClusters,4), "bases", nClusters);
  free(jointBases);

  remapedQualities = reduceQualities(nClusters, qualities, qmap1);
  free(qualities);
 
  jointQualities = join(nClusters, remapedQualities, 4);
  printArray(jointQualities, celingDev(nClusters,4), "rqualities",nClusters);
  free(jointQualities);
  
  interlevedQB = interleafe(remapedQualities,bases,nClusters);
  free(bases);
  free(remapedQualities);
  
  jointInterlevedQB = join(nClusters*2,interlevedQB,4);
  free(interlevedQB);

  printArray(jointInterlevedQB, celingDev(nClusters*2, 4),"basesandRqualities", nClusters); 
  free(jointInterlevedQB);
  
}

void analysis(void)
{
  unsigned int nClusters;
  unsigned char *basecalls, *joined;
  int i;
  nClusters = (readHeader(stdin)*151);
  printf("Number of Clusters: %d\n", nClusters);

  basecalls = malloc(nClusters);

  /* populate bases and qualities */
  readBaseCalls(nClusters, basecalls);

  joined = join(nClusters, basecalls, 4);
  free(basecalls);
  printArrayStdout(joined, celingDev(nClusters,4),nClusters/151);
  free(joined);
}


void sumfirst300(unsigned char array[])
{
  int i;
  int sum;
  sum = 0;
  for (i = 0; i < 300; ++i)
    sum += array[i];
  printf("%d\n",sum);
}
  

int celingDev(int dividend, int devisor)
{
  return ((dividend / devisor) + (dividend % devisor != 0));
}

/* read 32 bit little eden int from stdin */
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
void readBaseCalls(unsigned int nBases, unsigned char baseCalls[])
{
  int i;
  unsigned char c;
  for (i = 0; i < nBases; ++i){
    c = getchar();
    baseCalls[i] = c;
  }
}


void readAndSplitBaseCalls(unsigned int nBases, unsigned char bases[], unsigned char qualities[])
{
  int i;
  unsigned char c;
  for (i=0; i<nBases; ++i){
    c = getchar();
    bases[i] = (c % 4);
    qualities[i] = (c / 4);
  }
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

int getFilterClusters(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  int nPassed; /* number of passed bases */
  FILE *fp;
  fp = fopen(filename,"r");
  readHeader(fp);
  readHeader(fp);
  nClusters = readHeader(fp);
  nPassed = 0;
  fclose(fp);
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


void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[],unsigned char filter[])
{
  int j , i;
  j = 0;
  int nFilteredBases;
  unsigned char *filteredBaseCalls;
  filteredBaseCalls = malloc(nFilteredBases);
  for (i = 0; i < nBases; ++i){
    if (filter[i] == 1){
      filteredBaseCalls[j] = baseCalls[i];
      j += 1;
    }
  }
  free(filter);
  free(filteredBaseCalls);  
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


unsigned char * reduceQualities(unsigned int len, unsigned char qualities[],unsigned char qualityMap[8])
{
  unsigned char origonalQualities[8] = {0,7,11,22,27,32,37,42};
  unsigned char *reducedQualities;
  int i , j;
  reducedQualities = malloc(len);
  for (i = 0; i < len; ++i){
    for (j = 1; j < 8; ++j){
      if (qualities[i] == origonalQualities[j]){
	reducedQualities[i] = qualityMap[j];
      }
    }
  }
  return reducedQualities;
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


