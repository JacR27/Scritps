#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
/* parse, filter, demultiplex, reduse and split, bcl files */

int readHeader(FILE *input);

void printHist(int nUnique, char elements[], int count[]);

void charHist(unsigned char array[], int len);

void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[]);

void splitBaseCalls(unsigned char baseCalls[], unsigned int nBases, unsigned char bases[], unsigned char qualities[]);

void filterBaseCalls(unsigned int nPassed, unsigned int nBases, unsigned char baseCalls[],unsigned char filter[],unsigned char filteredBaseCalls[]);

unsigned char * join(unsigned int len, unsigned char array[],int bytesToJoin);

void demultiplex(unsigned int len, unsigned char baseCalls[]);

unsigned char * reduceQualities(unsigned int len, unsigned char qualities[], unsigned char qualityMap[]);

void printArray(unsigned char array[], unsigned int len, char *fileName, unsigned int header);

unsigned char * interleafe(unsigned char array1[], unsigned char array2[], unsigned int len);

void Firstoutputs(void);

void analysis(void);

void printArrayStdout(unsigned char array[], unsigned int len, unsigned int header);

int getFilterClusters(unsigned char filename[]);

int getFilterPasses(unsigned char filename[]);

int getFilterMask(unsigned char filter[], unsigned int len, unsigned char filename[]);

void stuff(void);

void printFileNames(char extention[]);

void frr(void);

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
  frr();
}

void frr(void)
{
  #define MAXFILENAMELEN 20
  int nSurfaces, nSwaths, nTiles, nFiles;
  int i, j, k , n;
  char lane[] = "6";
  
  nSurfaces = 2;
  nSwaths = 3;
  nTiles = 8;
  
  nFiles = (nSurfaces * nSwaths * nTiles);
  char filterNames[nFiles][MAXFILENAMELEN];
  char bclNames[nFiles][MAXFILENAMELEN];
  n = 0;
  for (i = 1; i <= nSurfaces ; ++i){
    for (j = 1; j <= nSwaths; ++j){
      for (k = 1; k <= nTiles; ++k){         
	sprintf(filterNames[n],"s_%s_%i%i%02d%s",lane,i,j,k,".filter");
	sprintf(bclNames[n],"s_%s_%i%i%02d%s",lane,i,j,k,".bcl");
	n++;
      }
    }
  }
  int nClusters;
  nClusters = getFilterClusters(filterNames[0]);
  unsigned char * filter;
  unsigned char * baseCalls, *bases, *qualities;
  
  filter = malloc(nClusters);
  baseCalls = malloc(nClusters);
  bases = malloc(nClusters);
  qualities = malloc(nClusters);

  int nPasses;
  
  int fi;
  for (fi=0; fi <1;++fi){
    nPasses = getFilterMask(filter, nClusters,filterNames[fi]);
    readBaseCalls(nClusters, baseCalls, bclNames[fi]);
    splitBaseCalls(baseCalls, nClusters, bases, qualities);
  }
  //printArray(qualities, nClusters,"qualities", nClusters);
  
  charHist(qualities,nClusters);
  free(baseCalls);
  free(filter);
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
  //printf("filename: %s, nClusters %d\n",filename,nClusters);
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

void charHist(unsigned char array[], int len)
{
  #define MAXUNIQUE 256
  int i;
  int nUnique;
  unsigned char c;
  unsigned int unique[MAXUNIQUE] = {0};

  nUnique = 0;
  for (i=0;i<len;++i){
    c = array[i];
    if (!unique[c]){
       nUnique += 1;
    }
    unique[c] +=1;
  }
  unsigned char elements[nUnique];
  unsigned int count[nUnique];
  int n;
  n = 0;
  for (i=0; i<MAXUNIQUE;i++){
      if (unique[i]){
        elements[n] = i;
        count[n] = unique[i];
	n++;
     }
  }
  printHist(nUnique, elements, count);
}

 void printHist(int nUnique, char elements[], int count[])
{
  int i;
  for (i=0;i<nUnique;i++){
     printf("element: %d, count %d \n", elements[i],count[i]);
  }
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
  //printf("Filter file:%s, version:%d number of clusters: %d, number of passed clusters %d (%f%%)\n",filename,filterVersion,nClusters,nPassed,(float)nPassed/nClusters *100);
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


void filterBaseCalls(unsigned int nPassed, unsigned int nBases, unsigned char baseCalls[],unsigned char filter[],unsigned char filteredBaseCalls[])
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

