#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
#include <ctype.h>
#include <unistd.h>



/*celing Devide*/
int celingDev(int dividend, int devisor)
{
  return ((dividend / devisor) + (dividend % devisor != 0));
}



/* read 32 bit little eden int from file */
int readHeader(gzFile input)
{
  #define BYTESININT32 4
  #define BITSINBYTE 8
  unsigned int nClusters;
  int c , i;
  nClusters = 0;
  for (i = 0; i<BYTESININT32; ++i){
    nClusters = nClusters + gzgetc(input) * pow(2,i*BITSINBYTE);
  }
  return nClusters;
}



/*interleafe two arrays of the same size */
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



/* read basecalls  */
void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[])
{
  #define BYTESININT32 4
  #define BITSINBYTE 8
  int i;
  int nClusters;
  unsigned char c;
  gzFile inputFile;
  inputFile = gzopen(filename,"r");
  nClusters = readHeader(inputFile);
  gzread(inputFile, baseCalls, nBases);  
  gzclose(inputFile);
  printf("filename: %s, nClusters %d\n",filename,nClusters);
}



/* read 4 bit basescall into 8 bit basecalls*/
void readReducedBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[])
{
  #define BYTESININT32 4
  #define BITSINBYTE 8
  char * buff;
  int i;
  unsigned int nBytes;
  unsigned char c;
  gzFile inputFile;
  unsigned char rawNClusters[BYTESININT32];
  
  nBytes =celingDev(nBases,2); 
  inputFile = gzopen(filename,"r");
  buff = malloc(nBytes);
  gzread(inputFile, rawNClusters, BYTESININT32);
  gzread(inputFile, buff, nBytes);
  
  for (i=0;i<nBases;i++){
    c = buff[i/2];
    if (i % 2)
      baseCalls[i] = c%16;
    else
      baseCalls[i] = c/16;
  }
  
  gzclose(inputFile);
  free(buff);
}



/* get bases and qualities for array of basecalls*/
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



/* make a historgram of chars from char array*/
void charHist(unsigned long unique[] , unsigned char array[], int len)
{
  #define MAXUNIQUE 256
  int i;
  unsigned char c;
  for (i=0;i<len;++i){
    c = array[i];
    unique[c] +=1;
  }
}



/* print char historgram*/
void printHist(unsigned long unique[])
{
  #define MAXUNIQUE 256
  unsigned long i, c, tot;
  tot = 0;
  for (i=0;i<MAXUNIQUE;i++){
    if (c = unique[i]){
      printf("number: %lu, count %lu \n", i, c);
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

  for (i = 0; i < len; ++i){
    reduced[i/bytesToJoin] += array[i] * pow(2,(maxMod-(i%bytesToJoin))*(BITSINBYTE/bytesToJoin));
    }

  return reduced;
}



/*rejoin bases and qualities*/
void rejoin(unsigned char *reduced, unsigned int len, unsigned char qualities[],unsigned char bases[]){
  int i;
  for (i = 0; i < len; ++i)
    reduced[i] = (qualities[i] << 2) + bases[i];
}



/*remap qualities*/
void reduceQualities(unsigned char reducedQualities[],unsigned char origonalQualities[], unsigned int len, unsigned char qualities[],unsigned char qualityMap[],int nQualities)
{
  int i , j;  
  for (i = 0; i < len; ++i){
    for (j = 0; j < nQualities; ++j){
      if (qualities[i] == origonalQualities[j]){
	reducedQualities[i] = qualityMap[j];
      }
    }
  }
}



/* outbut bcl to zlib file*/
void printArray(unsigned char array[], unsigned int len, char *fileName,unsigned int header)
{
  int i;
  gzFile outHandle;
  outHandle = gzopen(fileName,"w");
  gzwrite(outHandle, &header,sizeof(header));
  gzwrite(outHandle, array, len);
  gzclose(outHandle);
}

/* outbut bcl to stdout file*/
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
  FILE *fp; 
  fp = fopen(filename,"r");
  readLittleEInt(fp);
  filterVersion = readLittleEInt(fp);
  nClusters = readLittleEInt(fp);
  nPassed = 0;
  fread(filter,1,nClusters,fp);
  
  for (i = 0; i < nClusters; ++i){
    nPassed += filter[i];
  }

  printf("Filter file:%s, version:%d number of clusters: %d, number of passed clusters %d (%f%%)\n",filename,filterVersion,nClusters,nPassed,(float)nPassed/nClusters *100);
  fclose(fp);
  
  return nPassed;
}



/*get number of clusters from a bcl file*/
int getClusters(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  gzFile fp;
  fp = gzopen(filename,"r");
  nClusters = readHeader(fp);
  gzclose(fp);
  //nClusters = 69;
  return nClusters;
}


/*get number of clusters from a filter file*/
int getFilterClusters(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  FILE *fp;
  fp = fopen(filename,"r");
  readLittleEInt(fp);
  readLittleEInt(fp);
  nClusters = readLittleEInt(fp);
  fclose(fp);
  return nClusters;
}

/*get filter Passes*/
int getFilterPasses(unsigned char filename[])
{
  int nClusters; /* number contained in header of filter file */
  int nPassed; /* number of passed bases */
  int i; /* loop variable */
  unsigned char c; /*temp var*/
  FILE *fp;
  fp = fopen(filename,"r");
  readLittleEInt(fp);
  readLittleEInt(fp);
  nClusters = readLittleEInt(fp);
  nPassed = 0;
  for (i = 0; i < nClusters; ++i){
    c = fgetc(fp);
    nPassed += c;
  }
  fclose(fp);
  return nPassed;
}



/*filter basecalls*/
void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[], unsigned char filter[], unsigned char filteredBaseCalls[])
{
  int j , i, nPasses;
  int nFilteredBases;
  j = 0;
  nPasses = 0;
  printf("filtering\n");
  for (i = 0; i < nBases; ++i){
    if (filter[i] == 1){
      ++nPasses;
      filteredBaseCalls[j] = baseCalls[i];
      j++;
    }
  }
  printf("number of passed filtering: %d %d", nPasses, j); 
}



int readLittleEInt(FILE *input)
{
  #define BYTESININT32 4
  #define BITSINBYTE 8
  unsigned int nClusters;
  int c , i;
  nClusters = 0;
  for (i = 0; i<BYTESININT32; ++i){
    nClusters = nClusters + getc(input) * pow(2,i*BITSINBYTE);
  }
  return nClusters;
}



float score2prob(int score)
{
  return pow(10,(-score)/10.);
}



void qualities2probs(unsigned char qualities[],float probs[],unsigned int len)
{
  int i;
  for (i=0;i < len;++i){
    probs[i] = score2prob(qualities[i]);
  }
}


