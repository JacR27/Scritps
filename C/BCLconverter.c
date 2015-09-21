#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* parse, filter, demultiplex, reduse and split, bcl files */

int readHeader(void);
void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], unsigned char bases[], unsigned char qualities[]);
void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[]);
unsigned char * join(unsigned int len, unsigned char array[],int bytesToJoin);
void demultiplex(unsigned int len, unsigned char baseCalls[]);
void reduceQualities(unsigned int len, unsigned char qualities[]);
void printArray(unsigned char array[], unsigned int len, char *fileName, unsigned int header);

main()
{
  unsigned int nClusters;
  unsigned char *baseCalls, *bases, *qualities, *jointBases;

  nClusters = readHeader();
  bases = malloc(nClusters);
  baseCalls = malloc(nClusters);
  qualities = malloc(nClusters);

  readBaseCalls(nClusters, baseCalls, bases, qualities);
  
  jointBases = join(nClusters,bases, 4); 
  

  printArray(jointBases,1616708,"copy",nClusters);
  
  
  
  free(jointBases);
  free(baseCalls);
  free(bases);
  free(qualities);
}


/* read 32 bit little eden int from stdin */
int readHeader(void)
{
  #define BYTESININT32 4
  #define BITSINBYTE 8

  unsigned int nClusters;
  int c , i;
  
  nClusters = 0;

  for (i = 0; i < BYTESININT32; ++i){
    nClusters = nClusters + getchar() * pow(2,i*BITSINBYTE);
  }
  return nClusters;
}


/* read n bytes from standard in to array */
void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], unsigned char bases[], unsigned char qualities[])
{
  int i;
  unsigned char c;
  for (i = 0; i < nBases; ++i){
    c = getchar();
    baseCalls[i] = c;
    bases[i] = (c % 4);
    qualities[i] = (c / 4);
  }
}


/* remove element based on filter file */
void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[])
{
  unsigned char *filter;
  filter = malloc(nBases);
  FILE *fp;
  int i;
  unsigned int nFilteredBases;
  unsigned char c;
  nFilteredBases = 0;
  fp = fopen("s_4_1101.filter","r");
  for (i = 0; i < nBases; ++i){
    c = fgetc(fp);
    nFilteredBases = nFilteredBases + c;
    filter[i] = c;
  }
  fclose(fp);

  int j;
  j = 0;
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
  newLen = len / bytesToJoin + (len % bytesToJoin != 0);
  reduced = malloc(newLen);
  printf("%d",newLen);
  for (i = 0; i < len; ++i){
    reduced[i/bytesToJoin] = reduced[i/bytesToJoin] + array[i] * pow(2,(maxMod-(i%bytesToJoin))*(BITSINBYTE/bytesToJoin));
    //printf("base: %d index: %d newbase: %d\n",array[i],i/bytesToJoin,reduced[i/bytesToJoin]);
    }
  return reduced;
}


void reduceQualities(unsigned int len, unsigned char qualities[])
{
  unsigned char origonalQualities[9] = {0,7,11,22,27,32,37,42};
  unsigned char qualityMap[9] = {0,1,1,2,2,2,3,3};
  unsigned char *reducedQualities;
  int i , j;
  reducedQualities = malloc(len);
  for (i = 0; i < len; ++i){
    for (j = 1; j < 9; ++j){
      if (qualities[i] == origonalQualities[j]){
	reducedQualities[i] = qualityMap[j];
      }
    }
  }
  free(reducedQualities);
}

 void printArray(unsigned char array[], unsigned int len, char *fileName,unsigned int header)
{
  int i;
  FILE *outHandle;
  outHandle = fopen(fileName,"w");
  fwrite(&header,sizeof(header),1,outHandle);
  fwrite(array, 1, sizeof(array), outHandle);  
  fclose(outHandle);
}


