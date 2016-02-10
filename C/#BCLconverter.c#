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

float score2prob(int score);

void qualities2probs(unsigned char qualities[],float probs[],unsigned int len);

int readLittleEInt(FILE *input);

int getFilterClusters(unsigned char filename[]);
 
int readHeader(gzFile input);

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

void frr(void);

void printFileNames(void);

int mode, lvalue, tvalue, wvalue, svalue, cvalue, nFiles, cfvalue;
int fflag, rflag, mflag, jflag, iflag, dflag;
char bclFormat[100];
char printFile[50];
main( int argc, char *argv[])

    {
      extern int mode, lvalue, tvalue, wvalue, svalue, cvalue;
      extern char bclFormat[];
      
      lvalue = 1;
      tvalue = 310;
      wvalue = 1;
      svalue = 1;
      cvalue = 1;
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

      while ((c = getopt (argc, argv, "p:o:l:t:w:s:c:z:frmjid")) != -1)
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
	  case 'p':
	    strcpy(printFile, optarg);
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
      
      if (mode == 5)
	printBCL(printFile);
      else if (mode == 6)
        printcompressedBCL(printFile);
      else
	frr();
      
      return 0;
    }


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
    //printf("%d%c-%d\t%s",(i % 5) ? i : i , con[bases[i]], qualities[i], (i % 5) ? "" : "\n");
    cluster[0] = '\0';
  }
  free(qualities);
  free(bases);
}



void frr(void)
{
  #define MAXFOLDERNAMELEN 8
  printf("starting main\n");
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
    //sprintf(folderNames[t],"C%d.1/",f);
    sprintf(folderNames[t],"./");
  }


  printf("got folders\n");
#define MAXFILENAMELEN 150


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
	  snprintf(filterNames[n],MAXFILENAMELEN,"s_%d.%s",lvalue,"filter");
	  //snprintf(filterNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.%s",lvalue,i,j,k,"filter");
	//snprintf(bclNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.%s",lvalue,i,j,k,bclFormat);
	snprintf(bclNames[n],MAXFILENAMELEN,"%04d.%s",k,bclFormat);
	//snprintf(bclNames[n],MAXFILENAMELEN,"0d.%s",k,bclFormat);
	n++;
      }
    }
  }
  
  //printf("got file names\n");



  

  int nQualites = 8;

  //unsigned char qualityMap[] = {0,1,1,2,2,3,3,3};
  unsigned char qualityMap[] = {0,8,8,23,23,37,37,37};
  //unsigned char qualityMap[] = {0,30,30,30,30,30,30,30};
  //unsigned char qualityMap[] = {0,0,15,25,25,35,35,35};
  //unsigned char qualityMap[] = {0,35,35,35,35,35,35,35};
  //unsigned char qualityMap[] = {0,11,11,11,32,32,32,42};
  unsigned char origonalQualities[] = {0,7,14,22,27,32,37,42};
  //unsigned char qualityMap[] = {0,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,25,25,25,25,25,25,25,25,25,25,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35};
  //unsigned char qualityMap[] = {0,7,7,7,7,7,7,7,7,7,11,11,11,11,11,11,11,22,22,22,22,22,22,22,22,27,27,27,27,27,32,32,32,32,32,37,37,37,37,37,42,42,42,42,42,42,42,42,42,42,42};
  //unsigned char origonalQualities[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50};
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
      char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN];
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
        strcat(outpath, ".packed.gz");
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

      if(jflag){
	outpath[0] = '\0';
        strcpy(outpath, path);
        strcat(outpath, ".packed.gz");	
	joined = join(nPostFiltering, postFilteringBaseCalls, 2);
	printArray(joined, celingDev(nPostFiltering,2),outpath,nPostFiltering);
        free(joined);
      }

     
      if (rflag){
	bases = malloc(nPostFiltering);
	qualities = malloc(nPostFiltering);
	splitBaseCalls(postFilteringBaseCalls, nPostFiltering, bases, qualities);    
        outpath[0] = '\0';
	strcpy(outpath, path);	
        strcat(outpath, ".remapped.gz");
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
  gzFile inputFile;
  inputFile = gzopen(filename,"r");
  unsigned char rawNClusters[BYTESININT32];
  gzread(inputFile, rawNClusters, BYTESININT32);
  nClusters =0;
  for (i = 0; i<BYTESININT32; ++i){
    nClusters = nClusters + rawNClusters[i] * pow(2,i*BITSINBYTE);
  }
  printf("filename: %s, nClusters %d\n",filename,nClusters);
  gzread(inputFile, baseCalls, nBases);
  gzclose(inputFile);
}

void readReducedBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[])
{
  #define BYTESININT32 4
  #define BITSINBYTE 8
  char * buff;
  int i;
  unsigned int nBytes;
  nBytes =celingDev(nBases,2); 
  unsigned char c;
  unsigned char c1;
  unsigned char c2;
  gzFile inputFile;
  inputFile = gzopen(filename,"r");
  buff = malloc(nBytes);
  unsigned char rawNClusters[BYTESININT32];
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
    for (j = 0; j < nQualities; ++j){
      if (qualities[i] == origonalQualities[j]){
	reducedQualities[i] = qualityMap[j];
      }
    }
  }
}

 void printArray(unsigned char array[], unsigned int len, char *fileName,unsigned int header)
{
  int i;
  gzFile outHandle;
  outHandle = gzopen(fileName,"wb6h");
  gzwrite(outHandle, &header,sizeof(header));
  gzwrite(outHandle, array, len);
  gzclose(outHandle);
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
  readLittleEInt(fp);
  filterVersion = readLittleEInt(fp);
  nClusters = readLittleEInt(fp);
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
  gzFile fp;
  fp = gzopen(filename,"r");
  nClusters = readHeader(fp);
  gzclose(fp);
  //nClusters = 69;
  return nClusters;
}


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


