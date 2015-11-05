int celingDev(int dividend, int devisor);

int readHeader(gzFile input);

unsigned char * interleafe(unsigned char array1[], unsigned char array2[], unsigned int len);

void readBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[]);

void readReducedBaseCalls(unsigned int nBases, unsigned char baseCalls[], char filename[]);

void splitBaseCalls(unsigned char baseCalls[], unsigned int nBases, unsigned char bases[], unsigned char qualities[]);

void charHist(unsigned long unique[] , unsigned char array[], int len);

void printHist(unsigned long unique[]);

unsigned char * join(unsigned int len, unsigned char array[], int bytesToJoin);

void rejoin(unsigned char *reduced, unsigned int len, unsigned char qualities[],unsigned char bases[]);

void reduceQualities(unsigned char reducedQualities[],unsigned char origonalQualities[], unsigned int len, unsigned char qualities[],unsigned char qualityMap[],int nQualities);

void printArray(unsigned char array[], unsigned int len, char *fileName,unsigned int header);

void printArrayStdout(unsigned char array[], unsigned int len, unsigned int header);

int getFilterMask(unsigned char filter[], unsigned int len, unsigned char filename[]);

int getClusters(unsigned char filename[]);

int getFilterClusters(unsigned char filename[]);

int getFilterPasses(unsigned char filename[]);

void filterBaseCalls(unsigned int nBases, unsigned char baseCalls[], unsigned char filter[], unsigned char filteredBaseCalls[]);

int readLittleEInt(FILE *input);

float score2prob(int score);

void qualities2probs(unsigned char qualities[],float probs[],unsigned int len);

void printBCL(char filename[]);

void printcompressedBCL(char filename[]);
