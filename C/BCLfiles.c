#import <stdio.h>
#import <string.h>


#define MAXFOLDERNAMELEN 8

char *folderNames 


void getfolders(int start,int end, )
{
  int nCycles;
  nCycles = (end+ 1) - start;
  char folderNames[nCycles][MAXFOLDERNAMELEN];
  int f;
  int t;
  for (f = cfvalue, t = 0; f <= cvalue; ++f, ++t){
    printf("%d\n",f);
    sprintf(folderNames[t],"C%d.1/",f);
  }
}