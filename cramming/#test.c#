#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#define MAXNAMELEN 100
main()
{
  while (1){
    int i; 
    int32_t block_size;
    i = fread(&block_size, sizeof(int32_t), 1, stdin);
    //block_size = getchar();
    printf("%d\n%d\n%d\n",block_size,sizeof(int32_t),i);
    if(feof(stdin))
      break;
  }
}
