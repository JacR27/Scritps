#include <stdio.h>
#include <numa.h>
main()
{
  int numa;
  numa = numa_available();
  printf("%d\n",numa);
}
