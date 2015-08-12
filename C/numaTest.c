#include <stdio.h>
#include <numa.h>
main()
{
  int numa;
  numa = numa_max_node();
  printf("%d\n",numa);
}
