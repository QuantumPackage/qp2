#include <unistd.h>

void usleep_c(int s)
{
  usleep((useconds_t) s);
}
