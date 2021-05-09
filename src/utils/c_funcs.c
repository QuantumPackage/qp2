#include <unistd.h>
#include <stdio.h>
#include <string.h>

void usleep_c(int s)
{
  usleep((useconds_t) s);
}

void sscanf_ssds_c(const char* str, char* s1, char* s2, int* i, char* s3)
{
  sscanf(str, "%s %s %d %s", s1, s2, i, s3);
  s1[strlen(s1)] = ' ';
  s2[strlen(s2)] = ' ';
  s3[strlen(s3)] = ' ';
}

void sscanf_dd_c(const char* str, int* i1, int* i2)
{
  sscanf(str, "%d %d", i1, i2);
}

void sscanf_ddd_c(const char* str, int* i1, int* i2, int* i3)
{
  sscanf(str, "%d %d %d", i1, i2, i3);
}

void sscanf_ss_c(const char* str, char* s1, char* s2)
{
  sscanf(str, "%s %s", s1, s2);
  s1[strlen(s1)] = ' ';
  s2[strlen(s2)] = ' ';
}

void sscanf_sd_c(const char* str, char* s1, int* i)
{
  sscanf(str, "%s %d", s1, i);
  s1[strlen(s1)] = ' ';
}

