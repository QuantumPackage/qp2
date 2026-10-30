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

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>

/* Total physical memory of the host, in bytes (macOS). */
long long qp_phys_mem_bytes(void)
{
  long long size = 0;
  size_t len = sizeof(size);
  if (sysctlbyname("hw.memsize", &size, &len, NULL, 0) != 0)
    return 0;
  return size;
}

static mach_task_basic_info_data_t qp_task_basic(void)
{
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &count) != KERN_SUCCESS) {
    info.resident_size = 0;
    info.virtual_size = 0;
  }
  return info;
}

/* Resident set size of the current process, in bytes (macOS). */
long long qp_rss_bytes(void)
{
  return (long long) qp_task_basic().resident_size;
}

/* Virtual size of the current process, in bytes (macOS). */
long long qp_vsize_bytes(void)
{
  return (long long) qp_task_basic().virtual_size;
}
#endif

