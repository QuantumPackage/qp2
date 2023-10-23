#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>


void* mmap_fortran(char* filename, size_t bytes, int* file_descr, int read_only)
{
    int fd;
    int result;
    void* map;

    if (read_only == 1)
    {
        fd = open(filename, O_RDONLY, (mode_t)0600);
        if (fd == -1) {
            printf("%s:\n", filename);
            perror("Error opening mmap file for reading");
            exit(EXIT_FAILURE);
        }
        map = mmap(NULL, bytes, PROT_READ, MAP_SHARED, fd, 0);
    }
    else
    {
        fd = open(filename, O_RDWR | O_CREAT, (mode_t)0600);
        if (fd == -1) {
            printf("%s:\n", filename);
            perror("Error opening mmap file for writing");
            exit(EXIT_FAILURE);
        }

        result = lseek(fd, bytes+1, SEEK_SET);
        if (result == -1) {
            close(fd);
            printf("%s:\n", filename);
            perror("Error calling lseek() to stretch the file");
            exit(EXIT_FAILURE);
        }
        
        result = write(fd, "", 1);
        if (result != 1) {
            close(fd);
            printf("%s:\n", filename);
            perror("Error writing last byte of the file");
            exit(EXIT_FAILURE);
        }

        map = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    }

    if (map == MAP_FAILED) {
        close(fd);
        printf("%s: %lu\n", filename, bytes);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }

    *file_descr = fd;
    return map;
}

void munmap_fortran(size_t bytes, int fd, void* map)
{
    if (munmap(map, bytes) == -1) {
        perror("Error un-mmapping the file");
    }
    close(fd);
}


void msync_fortran(size_t bytes, int fd, void* map)
{
    if (msync(map, bytes, MS_SYNC) == -1) {
        perror("Error syncing the mmap file");
    }
}

