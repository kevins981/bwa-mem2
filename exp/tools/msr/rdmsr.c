#ident "$Id$"
/* ----------------------------------------------------------------------- *
 *   
 *   Copyright 2000 Transmeta Corporation - All Rights Reserved
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, Inc., 675 Mass Ave, Cambridge MA 02139,
 *   USA; either version 2 of the License, or (at your option) any later
 *   version; incorporated herein by reference.
 *
 * ----------------------------------------------------------------------- */

/*
 * rdmsr.c
 *
 * Utility to read an MSR.
 */

#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/types.h>

#include "version.h"

struct option long_options[] = {
  { "help",        0, 0, 'h' },
  { "version",     0, 0, 'v' },
  { "hexadecimal", 0, 0, 'x' },
  { "capital-hexadecimal", 0, 0, 'X' },
  { "decimal",     0, 0, 'd' },
  { "unsigned",    0, 0, 'u' },
  { "octal",       0, 0, 'o' },
  { "c-language",  0, 0, 'c' }, 
  { "zero-fill",   0, 0, '0' },
  { "raw",         0, 0, 'r' },
  { "processor",   1, 0, 'p' },
  { "cpu",         1, 0, 'p' },
  { 0, 0, 0, 0 }
};

#define mo_hex  0x01
#define mo_dec  0x02
#define mo_oct  0x03
#define mo_raw  0x04
#define mo_uns  0x05
#define mo_chx  0x06
#define mo_mask 0x0f
#define mo_fill 0x40
#define mo_c    0x80

const char *program;

void usage(void)
{
  fprintf(stderr, "Usage: %s [-p processor] [-hxXdoruc0] regno\n", program);
}

int main(int argc, char *argv[])
{
  uint32_t reg;
  uint64_t data;
  int c, fd;
  int mode = mo_hex;
  int cpu = 0;
  unsigned long arg;
  char *endarg;
  char *pat;
  char msr_file_name[64];

  program = argv[0];

  while ( (c = getopt_long(argc,argv,"hvxXdoruc0p:",long_options,NULL)) != -1 ) {
    switch ( c ) {
    case 'h':
      usage();
      exit(0);
    case 'v':
      fprintf(stderr, "%s: version %s\n", program, VERSION_STRING);
      exit(0);
    case 'x':
      mode = (mode & ~mo_mask) | mo_hex;
      break;
    case 'X':
      mode = (mode & ~mo_mask) | mo_chx;
      break;
    case 'o':
      mode = (mode & ~mo_mask) | mo_oct;
      break;
    case 'd':
      mode = (mode & ~mo_mask) | mo_dec;
      break;
    case 'r':
      mode = (mode & ~mo_mask) | mo_raw;
      break;
    case 'u':
      mode = (mode & ~mo_mask) | mo_uns;
      break;
    case 'c':
      mode |= mo_c;
      break;
    case '0':
      mode |= mo_fill;
      break;
    case 'p':
      arg = strtoul(optarg, &endarg, 0);
      if ( *endarg || arg > 255 ) {
	usage();
	exit(127);
      }
      cpu = (int)arg;
      break;
    default:
      usage();
      exit(127);
    }
  }

  if ( optind != argc-1 ) {
    /* Should have exactly one argument */
    usage();
    exit(127);
  }

  reg = strtoul(argv[optind], NULL, 0);

  sprintf(msr_file_name, "/dev/cpu/%d/msr", cpu);
  fd = open(msr_file_name, O_RDONLY);
  if ( fd < 0 ) {
    if ( errno == ENXIO ) {
      fprintf(stderr, "rdmsr: No CPU %d\n", cpu);
      exit(2);
    } else if ( errno == EIO ) {
      fprintf(stderr, "rdmsr: CPU %d doesn't support MSRs\n", cpu);
      exit(3);
    } else {
      perror("rdmsr:open");
      exit(127);
    }
  }
  
  if ( lseek(fd, reg, SEEK_SET) != reg ) {
    perror("rdmsr:seek");
    exit(127);
  }

  if ( read(fd, &data, sizeof data) != sizeof data ) {
    perror("rdmsr:read");
    exit(127);
  }

  close(fd);

  pat = NULL;

  switch(mode) {
  case mo_hex:
    pat = "%llx\n";
    break;
  case mo_chx:
    pat = "%llX\n";
    break;
  case mo_dec:
  case mo_dec|mo_c:
  case mo_dec|mo_fill|mo_c:
    pat = "%lld\n";
    break;
  case mo_uns:
    pat = "%llu\n";
    break;
  case mo_oct:
    pat = "%llo\n";
    break;
  case mo_hex|mo_c:
    pat = "0x%llx\n";
    break;
  case mo_chx|mo_c:
    pat = "0x%llX\n";
    break;
  case mo_oct|mo_c:
    pat = "0%llo\n";
    break;
  case mo_uns|mo_c:
  case mo_uns|mo_fill|mo_c:
    pat = "%lluU\n";
    break;
  case mo_hex|mo_fill:
    pat = "%016llx\n";
    break;
  case mo_chx|mo_fill:
    pat = "%016llX\n";
    break;
  case mo_dec|mo_fill:
    pat = "%020lld\n";
    break;
  case mo_uns|mo_fill:
    pat = "%020llu\n";
    break;
  case mo_oct|mo_fill:
    pat = "%022llo\n";
    break;
  case mo_hex|mo_fill|mo_c:
    pat = "0x%016llx\n";
    break;
  case mo_chx|mo_fill|mo_c:
    pat = "0x%016llX\n";
    break;
  case mo_oct|mo_fill|mo_c:
    pat = "0%022llo\n";
    break;
  case mo_raw:
  case mo_raw|mo_fill:
    fwrite(&data,sizeof data,1,stdout);
    break;
  case mo_raw|mo_c:
  case mo_raw|mo_fill|mo_c:
    {
      unsigned char *p = (unsigned char *)&data;
      int i;
      for ( i = 0 ; i < sizeof data ; i++ ) {
	printf("%s0x%02x", i?",":"{", (unsigned int)(*p++));
      }
      printf("}\n");
    }
  break;
  default:
    fprintf(stderr, "%s: Impossible case, line %d\n", program, __LINE__);
    exit(127);
  }

  if ( pat )
    printf(pat, data);

  exit(0);
}
