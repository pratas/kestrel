#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                         \n"
  "  -m <c>:<d>:<i>:<m/e>  reference context model (ex:-m 13:100:0:0/0),    \n"
  "  -m <c>:<d>:<i>:<m/e>  reference context model (ex:-m 20:1000:1:1/10),  \n"
  "  ...                                                                    \n"
  "                        templates use <c> for context-order size, <d> for\n"
  "                        alpha (1/<d>), <i> (0 or 1) to set the usage of  \n"
  "                        inverted repeats (set 1 to use) and <m> to the   \n"
  "                        maximum allowed mutation on the context without  \n"
  "                        being discarded (usefull in deep contexts), under\n"
  "                        the estimator <e>.                               \n"
  "                                                                         \n");
  }

void PrintMenu(void){
  fprintf(stderr,
  "Usage: KESTREL [OPTION]... [FILE1] [FILE2]                               \n"
  "A compression method to filter FASTQ reads by relative similarity.       \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                       give this help,                               \n"
  "  -F                       force mode (overwrites top file),             \n"
  "  -V                       display version number,                       \n"
  "  -v                       verbose mode (more information),              \n"
  "  -s                       show compression levels,                      \n"
  "  -l <level>               compression level [%u;%u],                    \n"
  "  -t <threshold>           get reads higher than threshold [0.0;1.0],    \n"
  "  -p <sample>              subsampling (default: %u),                    \n"
  "  -n <nThreads>            number of threads (default: %u),              \n"
  "  -o <FILE>                output filename,                              \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  [FILE1]                  reference filename (FASTA or FASTQ),          \n"
  "  [FILE2]                  raw reads filename (FASTQ).                   \n"
  "                                                                         \n"
  "Report issues to <{pratas,ap}@ua.pt>.                                    \n",
  MIN_LEV, MAX_LEV, DEFAULT_SAMPLE, DEFAULT_THREADS);
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                         \n"
  "                          ===================                            \n"
  "                          |   KESTREL %u.%u   |                          \n"
  "                          ===================                            \n"
  "                                                                         \n"
  "                     A compression method to filter                      \n"
  "                   FASTQ reads by relative similarity.                   \n"
  "                                                                         \n"
  "              Copyright (C) 2014-2017 University of Aveiro.              \n"
  "                                                                         \n"
  "                  This is a Free software, under GPLv3.                  \n"
  "                                                                         \n"
  "You may redistribute copies of it under the terms of the GNU - General   \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is NOT   \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by   \n"
  "Diogo Pratas and Armando J. Pinho.\n\n", VERSION, RELEASE);
  }

