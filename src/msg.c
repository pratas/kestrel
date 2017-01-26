#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                       \n"
  "  -rm <c>:<d>:<i>:<m/e>  reference context model (ex:-rm 13:100:0:0/0), \n"
  "  -rm <c>:<d>:<i>:<m/e>  reference context model (ex:-rm 18:1000:0:1/1000),\n"
  "  ...                                                                  \n"
  "  -tm <c>:<d>:<i>:<m/e>  target context model (ex:-tm 4:1:0:0/0),      \n"
  "  -tm <c>:<d>:<i>:<m/e>  target context model (ex:-tm 18:20:1:2/10),   \n"
  "  ...                                                                  \n"
  "                         target and reference templates use <c> for    \n"
  "                         context-order size, <d> for alpha (1/<d>),    \n"
  "                         <i> (0 or 1) to set the usage of inverted     \n"
  "                         repeats (1 to use) and <m> to the maximum     \n"
  "                         allowed mutation on the context without       \n"
  "                         being discarded (usefull in deep contexts),   \n"
  "                         under the estimator <e>,                      \n");
  } 

void PrintMenu(void){
  fprintf(stderr,
  "Usage: KESTREL [OPTION]... -r [FILE]  [FILE]:[...]                     \n"
  "Compression of genomic sequences from FASTQ files.                     \n"
  "                                                                       \n"
  "Non-mandatory arguments:                                               \n"
  "                                                                       \n"
  "  -h                     give this help,                               \n"
  "  -x                     show several running examples,                \n"
  "  -s                     show KESTREL compression levels,              \n"
  "  -v                     verbose mode (more information),              \n"
  "  -V                     display version number,                       \n"
  "  -f                     force overwrite of output,                    \n"
  "  -l <level>             level of compression [1;9] (lazy -tm setup),  \n"
  "  -g <gamma>             mixture decayment forgetting factor. It is    \n"
  "                         a real value in the interval [0;1),           \n"
  "  -c <cache>             maximum collisions for hash cache. Memory     \n"
  "                         values are higly dependent of the parameter   \n"
  "                         specification,                                \n");
  ModelsExplanation();
  fprintf(stderr,
  "                                                                       \n"
  "  -r <FILE>              reference file (\"-rm\" are loaded here),     \n"
  "                                                                       \n"
  "Mandatory arguments:                                                   \n"
  "                                                                       \n"
  "  <FILE>                 file to compress (last argument). For more    \n"
  "                         files use splitting \":\" characters.         \n"
  "                                                                       \n"
  "Report bugs to <{pratas,ap}@ua.pt>.                                    \n");
  }


void PrintVersion(void){
  fprintf(stderr,
  "                                                                       \n"
  "                            ===============                            \n"
  "                            | KESTREL %u.%u |                          \n"
  "                            ===============                            \n"
  "                                                                       \n"
  "[          compression of genomic sequences from FASTQ files         ] \n"
  "Copyright (C) 2016-2017 University of Aveiro. This is a Free software. \n"
  "You may redistribute copies of it under the terms of the GNU - General \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is NOT \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by \n"
  "Diogo Pratas and Armando J. Pinho.\n\n", VERSION, RELEASE);
  }

