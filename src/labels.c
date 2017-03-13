#include "labels.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SLABELS *CreateSLabels(void){
  uint32_t n;
  SLABELS *SL = (SLABELS *) Calloc(1, sizeof(SLABELS));
  SL->idx     = 0;
  SL->maxV    = SLCACHE;
  SL->maxH    = SLMAXSTR;
  SL->names   = (uint8_t **) Calloc(SL->maxV, sizeof(uint8_t *));
  for(n = 0 ; n < SL->maxV ; ++n)
    SL->names[n] = (uint8_t *) Calloc(SL->maxH+1, sizeof(uint8_t));
  return SL;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void AddSLabel(SLABELS *SL, char *str){
  uint32_t n;
  char c;
  for(n = 0 ; n < SL->maxH ; ++n){
    if((c = str[n]) == '\0')
      break;
    SL->names[SL->idx][n] = tolower(c);  
    }
  SL->names[SL->idx][n] = '\0';
  SL->idx++;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateSLabels(SLABELS *SL){
  uint32_t n;
  if(SL->idx == SL->maxV){
    SL->maxV += SLCACHE;
    SL->names = (uint8_t **) 
    Realloc(SL->names, SL->maxV * sizeof(uint8_t *), 0);
    for(n = SL->idx ; n < SL->maxV ; ++n)
      SL->names[n] = (uint8_t *) Calloc(SL->maxH+1, sizeof(uint8_t));
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int SearchSLabels(SLABELS *SL, char *str){
  uint32_t n;
  for(n = 0 ; n < SL->idx ; ++n)
    if(Strcasecmp(SL->names[n], str) == 0)
      return 1;
  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void DeleteSLabels(SLABELS *SL){
  uint32_t n;
  for(n = 0 ; n < SL->maxV ; ++n)
    Free(SL->names[n]);
  Free(SL->names);
  Free(SL);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
