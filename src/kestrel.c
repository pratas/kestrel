#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>

#include "mem.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "models.h"
#include "pmodels.h"
#include "stream.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - M O D E L S - - - - - - - - - - - - - - - -

CModel **Models;   // MEMORY SHARED BY THREADING


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - R E S E T   M O D E L S - - - - - - - - - - - -

void ResetModelsAndParam(CBUF *Buf, CModel **Shadow, CMWeight *CMW){
  uint32_t n;
  ResetCBuffer(Buf);
  for(n = 0 ; n < P->nModels ; ++n)
    ResetShadowModel(Shadow[n]);
  ResetWeightModel(CMW);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I O N - - - - - - - - - - - - - 

void CompressTarget(Threads T){
  FILE        *Reader = Fopen(P->base, "r");
  double      bits = 0;
  uint64_t    nBase = 0, r = 0, nSymbol, initNSymbol;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos, conName[MAX_NAME];
  PModel      **pModel, *MX;
  CModel      **Shadow; // SHADOWS FOR SUPPORTING MODELS WITH THREADING
  FloatPModel *PT;
  CMWeight    *CMW;
  int         action;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow      = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]); 
  pModel      = (PModel **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n] = CreatePModel(ALPHABET_SIZE);
  MX          = CreatePModel(ALPHABET_SIZE);
  PT          = CreateFloatPModel(ALPHABET_SIZE);
  CMW         = CreateWeightModel(totModels);

  FileType(PA, Reader);
  if(PA->type != 2){
    fprintf(stderr, "Error: input file must be in FASTQ format!\n");
    exit(1);
    }

  char name_o[4096];
  sprintf(name_o, "%s.%u", P->output, T.id);
  FILE *Writer = Fopen(name_o, "w");

  srand(T.id);
  initNSymbol = nSymbol = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      ++nSymbol;

      sym = readBuf[idxPos];

      if(sym == '@'){
        if((PA->nRead-1) % P->nThreads == T.id && PA->nRead > 1 && nBase > 1){
          //fprintf(stderr, "%g\n", bits / (2.0 * nBase));
          if(bits / (2.0 * nBase) > P->threshold)
            fprintf(Writer, "1"); // WRITE READ
          else
            fprintf(Writer, "0"); // IGNORE READ
          ResetModelsAndParam(symBuf, Shadow, CMW); // RESET MODELS
          r = nBase = bits = 0;
          }
        PA->nRead++;
        continue;
        }

      switch(PA->line){
        case 0: if(sym == '\n'){ PA->line = 1; PA->dna = 1; } break;
        case 1: if(sym == '\n'){ PA->line = 2; PA->dna = 0; } break;
        case 2: if(sym == '\n'){ PA->line = 3; PA->dna = 0; } break;
        case 3: if(sym == '\n'){ PA->line = 0; PA->dna = 0; } break;
        }
      if(PA->dna == 0 || sym == '\n') continue;

      if(PA->nRead % P->nThreads == T.id){

        if(sym == 'N') sym = rand() % 4;
        else           sym = DNASymToNum(sym);

        symBuf->buf[symBuf->idx] = sym;
        memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));
        n = 0;
        pos = &symBuf->buf[symBuf->idx-1];
        for(cModel = 0 ; cModel < P->nModels ; ++cModel){
          CModel *CM = Shadow[cModel];
          GetPModelIdx(pos, CM);
          ComputePModel(Models[cModel], pModel[n], CM->pModelIdx, CM->alphaDen);
          ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
          if(CM->edits != 0){
            ++n;
            CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
            CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+CM->SUBS.seq->idx
            -1, CM, CM->SUBS.idx);
            ComputePModel(Models[cModel], pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
            ComputeWeightedFreqs(CMW->weight[n], pModel[n], PT);
            }
          ++n;
          }

        ComputeMXProbs(PT, MX);
        bits += PModelSymbolLog(MX, sym);
        ++nBase;
        CalcDecayment(CMW, pModel, sym, P->gamma);
        RenormalizeWeights(CMW);
        CorrectXModels(Shadow, pModel, sym);
        UpdateCBuffer(symBuf);
        }
      }
        
  if(PA->nRead % P->nThreads == T.id){
    // ?
    }

  DeleteWeightModel(CMW);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  fclose(Writer);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReference(char *refName){
  FILE     *Reader = Fopen(refName, "r");
  uint32_t n;
  uint64_t idx = 0;
  uint64_t k, idxPos;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf  = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  *readBuf = Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t  sym, irSym = 0;
  FileType(PA, Reader);
  rewind(Reader);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1){ idx = 0; continue; }
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);
      for(n = 0 ; n < P->nModels ; ++n){
        CModel *CM = Models[n];
        GetPModelIdx(symBuf->buf+symBuf->idx-1, CM);
        if(CM->ir == 1) // INVERTED REPEATS
          irSym = GetPModelIdxIR(symBuf->buf+symBuf->idx, CM);
        if(++idx >= CM->ctx){
          UpdateCModelCounter(CM, sym, CM->pModelIdx);
          if(CM->ir == 1) // INVERTED REPEATS
            UpdateCModelCounter(CM, irSym, CM->pModelIdxIR);
          }
        }
      UpdateCBuffer(symBuf);
      }
 
  for(n = 0 ; n < P->nModels ; ++n)
    ResetCModelIdx(Models[n]);
  RemoveCBuffer(symBuf);
  Free(readBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - C O M P R E S S O R   M A I N - - - - - - - - - - - -

void CompressAction(Threads *T, char *refName, char *baseName){
  pthread_t t[P->nThreads];
  uint32_t n;
 
  Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Models[n] = CreateCModel(T[0].model[n].ctx, T[0].model[n].den, 
    T[0].model[n].ir, REFERENCE, P->col, T[0].model[n].edits, 
    T[0].model[n].eDen);
  fprintf(stderr, "  [+] Loading reference file ....... ");
  LoadReference(refName);
  fprintf(stderr, "Done!\n");

  fprintf(stderr, "  [+] Compressing FASTQ reads ...... ");
  for(n = 0 ; n < P->nThreads ; ++n)
    pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[n]));
  for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
    pthread_join(t[n+1], NULL);
  fprintf(stderr, "Done!\n");
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char     **p = *&argv, **xargv, *xpl = NULL;
  int32_t  xargc = 0;
  uint32_t n, k, col, ref;
  double   gamma;
  Threads  *T;
  
  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s")){
    PrintLevels(); 
    return EXIT_SUCCESS;
    }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v");
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F");
  P->sample   = ArgsNum    (DEFAULT_SAMPLE,  p, argc, "-p", MIN_SAP, MAX_SAP);
  P->level    = ArgsNum    (0,               p, argc, "-l", MIN_LEV, MAX_LEV);
  P->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-n", MIN_THREADS, 
  MAX_THREADS);

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-g") == 0) 
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-c") == 0) 
      col = atoi(xargv[n+1]);

  P->col       = ArgsNum    (col,   p, argc, "-c", 1, 253);
  P->gamma     = ArgsDouble (gamma, p, argc, "-g");
  P->threshold = fabs(ArgsDouble (0.9,   p, argc, "-t"));
  P->gamma     = ((int)(P->gamma * 65536)) / 65536.0;
  P->output    = ArgsFileGen(p, argc, "-o", "filtered", ".fq");

  FILE *OUTPUT = NULL;
  if(!P->force) 
    FAccessWPerm(P->output);
  OUTPUT = Fopen(P->output, "w");

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return EXIT_FAILURE;
    }

  // READ MODEL PARAMETERS FROM XARGS & ARGS
  T = (Threads *) Calloc(P->nThreads, sizeof(Threads));
  for(ref = 0 ; ref < P->nThreads ; ++ref){
    T[ref].model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));
    T[ref].id    = ref;
    k = 0;
    for(n = 1 ; n < argc ; ++n)
      if(strcmp(argv[n], "-m") == 0)
        T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
    if(P->level != 0){
      for(n = 1 ; n < xargc ; ++n)
        if(strcmp(xargv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
      }
    }

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P, T[0], argv[argc-2], argv[argc-1]);

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());
  P->base = argv[argc-1];
  CompressAction(T, argv[argc-2], P->base);

  fprintf(stderr, "  [+] Freeing compression models ... ");
  for(n = 0 ; n < P->nModels ; ++n)
    FreeCModel(Models[n]);
  Free(Models);
  fprintf(stderr, "Done!\n");

  StopTimeNDRM(Time, clock());
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveClock(Time);
  for(ref = 0 ; ref < P->nThreads ; ++ref){
    Free(T[ref].model);
    }
  Free(T);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
