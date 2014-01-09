/*==================================================================================================
 *  MergerRates:
 *
 *   - reads MergerTree's *_mtree files
 *   - counts 'mergers' in each *_mtree file defined via MERGER_RATIO below
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"


#define MERGER_RATIO   0.25           // writes output that readily allows to find mergers

//#define DEBUG


/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t  nprog;
  uint64_t *progid;
  uint64_t *nshared;
  uint64_t *npartprog;
} MTREE;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
uint64_t nhalos;  // number of haloes in *_mtree file
MTREE   *mtree;

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int      merger_rates     (char Filename[MAXSTRING]);
void     read_mtree       (char *infile);


/*==================================================================================================
 * main
 *==================================================================================================*/
int main()
{
  int      i, nFiles, isimu;
  uint64_t ihalo, ipart;
  char   **Filename;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  fprintf(stderr,"==========================================================\n");
  fprintf(stderr,"  count 1:%1d mergers in each of a list of *_mtree files\n",(int)(1./MERGER_RATIO+0.5));
  fprintf(stderr,"==========================================================\n");
  fprintf(stderr,"\nPlease give number of _mtree files:      ");
  scanf("%d", &nFiles);
  fprintf(stderr,"%d\n",nFiles);
  
  /* allocate memory for nFiles filenames, each of size MAXSTRING */
  Filename  = (char **) calloc(nFiles, sizeof(char *));
  for(i=0; i<nFiles; i++) {
    Filename[i]  = (char *) calloc(MAXSTRING, sizeof(char));
  }
  
  /* read input filenames from stdin */
  for(i=0; i<nFiles; i++) {
    fprintf(stderr,"Please give name of %5d. *_mtree file:            ",i+1);
    scanf("%s", Filename[i]);
    fprintf(stderr,"%s\n",Filename[i]);
  }
  fprintf(stderr,"\n");
  
  
  /*======================================================================*
   *  LOOP OVER ALL INPUT FILES                                           *
   *======================================================================*/
  for(i=0; i<nFiles; i++) {
    
    /* be verbose */
    fprintf(stderr,"Working on file '%s_mtree'\n",Filename[i]);

    /* read the next file into memory */
    merger_rates(Filename[i]);
    
  } // for(nFiles)
  
  
  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
  fprintf(stderr,"\nCleaning up ... ");
  
  /* remove filename storage */
  for(i=0; i<nFiles; i++) {
    if(Filename[i])  free(Filename[i]);
  }
  if(Filename)  free(Filename);
  
  printf("finished\n");
  return(1);
}


/*==================================================================================================
 * merger_rates:
 *
 *  get statistics for multiple progenitors (e.g. mass ratios, etc.)
 *
 *==================================================================================================*/
int merger_rates(char Filename[MAXSTRING])
{
  FILE    *fpout_progs, *fpout_merger;
  char     Filename_progs[MAXSTRING], Filename_merger[MAXSTRING];
  uint64_t ihalo;
  
  fprintf(stderr,"  o merger rates: ");
  
  // output files
  strcpy(Filename_progs, basename(Filename));
  strcat(Filename_progs, "_progs");
  strcpy(Filename_merger, basename(Filename));
  strcat(Filename_merger, "_merger");
  
  /* open output file */
  fpout_progs = fopen(Filename_progs,"w");
  if(fpout_progs == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",Filename_progs);
    exit(0);
  }
  
  /* open output file */
  fpout_merger = fopen(Filename_merger,"w");
  if(fpout_merger == NULL)  {
    fprintf(stderr,"could not open file %s\nexiting\n",Filename_merger);
    exit(0);
  }
  fprintf(fpout_progs, "#put a meaningful header here\n");
  fprintf(fpout_merger,"#put a meaningful header here\n");

  /*====================================================================
   *   read *_mtree file
   *   -> allocates memory for mtree[] array of structures from scratch
   *====================================================================*/
  read_mtree(Filename);
  
  
  /*=======================
   *     find mergers
   *=======================*/
#ifdef BE_CAREFUL_WITH_OPENMP
#pragma omp parallel for private(ihalo) shared(nhalos,mtree)
#endif
  for(ihalo=0; ihalo<nhalos; ihalo++) {
    //
    //   TODO
    //
  }
  
  
  /*====================================================================
   *   -> free the memory of the mtree[] array
   *====================================================================*/
  /* remove mtree[] from memory */
  for(ihalo=0; ihalo<nhalos; ihalo++) {
    if(mtree[ihalo].progid != NULL) free(mtree[ihalo].progid);
    if(mtree[ihalo].nshared != NULL) free(mtree[ihalo].nshared);
    if(mtree[ihalo].npartprog != NULL) free(mtree[ihalo].npartprog);
  }
  if(mtree != NULL) free(mtree);
  
  
  fclose(fpout_progs);
  fclose(fpout_merger);
  
  fprintf(stderr," done\n");
  
  return(1);
}

/*==================================================================================================
 * read_mtree:
 *
 *       simply reads in the *_mtree file and 
 *       puts it into the array of structures mtree[ihalo].XYZ
 *
 * Note: at this stage we just treat these entries as "lines" -> no connection to halos yet!
 *
 *==================================================================================================*/
void read_mtree(char *prefix)
{
  uint64_t ihalo, iprog;
  char     line[MAXSTRING], mtreename[MAXSTRING];
  FILE    *fpin;
  
  sprintf(mtreename,"%s_mtree",prefix);
  if((fpin = fopen(mtreename,"r")) == NULL) {
    fprintf(stderr,"cannot open  %s\nEXIT\n",mtreename);
    exit(0);
  }
  
  // ignore first two header lines
  fgets(line,MAXSTRING,fpin);
  fgets(line,MAXSTRING,fpin);
  
  // read rest of file into mtree[]
  mtree  = NULL;
  nhalos = 0;
  
  // read first halo line 
  fgets(line,MAXSTRING,fpin);
  while(!feof(fpin)) {
    // add halo to mtree[]
    nhalos++;
    mtree = (MTREEptr) realloc(mtree, nhalos*sizeof(MTREE));
    
    // scan halo properties
    sscanf(line,"%"SCNi64" %"SCNi64" %"SCNi64,
           &(mtree[nhalos-1].haloid),
           &(mtree[nhalos-1].npart),
           &(mtree[nhalos-1].nprog)  );
    
    // make room for progenitor properties
    mtree[nhalos-1].progid    = (uint64_t *) calloc(mtree[nhalos-1].nprog, sizeof(uint64_t));
    mtree[nhalos-1].nshared   = (uint64_t *) calloc(mtree[nhalos-1].nprog, sizeof(uint64_t));
    mtree[nhalos-1].npartprog = (uint64_t *) calloc(mtree[nhalos-1].nprog, sizeof(uint64_t));
    
    // read progenitor properties
    for(iprog=0; iprog<mtree[nhalos-1].nprog; iprog++) {
      // directly scan from file
      fscanf(fpin,"%"SCNi64" %"SCNi64" %"SCNi64,
             &(mtree[nhalos-1].progid[iprog]),
             &(mtree[nhalos-1].nshared[iprog]),
             &(mtree[nhalos-1].npartprog[iprog])  );
    }
    
    // get next halo
    fgets(line,MAXSTRING,fpin);
  }
  fclose(fpin);
  
  fprintf(stderr,"nhalos=%"PRIu64" ",nhalos);
}
