/*==================================================================================================
 *  MergerTree:   Merger Tree AHF_particles files
 *
 *
 *  input:    - how often to perform
 *            - 2x _particles files
 *
 *  output:   - 1x _mtree file
 *
 *
 * it is checked what halos in file2 make up the halos in file1, i.e.
 *
 *   file1   file2
 *
 *    0        0
 *    0       17
 *    0       31    -> halo #0 in file1 shares particles with halos #0,17,31 in file2
 *    1        2
 *    1       12
 *    1        4    -> halo #1 in file1 shares particles with halos #2,12,4  in file2
 *       etc.
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

#define MINCOMMON      10             // we only cross-correlate haloes if they at least share MINCOMMON particles
#define ONLY_USE_PTYPE 1              // restrict analysis to particles of this type (1 = dark matter)
//#define EXCLUSIVE_PARTICLES           // each particle is only allowed to belong to one object (i.e. the lowest mass one)
//#define MTREE_BOTH_WAYS               // make sure that every halo has only one descendant
//#define SUSSING2013                   // write _mtree in format used for Sussing Merger Trees 2013
//#define USE_LINENUMBER_AS_HALOID      // do not use the haloid as found in _particles
//#define MERGER_RATIO   0.25           // writes output that readily allows to find mergers

//#define DEBUG

#if (defined EXCLUSIVE_PARTICLES && defined WITH_OPENMP)
#define WITH_OPENMP2
#endif

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  uint64_t haloid[2];
  uint64_t id[2];
  uint64_t npart[2];
  uint64_t common;
} MTREE;

typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t *Pid;
  
  uint64_t  ncroco;
  MTREEptr  mtree;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  uint64_t  nhalos;
  uint64_t *Hid;
}PARTS;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

HALOptr halos[2];
PARTptr parts[2];
uint64_t    nHalos[2];
uint64_t    PidMax[2]={0,0};
uint64_t    PidMin=(1<<62);


#ifdef MERGER_RATIO
uint64_t nlines;  // number of lines in *_mtree file
#endif

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int      read_particles         (char filename[MAXSTRING], int isimu);
int      particle_halo_mapping  (int  isimu);
int      cross_correlation      (char OutFile[MAXSTRING]);
int      create_mtree           (uint64_t ihalo, int  isimu0,int isimu1);
int      clean_connection       (uint64_t ihalo, int isimu0, int isimu1);
int      write_mtree            (char OutFile[MAXSTRING]);
uint64_t max_merit              (uint64_t ihalo, int isimu);

#ifdef MERGER_RATIO
int  assign_progenitors     (char OutFile[MAXSTRING]);
void read_mtree             (char *infile);
#endif

/*==================================================================================================
 * main:
 *
 *       simply a wrapper for successive calls to create_mtree()
 *
 *==================================================================================================*/
int main()
{
  int      i, nFiles, isimu;
  uint64_t ihalo, ipart;
  char   **HaloFile;
  char   **OutFile;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"  construct a cross-correlation between consecutive *_particles files\n");
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"\nPlease give number of particles files (default=2):      ");
  scanf("%d", &nFiles);
  fprintf(stderr,"%d\n",nFiles);
  
  /* allocate memory for nFiles filenames, each of size MAXSTRING */
  HaloFile = (char **) calloc(nFiles, sizeof(char *));
  OutFile  = (char **) calloc(nFiles, sizeof(char *));
  for(i=0; i<nFiles; i++)
   {
    HaloFile[i] = (char *) calloc(MAXSTRING, sizeof(char));
    OutFile[i]  = (char *) calloc(MAXSTRING, sizeof(char));
   }
  
  /* read input filenames from stdin */
  for(i=0; i<nFiles; i++)
   {
    fprintf(stderr,"Please give name of %5d. *_particles file:            ",i+1);
    scanf("%s", HaloFile[i]);
    fprintf(stderr,"%s\n",HaloFile[i]);
   }
  
  /* read output filenames from stdin */
  for(i=0; i<nFiles-1; i++)
     {
      fprintf(stderr,"Please give prefix for %5d. output file:                 ",i+1);
      scanf("%s", OutFile[i]);
      fprintf(stderr,"%s\n",OutFile[i]);
     }
  fprintf(stderr,"\n");
  
  
  /*======================================================================*
   *  CREATE CROSS-CORRELATIONS                                           *
   *======================================================================*/
  
  /* read the first file into memory */
  fprintf(stderr,"Startup:\n");
  read_particles(HaloFile[0], 0);
  particle_halo_mapping(0);
  fprintf(stderr,"\n");
 
  for(i=0; i<nFiles-1; i++) {
    
    /* be verbose */
    fprintf(stderr,"Correlating '%s' to '%s'\n           -> writing to '%s'\n",
            HaloFile[i],HaloFile[i+1],OutFile[i]);

    /* read the next file into memory */
    read_particles(HaloFile[i+1], 1);
    particle_halo_mapping(1);
    
    /* cross correlate HaloFile[i] to HaloFile[i+1] */
    cross_correlation(OutFile[i]);
    
#ifdef MERGER_RATIO
    assign_progenitors(OutFile[i]); //dumps information about progenitors
#endif
    
    /* be verbose */
    fprintf(stderr,"  o making file 1 the new file 0 ...");

    /* remove HaloFile[0] from memory */
    for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
      if(halos[0][ihalo].Pid   != NULL) {
        free(halos[0][ihalo].Pid);
        halos[0][ihalo].Pid = NULL;
      }
      if(halos[0][ihalo].mtree != NULL) {
        free(halos[0][ihalo].mtree);
        halos[0][ihalo].mtree = NULL;
      }
    }
    for(ipart=0; ipart<PidMax[0]; ipart++) {
      if(parts[0][ipart].Hid != NULL) {
        free(parts[0][ipart].Hid);
        parts[0][ipart].Hid = NULL;
      }
    }
    if(halos[0] != NULL) {
      free(halos[0]);
      halos[0] = NULL;
    }
    if(parts[0] != NULL) {
      free(parts[0]);
      parts[0] = NULL;
    }
    
    /* make HaloFile[i+1] the new HaloFile[i] */
    nHalos[0] = nHalos[1];
    halos[0]  = halos[1];
    parts[0]  = parts[1];
    
    /* be verbose */
    fprintf(stderr," done\n");
  } // for(nFiles)
  
  
  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
  fprintf(stderr,"\nCleaning up ... ");

  /* remove HaloFile[0] from memory */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    if(halos[0][ihalo].Pid   != NULL) free(halos[0][ihalo].Pid);
#ifdef MTREE_BOTH_WAYS
    if(halos[0][ihalo].mtree != NULL) free(halos[0][ihalo].mtree);
#endif
  }
  for(ipart=0; ipart<PidMax[0]; ipart++) {
    if(parts[0][ipart].Hid != NULL) free(parts[0][ipart].Hid);
  }
  if(halos[0] != NULL) free(halos[0]);
  if(parts[0] != NULL) free(parts[0]);

  /* remove filename storage */
  for(i=0; i<nFiles; i++)
   {
    if(HaloFile[i]) free(HaloFile[i]);
    if(OutFile[i])  free(OutFile[i]);
   }
  if(HaloFile) free(HaloFile);
  if(OutFile)  free(OutFile);
  
  printf("finished\n");
  return(1);
}


/*==================================================================================================
 * read_particles:
 *
 * read the file storing the particle IDs for each halo
 *
 *      nHalos = number of halos found in file
 *      Pid    = id's of all those particles
 *
 *==================================================================================================*/
int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin;
  char      line[MAXSTRING];
  int64_t   ihalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, jpart, Pid, Ptype, numGoodHalos, haloid;
  uint64_t  PidMin_local=(1<<62);
  uint64_t  PidMax_local=0;
  time_t    elapsed = (time_t)0;

  elapsed -= time(NULL);
  fprintf(stderr,"  o reading file %s ...",filename);
  
  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
   }
  
  /* reset all variables */
  nHalos[isimu] = 0;
  ihalo         = -1;
  halos[isimu]  = NULL;
  
  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
  
  /* for AHF_particles files the first line is numGoodHalos which we can happily ignore */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
    fgets(line,MAXSTRING,fpin);  
  
  do {
    if(strncmp(line,"#",1) != 0)
     {
      /* has a haloid been written */
      if(sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
       {
        /* if not, just get the number of particles */
        sscanf(line,"%"SCNu64, &nPartInHalo);
        
        /* and use halo counter as id */
        haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
       }
#ifdef USE_LINENUMBER_AS_HALOID
      haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
#endif
      
      /* found yet another halo */
      ihalo++;
      nHalos[isimu] += 1;
      halos[isimu]   = (HALOptr) realloc(halos[isimu], nHalos[isimu]*sizeof(HALOS));
      
      /* store haloid */
      halos[isimu][ihalo].haloid = haloid;
      
      /* halos[][].Pid will be incrementally filled using realloc() */
      halos[isimu][ihalo].Pid = NULL;
      
      /* read all their id's */
      nPartInUse = 0;
      for(ipart=0; ipart<nPartInHalo; ipart++)
       {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
         }
        else if(Ptype < 0 || Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
         }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
        if(Ptype == ONLY_USE_PTYPE)
#endif
         {
          halos[isimu][ihalo].Pid             = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }
          halos[isimu][ihalo].Pid[nPartInUse] = Pid;
          
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] =       (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin =       (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local = (halos[isimu][ihalo].Pid[nPartInUse]);
          
          nPartInUse++;
         }
       }
      
      /* store number of particles in halo */
      halos[isimu][ihalo].npart = nPartInUse;      
     }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
  
  fclose(fpin);
  
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec. (full ID range = %"PRIu64" -> %"PRIu64", local ID range = %"PRIu64" -> %"PRIu64")\n",
          elapsed,PidMin,PidMax[isimu],PidMin_local,PidMax_local);
  
  return(1);
}


/*==================================================================================================
 * particle_halo_mapping:
 *
 *  for each particle remember to which halo(s) it belongs
 *
 *==================================================================================================*/
int particle_halo_mapping(int isimu)
{
  int64_t  ihalo;         // the downwards for-loop is running until ihalo=-1
  uint64_t ipart, jpart;
  time_t   elapsed = (time_t)0;

  elapsed -= time(NULL);
  fprintf(stderr,"  o creating particle<->halo mapping for file %d (PidMax=%"PRIu64") ...",isimu,PidMax[isimu]);
  
  parts[isimu] = (PARTptr) calloc(PidMax[isimu]+1, sizeof(PARTS));
  
  /* recording every halo it belongs to: running from low to high mass objects to allow for unique assignment! */
#ifdef WITH_OPENMP2
#pragma omp parallel for schedule(dynamic) private(ihalo,jpart,ipart) shared(nHalos,halos,parts,isimu)
#endif
  for(ihalo=nHalos[isimu]-1; ihalo>=0; ihalo--)
   {
    for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
     {
      ipart = halos[isimu][ihalo].Pid[jpart];
      
#ifdef EXCLUSIVE_PARTICLES
      if(parts[isimu][ipart].nhalos == 0)
#endif
       {
        parts[isimu][ipart].nhalos++;
        parts[isimu][ipart].Hid = (uint64_t *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(uint64_t));
        
        parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
       }
      
     }
   }
  
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);
  return(1);
}

/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(char OutFile[MAXSTRING])
{
  uint64_t  ihalo;
  time_t   elapsed = (time_t)0;

  
  /*---------------------------------------------------------
   * backwards correlation
   *---------------------------------------------------------*/
  elapsed -= time(NULL);
  fprintf(stderr,"  o generating cross-correlation 0->1 for %"PRIu64" haloes ...",nHalos[0]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* cross-correlation simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    create_mtree(ihalo, 0, 1);
  }
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);
  
  
  
  
  /*---------------------------------------------------------
   * forward correlation
   *---------------------------------------------------------*/
#ifdef MTREE_BOTH_WAYS
  elapsed -= time(NULL);
  fprintf(stderr,"  o generating cross-correlation 1->0 for %"PRIu64" haloes ...",nHalos[1]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* cross-correlation simu1->simu0 */
  for(ihalo=0; ihalo<nHalos[1]; ihalo++) {
    create_mtree(ihalo, 1, 0);
  }
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);
  
  
  
  elapsed -= time(NULL);
  fprintf(stderr,"  o removing network connections");
#ifdef WITH_OPENMP2
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* clean connections simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    clean_connection(ihalo, 0, 1);
  }
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);

#endif // MTREE_BOTH_WAYS

  
  
  write_mtree(OutFile);
  
  return(1);
}

/*==================================================================================================
 * clean_connection
 *==================================================================================================*/
int clean_connection(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t jhalo, icroco, ncroco_new;
  MTREE    *mtree;
#ifdef DEBUG
  uint64_t idesc;
#endif

  /* count number of new crocos */
  ncroco_new = 0;
  mtree      = NULL;
  
  /* loop over all cross-correlated haloes */
  for(icroco=0; icroco<halos[isimu0][ihalo].ncroco; icroco++) {
    jhalo = halos[isimu0][ihalo].mtree[icroco].id[1];
    
    /* check whether the present halo is the most likely descendant of this progenitor */
    if(max_merit(jhalo, isimu1) == ihalo) {
      // keep jhalo in mtree-list
      ncroco_new++;
      mtree = (MTREEptr) realloc(mtree, ncroco_new*sizeof(MTREE));
      
      mtree[ncroco_new-1].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
      mtree[ncroco_new-1].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
      mtree[ncroco_new-1].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
      mtree[ncroco_new-1].common    = halos[isimu0][ihalo].mtree[icroco].common;
      mtree[ncroco_new-1].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
      mtree[ncroco_new-1].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
      mtree[ncroco_new-1].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];
      
#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is     a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
#endif
    }
    else {
      // remove jhalo from mtree-list and henc do not add it to the new mtree[] list
#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is NOT a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
      for(idesc=0; idesc<halos[isimu1][jhalo].ncroco; idesc++) {
        fprintf(stderr,"    idesc=%ld haloidesc=%ld\n",idesc,halos[isimu1][jhalo].mtree[idesc].haloid[1]);
      }
#endif
    }
  } // for(icroco)
  
  /* replace halos[isimu0][ihalo].mtree[] with new structure array */
  if(halos[isimu0][ihalo].mtree) free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].ncroco = ncroco_new;
  halos[isimu0][ihalo].mtree  = mtree;
}

/*==================================================================================================
 * max_merit
 *==================================================================================================*/
uint64_t max_merit(uint64_t jhalo, int isimu)
{
  uint64_t ihalo;
  
  /* mtree[] is ordered by merit and hence we only need to check the first entry */
  if(halos[isimu][jhalo].ncroco > 0) {
    return(halos[isimu][jhalo].mtree[0].id[1]);
  }
  else {
#ifdef DEBUG
    fprintf(stderr,"jhalo=%ld in isimu=%d does not point to anywhere!?\n",jhalo,isimu);
#endif
    return(0);
  }
}

/*==================================================================================================
 * create_mtree
 *==================================================================================================*/
int create_mtree(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  uint64_t *common;
  MTREE    *mtree;
  
  double *merit;
  long unsigned *idx;
 
  /* temporary array pointers */
  halos[isimu0][ihalo].mtree = NULL;
  mtree  = NULL;
  common = NULL;
  merit  = NULL;
  idx    = NULL;
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
  }
    
  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }
  halos[isimu0][ihalo].ncroco = ncroco;
    
  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) {
    
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(ncroco, sizeof(MTREE));
    mtree  = (MTREEptr)        calloc(ncroco, sizeof(MTREE));
    idx    = (long unsigned *) calloc(ncroco, sizeof(long unsigned));
    merit  = (double *)        calloc(ncroco, sizeof(double));
    
    /* store cross-correlations temporarily for sorting by indexx() */
    icroco = 0;
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
      if(common[khalo] > MINCOMMON){
        mtree[icroco].id[0]     = ihalo;
        mtree[icroco].haloid[0] = halos[isimu0][ihalo].haloid;
        mtree[icroco].npart[0]  = halos[isimu0][ihalo].npart;
        mtree[icroco].common    = common[khalo];
        mtree[icroco].id[1]     = khalo;
        mtree[icroco].haloid[1] = halos[isimu1][khalo].haloid;
        mtree[icroco].npart[1]  = halos[isimu1][khalo].npart;
        
        merit[icroco] = pow2((double)mtree[icroco].common)/((double)mtree[icroco].npart[0]*(double)mtree[icroco].npart[1]);
        icroco++;
      }
    }
    
    /* order by merit function */
    indexx((long unsigned)ncroco, merit-1, idx-1);
    
    /* store in halos[isimu0][*].mtree structure (descending order!) */
    for(jcroco=0; jcroco<ncroco; jcroco++) {
      icroco = idx[ncroco-1-jcroco]-1;
      
      /* store mtree[] inside halos[][] structure in correct order */
      halos[isimu0][ihalo].mtree[jcroco].id[0]    = mtree[icroco].id[0];
      halos[isimu0][ihalo].mtree[jcroco].haloid[0]= mtree[icroco].haloid[0];
      halos[isimu0][ihalo].mtree[jcroco].npart[0] = mtree[icroco].npart[0];
      halos[isimu0][ihalo].mtree[jcroco].common   = mtree[icroco].common;
      halos[isimu0][ihalo].mtree[jcroco].id[1]    = mtree[icroco].id[1];
      halos[isimu0][ihalo].mtree[jcroco].haloid[1]= mtree[icroco].haloid[1];
      halos[isimu0][ihalo].mtree[jcroco].npart[1] = mtree[icroco].npart[1];
    }
    
    /* free temporary structures */
    if(mtree) {
      free(mtree);
      mtree = NULL;
    }
    if(idx) {
      free(idx);
      idx = NULL;
    }
    if(merit) {
      free(merit);
      merit = NULL;
    }
  } // if(ncroco)
  
  /* free temporary structures */
  if(common) {
    free(common);
    common = NULL;
  }
}

/*==================================================================================================
 * write_mtree:
 *==================================================================================================*/
int write_mtree(char OutFile[MAXSTRING])
{
  uint64_t  ihalo;
  int64_t   icroco;
  FILE *fpout, *fpout_idx;
  char outname[MAXSTRING], outname_idx[MAXSTRING];
  time_t   elapsed = (time_t)0;
  
  elapsed -= time(NULL);
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes ...",nHalos[0]);
  
  sprintf(outname,"%s_mtree",OutFile);
  strcpy(outname_idx, outname);
  strcat(outname_idx, "_idx");
  
  fpout = fopen(outname,"w");
  if(fpout == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname);
    exit(0);
  }
  
  fpout_idx = fopen(outname_idx,"w");
  if(fpout_idx == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname_idx);
    exit(0);
  }
  
  
#ifdef SUSSING2013
  fprintf(fpout,"%"PRIu64"\n",nHalos[0]);
#else // SUSSING2013
  fprintf(fpout,"#   HaloID(1)   HaloPart(2)  NumProgenitors(3)\n");
  fprintf(fpout,"#      SharedPart(1)    HaloID(2)   HaloPart(3)\n");
  fprintf(fpout_idx,"# HaloID(1) HaloID(2)\n");
#endif // SUSSING2013
  fflush(fpout);
  fflush(fpout_idx);

  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    if(halos[0][ihalo].ncroco > 0) {
      //      this is the old format where the haloid corresponds to the linenumber
      //      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].id[1]);
      
      // this is the case where we use the haloid as found in *_particles
      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
      
#ifdef SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);
#else // SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].npart,
              halos[0][ihalo].ncroco);
#endif // SUSSING2013
      fflush(fpout);
      
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
#ifdef SUSSING2013
        fprintf(fpout,"%"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].haloid[1]);
#else // SUSSING2013
        fprintf(fpout,"  %"PRIu64"  %"PRIu64"  %"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].common,
                halos[0][ihalo].mtree[icroco].haloid[1],
                halos[0][ihalo].mtree[icroco].npart[1]);
#endif // SUSSING2013
        fflush(fpout);
      }
    }
#ifdef SUSSING2013
    else {
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);      
    }
#endif // SUSSING2013
  }
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);

  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);
  return(1);
}

#ifdef MERGER_RATIO
/*==================================================================================================
 * assign_progenitors:
 *
 *  get statistics for multiple progenitors (e.g. mass ratios, etc.)
 *
 *  update 10/10/2007:
 *       each subhalo shares the most particles with its host :-(
 *       -> tried to fix this issue..
 *
 *==================================================================================================*/
int assign_progenitors(char OutFile[MAXSTRING])
{
  FILE   *fpin, *fpout, *fpout_merger;
  char    OutFile_mtree[MAXSTRING], OutFile_idx[MAXSTRING], OutFile_merger[MAXSTRING], line[MAXSTRING];
  long unsigned   *id1, *npart1, *common, *id2, *npart2, *idx, iline;
  double  xcommon, xnpart1, xnpart2, *ratio;
  long    prev_id1, nprog, iprog;
  
  fprintf(stderr,"  o assigning progenitors ...");
  
  strcpy(OutFile_mtree, OutFile);
  strcat(OutFile_mtree, "_mtree");
  strcpy(OutFile_idx, OutFile_mtree);
  strcat(OutFile_idx, "_progs");
  strcpy(OutFile_merger, OutFile_mtree);
  strcat(OutFile_merger, "_merger");
  
  /* open output file */
  fpout = fopen(OutFile_idx,"w");
  if(fpout == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_idx);
    exit(0);
   }
  fprintf(fpout,"#id(1) Np(2) iprog1(3) Np1(4) ncommon1(5) iprog2(6) Np2(7) ncommon2(8) iprog3(9) Np3(10) ncommon3(11)\n");
  
  /* open output file */
  fpout_merger = fopen(OutFile_merger,"w");
  if(fpout_merger == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_merger);
    exit(0);
   }
  fprintf(fpout_merger,"#id(1) iprog1(2) iprog2(3) common2/ncommon1(4)\n");

  /* read *_mtree file */
  read_mtree(OutFile);
  
  prev_id1  = mtree[0].id[0];  
  nprog     = 0;
  id1    = (long *) calloc(1, sizeof(long));
  id2    = (long *) calloc(1, sizeof(long));
  npart1 = (long *) calloc(1, sizeof(long));
  npart2 = (long *) calloc(1, sizeof(long));
  common = (long *) calloc(1, sizeof(long));
  
  for(iline=0; iline<nlines; iline++)
   {
    if(mtree[iline].id[0] == prev_id1)
     {
      /* make room for one more progenitor */
      nprog++;
      id1    = (long *) realloc(id1,    (nprog)*sizeof(long));
      id2    = (long *) realloc(id2,    (nprog)*sizeof(long));
      npart1 = (long *) realloc(npart1, (nprog)*sizeof(long));
      npart2 = (long *) realloc(npart2, (nprog)*sizeof(long));
      common = (long *) realloc(common, (nprog)*sizeof(long));
      
      /* copy information from mtree[] */
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
     }
    else
     {
      ratio = (double *)        calloc(nprog+1, sizeof(double));
      idx   = (long unsigned *) calloc(nprog+1, sizeof(long unsigned));
      
      for(iprog=0; iprog<nprog; iprog++)
       {
        /* calculate progenitor criterion */
        xcommon          = (double)common[iprog];
        xnpart1          = (double)npart1[iprog];
        xnpart2          = (double)npart2[iprog];
        ratio[iprog]     = pow2(xcommon)/(xnpart1*xnpart2);
       }
      /* sort all progenitor by merit function */
      indexx(nprog, ratio-1, idx-1);
      
      /*-----------------------------------------------------
       * dump information to file
       *-----------------------------------------------------*/
      if(nprog > 2)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1],
                id2[idx[nprog-3]-1],npart2[idx[nprog-3]-1],common[idx[nprog-3]-1]);
        
        /* iprog2 is most credible second progenitor */
        if(id2[idx[nprog-1]-1]<id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* iprog3 is most credible second progenitor */
        else if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-3]-1],
                    (double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* nothing else to do as only one real progenitor exists */
       }
      
      else if (nprog > 1)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1]);            

        /* this indicates that the halo itself is a subhalo */
        if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1])
         {
          /* nothing to do as only one real progenitor exists */
         }
        else
         {
          /* iprog2 is most credible second progenitor */
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
       }
      else
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       -1 -1 -1       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1]);
        
        /* nothing else to do as there are no multiple progenitors */
       }
      
      
      free(ratio);
      free(idx);
      
      /* start a new progenitor list */
      free(id1);
      free(id2);
      free(npart1);
      free(npart2);
      free(common);
      
      nprog  = 1;
      id1    = (long *) calloc(nprog, sizeof(long));
      id2    = (long *) calloc(nprog, sizeof(long));
      npart1 = (long *) calloc(nprog, sizeof(long));
      npart2 = (long *) calloc(nprog, sizeof(long));
      common = (long *) calloc(nprog, sizeof(long));
      
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
      
      prev_id1        = id1[nprog-1];
      
     }
   }
  
  if(id1 != NULL) free(id1);
  if(id2 != NULL) free(id2);
  if(npart1 != NULL) free(npart1);
  if(npart2 != NULL) free(npart2);
  if(common != NULL) free(common);
  
  fclose(fpout);
  fclose(fpout_merger);
  
  fprintf(stderr," done\n");
  
  return(1);
}

/*==================================================================================================
 * read_mtree:
 *
 *       simply reads in the *_mtree file and 
 *       puts it into the array of structures mtree[iline].XYZ
 *
 * Note: at this stage we just treat these entries as "lines" -> no connection to halos yet!
 *
 *==================================================================================================*/
void read_mtree(char *prefix)
{
  long unsigned iline;
  char          line[MAXSTRING], outname[MAXSTRING];
  FILE         *fpin;
  
  sprintf(outname,"%s_mtree",prefix);
  if((fpin = fopen(outname,"r")) == NULL)
   {
    fprintf(stderr,"cannot open  %s\nEXIT\n",outname);
    exit(0);
   }
  
  // count number of lines
  nlines = 0;
  while(!feof(fpin))
   {
    nlines++;
    fgets(line,MAXSTRING,fpin);
   }
  
  //fprintf(stderr,"  o found %ld lines in:  %s\n",nlines,infile);
  
  // allocate memory
  mtree = (MTREEptr) realloc((MTREEptr)mtree, nlines*sizeof(MTREE));
  
  // actually read the file
  rewind(fpin);
  for(iline=0; iline<nlines; iline++)
   {
    // read next line from file
    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64,
           &(mtree[iline].id[0]), 
           &(mtree[iline].npart[0]), 
           &(mtree[iline].common), 
           &(mtree[iline].id[1]), 
           &(mtree[iline].npart[1]));
    //fprintf(stderr,"iline=%ld id[0]=%ld\n",iline,mtree[iline].id[0]);
   }
  
  fclose(fpin);
}

#endif // MERGER_RATIO
